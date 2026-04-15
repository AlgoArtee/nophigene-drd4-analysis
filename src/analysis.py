#!/usr/bin/env python3
"""
DRD4 Gene Analysis Pipeline

Usage examples:
    python src/analysis.py \
        --vcf data/GFXC926398.filtered.snp.vcf.gz \
        --idat data/202277800037_R01C01 \
        --out results/drd4_report.html \
        --region 11:637269-640706
"""

from __future__ import annotations

import argparse
import html
import json
import logging
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import allel
import pandas as pd
from methylprep import run_pipeline

try:
    from .helper_functions.filter_manifest_region import (
        sanitize_gene_name_for_filename,
        save_filtered_manifest,
    )
except ImportError:
    from helper_functions.filter_manifest_region import sanitize_gene_name_for_filename, save_filtered_manifest

DEFAULT_REGION = "11:637269-640706"
DEFAULT_REPORT_NAME = "drd4_report.html"
DEFAULT_GENE_NAME = "DRD4"
GENE_DATA_DIR = Path(__file__).resolve().parent / "gene_data"
INTERPRETATION_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "drd4_interpretation_db.json"
POPULATION_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "drd4_population_db.json"

# Configure the root logger once so both CLI and web runs stream progress.
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def _candidate_gene_database_paths(gene_name: str, suffix: str) -> list[Path]:
    """Return likely per-gene knowledge-base paths in priority order."""
    sanitized_gene_name = sanitize_gene_name_for_filename(gene_name)
    candidates = [
        GENE_DATA_DIR / f"{sanitized_gene_name.lower()}_{suffix}",
        GENE_DATA_DIR / f"{sanitized_gene_name}_{suffix}",
        GENE_DATA_DIR / f"{sanitized_gene_name.upper()}_{suffix}",
    ]

    unique_candidates: list[Path] = []
    seen_paths: set[Path] = set()
    for path in candidates:
        if path in seen_paths:
            continue
        seen_paths.add(path)
        unique_candidates.append(path)
    return unique_candidates


def find_gene_database_path(gene_name: str, suffix: str) -> Path | None:
    """Locate a bundled per-gene database when one is available."""
    for candidate in _candidate_gene_database_paths(gene_name, suffix):
        if candidate.exists():
            return candidate
    return None


class AnalysisError(RuntimeError):
    """Raised when the DRD4 workflow cannot complete successfully."""


@dataclass
class AnalysisResult:
    """Structured result returned by :func:`run_analysis`.

    Attributes
    ----------
    variants : pd.DataFrame
        Regional PASS variants loaded from the source VCF.
    methylation : pd.DataFrame
        Probe-level methylation table after joining with the curated DRD4
        manifest subset.
    popstats : Any | None
        Optional population statistics payload loaded from a user-supplied CSV
        or JSON file.
    report_path : Path
        Path to the generated output report.
    methylation_output_path : Path
        Path to the exported methylation CSV companion file.
    region : str
        Genomic interval used for the run.
    vcf_path : Path
        Input VCF path used during execution.
    idat_base : Path
        Input IDAT prefix used during execution.
    variant_interpretations : dict[str, Any]
        Curated DRD4 variant interpretations matched against the loaded PASS
        variants using the local interpretation database.
    methylation_insights : dict[str, Any]
        Gene-level methylation interpretation assembled from the current probe
        table plus the local interpretation database.
    knowledge_base : dict[str, Any]
        Parsed local interpretation database used to build the biological and
        clinical insights shown in the UI.
    population_insights : dict[str, Any]
        Population-frequency and geography summaries assembled from the local
        DRD4 population database.
    population_database : dict[str, Any]
        Parsed local population database used to add location-based frequency
        context for common DRD4 variants.
    """

    variants: pd.DataFrame
    methylation: pd.DataFrame
    popstats: Any | None
    report_path: Path
    methylation_output_path: Path
    region: str
    vcf_path: Path
    idat_base: Path
    variant_interpretations: dict[str, Any]
    methylation_insights: dict[str, Any]
    knowledge_base: dict[str, Any]
    population_insights: dict[str, Any]
    population_database: dict[str, Any]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the DRD4 analysis workflow.

    Parameters
    ----------
    argv : list[str] | None, optional
        Argument list to parse. When omitted, argparse reads directly from
        ``sys.argv``. Accepting an explicit list lets the Docker launcher reuse
        the CLI parser without shelling out to a subprocess.

    Returns
    -------
    argparse.Namespace
        Parsed arguments ready to be consumed by :func:`main`.
    """
    parser = argparse.ArgumentParser(
        description="DRD4 gene analysis: variants, methylation, population stats, and report generation."
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Tabix-indexed VCF file with filtered SNPs (for example, *.vcf.gz).",
    )
    parser.add_argument(
        "--idat",
        required=True,
        help="Path prefix to Illumina IDATs, without the _Grn/_Red suffixes.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Path for the output report (HTML, CSV, or JSON).",
    )
    parser.add_argument(
        "--region",
        default=DEFAULT_REGION,
        help="Genomic region in chr:start-end format. Defaults to the DRD4 GRCh37 interval.",
    )
    parser.add_argument(
        "--popstats",
        default=None,
        help="Optional JSON or CSV with population-frequency data.",
    )
    parser.add_argument(
        "--manifest-file",
        default=None,
        help="Optional manifest file passed through to methylprep.",
    )

    return parser.parse_args(argv)


def _serialize_popstats(popstats: Any | None) -> Any:
    """Convert population statistics to a JSON-friendly representation."""
    if popstats is None:
        return None
    if isinstance(popstats, pd.DataFrame):
        return popstats.to_dict(orient="records")
    return popstats


def _derive_methylation_output_path(output_path: str | Path) -> Path:
    """Derive the companion methylation CSV path from the requested report path."""
    report_path = Path(output_path)
    stem = report_path.stem if report_path.suffix else report_path.name
    return report_path.with_name(f"{stem}_methylation.csv")


def load_interpretation_database(database_path: str | Path = INTERPRETATION_DB_PATH) -> dict[str, Any]:
    """Load the curated local DRD4 interpretation database.

    Parameters
    ----------
    database_path : str | Path, optional
        Path to the JSON-backed local knowledge base that stores DRD4 variant
        and methylation interpretations.

    Returns
    -------
    dict[str, Any]
        Parsed JSON payload describing gene context plus curated variant
        records.

    Raises
    ------
    AnalysisError
        Raised when the JSON file cannot be found or parsed.
    """
    payload_path = Path(database_path)
    if not payload_path.exists():
        raise AnalysisError(f"Interpretation database not found: {payload_path}")

    try:
        return json.loads(payload_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AnalysisError(f"Interpretation database is not valid JSON: {payload_path}") from exc


def load_population_database(database_path: str | Path = POPULATION_DB_PATH) -> dict[str, Any]:
    """Load the curated local DRD4 population database."""
    payload_path = Path(database_path)
    if not payload_path.exists():
        raise AnalysisError(f"Population database not found: {payload_path}")

    try:
        return json.loads(payload_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AnalysisError(f"Population database is not valid JSON: {payload_path}") from exc


def load_gene_interpretation_database(gene_name: str) -> dict[str, Any] | None:
    """Load a bundled gene-specific interpretation database when one exists."""
    database_path = find_gene_database_path(gene_name, "interpretation_db.json")
    if database_path is None:
        return None
    return load_interpretation_database(database_path)


def load_gene_population_database(gene_name: str) -> dict[str, Any] | None:
    """Load a bundled gene-specific population database when one exists."""
    database_path = find_gene_database_path(gene_name, "population_db.json")
    if database_path is None:
        return None
    return load_population_database(database_path)


def _decode_scalar(value: Any) -> Any:
    """Decode byte-valued VCF fields into plain Python strings when needed."""
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8")
    return value


def _normalize_lookup_key(value: str) -> str:
    """Canonicalize lookup keys so IDs and coordinate aliases can match reliably."""
    normalized = value.strip().lower().replace(" ", "")
    if normalized.startswith("chr"):
        normalized = normalized[3:]
    return normalized


def _parse_region_string(region: str) -> dict[str, Any]:
    """Parse a ``chr:start-end`` style region string into normalized components."""
    cleaned_region = region.strip().replace(",", "")
    match = re.fullmatch(r"(?:chr)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)", cleaned_region)
    if match is None:
        raise AnalysisError(
            f"Unsupported region format '{region}'. Use chr:start-end, for example 11:637269-640706."
        )

    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        raise AnalysisError(f"Invalid region '{region}': start must be <= end.")

    return {
        "chrom": match.group("chrom"),
        "start": start,
        "end": end,
    }


def _format_interval(chrom: str, start: int, end: int) -> str:
    """Format a genomic interval for human-readable summaries."""
    return f"chr{chrom}:{start:,}-{end:,}"


def _format_point(chrom: str, position: int) -> str:
    """Format a single genomic coordinate for human-readable summaries."""
    return f"chr{chrom}:{position:,}"


def _format_frequency(value: float | None) -> str:
    """Render an allele frequency as a percentage string for the UI."""
    if value is None:
        return "n/a"
    return f"{value * 100:.1f}%"


def _intervals_overlap(first: dict[str, Any], second: dict[str, Any]) -> bool:
    """Return whether two intervals on the same chromosome overlap."""
    first_chrom = str(first["chrom"]).removeprefix("chr")
    second_chrom = str(second["chrom"]).removeprefix("chr")
    if first_chrom != second_chrom:
        return False
    return first["start"] <= second["end"] and second["start"] <= first["end"]


def _build_interval_record(label: str, chrom: str, start: int, end: int, definition: str) -> dict[str, Any]:
    """Create a normalized interval record used in the UI summaries."""
    return {
        "label": label,
        "chrom": chrom,
        "start": start,
        "end": end,
        "length_bp": end - start + 1,
        "display": _format_interval(chrom, start, end),
        "definition": definition,
    }


def _build_variant_lookup_keys(row: pd.Series) -> set[str]:
    """Create rsID and coordinate aliases for a single variant row."""
    keys: set[str] = set()

    raw_id = row.get("id")
    if pd.notna(raw_id) and str(raw_id).strip() and str(raw_id).strip() != ".":
        keys.update(_normalize_lookup_key(part) for part in str(raw_id).split(";") if part.strip())

    chrom = str(row.get("chrom", "")).strip()
    pos = row.get("pos")
    ref = str(row.get("ref", "") or "").strip().upper()
    alt = str(row.get("alt", "") or "").strip().upper()

    if chrom and pd.notna(pos):
        chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
        pos_text = str(int(pos))
        keys.add(_normalize_lookup_key(f"{chrom}:{pos_text}"))

        if ref and alt and ref not in {"NAN", "NONE"} and alt not in {"NAN", "NONE"}:
            keys.add(_normalize_lookup_key(f"{chrom}:{pos_text}:{ref}>{alt}"))

    return keys


def _format_variant_display(row: pd.Series) -> str:
    """Render a human-readable label for a variant row in the UI."""
    raw_id = row.get("id")
    if pd.notna(raw_id) and str(raw_id).strip() and str(raw_id).strip() != ".":
        return str(raw_id)
    return f"{row['chrom']}:{int(row['pos'])} {row['ref']}>{row['alt']}"


def _match_variant_record(row: pd.Series, variant_records: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Return the first curated record that matches the observed variant row."""
    lookup_keys = _build_variant_lookup_keys(row)
    for record in variant_records:
        record_keys = {
            _normalize_lookup_key(candidate)
            for candidate in record.get("lookup_keys", [])
            if candidate
        }
        if lookup_keys & record_keys:
            return record
    return None


def _build_known_variant_summary(record: dict[str, Any]) -> dict[str, Any]:
    """Project a curated database record into a compact UI-friendly summary."""
    position = record.get("position")
    chrom = record.get("chromosome", "11")
    assay_note = None
    if not record.get("is_assayable_in_snp_vcf", True):
        assay_note = "Not directly called by the current SNP-oriented VCF preview."

    return {
        "variant": record.get("display_name", record["variant"]),
        "common_name": record.get("common_name"),
        "position": _format_point(chrom, int(position)) if position is not None else "Repeat / structural locus",
        "clinical_significance": record.get("clinical_significance", "Clinical significance not specified."),
        "summary": record.get("clinical_interpretation", ""),
        "interpretation_scope": record.get("interpretation_scope", "Research context"),
        "functional_effects": record.get("functional_effects", []),
        "associated_conditions": record.get("associated_conditions", []),
        "research_context": record.get("research_context", []),
        "literature_findings": _build_literature_findings(record),
        "assay_note": assay_note,
        "evidence": record.get("evidence", []),
    }


def _build_observed_variant_summary(
    row: pd.Series,
    matched_record: dict[str, Any] | None,
    *,
    gene_name: str = DEFAULT_GENE_NAME,
) -> dict[str, Any]:
    """Render an observed variant plus whatever curated significance is available."""
    chrom = str(row.get("chrom", "")).removeprefix("chr")
    position = int(row["pos"])
    quality = row.get("qual")
    quality_display = f"{float(quality):.2f}" if pd.notna(quality) else "n/a"

    if matched_record is None:
        return {
            "display": _format_variant_display(row),
            "position": _format_point(chrom, position),
            "quality": quality_display,
            "matched_variant": None,
            "clinical_significance": (
                f"No curated clinical significance available in the local {gene_name} database for this observed site."
            ),
            "summary": (
                f"Observed inside the selected {gene_name} search window, but not one of the seeded "
                f"common {gene_name} promoter or gene variants."
            ),
            "functional_effects": [],
            "associated_conditions": [],
            "research_context": [],
            "literature_findings": [],
            "evidence": [],
        }

    return {
        "display": _format_variant_display(row),
        "position": _format_point(chrom, position),
        "quality": quality_display,
        "matched_variant": matched_record.get("display_name", matched_record["variant"]),
        "clinical_significance": matched_record.get(
            "clinical_significance", "Clinical significance not specified."
        ),
        "summary": matched_record.get("clinical_interpretation", ""),
        "functional_effects": matched_record.get("functional_effects", []),
        "associated_conditions": matched_record.get("associated_conditions", []),
        "research_context": matched_record.get("research_context", []),
        "literature_findings": _build_literature_findings(matched_record),
        "evidence": matched_record.get("evidence", []),
    }


def _build_literature_findings(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Normalize literature findings so the UI can render paper-level bullets consistently."""
    findings: list[dict[str, Any]] = []
    for finding in record.get("literature_findings", []):
        paper = str(finding.get("paper", "")).strip()
        phenotype = str(finding.get("phenotype", "")).strip()
        finding_text = str(finding.get("finding", "")).strip()
        if not paper or not finding_text:
            continue
        findings.append(
            {
                "paper": paper,
                "genotypes": str(finding.get("genotypes", "")).strip(),
                "phenotype": phenotype,
                "finding": finding_text,
                "url": str(finding.get("url", "")).strip(),
            }
        )
    return findings


def _summarize_observed_region_variants(
    variants: pd.DataFrame,
    variant_records: list[dict[str, Any]],
    *,
    chrom: str,
    start: int,
    end: int,
    gene_name: str,
    limit: int = 12,
) -> list[dict[str, Any]]:
    """Summarize observed PASS variants in a given region for UI display."""
    observed = variants[
        (variants["chrom"].astype(str).str.removeprefix("chr") == chrom)
        & (variants["pos"] >= start)
        & (variants["pos"] <= end)
    ].sort_values("pos")

    summaries: list[dict[str, Any]] = []
    for _, row in observed.head(limit).iterrows():
        summaries.append(
            _build_observed_variant_summary(
                row,
                _match_variant_record(row, variant_records),
                gene_name=gene_name,
            )
        )
    return summaries


def _build_region_variant_analysis(
    *,
    region_record: dict[str, Any],
    search_region: dict[str, Any],
    variants: pd.DataFrame,
    curated_records: list[dict[str, Any]],
    gene_name: str,
    inclusion_hint: str,
) -> dict[str, Any]:
    """Build a structured analysis block for promoter or gene-body coverage."""
    included = _intervals_overlap(search_region, region_record)
    found_variants = (
        _summarize_observed_region_variants(
            variants,
            curated_records,
            chrom=region_record["chrom"],
            start=region_record["start"],
            end=region_record["end"],
            gene_name=gene_name,
        )
        if included
        else []
    )

    if included and found_variants:
        analysis_note = (
            f"The current search window overlaps {region_record['label'].lower()} "
            f"and found {len(found_variants)} PASS variant(s) in the first preview slice."
        )
    elif included:
        analysis_note = (
            f"The current search window overlaps {region_record['label'].lower()}, "
            "but no PASS variants from that region were retained in the current preview."
        )
    else:
        analysis_note = (
            f"The current search window does not overlap {region_record['label'].lower()}. "
            f"{inclusion_hint}"
        )

    return {
        "label": region_record["label"],
        "window": region_record["display"],
        "length_bp": region_record["length_bp"],
        "definition": region_record["definition"],
        "included": included,
        "analysis_note": analysis_note,
        "found_variant_count": len(found_variants),
        "found_variants": found_variants,
        "known_variants": [_build_known_variant_summary(record) for record in curated_records],
    }


def _build_population_frequency_rows(
    entries: list[dict[str, Any]],
    *,
    focus_alleles: list[str],
    effect_allele: str | None,
) -> list[dict[str, Any]]:
    """Normalize allele-frequency rows so the template can render them directly."""
    normalized_rows: list[dict[str, Any]] = []
    for entry in entries:
        allele_frequencies = entry.get("allele_frequencies", {})
        normalized_rows.append(
            {
                "population_code": entry.get("population_code"),
                "location_group": entry.get("location_group", "Unspecified"),
                "label": entry.get("label", entry.get("population_code", "Population")),
                "granularity": entry.get("granularity", "reference"),
                "effect_allele": effect_allele,
                "effect_allele_frequency": allele_frequencies.get(effect_allele) if effect_allele else None,
                "effect_allele_display": (
                    _format_frequency(allele_frequencies.get(effect_allele))
                    if effect_allele
                    else "n/a"
                ),
                "allele_frequencies": [
                    {
                        "allele": allele,
                        "frequency": allele_frequencies.get(allele),
                        "display": _format_frequency(allele_frequencies.get(allele)),
                    }
                    for allele in focus_alleles
                    if allele in allele_frequencies
                ],
            }
        )
    return normalized_rows


def _group_population_rows_by_location(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group population rows by geography for collapsible UI sections."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["location_group"], []).append(row)

    return [
        {
            "location_group": location_group,
            "entries": sorted(entries, key=lambda item: item["label"]),
        }
        for location_group, entries in sorted(grouped.items(), key=lambda item: item[0])
    ]


def _summarize_population_extremes(
    rows: list[dict[str, Any]], effect_allele: str | None
) -> dict[str, Any] | None:
    """Find the highest and lowest superpopulation frequencies for a focal allele."""
    if not effect_allele:
        return None

    comparable_rows = [
        row
        for row in rows
        if row.get("granularity") == "superpopulation" and row.get("effect_allele_frequency") is not None
    ]
    if not comparable_rows:
        return None

    highest = max(comparable_rows, key=lambda item: item["effect_allele_frequency"])
    lowest = min(comparable_rows, key=lambda item: item["effect_allele_frequency"])
    return {
        "effect_allele": effect_allele,
        "highest": highest,
        "lowest": lowest,
        "summary": (
            f"The {effect_allele} allele is most frequent in {highest['location_group']} "
            f"({highest['effect_allele_display']}) and least frequent in {lowest['location_group']} "
            f"({lowest['effect_allele_display']}) across the 1000 Genomes superpopulation layer."
        ),
    }


def _build_sample_variant_highlights(
    *,
    matched_records: list[dict[str, Any]],
    promoter_analysis: dict[str, Any],
    gene_analysis: dict[str, Any],
    gene_name: str,
) -> dict[str, Any]:
    """Summarize the most visible sample-specific variant findings first."""
    if matched_records:
        summary = (
            f"This sample matched {len(matched_records)} curated {gene_name} research variant(s). "
            "Those direct sample hits are shown first below."
        )
        items = [
            {
                "title": record["variant"],
                "observed_variant": record["observed_variant"],
                "category": record.get("interpretation_scope", "Research context"),
                "description": record.get("clinical_interpretation", ""),
                "conditions": record.get("associated_conditions", []),
                "literature_findings": record.get("literature_findings", []),
            }
            for record in matched_records
        ]
        return {
            "summary": summary,
            "highlight_items": items,
        }

    found_items: list[dict[str, Any]] = []
    for region_label, region_analysis in (
        ("Promoter review", promoter_analysis),
        ("Gene body review", gene_analysis),
    ):
        for record in region_analysis.get("found_variants", []):
            found_items.append(
                {
                    "title": record["display"],
                    "observed_variant": record["position"],
                    "category": region_label,
                    "description": record.get("summary", ""),
                    "conditions": record.get("associated_conditions", []),
                    "literature_findings": record.get("literature_findings", []),
                }
            )

    if found_items:
        summary = (
            f"This sample did not hit one of the curated named {gene_name} markers, but it did contain "
            f"{len(found_items)} observed PASS variant(s) inside the reviewed promoter or gene intervals."
        )
    else:
        summary = (
            f"This sample did not yield a visible {gene_name} promoter or gene-body PASS variant in the current preview slice."
        )

    return {
        "summary": summary,
        "highlight_items": found_items,
    }


def _build_region_recommendations(
    promoter_region_record: dict[str, Any],
    gene_region_record: dict[str, Any],
    combined_region: str,
    *,
    gene_name: str,
) -> list[dict[str, str]]:
    """Return practical region-span recommendations for common DRD4 review goals."""
    return [
        {
            "title": "Promoter only",
            "region": promoter_region_record["display"],
            "purpose": (
                "Use this when you want to focus on the upstream promoter-review window "
                f"and the classically studied {gene_name} promoter hotspot without loading the transcribed gene body."
            ),
        },
        {
            "title": "Gene body only",
            "region": gene_region_record["display"],
            "purpose": (
                f"Use this when you want the canonical {gene_name} transcribed interval but do not need the upstream promoter review window."
            ),
        },
        {
            "title": "Promoter plus gene body",
            "region": combined_region,
            "purpose": (
                f"Use this when you want the full audit: upstream promoter context plus the canonical {gene_name} gene interval in one search."
            ),
        },
    ]


def _categorize_beta(mean_beta: float | None) -> str:
    """Map average beta values to a coarse descriptive band for UI summaries."""
    if mean_beta is None or pd.isna(mean_beta):
        return "unavailable"
    if mean_beta < 0.20:
        return "low"
    if mean_beta < 0.60:
        return "intermediate"
    if mean_beta < 0.80:
        return "moderately high"
    return "high"


def build_variant_interpretations(
    variants: pd.DataFrame, knowledge_base: dict[str, Any], *, region: str
) -> dict[str, Any]:
    """Build a structured DRD4 locus audit for the observed PASS variants."""
    variant_records = knowledge_base.get("variant_records", [])
    gene_context = knowledge_base.get("gene_context", {})
    gene_name = str(gene_context.get("gene_name", DEFAULT_GENE_NAME))
    chrom = str(gene_context.get("chromosome", "11")).removeprefix("chr")
    gene_region_source = gene_context.get("gene_region", {})
    promoter_region_source = gene_context.get("promoter_review_region", {})
    promoter_hotspot_source = gene_context.get("promoter_hotspot_region", {})

    search_region = _parse_region_string(region)
    search_region_record = _build_interval_record(
        "Current search interval",
        str(search_region["chrom"]).removeprefix("chr"),
        int(search_region["start"]),
        int(search_region["end"]),
        "Exact interval used to pull PASS variants from the selected VCF.",
    )
    gene_region_record = _build_interval_record(
        gene_region_source.get("label", f"{gene_name} transcribed interval"),
        chrom,
        int(gene_region_source["start"]),
        int(gene_region_source["end"]),
        gene_region_source.get("definition", ""),
    )
    promoter_region_record = _build_interval_record(
        promoter_region_source.get("label", "Operational promoter review window"),
        chrom,
        int(promoter_region_source["start"]),
        int(promoter_region_source["end"]),
        promoter_region_source.get("definition", ""),
    )
    promoter_hotspot_record = _build_interval_record(
        promoter_hotspot_source.get("label", "Promoter polymorphism hotspot"),
        chrom,
        int(promoter_hotspot_source["start"]),
        int(promoter_hotspot_source["end"]),
        promoter_hotspot_source.get("definition", ""),
    )

    promoter_records = [
        record
        for record in variant_records
        if record.get("region_class") in {"promoter", "promoter_structural", "upstream_regulatory"}
    ]
    gene_records = [
        record
        for record in variant_records
        if record.get("region_class") in {"coding_repeat", "gene_body"}
    ]

    matched_records: list[dict[str, Any]] = []
    seen_matches: set[tuple[str, str]] = set()
    for _, row in variants.iterrows():
        matched_record = _match_variant_record(row, variant_records)
        if matched_record is None:
            continue

        observed_label = _format_variant_display(row)
        dedupe_key = (matched_record["variant"], observed_label)
        if dedupe_key in seen_matches:
            continue

        matched_records.append(
            {
                "observed_variant": observed_label,
                "variant": matched_record.get("display_name", matched_record["variant"]),
                "interpretation_scope": matched_record.get("interpretation_scope", "Research context"),
                "clinical_interpretation": matched_record.get("clinical_interpretation", ""),
                "clinical_significance": matched_record.get(
                    "clinical_significance", "Clinical significance not specified."
                ),
                "methylation_interpretation": matched_record.get("methylation_interpretation", ""),
                "functional_effects": matched_record.get("functional_effects", []),
                "associated_conditions": matched_record.get("associated_conditions", []),
                "research_context": matched_record.get("research_context", []),
                "literature_findings": _build_literature_findings(matched_record),
                "relevant_probe_ids": matched_record.get("relevant_methylation_probe_ids", []),
                "evidence": matched_record.get("evidence", []),
            }
        )
        seen_matches.add(dedupe_key)

    promoter_analysis = _build_region_variant_analysis(
        region_record=promoter_region_record,
        search_region=search_region,
        variants=variants,
        curated_records=promoter_records,
        gene_name=gene_name,
        inclusion_hint=(
            "Use the promoter-plus-gene interval "
            f"{gene_context.get('recommended_promoter_plus_gene_region', _format_interval(chrom, promoter_region_record['start'], gene_region_record['end']))} "
            f"if you want the upstream {gene_name} promoter reviewed alongside the gene."
        ),
    )
    gene_analysis = _build_region_variant_analysis(
        region_record=gene_region_record,
        search_region=search_region,
        variants=variants,
        curated_records=gene_records,
        gene_name=gene_name,
        inclusion_hint=f"Choose a region that overlaps the canonical {gene_name} gene interval if you want gene-body variants interpreted.",
    )

    promoter_phrase = "overlaps" if promoter_analysis["included"] else "does not overlap"
    gene_phrase = "overlaps" if gene_analysis["included"] else "does not overlap"
    combined_region = gene_context.get(
        "recommended_promoter_plus_gene_region",
        _format_interval(chrom, promoter_region_record["start"], gene_region_record["end"]),
    )
    summary = (
        f"{gene_name} is located at {gene_region_record['display']} on {gene_context.get('cytoband', 'the reported cytoband')} "
        f"and spans {gene_region_record['length_bp']:,} bp on the {gene_context.get('assembly', 'GRCh37 / hg19')} assembly. "
        f"The current search interval {search_region_record['display']} {promoter_phrase} the operational promoter review window "
        f"{promoter_region_record['display']} and {gene_phrase} the {gene_name} transcribed interval {gene_region_record['display']}. "
        f"In the current preview, {promoter_analysis['found_variant_count']} promoter-window variant(s) and "
        f"{gene_analysis['found_variant_count']} gene-interval variant(s) were surfaced."
    )

    return {
        "summary": summary,
        "matched_records": matched_records,
        "unclassified_variant_count": max(len(variants) - len(matched_records), 0),
        "gene_summary": gene_context.get("gene_summary", ""),
        "database_name": knowledge_base.get("database_name", f"Local {gene_name} interpretation database"),
        "gene_name": gene_name,
        "clinical_context": gene_context.get("clinical_context", ""),
        "variant_effect_overview": gene_context.get("variant_effect_overview", []),
        "condition_research_overview": gene_context.get("condition_research_overview", []),
        "sample_highlights": _build_sample_variant_highlights(
            matched_records=matched_records,
            promoter_analysis=promoter_analysis,
            gene_analysis=gene_analysis,
            gene_name=gene_name,
        ),
        "region_recommendations": _build_region_recommendations(
            promoter_region_record,
            gene_region_record,
            combined_region,
            gene_name=gene_name,
        ),
        "gene_region": gene_region_record,
        "promoter_region": promoter_region_record,
        "promoter_hotspot_region": promoter_hotspot_record,
        "search_region": search_region_record,
        "promoter_analysis": promoter_analysis,
        "gene_analysis": gene_analysis,
        "recommended_promoter_plus_gene_region": combined_region,
    }


def build_population_insights(
    variants: pd.DataFrame,
    knowledge_base: dict[str, Any],
    population_database: dict[str, Any],
) -> dict[str, Any]:
    """Build geography-aware population summaries for common curated gene variants."""
    variant_records = knowledge_base.get("variant_records", [])
    gene_name = str(knowledge_base.get("gene_context", {}).get("gene_name", DEFAULT_GENE_NAME))
    curated_record_map = {record["variant"]: record for record in variant_records}

    observed_variant_map: dict[str, list[str]] = {}
    for _, row in variants.iterrows():
        matched_record = _match_variant_record(row, variant_records)
        if matched_record is None:
            continue
        observed_variant_map.setdefault(matched_record["variant"], []).append(_format_variant_display(row))

    population_variant_records: list[dict[str, Any]] = []
    location_groups: set[str] = set()
    for record in population_database.get("variant_population_records", []):
        knowledge_record = curated_record_map.get(record["variant"], {})
        focus_alleles = list(record.get("focus_alleles", []))
        effect_allele = record.get("effect_allele")
        top_level_rows = _build_population_frequency_rows(
            record.get("top_level_location_frequencies", []),
            focus_alleles=focus_alleles,
            effect_allele=effect_allele,
        )
        detailed_rows = _build_population_frequency_rows(
            record.get("detailed_population_frequencies", []),
            focus_alleles=focus_alleles,
            effect_allele=effect_allele,
        )

        for row in top_level_rows:
            location_groups.add(row["location_group"])

        population_variant_records.append(
            {
                "variant": record["variant"],
                "display_name": record.get("display_name", record["variant"]),
                "common_name": record.get("common_name"),
                "effect_allele": effect_allele,
                "effect_summary": record.get("effect_summary", ""),
                "functional_effects": knowledge_record.get("functional_effects", []),
                "associated_conditions": knowledge_record.get("associated_conditions", []),
                "research_context": knowledge_record.get("research_context", []),
                "observed_in_run": record["variant"] in observed_variant_map,
                "observed_variants": observed_variant_map.get(record["variant"], []),
                "top_level_location_frequencies": top_level_rows,
                "detailed_population_groups": _group_population_rows_by_location(detailed_rows),
                "population_extremes": _summarize_population_extremes(top_level_rows, effect_allele),
                "source_url": record.get("source_url"),
            }
        )

    matched_population_variants = [
        record["display_name"] for record in population_variant_records if record["observed_in_run"]
    ]
    if matched_population_variants:
        overlap_note = (
            "The current run directly overlaps curated population-backed variants: "
            + ", ".join(matched_population_variants)
            + "."
        )
    else:
        overlap_note = (
            "The current run did not directly hit one of the curated SNPs with built-in population frequencies, "
            f"so this section provides reference context for commonly studied {gene_name} variants."
        )

    if population_variant_records:
        location_summary = ", ".join(sorted(location_groups)) if location_groups else "available reference panels"
        summary = (
            f"Population reference data are available for {len(population_variant_records)} curated {gene_name} SNPs "
            f"across {location_summary} panels. "
            f"{overlap_note}"
        )
    elif population_database.get("gene_population_patterns"):
        summary = (
            f"No embedded allele-frequency panel is bundled for {gene_name} yet, but literature-backed gene-level "
            f"population notes are available. {overlap_note}"
        )
    else:
        summary = (
            f"No curated population reference panel is bundled for {gene_name} yet. "
            "The app is therefore showing raw variant and methylation results without population-frequency overlays."
        )

    return {
        "database_name": population_database.get("database_name", f"Local {gene_name} population database"),
        "summary": summary,
        "location_groups": sorted(location_groups),
        "variant_population_records": population_variant_records,
        "gene_population_patterns": population_database.get("gene_population_patterns", []),
        "gene_population_patterns_intro": population_database.get(
            "gene_population_patterns_intro",
            f"Broader population patterns curated from the {gene_name} literature.",
        ),
        "sources": population_database.get("sources", []),
    }


def build_methylation_insights(
    methylation: pd.DataFrame, knowledge_base: dict[str, Any]
) -> dict[str, Any]:
    """Build a gene-level methylation interpretation from the current probe table."""
    gene_context = knowledge_base.get("gene_context", {})
    gene_name = str(gene_context.get("gene_name", DEFAULT_GENE_NAME))
    relevant_probe_ids = list(gene_context.get("relevant_methylation_probe_ids", []))
    relevant_probe_lookup = set(relevant_probe_ids)

    observed_relevant = methylation[methylation["probe_id"].isin(relevant_probe_lookup)].copy()
    if observed_relevant.empty:
        observed_relevant = methylation.copy()

    mean_beta = float(observed_relevant["beta"].mean()) if not observed_relevant.empty else None
    beta_band = _categorize_beta(mean_beta)

    if "UCSC_RefGene_Group" in observed_relevant.columns:
        group_breakdown = (
            observed_relevant["UCSC_RefGene_Group"]
            .fillna("Unannotated")
            .value_counts()
            .to_dict()
        )
    else:
        group_breakdown = {}

    group_summary = ", ".join(f"{count} {group}" for group, count in group_breakdown.items())
    if not group_summary:
        group_summary = "annotation breakdown unavailable"

    summary_prefix = (
        f"The current run captured {len(observed_relevant)} curated {gene_name} probe(s) "
        f"with a mean beta of {mean_beta:.2f}, which sits in the {beta_band} range. "
        if mean_beta is not None
        else f"The current run did not yield a mean beta estimate for the curated {gene_name} probes. "
    )

    summary = (
        f"{summary_prefix}The bundled {gene_name} array subset is dominated by {group_summary}. "
        f"{gene_context.get('methylation_interpretation', '')}"
    ).strip()

    preview_columns = [
        "probe_id",
        "beta",
        "UCSC_RefGene_Group",
        "Relation_to_UCSC_CpG_Island",
        "UCSC_CpG_Islands_Name",
    ]
    available_preview_columns = [
        column for column in preview_columns if column in observed_relevant.columns
    ]

    return {
        "gene_name": gene_name,
        "clinical_context": gene_context.get("clinical_context", ""),
        "summary": summary,
        "mean_beta": round(mean_beta, 3) if mean_beta is not None else None,
        "beta_band": beta_band,
        "observed_probe_count": int(len(observed_relevant)),
        "curated_probe_count": len(relevant_probe_ids),
        "probe_ids": observed_relevant["probe_id"].tolist(),
        "group_breakdown": group_breakdown,
        "methylation_effects": gene_context.get("methylation_effects", []),
        "methylation_condition_research": gene_context.get("methylation_condition_research", []),
        "evidence": gene_context.get("evidence", []),
        "probe_preview": observed_relevant[available_preview_columns].copy(),
    }


def build_generic_variant_interpretations(
    variants: pd.DataFrame,
    *,
    region: str,
    gene_name: str,
) -> dict[str, Any]:
    """Build a generic region-based interpretation payload for genes without a curated database."""
    search_region = _parse_region_string(region)
    chrom = str(search_region["chrom"]).removeprefix("chr")
    search_region_record = _build_interval_record(
        f"{gene_name} selected interval",
        chrom,
        int(search_region["start"]),
        int(search_region["end"]),
        "Interval selected during preprocessing or manual entry for the current gene.",
    )

    observed_summaries = [
        _build_observed_variant_summary(row, None, gene_name=gene_name)
        for _, row in variants.sort_values("pos").head(12).iterrows()
    ]
    highlight_items = [
        {
            "title": item["display"],
            "observed_variant": item["position"],
            "category": "Observed PASS variant",
            "description": item["summary"],
            "conditions": [],
            "literature_findings": [],
        }
        for item in observed_summaries
    ]

    no_model_region = {
        "label": "Promoter model unavailable",
        "window": "Manual promoter interval required",
        "length_bp": 0,
        "definition": (
            f"No curated promoter interval is bundled for {gene_name} yet. "
            "Use a gene-specific promoter window if you need promoter-only review."
        ),
        "included": False,
        "analysis_note": (
            f"Promoter-specific interpretation is not currently modeled for {gene_name}. "
            "The app is using the selected interval as a generic gene-level review window."
        ),
        "found_variant_count": 0,
        "found_variants": [],
        "known_variants": [],
    }

    return {
        "summary": (
            f"{gene_name} is being analyzed with a generic region-based workflow. "
            f"No curated interpretation database is bundled for this gene yet, so the app is prioritizing "
            f"observed PASS variants and raw previews inside {search_region_record['display']}."
        ),
        "matched_records": [],
        "unclassified_variant_count": int(len(variants)),
        "gene_summary": (
            f"The current run treats {search_region_record['display']} as the active {gene_name} review interval."
        ),
        "database_name": f"No curated {gene_name} interpretation database loaded",
        "gene_name": gene_name,
        "clinical_context": (
            f"Variant interpretation for {gene_name} is currently generic. "
            "Observed calls are shown directly, but no bundled disease- or gene-specific assertions are being applied."
        ),
        "variant_effect_overview": [],
        "condition_research_overview": [],
        "sample_highlights": {
            "summary": (
                f"This sample yielded {len(observed_summaries)} visible PASS variant(s) in the current {gene_name} interval preview."
            ),
            "highlight_items": highlight_items,
        },
        "region_recommendations": [
            {
                "title": "Promoter only",
                "region": "Set manually for this gene",
                "purpose": (
                    f"No curated promoter span is bundled for {gene_name} yet. "
                    "Enter a promoter-focused coordinate window manually if you need upstream review."
                ),
            },
            {
                "title": "Gene body only",
                "region": search_region_record["display"],
                "purpose": (
                    f"Use the selected {gene_name} interval when you want a direct gene-body review with the current generic workflow."
                ),
            },
            {
                "title": "Promoter plus gene body",
                "region": search_region_record["display"],
                "purpose": (
                    f"Start from the selected {gene_name} interval and extend it upstream manually if you want promoter-plus-gene coverage."
                ),
            },
        ],
        "gene_region": search_region_record,
        "promoter_region": {
            "label": "Promoter model unavailable",
            "chrom": chrom,
            "start": int(search_region["start"]),
            "end": int(search_region["start"]),
            "length_bp": 0,
            "display": "Manual promoter interval required",
            "definition": no_model_region["definition"],
        },
        "promoter_hotspot_region": {
            "label": "Promoter hotspot unavailable",
            "chrom": chrom,
            "start": int(search_region["start"]),
            "end": int(search_region["start"]),
            "length_bp": 0,
            "display": "No curated hotspot available",
            "definition": f"No curated promoter hotspot interval is bundled for {gene_name}.",
        },
        "search_region": search_region_record,
        "promoter_analysis": no_model_region,
        "gene_analysis": {
            "label": f"{gene_name} selected interval",
            "window": search_region_record["display"],
            "length_bp": search_region_record["length_bp"],
            "definition": search_region_record["definition"],
            "included": True,
            "analysis_note": (
                f"The current run found {len(observed_summaries)} PASS variant(s) in the first preview slice for the selected {gene_name} interval."
            ),
            "found_variant_count": len(observed_summaries),
            "found_variants": observed_summaries,
            "known_variants": [],
        },
        "recommended_promoter_plus_gene_region": search_region_record["display"],
    }


def build_generic_methylation_insights(methylation: pd.DataFrame, *, gene_name: str) -> dict[str, Any]:
    """Build a generic methylation payload for genes without a curated interpretation database."""
    mean_beta = float(methylation["beta"].mean()) if "beta" in methylation.columns and not methylation.empty else None
    beta_band = _categorize_beta(mean_beta)

    if "UCSC_RefGene_Group" in methylation.columns:
        group_breakdown = (
            methylation["UCSC_RefGene_Group"].fillna("Unannotated").value_counts().to_dict()
        )
    else:
        group_breakdown = {}

    preview_columns = [
        "probe_id",
        "beta",
        "chrom",
        "pos",
        "UCSC_RefGene_Group",
        "Relation_to_UCSC_CpG_Island",
        "UCSC_CpG_Islands_Name",
    ]
    available_preview_columns = [column for column in preview_columns if column in methylation.columns]

    summary_prefix = (
        f"The current run captured {len(methylation)} probe(s) from the selected {gene_name} manifest subset "
        f"with a mean beta of {mean_beta:.2f}. "
        if mean_beta is not None
        else f"The current run captured {len(methylation)} probe(s) from the selected {gene_name} manifest subset. "
    )

    return {
        "gene_name": gene_name,
        "clinical_context": (
            f"No curated methylation interpretation database is bundled for {gene_name} yet."
        ),
        "summary": (
            summary_prefix
            + "These values are shown as gene-focused methylation context rather than as a curated clinical interpretation."
        ),
        "mean_beta": round(mean_beta, 3) if mean_beta is not None else None,
        "beta_band": beta_band,
        "observed_probe_count": int(len(methylation)),
        "curated_probe_count": int(len(methylation)),
        "probe_ids": methylation["probe_id"].tolist() if "probe_id" in methylation.columns else [],
        "group_breakdown": group_breakdown,
        "methylation_effects": [
            f"The current methylation summary for {gene_name} is based on the filtered manifest subset selected during preprocessing.",
            "Beta values are being shown as region-level context until a curated gene-specific methylation knowledge base is added.",
        ],
        "methylation_condition_research": [],
        "evidence": [],
        "probe_preview": methylation[available_preview_columns].copy() if available_preview_columns else methylation.head(12).copy(),
    }


def build_empty_population_insights(*, gene_name: str) -> dict[str, Any]:
    """Return a placeholder population-insight payload for genes without curated reference data."""
    return {
        "summary": (
            f"No curated population reference database is bundled for {gene_name} yet. "
            "The app is therefore showing raw variant and methylation results without population-frequency overlays."
        ),
        "location_groups": [],
        "sources": [],
        "variant_population_records": [],
        "gene_population_patterns": [],
        "gene_population_patterns_intro": "",
        "database_name": f"No curated {gene_name} population database loaded",
        "database_version": "generic",
    }


def load_variants(vcf_path: str, region: str) -> pd.DataFrame:
    """Load PASS variants for the requested region from a tabix-indexed VCF.

    Parameters
    ----------
    vcf_path : str
        Filesystem path to a bgzip-compressed and tabix-indexed VCF.
    region : str
        Region string understood by ``scikit-allel``.

    Returns
    -------
    pd.DataFrame
        Variant table with one row per PASS site and the columns ``chrom``,
        ``pos``, ``ref``, ``alt``, ``qual``, and ``filter_pass``.

    Raises
    ------
    AnalysisError
        Raised when the VCF cannot be found, parsed, or filtered to any PASS
        variants in the requested region.
    """
    if not os.path.exists(vcf_path):
        raise AnalysisError(f"VCF not found: {vcf_path}")

    try:
        callset = allel.read_vcf(
            vcf_path,
            region=region,
            fields=[
                "variants/CHROM",
                "variants/ID",
                "variants/POS",
                "variants/REF",
                "variants/ALT",
                "variants/QUAL",
                "variants/FILTER_PASS",
            ],
        )
    except Exception as exc:
        raise AnalysisError(f"Failed to read VCF '{vcf_path}' for region '{region}': {exc}") from exc

    if not callset or "variants/CHROM" not in callset:
        raise AnalysisError(f"No variants could be read from {vcf_path} in region {region}")

    df = pd.DataFrame(
        {
            "chrom": [_decode_scalar(value) for value in callset["variants/CHROM"]],
            "id": [
                None
                if _decode_scalar(value) in {None, "."}
                else _decode_scalar(value)
                for value in callset["variants/ID"]
            ],
            "pos": callset["variants/POS"],
            "ref": [_decode_scalar(value) for value in callset["variants/REF"]],
            "alt": [
                _decode_scalar(alts[0]) if len(alts) > 0 else None
                for alts in callset["variants/ALT"]
            ],
            "qual": callset["variants/QUAL"],
            "filter_pass": callset["variants/FILTER_PASS"],
        }
    )

    # Keep only PASS variants so downstream summaries and the visible output
    # consistently reflect the analysis-ready subset.
    df = df[df["filter_pass"].fillna(False)].reset_index(drop=True)
    if df.empty:
        raise AnalysisError(f"No PASS variants found in {region} for {vcf_path}")

    logger.info("Loaded %d PASS variants from %s in region %s", len(df), vcf_path, region)
    return df


def _gene_manifest_filename(gene_name: str, genome_build: str = "hg19") -> str:
    """Build the conventional gene-specific manifest subset filename."""
    return f"{sanitize_gene_name_for_filename(gene_name)}_epigenetics_{genome_build}.csv"


def _prepare_gene_manifest_subset(
    *,
    manifest_filepath: str,
    gene_name: str,
    region: str,
    genome_build: str = "hg19",
) -> Path:
    """Save the selected gene-region manifest subset into `src/gene_data`.

    The full manifest path is still passed through to methylprep, but the
    analysis also keeps a smaller gene-specific CSV next to the bundled
    knowledge-base files so the post-pipeline annotation join can reuse it.
    """
    try:
        selection = save_filtered_manifest(
            gene_name=gene_name,
            manifest_path=manifest_filepath,
            region=region,
            genome_build=genome_build,
            output_dir=GENE_DATA_DIR,
        )
    except Exception as exc:
        raise AnalysisError(
            f"Failed to prepare the {gene_name} methylation subset from '{manifest_filepath}': {exc}"
        ) from exc

    output_path = Path(selection["output_path"])
    logger.info("Saved %s gene manifest subset to %s", gene_name, output_path)
    return output_path


def _run_methylprep_pipeline(data_dir: str, *, manifest_filepath: str | None) -> pd.DataFrame:
    """Run methylprep with the requested manifest override, if any."""
    return run_pipeline(
        data_dir,
        export=True,
        betas=True,
        manifest_filepath=manifest_filepath,
    )


def load_methylation(
    idat_base: str,
    manifest_filepath: str | None = None,
    *,
    gene_name: str = DEFAULT_GENE_NAME,
    region: str = DEFAULT_REGION,
) -> pd.DataFrame:
    """Load methylation beta values and annotate them with a selected manifest subset.

    Parameters
    ----------
    idat_base : str
        Path prefix shared by the two Illumina IDAT files. For a sample stored as
        ``data/R01C01_Grn.idat`` and ``data/R01C01_Red.idat``, pass
        ``data/R01C01``.
    manifest_filepath : str | None, optional
        Optional path to the full Illumina manifest. When provided, the
        function saves a gene-specific subset CSV to ``src/gene_data`` before
        running methylprep, then reuses that subset for the final probe
        annotation join.
    gene_name : str, optional
        Gene symbol used to derive the gene-specific manifest subset filename.
    region : str, optional
        Genomic interval used to select probe rows for the gene-specific
        manifest subset.

    Returns
    -------
    pd.DataFrame
        Probe-level methylation table limited to probes present in the selected
        region-specific manifest subset.

    Raises
    ------
    AnalysisError
        Raised when the IDAT pair is incomplete, methylprep fails, or the local
        manifest subset cannot be loaded.
    """
    logger.info("Starting methylation loading for sample")
    normalized_gene_name = gene_name.strip().upper() or DEFAULT_GENE_NAME

    data_dir = os.path.dirname(idat_base) or "."
    sample_name = os.path.basename(idat_base)

    logger.debug("Data directory: %s", data_dir)
    logger.debug("Sample name: %s", sample_name)

    # Validate both channels first so failures are short and obvious.
    for suffix in ("_Grn.idat", "_Red.idat"):
        path = os.path.join(data_dir, sample_name + suffix)
        logger.debug("Checking IDAT: %s", path)
        if not os.path.isfile(path):
            raise AnalysisError(f"Missing IDAT file: {path}")

    region_manifest_file = (
        Path(__file__).resolve().parent
        / "gene_data"
        / _gene_manifest_filename(normalized_gene_name)
    )
    pipeline_manifest_path = manifest_filepath
    if manifest_filepath:
        region_manifest_file = _prepare_gene_manifest_subset(
            manifest_filepath=manifest_filepath,
            gene_name=normalized_gene_name,
            region=region,
        )

    try:
        logger.info("Running methylprep pipeline with betas=True")
        beta_values = _run_methylprep_pipeline(
            data_dir,
            manifest_filepath=pipeline_manifest_path,
        )
        logger.debug("Beta-values DataFrame loaded, shape = %s", beta_values.shape)
    except Exception as exc:
        if pipeline_manifest_path:
            logger.warning(
                "methylprep rejected custom manifest '%s'; retrying with methylprep's default manifest.",
                pipeline_manifest_path,
                exc_info=True,
            )
            try:
                beta_values = _run_methylprep_pipeline(data_dir, manifest_filepath=None)
                logger.debug("Beta-values DataFrame loaded, shape = %s", beta_values.shape)
            except Exception as retry_exc:
                logger.exception("run_pipeline failed even after retrying without a custom manifest")
                raise AnalysisError(
                    f"methylprep failed for sample '{sample_name}' with custom manifest "
                    f"'{pipeline_manifest_path}', and the retry without that manifest also failed: {retry_exc}"
                ) from retry_exc
        else:
            logger.exception("run_pipeline failed")
            raise AnalysisError(f"methylprep failed for sample '{sample_name}': {exc}") from exc

    if sample_name not in beta_values.columns:
        raise AnalysisError(f"No beta column named '{sample_name}' in methylprep output")

    # Reset the index so probe IDs become a merge key instead of staying hidden
    # in the wide matrix index produced by methylprep.
    beta_df = beta_values[sample_name].rename("beta").reset_index()
    beta_df = beta_df.rename(columns={beta_df.columns[0]: "probe_id"})
    logger.debug("First 5 rows of beta-values DataFrame after renaming:\n%s", beta_df.head(5))
    try:
        manifest_region = pd.read_csv(region_manifest_file)
    except Exception as exc:
        raise AnalysisError(f"Failed to read region manifest file '{region_manifest_file}': {exc}") from exc

    rename_cols = {
        "IlmnID": "probe_id",
        "CHR": "chrom",
        "MAPINFO": "pos",
        "UCSC_RefGene_Name": "gene",
    }
    manifest_region = manifest_region.rename(columns=rename_cols)

    required = {"probe_id", "chrom", "pos", "gene"}
    missing = required - set(manifest_region.columns)
    if missing:
        raise AnalysisError(f"Region manifest is missing required columns: {sorted(missing)}")

    # An inner join keeps only probes observed in both the sample output and the
    # curated DRD4 manifest subset.
    merged = pd.merge(beta_df, manifest_region, on="probe_id", how="inner")
    logger.info("After merging, %d probes remain", len(merged))

    keep_columns = [
        "probe_id",
        "beta",
        "chrom",
        "pos",
        "GencodeBasicV12_NAME",
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Accession",
        "UCSC_RefGene_Group",
        "UCSC_CpG_Islands_Name",
        "Relation_to_UCSC_CpG_Island",
        "Phantom4_Enhancers",
        "Phantom5_Enhancers,DMR,450k_Enhancer,HMM_Island",
        "DNase_Hypersensitivity_NAME",
        "DNase_Hypersensitivity_Evidence_Count",
    ]

    available = set(merged.columns)
    missing_columns = [column for column in keep_columns if column not in available]
    if missing_columns:
        logger.warning("Some expected columns are missing: %s", missing_columns)

    final_df = merged[[column for column in keep_columns if column in available]]
    logger.debug("First 5 rows of filtered epigenetic data table:\n%s", final_df.head(5))
    return final_df


def fetch_population_stats(popstats_source: str, variants: pd.DataFrame) -> Any:
    """Load optional population statistics from a CSV or JSON sidecar file.

    Parameters
    ----------
    popstats_source : str
        Path to a CSV or JSON file containing population-level annotations.
    variants : pd.DataFrame
        Variant table returned by :func:`load_variants`. The current
        implementation does not merge by variant yet, but the argument is kept so
        the function signature already matches the future enrichment step.

    Returns
    -------
    Any
        Parsed CSV as a DataFrame or JSON as a Python object.

    Raises
    ------
    AnalysisError
        Raised when the file does not exist or uses an unsupported extension.
    """
    _ = variants
    popstats_path = Path(popstats_source)
    if not popstats_path.exists():
        raise AnalysisError(f"Population statistics file not found: {popstats_source}")

    suffix = popstats_path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(popstats_path)
    if suffix == ".json":
        return json.loads(popstats_path.read_text(encoding="utf-8"))

    raise AnalysisError(
        f"Unsupported population statistics format '{popstats_path.suffix}'. Use CSV or JSON."
    )


def _render_section_table(df: pd.DataFrame, title: str, rows: int = 20) -> str:
    """Render a DataFrame preview section for the HTML report."""
    preview = df.head(rows).to_html(index=False, classes="data-table", border=0)
    return f"<section><h2>{html.escape(title)}</h2>{preview}</section>"


def generate_report(
    variants: pd.DataFrame,
    methylation: pd.DataFrame,
    popstats: Any,
    output_path: str,
    *,
    gene_name: str = DEFAULT_GENE_NAME,
    region: str,
    methylation_output_path: Path | None = None,
) -> Path:
    """Generate a report artifact from the assembled analysis tables.

    Parameters
    ----------
    variants : pd.DataFrame
        PASS variant table for the requested gene interval.
    methylation : pd.DataFrame
        Annotated probe-level beta-value table returned by
        :func:`load_methylation`.
    popstats : Any
        Optional population statistics payload loaded from CSV or JSON.
    output_path : str
        Destination path for the final report artifact.
    region : str
        Genomic interval used during the run. It is surfaced in the report
        summary so the output remains self-describing.
    gene_name : str, optional
        Gene name displayed in the report heading and summary copy.
    methylation_output_path : Path | None, optional
        Path to the exported methylation CSV, shown in the report when provided.

    Returns
    -------
    Path
        Final report path written to disk.

    Raises
    ------
    AnalysisError
        Raised when the requested report extension is unsupported.
    """
    report_path = Path(output_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = report_path.suffix.lower() or ".html"

    popstats_section = ""
    if isinstance(popstats, pd.DataFrame):
        popstats_section = _render_section_table(popstats, "Population Statistics Preview")
    elif popstats is not None:
        payload = html.escape(json.dumps(popstats, indent=2))
        popstats_section = f"<section><h2>Population Statistics Preview</h2><pre>{payload}</pre></section>"

    if suffix == ".html":
        methylation_path_markup = ""
        if methylation_output_path is not None:
            methylation_path_markup = (
                "<p><strong>Methylation CSV:</strong> "
                f"{html.escape(str(methylation_output_path))}</p>"
            )

        report_html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(gene_name)} Analysis Report</title>
  <style>
    :root {{
      --bg: #f6efe3;
      --panel: rgba(255, 252, 245, 0.92);
      --ink: #1f2a2e;
      --muted: #51666a;
      --accent: #0f766e;
      --accent-2: #c26a3d;
      --line: rgba(31, 42, 46, 0.14);
      --shadow: 0 24px 70px rgba(31, 42, 46, 0.12);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(194, 106, 61, 0.20), transparent 30rem),
        radial-gradient(circle at top right, rgba(15, 118, 110, 0.20), transparent 28rem),
        linear-gradient(180deg, #fbf6ee 0%, var(--bg) 100%);
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 48px 20px 64px;
    }}
    .hero {{
      padding: 28px 30px;
      border-radius: 28px;
      background: var(--panel);
      border: 1px solid var(--line);
      box-shadow: var(--shadow);
      backdrop-filter: blur(14px);
    }}
    .hero h1 {{
      margin: 0 0 12px;
      font-size: clamp(2rem, 5vw, 3.6rem);
      line-height: 1;
      letter-spacing: -0.04em;
    }}
    .hero p {{
      margin: 8px 0;
      color: var(--muted);
      max-width: 52rem;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 14px;
      margin-top: 24px;
    }}
    .metric {{
      padding: 18px 20px;
      border-radius: 20px;
      background: rgba(255, 255, 255, 0.7);
      border: 1px solid var(--line);
    }}
    .metric span {{
      display: block;
      color: var(--muted);
      font-size: 0.9rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 1.8rem;
    }}
    section {{
      margin-top: 24px;
      padding: 24px;
      border-radius: 24px;
      background: var(--panel);
      border: 1px solid var(--line);
      box-shadow: var(--shadow);
    }}
    h2 {{
      margin-top: 0;
      font-size: 1.3rem;
      letter-spacing: -0.02em;
    }}
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.92rem;
    }}
    .data-table th,
    .data-table td {{
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
    }}
    .data-table th {{
      background: rgba(15, 118, 110, 0.08);
    }}
    pre {{
      padding: 16px;
      overflow-x: auto;
      border-radius: 16px;
      background: #172023;
      color: #f4f0e8;
    }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>{html.escape(gene_name)} Analysis Report</h1>
      <p>This report summarizes the current {html.escape(gene_name)} variant and methylation analysis run.</p>
      <p><strong>Region:</strong> {html.escape(region)}</p>
      <p><strong>Report path:</strong> {html.escape(str(report_path))}</p>
      {methylation_path_markup}
      <div class="metrics">
        <article class="metric">
          <span>PASS variants</span>
          <strong>{len(variants)}</strong>
        </article>
        <article class="metric">
          <span>Methylation probes</span>
          <strong>{len(methylation)}</strong>
        </article>
        <article class="metric">
          <span>Population stats</span>
          <strong>{"Yes" if popstats is not None else "No"}</strong>
        </article>
      </div>
    </section>
    {_render_section_table(variants, "Variant Preview")}
    {_render_section_table(methylation, "Methylation Preview")}
    {popstats_section}
  </main>
</body>
</html>
"""
        report_path.write_text(report_html, encoding="utf-8")
        return report_path

    if suffix == ".json":
        payload = {
            "region": region,
            "variants": variants.to_dict(orient="records"),
            "methylation": methylation.to_dict(orient="records"),
            "population_statistics": _serialize_popstats(popstats),
            "methylation_output_path": str(methylation_output_path) if methylation_output_path else None,
        }
        report_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return report_path

    if suffix == ".csv":
        summary = pd.DataFrame(
            [
                {"metric": "region", "value": region},
                {"metric": "variant_count", "value": len(variants)},
                {"metric": "methylation_probe_count", "value": len(methylation)},
                {"metric": "has_population_stats", "value": popstats is not None},
                {
                    "metric": "methylation_output_path",
                    "value": str(methylation_output_path) if methylation_output_path else "",
                },
            ]
        )
        summary.to_csv(report_path, index=False)
        return report_path

    raise AnalysisError(
        f"Unsupported output format '{report_path.suffix}'. Use .html, .json, or .csv."
    )


def run_analysis(
    *,
    vcf_path: str,
    idat_base: str,
    output_path: str,
    gene_name: str = DEFAULT_GENE_NAME,
    region: str = DEFAULT_REGION,
    popstats_source: str | None = None,
    manifest_filepath: str | None = None,
) -> AnalysisResult:
    """Run the end-to-end gene analysis workflow.

    Parameters
    ----------
    vcf_path : str
        Input VCF path containing DRD4-region variants.
    idat_base : str
        Input IDAT prefix without the color suffix.
    output_path : str
        Output report path.
    gene_name : str, optional
        Gene symbol associated with the current run.
    region : str, optional
        Genomic region to inspect.
    popstats_source : str | None, optional
        Optional population statistics sidecar file path.
    manifest_filepath : str | None, optional
        Optional full manifest file path passed through to methylprep and used
        to refresh the gene-specific subset stored in ``src/gene_data``.

    Returns
    -------
    AnalysisResult
        Structured result containing the in-memory tables plus the generated
        output paths.
    """
    normalized_gene_name = gene_name.strip().upper() or DEFAULT_GENE_NAME

    variants = load_variants(vcf_path, region)
    methylation = load_methylation(
        idat_base,
        manifest_filepath=manifest_filepath,
        gene_name=normalized_gene_name,
        region=region,
    )

    knowledge_base = load_gene_interpretation_database(normalized_gene_name)
    if knowledge_base is not None:
        variant_interpretations = build_variant_interpretations(variants, knowledge_base, region=region)
        methylation_insights = build_methylation_insights(methylation, knowledge_base)
    else:
        knowledge_base = {
            "database_name": f"No curated {normalized_gene_name} interpretation database loaded",
            "version": "generic",
            "gene_context": {"gene_name": normalized_gene_name},
        }
        population_database = {
            "database_name": f"No curated {normalized_gene_name} population database loaded",
            "version": "generic",
        }
        variant_interpretations = build_generic_variant_interpretations(
            variants,
            region=region,
            gene_name=normalized_gene_name,
        )
        methylation_insights = build_generic_methylation_insights(
            methylation,
            gene_name=normalized_gene_name,
        )

    population_database = load_gene_population_database(normalized_gene_name)
    if population_database is not None and knowledge_base.get("version") != "generic":
        population_insights = build_population_insights(variants, knowledge_base, population_database)
    else:
        population_database = {
            "database_name": f"No curated {normalized_gene_name} population database loaded",
            "version": "generic",
        }
        population_insights = build_empty_population_insights(gene_name=normalized_gene_name)

    report_path = Path(output_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    methylation_output_path = _derive_methylation_output_path(report_path)
    methylation_output_path.parent.mkdir(parents=True, exist_ok=True)
    methylation.to_csv(methylation_output_path, index=False)
    logger.info("Saved methylation data to %s", methylation_output_path)

    popstats = None
    if popstats_source:
        popstats = fetch_population_stats(popstats_source, variants)

    final_report_path = generate_report(
        variants,
        methylation,
        popstats,
        str(report_path),
        gene_name=normalized_gene_name,
        region=region,
        methylation_output_path=methylation_output_path,
    )

    return AnalysisResult(
        variants=variants,
        methylation=methylation,
        popstats=popstats,
        report_path=final_report_path,
        methylation_output_path=methylation_output_path,
        region=region,
        vcf_path=Path(vcf_path),
        idat_base=Path(idat_base),
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        knowledge_base=knowledge_base,
        population_insights=population_insights,
        population_database=population_database,
    )


def main(argv: list[str] | None = None) -> int:
    """Run the DRD4 workflow from the command line and return an exit code."""
    try:
        args = parse_args(argv)
        result = run_analysis(
            vcf_path=args.vcf,
            idat_base=args.idat,
            output_path=args.out,
            gene_name=DEFAULT_GENE_NAME,
            region=args.region,
            popstats_source=args.popstats,
            manifest_filepath=args.manifest_file,
        )
    except AnalysisError as exc:
        logger.error("%s", exc)
        return 1

    print(f"Saved report to {result.report_path}")
    print(f"Saved methylation CSV to {result.methylation_output_path}")
    return 0


if __name__ == "__main__":
    print(__doc__)
    raise SystemExit(main())
