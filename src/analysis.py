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

DEFAULT_REGION = "11:637269-640706"
DEFAULT_REPORT_NAME = "drd4_report.html"
INTERPRETATION_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "drd4_interpretation_db.json"
POPULATION_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "drd4_population_db.json"

# Configure the root logger once so both CLI and web runs stream progress.
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


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
        "assay_note": assay_note,
        "evidence": record.get("evidence", []),
    }


def _build_observed_variant_summary(
    row: pd.Series,
    matched_record: dict[str, Any] | None,
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
            "clinical_significance": "No curated clinical significance available in the local DRD4 database for this observed site.",
            "summary": "Observed inside the selected DRD4 search window, but not one of the seeded common DRD4 promoter or gene variants.",
            "functional_effects": [],
            "associated_conditions": [],
            "research_context": [],
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
        "evidence": matched_record.get("evidence", []),
    }


def _summarize_observed_region_variants(
    variants: pd.DataFrame,
    variant_records: list[dict[str, Any]],
    *,
    chrom: str,
    start: int,
    end: int,
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
        summaries.append(_build_observed_variant_summary(row, _match_variant_record(row, variant_records)))
    return summaries


def _build_region_variant_analysis(
    *,
    region_record: dict[str, Any],
    search_region: dict[str, Any],
    variants: pd.DataFrame,
    curated_records: list[dict[str, Any]],
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
) -> dict[str, Any]:
    """Summarize the most visible sample-specific variant findings first."""
    if matched_records:
        summary = (
            f"This sample matched {len(matched_records)} curated DRD4 research variant(s). "
            "Those direct sample hits are shown first below."
        )
        items = [
            {
                "title": record["variant"],
                "observed_variant": record["observed_variant"],
                "category": record.get("interpretation_scope", "Research context"),
                "description": record.get("clinical_interpretation", ""),
                "conditions": record.get("associated_conditions", []),
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
                }
            )

    if found_items:
        summary = (
            "This sample did not hit one of the curated named DRD4 markers, but it did contain "
            f"{len(found_items)} observed PASS variant(s) inside the reviewed promoter or gene intervals."
        )
    else:
        summary = (
            "This sample did not yield a visible DRD4 promoter or gene-body PASS variant in the current preview slice."
        )

    return {
        "summary": summary,
        "highlight_items": found_items,
    }


def _build_region_recommendations(
    promoter_region_record: dict[str, Any],
    gene_region_record: dict[str, Any],
    combined_region: str,
) -> list[dict[str, str]]:
    """Return practical region-span recommendations for common DRD4 review goals."""
    return [
        {
            "title": "Promoter only",
            "region": promoter_region_record["display"],
            "purpose": "Use this when you want to focus on the upstream promoter-review window and the classic DRD4 promoter polymorphism hotspot without loading the transcribed gene body.",
        },
        {
            "title": "Gene body only",
            "region": gene_region_record["display"],
            "purpose": "Use this when you want the canonical DRD4 transcribed interval but do not need the upstream promoter review window.",
        },
        {
            "title": "Promoter plus gene body",
            "region": combined_region,
            "purpose": "Use this when you want the full audit: upstream promoter context plus the canonical DRD4 gene interval in one search.",
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
        gene_region_source.get("label", "DRD4 transcribed interval"),
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
        inclusion_hint=(
            "Use the promoter-plus-gene interval "
            f"{gene_context.get('recommended_promoter_plus_gene_region', _format_interval(chrom, promoter_region_record['start'], gene_region_record['end']))} "
            "if you want the upstream DRD4 promoter reviewed alongside the gene."
        ),
    )
    gene_analysis = _build_region_variant_analysis(
        region_record=gene_region_record,
        search_region=search_region,
        variants=variants,
        curated_records=gene_records,
        inclusion_hint="Choose a region that overlaps the canonical DRD4 gene interval if you want gene-body variants interpreted.",
    )

    promoter_phrase = "overlaps" if promoter_analysis["included"] else "does not overlap"
    gene_phrase = "overlaps" if gene_analysis["included"] else "does not overlap"
    combined_region = gene_context.get(
        "recommended_promoter_plus_gene_region",
        _format_interval(chrom, promoter_region_record["start"], gene_region_record["end"]),
    )
    summary = (
        f"DRD4 is located at {gene_region_record['display']} on {gene_context.get('cytoband', '11p15.5')} "
        f"and spans {gene_region_record['length_bp']:,} bp on the {gene_context.get('assembly', 'GRCh37 / hg19')} assembly. "
        f"The current search interval {search_region_record['display']} {promoter_phrase} the operational promoter review window "
        f"{promoter_region_record['display']} and {gene_phrase} the DRD4 transcribed interval {gene_region_record['display']}. "
        f"In the current preview, {promoter_analysis['found_variant_count']} promoter-window variant(s) and "
        f"{gene_analysis['found_variant_count']} gene-interval variant(s) were surfaced."
    )

    return {
        "summary": summary,
        "matched_records": matched_records,
        "unclassified_variant_count": max(len(variants) - len(matched_records), 0),
        "gene_summary": gene_context.get("gene_summary", ""),
        "database_name": knowledge_base.get("database_name", "Local DRD4 interpretation database"),
        "gene_name": gene_context.get("gene_name", "DRD4"),
        "clinical_context": gene_context.get("clinical_context", ""),
        "variant_effect_overview": gene_context.get("variant_effect_overview", []),
        "condition_research_overview": gene_context.get("condition_research_overview", []),
        "sample_highlights": _build_sample_variant_highlights(
            matched_records=matched_records,
            promoter_analysis=promoter_analysis,
            gene_analysis=gene_analysis,
        ),
        "region_recommendations": _build_region_recommendations(
            promoter_region_record,
            gene_region_record,
            combined_region,
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
    """Build geography-aware population summaries for common curated DRD4 variants."""
    variant_records = knowledge_base.get("variant_records", [])
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
            "so this section provides reference context for commonly studied DRD4 variants."
        )

    summary = (
        f"Population reference data are available for {len(population_variant_records)} curated DRD4 SNPs "
        f"across {', '.join(sorted(location_groups))} panels, plus literature-based geographic patterns for the "
        "exon III VNTR. "
        f"{overlap_note}"
    )

    return {
        "database_name": population_database.get("database_name", "Local DRD4 population database"),
        "summary": summary,
        "location_groups": sorted(location_groups),
        "variant_population_records": population_variant_records,
        "gene_population_patterns": population_database.get("gene_population_patterns", []),
        "sources": population_database.get("sources", []),
    }


def build_methylation_insights(
    methylation: pd.DataFrame, knowledge_base: dict[str, Any]
) -> dict[str, Any]:
    """Build a gene-level methylation interpretation from the current probe table."""
    gene_context = knowledge_base.get("gene_context", {})
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
        f"The current run captured {len(observed_relevant)} curated DRD4 probe(s) "
        f"with a mean beta of {mean_beta:.2f}, which sits in the {beta_band} range. "
        if mean_beta is not None
        else "The current run did not yield a mean beta estimate for the curated DRD4 probes. "
    )

    summary = (
        f"{summary_prefix}The bundled DRD4 array subset is dominated by {group_summary}. "
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
        "gene_name": gene_context.get("gene_name", "DRD4"),
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


def load_methylation(idat_base: str, manifest_filepath: str | None = None) -> pd.DataFrame:
    """Load methylation beta values and annotate them with DRD4 region metadata.

    Parameters
    ----------
    idat_base : str
        Path prefix shared by the two Illumina IDAT files. For a sample stored as
        ``data/R01C01_Grn.idat`` and ``data/R01C01_Red.idat``, pass
        ``data/R01C01``.
    manifest_filepath : str | None, optional
        Optional path to a custom manifest file consumed directly by
        ``methylprep.run_pipeline``.

    Returns
    -------
    pd.DataFrame
        Probe-level methylation table limited to probes present in the curated
        DRD4 region manifest.

    Raises
    ------
    AnalysisError
        Raised when the IDAT pair is incomplete, methylprep fails, or the local
        DRD4 manifest subset cannot be loaded.
    """
    logger.info("Starting methylation loading for sample")

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

    try:
        logger.info("Running methylprep pipeline with betas=True")
        beta_values = run_pipeline(
            data_dir,
            export=True,
            betas=True,
            manifest_filepath=manifest_filepath,
        )
        logger.debug("Beta-values DataFrame loaded, shape = %s", beta_values.shape)
    except Exception as exc:
        logger.exception("run_pipeline failed")
        raise AnalysisError(f"methylprep failed for sample '{sample_name}': {exc}") from exc

    if sample_name not in beta_values.columns:
        raise AnalysisError(f"No beta column named '{sample_name}' in methylprep output")

    # Reset the index so probe IDs become a merge key instead of staying hidden
    # in the wide matrix index produced by methylprep.
    beta_df = beta_values[sample_name].rename("beta").reset_index()
    beta_df = beta_df.rename(columns={beta_df.columns[0]: "probe_id"})
    logger.debug("First 5 rows of beta-values DataFrame after renaming:\n%s", beta_df.head(5))

    region_manifest_file = Path(__file__).resolve().parent / "gene_data" / "drd4_epigenetics_hg19.csv"
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
    region: str,
    methylation_output_path: Path | None = None,
) -> Path:
    """Generate a report artifact from the assembled analysis tables.

    Parameters
    ----------
    variants : pd.DataFrame
        PASS variant table for the requested DRD4 interval.
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
  <title>DRD4 Analysis Report</title>
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
      <h1>DRD4 Analysis Report</h1>
      <p>This report summarizes the current DRD4 variant and methylation analysis run.</p>
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
    region: str = DEFAULT_REGION,
    popstats_source: str | None = None,
    manifest_filepath: str | None = None,
) -> AnalysisResult:
    """Run the end-to-end DRD4 analysis workflow.

    Parameters
    ----------
    vcf_path : str
        Input VCF path containing DRD4-region variants.
    idat_base : str
        Input IDAT prefix without the color suffix.
    output_path : str
        Output report path.
    region : str, optional
        Genomic region to inspect.
    popstats_source : str | None, optional
        Optional population statistics sidecar file path.
    manifest_filepath : str | None, optional
        Optional manifest file path passed through to methylprep.

    Returns
    -------
    AnalysisResult
        Structured result containing the in-memory tables plus the generated
        output paths.
    """
    variants = load_variants(vcf_path, region)
    methylation = load_methylation(idat_base, manifest_filepath=manifest_filepath)
    knowledge_base = load_interpretation_database()
    population_database = load_population_database()
    variant_interpretations = build_variant_interpretations(variants, knowledge_base, region=region)
    methylation_insights = build_methylation_insights(methylation, knowledge_base)
    population_insights = build_population_insights(variants, knowledge_base, population_database)

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
