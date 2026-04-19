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
PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENE_DATA_DIR = Path(__file__).resolve().parent / "gene_data"
INTERPRETATION_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "drd4_interpretation_db.json"
POPULATION_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "drd4_population_db.json"
SYNTHESIS_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "drd4_synthesis.json"
GENERAL_ANALYSIS_DATABASE_PATH = PROJECT_ROOT / "results" / "general_gene_analysis_database.csv"
GENERAL_ANALYSIS_DATABASE_COLUMNS = [
    "gene",
    "variant key",
    "observed gene variant",
    "gene variant label",
    "change",
    "chromosome",
    "position",
    "variant location",
    "gene location",
    "source",
    "(VCF) quality (qual)",
    "matched curated marker",
    "variant interpretation scope",
    "curated biological significance",
    "functional effects",
    "associated conditions",
    "methylation-linked probes",
    "mean beta whitelist",
    "mean beta related to gene",
    "mean beta on found probes in the area (numerical rows)",
]

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
    predictive_theses : dict[str, Any]
        Gene-level predictive synthesis payload assembled from the local
        synthesis database plus the current sample's variant and methylation
        results.
    general_database_path : Path
        Path to the central variant-level analysis database.
    general_database_status : str
        Human-readable status describing whether this run added, skipped, or
        overwrote the central database row.
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
    predictive_theses: dict[str, Any]
    general_database_path: Path
    general_database_status: str


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


def _load_json_database(
    database_path: str | Path,
    *,
    missing_label: str,
    invalid_label: str,
) -> dict[str, Any]:
    """Load a JSON-backed local database with a consistent error surface."""
    payload_path = Path(database_path)
    if not payload_path.exists():
        raise AnalysisError(f"{missing_label}: {payload_path}")

    try:
        return json.loads(payload_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AnalysisError(f"{invalid_label}: {payload_path}") from exc


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
    return _load_json_database(
        database_path,
        missing_label="Interpretation database not found",
        invalid_label="Interpretation database is not valid JSON",
    )


def load_population_database(database_path: str | Path = POPULATION_DB_PATH) -> dict[str, Any]:
    """Load the curated local DRD4 population database."""
    return _load_json_database(
        database_path,
        missing_label="Population database not found",
        invalid_label="Population database is not valid JSON",
    )


def load_synthesis_database(database_path: str | Path = SYNTHESIS_DB_PATH) -> dict[str, Any]:
    """Load the curated local predictive synthesis database."""
    return _load_json_database(
        database_path,
        missing_label="Synthesis database not found",
        invalid_label="Synthesis database is not valid JSON",
    )


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


def load_gene_synthesis_database(gene_name: str) -> dict[str, Any] | None:
    """Load a bundled gene-specific predictive synthesis database when one exists."""
    database_path = find_gene_database_path(gene_name, "synthesis.json")
    if database_path is None:
        return None
    return load_synthesis_database(database_path)


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


def _normalize_allele_change(value: Any) -> str:
    """Canonicalize allele-change strings such as ``A -> G`` and ``A>G``."""
    text = str(value or "").strip().upper().replace(" ", "")
    return text.replace("->", ">")


def _extract_alt_allele_from_change(value: Any) -> str:
    """Return the ALT allele from a display change such as ``A -> G`` when available."""
    normalized = _normalize_allele_change(value)
    if ">" not in normalized:
        return ""
    return normalized.rsplit(">", 1)[-1].strip()


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


def _format_variant_label(value: Any) -> str:
    """Render a compact variant label while preserving explicit missing IDs as ``None``."""
    if pd.isna(value):
        return "None"
    text = str(value).strip()
    if not text or text == ".":
        return "None"
    return text


def _format_variant_change(ref: Any, alt: Any) -> str:
    """Render a concise allele-change display for sample result tables."""
    ref_text = "" if pd.isna(ref) else str(ref).strip()
    alt_text = "" if pd.isna(alt) else str(alt).strip()
    if not ref_text and not alt_text:
        return "Unavailable"
    if not ref_text or not alt_text:
        return f"{ref_text or '?'} -> {alt_text or '?'}"
    return f"{ref_text} -> {alt_text}"


def _dedupe_text_items(items: list[str]) -> list[str]:
    """Return non-empty text items while preserving the first occurrence order."""
    seen: set[str] = set()
    deduped: list[str] = []
    for item in items:
        text = str(item).strip()
        if not text:
            continue
        normalized = text.casefold()
        if normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(text)
    return deduped


def _build_specific_variant_link_summary(
    record: dict[str, Any] | None,
    *,
    fallback: str,
) -> str:
    """Return a compact, variant-specific association summary for sample tables."""
    if not record:
        return fallback

    link_items: list[str] = []
    usual_variant_note = str(record.get("usual_variant_note", "")).strip()
    common_name = str(record.get("common_name", "")).strip()
    if usual_variant_note:
        link_items.append(usual_variant_note)
    if common_name:
        link_items.append(common_name)

    associated_conditions = [str(item).strip() for item in record.get("associated_conditions", []) if str(item).strip()]
    if associated_conditions:
        link_items.append("Studied in " + ", ".join(associated_conditions))

    literature_findings = record.get("literature_findings", [])
    if not associated_conditions and literature_findings:
        first_finding = literature_findings[0]
        phenotype = str(first_finding.get("phenotype", "")).strip()
        if phenotype:
            link_items.append(phenotype)

    functional_effects = [str(item).strip() for item in record.get("functional_effects", []) if str(item).strip()]
    if not associated_conditions and not literature_findings and functional_effects:
        link_items.append(functional_effects[0])

    deduped_items = _dedupe_text_items(link_items)
    if deduped_items:
        return "; ".join(deduped_items[:2])

    clinical_significance = str(record.get("clinical_significance", "")).strip()
    if clinical_significance:
        return clinical_significance

    clinical_interpretation = str(record.get("clinical_interpretation", "")).strip()
    if clinical_interpretation:
        return clinical_interpretation

    return fallback


def annotate_known_variant_ids(
    variants: pd.DataFrame,
    knowledge_base: dict[str, Any] | None,
) -> pd.DataFrame:
    """Fill display-friendly IDs when the curated bundle can name an unlabeled variant."""
    labeled = variants.copy()
    if labeled.empty:
        if "id_source" not in labeled.columns:
            labeled["id_source"] = pd.Series(dtype="object")
        return labeled

    variant_records = knowledge_base.get("variant_records", []) if knowledge_base else []
    resolved_ids: list[Any] = []
    id_sources: list[str] = []

    for _, row in labeled.iterrows():
        raw_id = row.get("id")
        raw_id_text = "" if pd.isna(raw_id) else str(raw_id).strip()
        matched_record = _match_variant_record(row, variant_records) if variant_records else None

        if raw_id_text and raw_id_text != ".":
            resolved_ids.append(raw_id_text)
            id_sources.append("Source VCF")
            continue

        if matched_record is not None:
            resolved_ids.append(
                str(matched_record.get("display_name", matched_record["variant"])).strip()
            )
            id_sources.append("Knowledge base match")
            continue

        resolved_ids.append(None)
        id_sources.append("Unlabeled in source VCF")

    labeled["id"] = resolved_ids
    labeled["id_source"] = id_sources
    return labeled


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
        "region_class": record.get("region_class", "research_marker"),
        "functional_effects": record.get("functional_effects", []),
        "associated_conditions": record.get("associated_conditions", []),
        "research_context": record.get("research_context", []),
        "literature_findings": _build_literature_findings(record),
        "assay_note": assay_note,
        "evidence": record.get("evidence", []),
    }


def _build_curated_named_marker_catalog(
    variants: pd.DataFrame,
    variant_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Return all curated named markers with run-specific observation flags."""
    observed_marker_map: dict[str, list[str]] = {}
    for _, row in variants.iterrows():
        matched_record = _match_variant_record(row, variant_records)
        if matched_record is None:
            continue
        observed_marker_map.setdefault(matched_record["variant"], []).append(_format_variant_display(row))

    marker_catalog: list[dict[str, Any]] = []
    for record in variant_records:
        summary = _build_known_variant_summary(record)
        observed_variants = observed_marker_map.get(record["variant"], [])
        marker_catalog.append(
            {
                **summary,
                "observed_in_run": bool(observed_variants),
                "observed_variants": observed_variants,
            }
        )
    return marker_catalog


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
            "variant_label": _format_variant_label(row.get("id")),
            "change": _format_variant_change(row.get("ref"), row.get("alt")),
            "linked_to": f"No curated local {gene_name} link is bundled for this PASS variant yet.",
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
        "variant_label": _format_variant_label(row.get("id")),
        "change": _format_variant_change(row.get("ref"), row.get("alt")),
        "linked_to": _build_specific_variant_link_summary(
            matched_record,
            fallback=f"No curated local {gene_name} link is bundled for this PASS variant yet.",
        ),
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
                "change": record.get("change", "Unavailable"),
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
            "result_table_rows": [
                {
                    "variant_label": record.get("variant_label", "None"),
                    "change": record.get("change", "Unavailable"),
                    "linked_to": record.get("linked_to", ""),
                }
                for record in matched_records
            ],
        }

    found_items: list[dict[str, Any]] = []
    result_table_rows: list[dict[str, str]] = []
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
            result_table_rows.append(
                {
                    "variant_label": record.get("variant_label", "None"),
                    "change": record.get("change", "Unavailable"),
                    "linked_to": record.get("linked_to", ""),
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
        "result_table_rows": result_table_rows,
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


def _coerce_numeric_beta_values(df: pd.DataFrame) -> pd.Series:
    """Return non-null numeric beta values from a methylation table."""
    if "beta" not in df.columns or df.empty:
        return pd.Series(dtype="float64")
    return pd.to_numeric(df["beta"], errors="coerce").dropna()


def _mean_beta_or_none(beta_values: pd.Series) -> float | None:
    """Return the mean beta value for a numeric series or ``None`` when empty."""
    if beta_values.empty:
        return None
    return float(beta_values.mean())


def _round_beta(mean_beta: float | None) -> float | None:
    """Round a mean beta value for UI display while preserving ``None``."""
    return round(mean_beta, 3) if mean_beta is not None else None


def _build_methylation_metric_rows(
    *,
    gene_name: str,
    whitelist_mean_beta: float | None,
    whitelist_count: int,
    gene_name_mean_beta: float | None,
    gene_name_count: int,
    raw_mean_beta: float | None,
    raw_count: int,
) -> list[dict[str, Any]]:
    """Build table-ready methylation summary rows for the UI and generated reports."""
    return [
        {
            "metric": "Whitelist mean beta",
            "mean_beta": _round_beta(whitelist_mean_beta),
            "mean_beta_display": str(_round_beta(whitelist_mean_beta)) if whitelist_mean_beta is not None else "Unavailable",
            "numeric_values": int(whitelist_count),
            "summary": (
                f"Observed whitelist probe mean across {whitelist_count} numeric value(s)."
                if whitelist_mean_beta is not None
                else "No numeric observed whitelist probe was available in this run."
            ),
        },
        {
            "metric": f"{gene_name}-named row mean beta",
            "mean_beta": _round_beta(gene_name_mean_beta),
            "mean_beta_display": str(_round_beta(gene_name_mean_beta)) if gene_name_mean_beta is not None else "Unavailable",
            "numeric_values": int(gene_name_count),
            "summary": (
                f"Rows where {gene_name} was explicitly annotated in the gene-name columns."
                if gene_name_mean_beta is not None
                else f"No numeric row explicitly annotated {gene_name} in the available gene-name columns."
            ),
        },
        {
            "metric": "All numeric-row mean beta",
            "mean_beta": _round_beta(raw_mean_beta),
            "mean_beta_display": str(_round_beta(raw_mean_beta)) if raw_mean_beta is not None else "Unavailable",
            "numeric_values": int(raw_count),
            "summary": (
                "All numeric beta values across the full raw methylation table."
                if raw_mean_beta is not None
                else "No numeric beta values were available in the raw methylation table."
            ),
        },
    ]


def _select_gene_named_methylation_rows(
    methylation: pd.DataFrame, gene_name: str
) -> tuple[pd.DataFrame, list[str]]:
    """Return rows whose gene-annotation columns explicitly mention ``gene_name``."""
    if methylation.empty:
        return methylation.iloc[0:0].copy(), []

    normalized_gene_name = gene_name.strip().upper()
    if not normalized_gene_name:
        return methylation.iloc[0:0].copy(), []

    annotation_columns = [
        column
        for column in ("GencodeBasicV12_NAME", "gene", "UCSC_RefGene_Name")
        if column in methylation.columns
    ]
    if not annotation_columns:
        return methylation.iloc[0:0].copy(), []

    match_pattern = rf"(?:^|;)\s*{re.escape(normalized_gene_name)}\s*(?:;|$)"
    combined_mask = pd.Series(False, index=methylation.index, dtype="bool")
    matched_columns: list[str] = []

    for column in annotation_columns:
        mask = (
            methylation[column]
            .fillna("")
            .astype(str)
            .str.upper()
            .str.contains(match_pattern, regex=True, na=False)
        )
        if bool(mask.any()):
            matched_columns.append(column)
            combined_mask = combined_mask | mask

    return methylation.loc[combined_mask].copy(), matched_columns


def _build_whitelist_explanation(gene_name: str, relevant_probe_ids: list[str]) -> str:
    """Explain how the curated methylation whitelist is assembled."""
    if not relevant_probe_ids:
        return (
            f"No curated {gene_name} methylation whitelist is bundled yet, so whitelist-only "
            "statistics are unavailable for this gene."
        )
    return (
        f"The curated {gene_name} methylation whitelist is bundled manually in the local "
        "interpretation database under `relevant_methylation_probe_ids`. It is a literature-guided "
        "hotspot subset used for interpretation, not an automatic list of every probe row in the "
        "current manifest slice."
    )


def _split_semicolon_tokens(value: Any) -> list[str]:
    """Split a semicolon-delimited annotation field into unique non-empty tokens."""
    if value is None or pd.isna(value):
        return []
    tokens = [token.strip() for token in str(value).split(";")]
    unique_tokens: list[str] = []
    seen_tokens: set[str] = set()
    for token in tokens:
        if not token or token in seen_tokens:
            continue
        seen_tokens.add(token)
        unique_tokens.append(token)
    return unique_tokens


def _load_gene_manifest_probe_lookup(
    gene_name: str, probe_ids: list[str]
) -> dict[str, dict[str, Any]]:
    """Load bundled manifest annotations for the requested whitelist probes."""
    if not probe_ids:
        return {}

    manifest_path = GENE_DATA_DIR / _gene_manifest_filename(gene_name)
    if not manifest_path.exists():
        logger.warning("Bundled manifest subset is missing for %s: %s", gene_name, manifest_path)
        return {}

    try:
        manifest = pd.read_csv(manifest_path)
    except Exception:
        logger.exception("Failed to read bundled manifest subset for %s", gene_name)
        return {}

    manifest = manifest.rename(columns={"IlmnID": "probe_id", "CHR": "chrom", "MAPINFO": "pos"})
    if "probe_id" not in manifest.columns:
        return {}

    subset = manifest[manifest["probe_id"].isin(set(probe_ids))].copy()
    if subset.empty:
        return {}

    subset = subset.drop_duplicates(subset=["probe_id"], keep="first")
    return {str(row["probe_id"]): row.to_dict() for _, row in subset.iterrows()}


def _format_probe_locus(annotation: dict[str, Any]) -> str | None:
    """Format a probe genomic locus from manifest or observed-row annotations."""
    chrom = annotation.get("chrom")
    pos = annotation.get("pos")
    if chrom is None or pos is None or pd.isna(chrom) or pd.isna(pos):
        return None
    chrom_text = str(chrom).removeprefix("chr")
    try:
        pos_value = int(float(pos))
    except (TypeError, ValueError):
        return None
    return f"chr{chrom_text}:{pos_value:,}"


def _build_nearby_manifest_variant_rows(
    annotation: dict[str, Any],
    *,
    allowed_variant_ids: set[str] | None = None,
) -> list[dict[str, Any]]:
    """Return nearby manifest SNP annotations paired with their reported distances."""
    variant_ids = _split_semicolon_tokens(annotation.get("SNP_ID"))
    distance_tokens = _split_semicolon_tokens(annotation.get("SNP_DISTANCE"))
    allowed_lookup = (
        {_normalize_lookup_key(variant_id) for variant_id in allowed_variant_ids if variant_id}
        if allowed_variant_ids is not None
        else None
    )
    nearby_rows: list[dict[str, Any]] = []
    for index, variant_id in enumerate(variant_ids):
        if allowed_lookup is not None and _normalize_lookup_key(variant_id) not in allowed_lookup:
            continue
        distance = distance_tokens[index] if index < len(distance_tokens) else ""
        nearby_rows.append(
            {
                "variant": variant_id,
                "distance": distance,
            }
        )
    return nearby_rows


def _collect_variant_record_papers(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Collect deduplicated papers and evidence links for one curated variant record."""
    papers: list[dict[str, Any]] = []
    seen_keys: set[tuple[str, str]] = set()

    for finding in _build_literature_findings(record):
        label = str(finding.get("paper", "")).strip()
        url = str(finding.get("url", "")).strip()
        if not label and not url:
            continue
        dedupe_key = (label, url)
        if dedupe_key in seen_keys:
            continue
        seen_keys.add(dedupe_key)
        papers.append(
            {
                "label": label or url,
                "url": url,
                "source_variant": record.get("display_name", record.get("variant", "")),
            }
        )

    for source in record.get("evidence", []):
        label = str(source.get("label", "")).strip()
        url = str(source.get("url", "")).strip()
        if not label and not url:
            continue
        dedupe_key = (label, url)
        if dedupe_key in seen_keys:
            continue
        seen_keys.add(dedupe_key)
        papers.append(
            {
                "label": label or url,
                "url": url,
                "source_variant": record.get("display_name", record.get("variant", "")),
            }
        )

    return papers


def _build_whitelist_probe_reference_rows(
    methylation: pd.DataFrame,
    knowledge_base: dict[str, Any],
    *,
    matched_variant_ids: set[str] | None = None,
) -> list[dict[str, Any]]:
    """Map each whitelist probe to curated variants, nearby loci, and supporting papers."""
    gene_context = knowledge_base.get("gene_context", {})
    gene_name = str(gene_context.get("gene_name", DEFAULT_GENE_NAME))
    relevant_probe_ids = list(gene_context.get("relevant_methylation_probe_ids", []))
    if not relevant_probe_ids:
        return []

    observed_rows = (
        methylation.drop_duplicates(subset=["probe_id"], keep="first")
        if "probe_id" in methylation.columns
        else methylation.iloc[0:0].copy()
    )
    observed_lookup = (
        {str(row["probe_id"]): row.to_dict() for _, row in observed_rows.iterrows()}
        if "probe_id" in observed_rows.columns
        else {}
    )
    manifest_lookup = _load_gene_manifest_probe_lookup(gene_name, relevant_probe_ids)
    variant_records = knowledge_base.get("variant_records", [])
    matched_variant_lookup = (
        {_normalize_lookup_key(variant_id) for variant_id in matched_variant_ids if variant_id}
        if matched_variant_ids is not None
        else None
    )

    reference_rows: list[dict[str, Any]] = []
    for probe_id in relevant_probe_ids:
        observed_annotation = observed_lookup.get(probe_id, {})
        manifest_annotation = manifest_lookup.get(probe_id, {})
        merged_annotation = {**manifest_annotation, **observed_annotation}
        observed_beta = pd.to_numeric(observed_annotation.get("beta"), errors="coerce")

        linked_variant_records = [
            record
            for record in variant_records
            if probe_id in record.get("relevant_methylation_probe_ids", [])
            and (
                matched_variant_lookup is None
                or _normalize_lookup_key(str(record.get("variant", ""))) in matched_variant_lookup
            )
        ]
        linked_variants = []
        for record in linked_variant_records:
            label = str(record.get("display_name", record.get("variant", probe_id))).strip()
            common_name = str(record.get("common_name", "")).strip()
            locus = None
            position = record.get("position")
            chromosome = record.get("chromosome")
            if chromosome is not None and position is not None and not pd.isna(position):
                try:
                    locus = f"chr{str(chromosome).removeprefix('chr')}:{int(position):,}"
                except (TypeError, ValueError):
                    locus = None
            linked_variants.append(
                {
                    "label": label,
                    "common_name": common_name,
                    "locus": locus,
                }
            )

        papers: list[dict[str, Any]] = []
        seen_papers: set[tuple[str, str, str]] = set()
        for record in linked_variant_records:
            for paper in _collect_variant_record_papers(record):
                dedupe_key = (
                    paper.get("label", ""),
                    paper.get("url", ""),
                    paper.get("source_variant", ""),
                )
                if dedupe_key in seen_papers:
                    continue
                seen_papers.add(dedupe_key)
                papers.append(paper)

        nearby_manifest_variants = _build_nearby_manifest_variant_rows(
            merged_annotation,
            allowed_variant_ids=matched_variant_ids,
        )
        if matched_variant_lookup is not None and not (
            linked_variants or nearby_manifest_variants or papers
        ):
            continue

        reference_rows.append(
            {
                "probe_id": probe_id,
                "observed_in_run": probe_id in observed_lookup,
                "beta": _round_beta(float(observed_beta)) if pd.notna(observed_beta) else None,
                "probe_locus": _format_probe_locus(merged_annotation),
                "linked_variants": linked_variants,
                "nearby_manifest_variants": nearby_manifest_variants,
                "papers": papers,
            }
        )

    return reference_rows


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
                "variant_label": _format_variant_label(row.get("id")),
                "change": _format_variant_change(row.get("ref"), row.get("alt")),
                "reference_allele": "" if pd.isna(row.get("ref")) else str(row.get("ref")).strip(),
                "alternate_allele": "" if pd.isna(row.get("alt")) else str(row.get("alt")).strip(),
                "linked_to": _build_specific_variant_link_summary(
                    matched_record,
                    fallback=f"No curated local {gene_name} link is bundled for this PASS variant yet.",
                ),
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
    curated_named_markers = _build_curated_named_marker_catalog(variants, variant_records)
    observed_named_marker_count = sum(1 for item in curated_named_markers if item["observed_in_run"])
    if curated_named_markers:
        curated_named_markers_summary = (
            f"The local {gene_name} bundle seeds {len(curated_named_markers)} curated named marker(s). "
            f"The current run directly matched {observed_named_marker_count} of them."
        )
    else:
        curated_named_markers_summary = (
            f"No curated named markers are bundled for {gene_name} yet."
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
        "curated_named_markers": curated_named_markers,
        "curated_named_markers_summary": curated_named_markers_summary,
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
    methylation: pd.DataFrame,
    knowledge_base: dict[str, Any],
    *,
    matched_variant_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Build a gene-level methylation interpretation from the current probe table."""
    gene_context = knowledge_base.get("gene_context", {})
    gene_name = str(gene_context.get("gene_name", DEFAULT_GENE_NAME))
    relevant_probe_ids = list(gene_context.get("relevant_methylation_probe_ids", []))
    relevant_probe_lookup = set(relevant_probe_ids)

    observed_relevant = (
        methylation[methylation["probe_id"].isin(relevant_probe_lookup)].copy()
        if relevant_probe_lookup
        else methylation.iloc[0:0].copy()
    )
    gene_named_rows, gene_named_match_columns = _select_gene_named_methylation_rows(
        methylation, gene_name
    )

    whitelist_beta_values = _coerce_numeric_beta_values(observed_relevant)
    gene_named_beta_values = _coerce_numeric_beta_values(gene_named_rows)
    full_table_beta_values = _coerce_numeric_beta_values(methylation)

    whitelist_mean_beta = _mean_beta_or_none(whitelist_beta_values)
    gene_name_mean_beta = _mean_beta_or_none(gene_named_beta_values)
    raw_mean_beta = _mean_beta_or_none(full_table_beta_values)

    primary_mean_beta = (
        whitelist_mean_beta
        if whitelist_mean_beta is not None
        else gene_name_mean_beta
        if gene_name_mean_beta is not None
        else raw_mean_beta
    )
    beta_band = _categorize_beta(primary_mean_beta)
    beta_band_source_label = (
        "Whitelist mean beta"
        if whitelist_mean_beta is not None
        else f"{gene_name}-named row mean beta"
        if gene_name_mean_beta is not None
        else "All numeric-row mean beta"
    )

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

    whitelist_summary = (
        f"the whitelist mean beta is {whitelist_mean_beta:.2f} from {len(whitelist_beta_values)} numeric whitelist probe(s)"
        if whitelist_mean_beta is not None
        else "the whitelist mean beta is unavailable because no observed whitelist probe carried a numeric beta value"
    )
    gene_named_summary = (
        f"the {gene_name}-named row mean beta is {gene_name_mean_beta:.2f} from {len(gene_named_beta_values)} numeric row(s)"
        if gene_name_mean_beta is not None
        else f"the {gene_name}-named row mean beta is unavailable because no numeric row explicitly mentioned {gene_name}"
    )
    raw_summary = (
        f"the all-numeric-row mean beta is {raw_mean_beta:.2f} from {len(full_table_beta_values)} numeric row(s)"
        if raw_mean_beta is not None
        else "the all-numeric-row mean beta is unavailable because the table did not contain any numeric beta values"
    )

    summary = (
        f"The current run observed {len(observed_relevant)} of {len(relevant_probe_ids)} curated whitelist probe(s), "
        f"{len(gene_named_rows)} row(s) whose gene annotation explicitly mentions {gene_name}, "
        f"and {len(methylation)} row(s) in the full methylation table. "
        f"Across those views, {whitelist_summary}; {gene_named_summary}; and {raw_summary}. "
        f"The observed whitelist subset is dominated by {group_summary}. "
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
    observed_relevant_lookup = set(observed_relevant["probe_id"].tolist()) if "probe_id" in observed_relevant.columns else set()
    gene_name_match_rule = (
        f"Rows count toward the {gene_name}-named mean when {gene_name} appears as a semicolon-delimited token in "
        + ", ".join(gene_named_match_columns)
        + "."
        if gene_named_match_columns
        else f"No gene-annotation column in the current methylation table explicitly mentioned {gene_name}."
    )
    whitelist_probe_reference_rows = _build_whitelist_probe_reference_rows(
        methylation,
        knowledge_base,
        matched_variant_ids=matched_variant_ids,
    )
    summary_metric_rows = _build_methylation_metric_rows(
        gene_name=gene_name,
        whitelist_mean_beta=whitelist_mean_beta,
        whitelist_count=len(whitelist_beta_values),
        gene_name_mean_beta=gene_name_mean_beta,
        gene_name_count=len(gene_named_beta_values),
        raw_mean_beta=raw_mean_beta,
        raw_count=len(full_table_beta_values),
    )
    if matched_variant_ids is None:
        whitelist_probe_reference_summary = (
            "Each whitelist probe is cross-referenced to any curated variant records that explicitly cite it, "
            "plus nearby manifest SNP annotations and the bundled papers behind those variant interpretations."
        )
    elif matched_variant_ids:
        whitelist_probe_reference_summary = (
            "Each whitelist probe is cross-referenced only to curated variant records observed in the current run, "
            "plus any matching manifest SNP annotations and bundled papers for those same observed variants."
        )
    else:
        whitelist_probe_reference_summary = (
            "No curated variant observed in the current run has a probe-specific reference map, so the probe-to-variant table is hidden."
        )

    return {
        "gene_name": gene_name,
        "clinical_context": gene_context.get("clinical_context", ""),
        "summary": summary,
        "mean_beta": _round_beta(whitelist_mean_beta),
        "mean_beta_label": "Whitelist mean beta",
        "mean_beta_probe_count": int(len(whitelist_beta_values)),
        "whitelist_mean_beta": _round_beta(whitelist_mean_beta),
        "whitelist_mean_beta_label": "Whitelist mean beta",
        "whitelist_mean_beta_probe_count": int(len(whitelist_beta_values)),
        "whitelist_probe_count": len(relevant_probe_ids),
        "whitelist_observed_probe_count": int(len(observed_relevant)),
        "whitelist_explanation": _build_whitelist_explanation(gene_name, relevant_probe_ids),
        "whitelist_literature_context": gene_context.get("methylation_interpretation", ""),
        "whitelist_probe_statuses": [
            {
                "probe_id": probe_id,
                "observed_in_run": probe_id in observed_relevant_lookup,
            }
            for probe_id in relevant_probe_ids
        ],
        "gene_name_mean_beta": _round_beta(gene_name_mean_beta),
        "gene_name_mean_beta_label": f"{gene_name}-named row mean beta",
        "gene_name_mean_beta_probe_count": int(len(gene_named_beta_values)),
        "gene_name_row_count": int(len(gene_named_rows)),
        "gene_name_match_columns": gene_named_match_columns,
        "gene_name_match_rule": gene_name_match_rule,
        "raw_mean_beta": _round_beta(raw_mean_beta),
        "raw_mean_beta_label": "All numeric-row mean beta",
        "raw_probe_count": int(len(methylation)),
        "raw_mean_beta_probe_count": int(len(full_table_beta_values)),
        "all_numeric_mean_beta": _round_beta(raw_mean_beta),
        "all_numeric_mean_beta_label": "All numeric-row mean beta",
        "all_numeric_mean_beta_probe_count": int(len(full_table_beta_values)),
        "summary_metric_rows": summary_metric_rows,
        "beta_band": beta_band,
        "beta_band_source_label": beta_band_source_label,
        "observed_probe_count": int(len(observed_relevant)),
        "curated_probe_count": len(relevant_probe_ids),
        "probe_ids": observed_relevant["probe_id"].tolist(),
        "group_breakdown": group_breakdown,
        "methylation_effects": gene_context.get("methylation_effects", []),
        "methylation_condition_research": gene_context.get("methylation_condition_research", []),
        "evidence": gene_context.get("evidence", []),
        "probe_preview": observed_relevant[available_preview_columns].copy(),
        "whitelist_probe_reference_rows": whitelist_probe_reference_rows,
        "whitelist_probe_reference_summary": whitelist_probe_reference_summary,
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
            "result_table_rows": [
                {
                    "variant_label": item.get("variant_label", "None"),
                    "change": item.get("change", "Unavailable"),
                    "linked_to": item.get("linked_to", ""),
                }
                for item in observed_summaries
            ],
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
        "curated_named_markers": [],
        "curated_named_markers_summary": (
            f"No curated named markers are bundled for {gene_name} yet."
        ),
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
    gene_named_rows, gene_named_match_columns = _select_gene_named_methylation_rows(
        methylation, gene_name
    )
    gene_named_beta_values = _coerce_numeric_beta_values(gene_named_rows)
    beta_values = _coerce_numeric_beta_values(methylation)
    gene_name_mean_beta = _mean_beta_or_none(gene_named_beta_values)
    mean_beta = _mean_beta_or_none(beta_values)
    primary_mean_beta = gene_name_mean_beta if gene_name_mean_beta is not None else mean_beta
    beta_band = _categorize_beta(primary_mean_beta)
    beta_band_source_label = (
        f"{gene_name}-named row mean beta"
        if gene_name_mean_beta is not None
        else "All numeric-row mean beta"
    )

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
        f"The current run captured {len(methylation)} probe(s) from the selected {gene_name} manifest subset, "
        f"with a {gene_name}-named row mean beta of {gene_name_mean_beta:.2f} from {len(gene_named_beta_values)} numeric row(s) "
        f"and an all-numeric-row mean beta of {mean_beta:.2f} from {len(beta_values)} numeric row(s). "
        if mean_beta is not None and gene_name_mean_beta is not None
        else f"The current run captured {len(methylation)} probe(s) from the selected {gene_name} manifest subset "
        f"with an all-numeric-row mean beta of {mean_beta:.2f}. "
        if mean_beta is not None
        else f"The current run captured {len(methylation)} probe(s) from the selected {gene_name} manifest subset. "
    )
    gene_name_match_rule = (
        f"Rows count toward the {gene_name}-named mean when {gene_name} appears as a semicolon-delimited token in "
        + ", ".join(gene_named_match_columns)
        + "."
        if gene_named_match_columns
        else f"No gene-annotation column in the current methylation table explicitly mentioned {gene_name}."
    )
    summary_metric_rows = _build_methylation_metric_rows(
        gene_name=gene_name,
        whitelist_mean_beta=None,
        whitelist_count=0,
        gene_name_mean_beta=gene_name_mean_beta,
        gene_name_count=len(gene_named_beta_values),
        raw_mean_beta=mean_beta,
        raw_count=len(beta_values),
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
        "mean_beta": _round_beta(gene_name_mean_beta if gene_name_mean_beta is not None else mean_beta),
        "mean_beta_label": (
            f"{gene_name}-named row mean beta"
            if gene_name_mean_beta is not None
            else "All numeric-row mean beta"
        ),
        "mean_beta_probe_count": int(
            len(gene_named_beta_values) if gene_name_mean_beta is not None else len(beta_values)
        ),
        "whitelist_mean_beta": None,
        "whitelist_mean_beta_label": "Whitelist mean beta",
        "whitelist_mean_beta_probe_count": 0,
        "whitelist_probe_count": 0,
        "whitelist_observed_probe_count": 0,
        "whitelist_explanation": _build_whitelist_explanation(gene_name, []),
        "whitelist_literature_context": "",
        "whitelist_probe_statuses": [],
        "gene_name_mean_beta": _round_beta(gene_name_mean_beta),
        "gene_name_mean_beta_label": f"{gene_name}-named row mean beta",
        "gene_name_mean_beta_probe_count": int(len(gene_named_beta_values)),
        "gene_name_row_count": int(len(gene_named_rows)),
        "gene_name_match_columns": gene_named_match_columns,
        "gene_name_match_rule": gene_name_match_rule,
        "raw_mean_beta": _round_beta(mean_beta),
        "raw_mean_beta_label": "All numeric-row mean beta",
        "raw_probe_count": int(len(methylation)),
        "raw_mean_beta_probe_count": int(len(beta_values)),
        "all_numeric_mean_beta": _round_beta(mean_beta),
        "all_numeric_mean_beta_label": "All numeric-row mean beta",
        "all_numeric_mean_beta_probe_count": int(len(beta_values)),
        "summary_metric_rows": summary_metric_rows,
        "beta_band": beta_band,
        "beta_band_source_label": beta_band_source_label,
        "observed_probe_count": int(len(methylation)),
        "curated_probe_count": 0,
        "probe_ids": methylation["probe_id"].tolist() if "probe_id" in methylation.columns else [],
        "group_breakdown": group_breakdown,
        "methylation_effects": [
            f"The current methylation summary for {gene_name} is based on the filtered manifest subset selected during preprocessing.",
            "Beta values are being shown as region-level context until a curated gene-specific methylation knowledge base is added.",
        ],
        "methylation_condition_research": [],
        "evidence": [],
        "probe_preview": methylation[available_preview_columns].copy() if available_preview_columns else methylation.head(12).copy(),
        "whitelist_probe_reference_rows": [],
        "whitelist_probe_reference_summary": (
            f"No curated whitelist-probe reference map is bundled for {gene_name} yet."
        ),
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


def _categorize_predictive_beta_band(mean_beta: float | None) -> str:
    """Collapse UI methylation bands into the three predictive case buckets."""
    descriptive_band = _categorize_beta(mean_beta)
    if descriptive_band == "unavailable":
        return "unavailable"
    if descriptive_band == "low":
        return "low"
    if descriptive_band == "high":
        return "high"
    return "medium"


def _format_predictive_beta_display(mean_beta: float | None) -> str:
    """Render predictive beta values consistently for the UI."""
    rounded = _round_beta(mean_beta)
    return str(rounded) if rounded is not None else "Unavailable"


def _summarize_predictive_observed_variants(
    variant_interpretations: dict[str, Any],
) -> list[str]:
    """Return deduplicated human-readable labels for the current sample's observed variants."""
    observed_items: list[str] = []
    for item in variant_interpretations.get("sample_highlights", {}).get("highlight_items", []):
        title = str(item.get("title", "")).strip()
        observed_variant = str(item.get("observed_variant", "")).strip()
        change = str(item.get("change", "")).strip()
        has_change = bool(change and change.casefold() != "unavailable")
        if title and has_change:
            observed_items.append(f"{title} {change}")
        elif title and observed_variant and observed_variant != title:
            observed_items.append(f"{title} ({observed_variant})")
        elif title:
            observed_items.append(title)
        elif observed_variant:
            observed_items.append(observed_variant)
    return _dedupe_text_items(observed_items)


def _build_synthesis_variant_prediction_rule_lookup(
    synthesis_database: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Index concrete variant prediction rules by the labels seen in analysis output."""
    rule_lookup: dict[str, dict[str, Any]] = {}
    for rule in synthesis_database.get("variant_prediction_rules", []):
        lookup_candidates = [
            rule.get("variant"),
            rule.get("display_name"),
            rule.get("common_name"),
            *rule.get("lookup_keys", []),
        ]
        for candidate in lookup_candidates:
            candidate_text = str(candidate or "").strip()
            if not candidate_text:
                continue
            rule_lookup[_normalize_lookup_key(candidate_text)] = rule
    return rule_lookup


def _find_synthesis_variant_prediction_rule(
    record: dict[str, Any],
    rule_lookup: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    """Return the concrete prediction rule for a matched variant record when available."""
    lookup_candidates = [
        record.get("variant"),
        record.get("variant_label"),
        record.get("observed_variant"),
    ]
    for candidate in lookup_candidates:
        candidate_text = str(candidate or "").strip()
        if not candidate_text:
            continue
        matched_rule = rule_lookup.get(_normalize_lookup_key(candidate_text))
        if matched_rule is not None:
            return matched_rule
    return None


def _format_predictive_observed_signal(record: dict[str, Any]) -> str:
    """Render the observed variant and its actual REF -> ALT change together."""
    observed_signal = str(record.get("observed_variant", record.get("variant", ""))).strip()
    change = str(record.get("change", "")).strip()
    if change and change.casefold() != "unavailable" and change not in observed_signal:
        return f"{observed_signal} ({change})" if observed_signal else change
    return observed_signal


def _render_sample_change_template(
    template: str,
    *,
    record: dict[str, Any],
    concrete_rule: dict[str, Any],
) -> str:
    """Fill a controlled sample-change template from the synthesis database."""
    replacements = {
        "{change}": str(record.get("change", "Unavailable")).strip() or "Unavailable",
        "{variant}": str(record.get("variant", concrete_rule.get("variant", ""))).strip(),
        "{display_name}": str(concrete_rule.get("display_name", record.get("variant", ""))).strip(),
        "{observed_variant}": str(record.get("observed_variant", "")).strip(),
        "{alt_allele}": _extract_alt_allele_from_change(record.get("change")),
    }
    rendered = str(template or "")
    for placeholder, value in replacements.items():
        rendered = rendered.replace(placeholder, value)
    return rendered.strip()


def _select_synthesis_prediction_for_record(
    record: dict[str, Any],
    concrete_rule: dict[str, Any] | None,
) -> dict[str, str]:
    """Choose the most sample-specific prediction text for one matched variant."""
    if concrete_rule is None:
        return {
            "prediction": str(record.get("clinical_interpretation", "")).strip(),
            "research_focus": "; ".join(
                _dedupe_text_items(record.get("associated_conditions", []))[:3]
            ),
            "source": str(record.get("interpretation_scope", "Curated marker")).strip(),
        }

    change = str(record.get("change", "")).strip()
    normalized_change = _normalize_allele_change(change)
    observed_alt_allele = _extract_alt_allele_from_change(change)

    for allele_rule in concrete_rule.get("allele_change_rules", []):
        rule_change = _normalize_allele_change(allele_rule.get("change"))
        rule_alt_allele = str(allele_rule.get("alt_allele", "")).strip().upper()
        change_matches = bool(rule_change and normalized_change and rule_change == normalized_change)
        alt_matches = bool(rule_alt_allele and observed_alt_allele and rule_alt_allele == observed_alt_allele)
        if not change_matches and not alt_matches:
            continue

        prediction = str(allele_rule.get("prediction", "")).strip()
        if not prediction:
            prediction = str(concrete_rule.get("prediction", "")).strip()
        research_focus = str(allele_rule.get("basis", "")).strip()
        if not research_focus:
            research_focus = str(concrete_rule.get("basis", "")).strip()
        return {
            "prediction": prediction,
            "research_focus": research_focus,
            "source": "Sample allele-change thesis",
        }

    prediction = str(concrete_rule.get("prediction", "")).strip()
    research_focus = str(concrete_rule.get("basis", "")).strip()
    sample_change_template = str(concrete_rule.get("sample_change_template", "")).strip()
    if change and change.casefold() != "unavailable" and sample_change_template:
        change_anchor = _render_sample_change_template(
            sample_change_template,
            record=record,
            concrete_rule=concrete_rule,
        )
        if change_anchor and prediction:
            prediction = f"{change_anchor} {prediction}"
        elif change_anchor:
            prediction = change_anchor
        return {
            "prediction": prediction,
            "research_focus": research_focus,
            "source": "Sample change-anchored thesis",
        }

    return {
        "prediction": prediction,
        "research_focus": research_focus,
        "source": "Concrete variant thesis",
    }


def build_predictive_theses(
    *,
    variant_interpretations: dict[str, Any],
    methylation_insights: dict[str, Any],
    knowledge_base: dict[str, Any] | None = None,
    synthesis_database: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the predictive-thesis payload shown after a completed analysis run."""
    knowledge_base = knowledge_base or {}
    synthesis_database = synthesis_database or {}
    gene_context = knowledge_base.get("gene_context", {})
    gene_name = str(
        synthesis_database.get("gene_name")
        or variant_interpretations.get("gene_name")
        or methylation_insights.get("gene_name")
        or gene_context.get("gene_name")
        or DEFAULT_GENE_NAME
    ).strip() or DEFAULT_GENE_NAME

    case_catalog = [
        case
        for case in synthesis_database.get("cases", [])
        if str(case.get("case_id", "")).strip()
    ]
    case_lookup = {
        str(case["case_id"]).strip(): case
        for case in case_catalog
    }
    variant_prediction_rule_lookup = _build_synthesis_variant_prediction_rule_lookup(
        synthesis_database
    )

    matched_records = variant_interpretations.get("matched_records", [])
    highlight_items = variant_interpretations.get("sample_highlights", {}).get("highlight_items", [])
    sample_variant_rows = variant_interpretations.get("sample_highlights", {}).get("result_table_rows", [])
    observed_variant_labels = _summarize_predictive_observed_variants(variant_interpretations)
    variant_found = bool(sample_variant_rows)

    variant_case = case_lookup.get("gene_variant_found")
    variant_summary = (
        str(variant_case.get("prediction", "")).strip()
        if variant_found and variant_case is not None
        else ""
    )
    if not variant_summary and variant_found:
        variant_summary = (
            f"Observed {gene_name} variation is present in this sample, so the gene-level predictive thesis "
            "should be read as locus-specific research context rather than as a stand-alone diagnosis."
        )
    if not variant_found:
        variant_summary = (
            f"No promoter or gene-body {gene_name} variant was visible in the current preview, so the "
            "variant-gated predictive thesis cases did not match this sample."
        )
    if observed_variant_labels:
        variant_summary = (
            f"{variant_summary} Observed sample signal: {', '.join(observed_variant_labels[:4])}."
        ).strip()

    variant_prediction_rows: list[dict[str, str]] = []
    if variant_found and variant_case is not None:
        variant_prediction_rows.append(
            {
                "observed_signal": (
                    ", ".join(observed_variant_labels[:3])
                    if observed_variant_labels
                    else f"{gene_name} interval variant observed"
                ),
                "source": "Gene-level thesis",
                "prediction": str(variant_case.get("prediction", "")).strip(),
                "research_focus": "; ".join(
                    _dedupe_text_items(variant_case.get("research_focus", []))[:3]
                ),
            }
        )

    for record in matched_records:
        concrete_rule = _find_synthesis_variant_prediction_rule(
            record,
            variant_prediction_rule_lookup,
        )
        selected_prediction = _select_synthesis_prediction_for_record(record, concrete_rule)

        variant_prediction_rows.append(
            {
                "observed_signal": _format_predictive_observed_signal(record),
                "source": selected_prediction["source"],
                "prediction": selected_prediction["prediction"],
                "research_focus": selected_prediction["research_focus"],
            }
        )

    if not matched_records:
        for item in highlight_items:
            variant_prediction_rows.append(
                {
                    "observed_signal": str(item.get("title", item.get("observed_variant", ""))).strip(),
                    "source": str(item.get("category", "Interval variant")).strip(),
                    "prediction": str(item.get("description", "")).strip(),
                    "research_focus": "; ".join(
                        _dedupe_text_items(item.get("conditions", []))[:3]
                    ),
                }
            )

    methylation_source_rows = [
        {
            "metric_key": "whitelist",
            "label": str(methylation_insights.get("whitelist_mean_beta_label", "Whitelist mean beta")).strip(),
            "mean_beta": methylation_insights.get("whitelist_mean_beta"),
            "probe_count": int(methylation_insights.get("whitelist_mean_beta_probe_count", 0) or 0),
        },
        {
            "metric_key": "gene_name_related",
            "label": str(
                methylation_insights.get("gene_name_mean_beta_label", f"{gene_name}-named row mean beta")
            ).strip(),
            "mean_beta": methylation_insights.get("gene_name_mean_beta"),
            "probe_count": int(methylation_insights.get("gene_name_mean_beta_probe_count", 0) or 0),
        },
        {
            "metric_key": "all_numeric",
            "label": str(
                methylation_insights.get("all_numeric_mean_beta_label", "All numeric-row mean beta")
            ).strip(),
            "mean_beta": methylation_insights.get("all_numeric_mean_beta"),
            "probe_count": int(methylation_insights.get("all_numeric_mean_beta_probe_count", 0) or 0),
        },
    ]

    methylation_prediction_rows: list[dict[str, Any]] = []
    matched_cases: list[dict[str, str]] = []

    if variant_found and variant_case is not None:
        matched_cases.append(
            {
                "case_label": str(variant_case.get("label", "Gene variant found")).strip(),
                "trigger": "Observed promoter or gene-body variant",
                "source": "Variant-only synthesis",
                "mean_beta_display": "n/a",
                "band": "n/a",
                "prediction": str(variant_case.get("prediction", "")).strip(),
                "research_focus": "; ".join(
                    _dedupe_text_items(variant_case.get("research_focus", []))[:3]
                ),
            }
        )

    for source_row in methylation_source_rows:
        mean_beta = source_row["mean_beta"]
        band = _categorize_predictive_beta_band(mean_beta)
        case_id = (
            f"gene_variant_found__{source_row['metric_key']}__{band}"
            if band != "unavailable"
            else ""
        )
        case = case_lookup.get(case_id) if case_id else None
        matched = variant_found and case is not None

        if mean_beta is None:
            prediction = f"No numeric beta value was available for {source_row['label'].lower()}."
        elif not variant_found:
            prediction = (
                f"{source_row['label']} was computed, but the predictive thesis matrix only matches after "
                f"a {gene_name} variant is observed."
            )
        elif case is not None:
            prediction = str(case.get("prediction", "")).strip()
        else:
            prediction = (
                f"No bundled predictive thesis case is available for {source_row['label']} with a {band} methylation band."
            )

        research_focus_items = (
            _dedupe_text_items(case.get("research_focus", []))[:3]
            if case is not None
            else []
        )

        methylation_prediction_rows.append(
            {
                "metric_key": source_row["metric_key"],
                "metric_label": source_row["label"],
                "mean_beta": _round_beta(mean_beta),
                "mean_beta_display": _format_predictive_beta_display(mean_beta),
                "probe_count": source_row["probe_count"],
                "band": band,
                "band_display": band.title() if band != "unavailable" else "Unavailable",
                "prediction": prediction,
                "matched": matched,
                "matched_case_label": str(case.get("label", "")).strip() if case is not None else "",
                "research_focus": "; ".join(research_focus_items),
            }
        )

        if matched and case is not None:
            matched_cases.append(
                {
                    "case_label": str(case.get("label", "")).strip(),
                    "trigger": f"Variant found + {source_row['label']}",
                    "source": source_row["label"],
                    "mean_beta_display": _format_predictive_beta_display(mean_beta),
                    "band": band.title(),
                    "prediction": str(case.get("prediction", "")).strip(),
                    "research_focus": "; ".join(research_focus_items),
                }
            )

    if matched_cases:
        summary = (
            f"{gene_name} matched {len(matched_cases)} predictive thesis case(s) in this run: "
            f"the base variant case plus {max(len(matched_cases) - 1, 0)} methylation-linked case(s)."
        )
    elif variant_found:
        summary = (
            f"{gene_name} variation was observed, but none of the bundled predictive thesis cases could be matched "
            "to the available methylation values."
        )
    else:
        summary = (
            f"No predictive thesis case matched because the current {gene_name} run did not surface a promoter or gene-body variant."
        )

    return {
        "gene_name": gene_name,
        "database_name": synthesis_database.get(
            "database_name",
            f"No curated {gene_name} predictive synthesis database loaded",
        ),
        "database_version": synthesis_database.get("version", "generic"),
        "matching_rule": synthesis_database.get(
            "matching_rule",
            "Cases match only when a gene variant is present, plus the requested methylation source resolves to a low, medium, or high beta band.",
        ),
        "disclaimer": synthesis_database.get(
            "disclaimer",
            "Predictive theses in this app are literature-guided research summaries, not diagnostic claims.",
        ),
        "seeded_markers": synthesis_database.get("seeded_markers", []),
        "concrete_variant_prediction": synthesis_database.get("concrete_variant_prediction", ""),
        "variant_found": variant_found,
        "variant_found_label": "Yes" if variant_found else "No",
        "variant_summary": variant_summary,
        "variant_prediction_rows": variant_prediction_rows,
        "methylation_prediction_rows": methylation_prediction_rows,
        "matched_cases": matched_cases,
        "matched_case_count": len(matched_cases),
        "case_catalog_size": len(case_catalog),
        "summary": summary,
    }


def _join_unique_database_values(values: list[Any]) -> str:
    """Join compact central-database values while preserving order."""
    cleaned_values: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = "" if value is None else str(value).strip()
        if not text or text.casefold() in seen:
            continue
        seen.add(text.casefold())
        cleaned_values.append(text)
    return "; ".join(cleaned_values)


def _format_database_beta(value: Any) -> float | str:
    """Return beta values as stable numeric CSV fields when available."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
        return round(float(value), 6)
    except (TypeError, ValueError):
        return str(value).strip()


def _format_database_quality(value: Any) -> str:
    """Render VCF quality values for the central variant-level database."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return str(value).strip()


def _format_database_position(value: Any) -> str:
    """Render a stable genomic coordinate for a central database row."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
        return str(int(value))
    except (TypeError, ValueError):
        return str(value).strip()


def _build_observed_variant_key(row: pd.Series) -> str:
    """Return a biologically stable key for one observed variant allele."""
    chrom = str(row.get("chrom", "")).strip().removeprefix("chr")
    pos = _format_database_position(row.get("pos"))
    ref = "" if pd.isna(row.get("ref")) else str(row.get("ref")).strip()
    alt = "" if pd.isna(row.get("alt")) else str(row.get("alt")).strip()
    if chrom and pos and (ref or alt):
        return f"chr{chrom}:{pos}:{ref}>{alt}"
    return _format_variant_display(row)


def _format_variant_location_for_database(row: pd.Series) -> str:
    """Render the point location used for variant-level central database rows."""
    chrom = str(row.get("chrom", "")).strip().removeprefix("chr")
    pos = _format_database_position(row.get("pos"))
    if chrom and pos:
        return f"chr{chrom}:{int(pos):,}" if pos.isdigit() else f"chr{chrom}:{pos}"
    return ""


def _build_general_analysis_database_rows(
    *,
    gene_name: str,
    variants: pd.DataFrame,
    variant_interpretations: dict[str, Any],
    methylation_insights: dict[str, Any],
) -> list[dict[str, Any]]:
    """Build central database rows with one entry per observed variant."""
    normalized_gene_name = gene_name.strip().upper() or DEFAULT_GENE_NAME

    gene_region = variant_interpretations.get("gene_region", {})
    search_region = variant_interpretations.get("search_region", {})
    gene_location = str(
        gene_region.get("display")
        or search_region.get("display")
        or variant_interpretations.get("recommended_promoter_plus_gene_region")
        or ""
    ).strip()

    matched_records_by_variant = {
        str(record.get("observed_variant", "")).strip(): record
        for record in variant_interpretations.get("matched_records", [])
        if str(record.get("observed_variant", "")).strip()
    }

    output_rows: list[dict[str, Any]] = []
    for _, row in variants.iterrows():
        observed_variant = _format_variant_display(row)
        matched_record = matched_records_by_variant.get(observed_variant, {})
        variant_label = _format_variant_label(row.get("id"))
        chrom = str(row.get("chrom", "")).strip().removeprefix("chr")
        position = _format_database_position(row.get("pos"))
        output_rows.append(
            {
                "gene": normalized_gene_name,
                "variant key": _build_observed_variant_key(row),
                "observed gene variant": observed_variant,
                "gene variant label": variant_label,
                "change": _format_variant_change(row.get("ref"), row.get("alt")),
                "chromosome": f"chr{chrom}" if chrom else "",
                "position": position,
                "variant location": _format_variant_location_for_database(row),
                "gene location": gene_location,
                "source": "VCF",
                "(VCF) quality (qual)": _format_database_quality(row.get("qual")),
                "matched curated marker": str(matched_record.get("variant", "")).strip(),
                "variant interpretation scope": str(
                    matched_record.get("interpretation_scope", "Unclassified observed variant")
                ).strip(),
                "curated biological significance": str(
                    matched_record.get("clinical_significance")
                    or matched_record.get("clinical_interpretation")
                    or f"No curated local {normalized_gene_name} significance is bundled for this observed variant."
                ).strip(),
                "functional effects": _join_unique_database_values(
                    list(matched_record.get("functional_effects") or [])
                ),
                "associated conditions": _join_unique_database_values(
                    list(matched_record.get("associated_conditions") or [])
                ),
                "methylation-linked probes": _join_unique_database_values(
                    list(matched_record.get("relevant_probe_ids") or [])
                ),
                "mean beta whitelist": _format_database_beta(methylation_insights.get("whitelist_mean_beta")),
                "mean beta related to gene": _format_database_beta(methylation_insights.get("gene_name_mean_beta")),
                "mean beta on found probes in the area (numerical rows)": _format_database_beta(
                    methylation_insights.get("all_numeric_mean_beta")
                ),
            }
        )
    return output_rows


def update_general_analysis_database(
    *,
    gene_name: str,
    variants: pd.DataFrame,
    variant_interpretations: dict[str, Any],
    methylation_insights: dict[str, Any],
    overwrite: bool = False,
    database_path: str | Path = GENERAL_ANALYSIS_DATABASE_PATH,
) -> dict[str, Any]:
    """Add or optionally replace observed variant rows in the central database."""
    output_path = Path(database_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    normalized_gene_name = gene_name.strip().upper() or DEFAULT_GENE_NAME
    new_rows = _build_general_analysis_database_rows(
        gene_name=normalized_gene_name,
        variants=variants,
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
    )
    if not new_rows:
        return {
            "action": "skipped_empty",
            "path": output_path,
            "message": f"Central database did not add {normalized_gene_name}; no observed variant rows were available.",
        }

    if output_path.exists():
        try:
            database = pd.read_csv(output_path, dtype=object)
        except pd.errors.EmptyDataError:
            database = pd.DataFrame(columns=GENERAL_ANALYSIS_DATABASE_COLUMNS)
    else:
        database = pd.DataFrame(columns=GENERAL_ANALYSIS_DATABASE_COLUMNS)

    for column in GENERAL_ANALYSIS_DATABASE_COLUMNS:
        if column not in database.columns:
            database[column] = ""

    primary_columns = list(GENERAL_ANALYSIS_DATABASE_COLUMNS)
    extra_columns = [column for column in database.columns if column not in primary_columns]
    database = database[primary_columns + extra_columns].fillna("")

    existing_mask = database["gene"].astype(str).str.strip().str.upper() == normalized_gene_name
    if overwrite and existing_mask.any():
        database = database.loc[~existing_mask].copy()
        database = pd.concat([database, pd.DataFrame(new_rows)], ignore_index=True)
        database = database[primary_columns + extra_columns].fillna("")
        database.to_csv(output_path, index=False)
        return {
            "action": "overwritten",
            "path": output_path,
            "message": (
                f"Overwrote {len(new_rows)} observed {normalized_gene_name} variant row(s) "
                "in the central analysis database."
            ),
        }

    legacy_gene_mask = existing_mask & (
        database["variant key"].astype(str).str.strip() == ""
    )
    if legacy_gene_mask.any():
        database = database.loc[~legacy_gene_mask].copy()
        existing_mask = database["gene"].astype(str).str.strip().str.upper() == normalized_gene_name

    existing_variant_keys = set(
        database.loc[existing_mask, "variant key"].astype(str).str.strip().str.casefold()
    )
    rows_to_add = [
        row
        for row in new_rows
        if str(row.get("variant key", "")).strip().casefold() not in existing_variant_keys
    ]

    if not rows_to_add:
        return {
            "action": "skipped_existing",
            "path": output_path,
            "message": (
                f"Central database already contains all {len(new_rows)} observed {normalized_gene_name} "
                "variant row(s); existing entries were kept. Check overwrite in Run Analysis to replace them."
            ),
        }

    database = pd.concat([database, pd.DataFrame(rows_to_add)], ignore_index=True)
    database = database[primary_columns + extra_columns].fillna("")
    database.to_csv(output_path, index=False)

    skipped_count = len(new_rows) - len(rows_to_add)
    skipped_note = f" Kept {skipped_count} existing row(s)." if skipped_count else ""
    return {
        "action": "added",
        "path": output_path,
        "message": (
            f"Added {len(rows_to_add)} observed {normalized_gene_name} variant row(s) "
            f"to the central analysis database.{skipped_note}"
        ),
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

    named_id_count = int(df["id"].notna().sum())
    logger.info(
        "Loaded %d PASS variants from %s in region %s (%d named IDs, %d unlabeled in source VCF)",
        len(df),
        vcf_path,
        region,
        named_id_count,
        len(df) - named_id_count,
    )
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


def _render_section_table(df: pd.DataFrame, title: str, rows: int | None = 20) -> str:
    """Render a DataFrame preview section for the HTML report."""
    preview_df = df if rows is None else df.head(rows)
    preview = preview_df.to_html(index=False, classes="data-table", border=0)
    return (
        f"<section><h2>{html.escape(title)}</h2>"
        f'<div class="report-table-shell">{preview}</div>'
        "</section>"
    )


def _with_preferred_column_order(df: pd.DataFrame, preferred_columns: list[str]) -> pd.DataFrame:
    """Move preferred columns to the front while preserving the remaining order."""
    ordered_columns = [column for column in preferred_columns if column in df.columns]
    ordered_columns.extend(column for column in df.columns if column not in ordered_columns)
    return df.loc[:, ordered_columns]


def _prepare_variant_table_for_output(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize missing IDs so exported reports stay readable."""
    preview_df = df.copy()
    if "id" in preview_df.columns:
        preview_df["id"] = preview_df["id"].where(preview_df["id"].notna(), "Unlabeled in source VCF")
        preview_df["id"] = preview_df["id"].replace({"": "Unlabeled in source VCF", ".": "Unlabeled in source VCF"})
    return _with_preferred_column_order(
        preview_df,
        ["chrom", "id", "pos", "ref", "alt", "qual", "filter_pass", "id_source"],
    )


def _prepare_methylation_table_for_output(df: pd.DataFrame) -> pd.DataFrame:
    """Apply stable column ordering to methylation tables for the UI and reports."""
    return _with_preferred_column_order(
        df.copy(),
        [
            "probe_id",
            "beta",
            "chrom",
            "pos",
            "ref",
            "alt",
            "GencodeBasicV12_NAME",
            "UCSC_RefGene_Name",
            "UCSC_RefGene_Group",
            "UCSC_CpG_Islands_Name",
            "Relation_to_UCSC_CpG_Island",
        ],
    )


def _flatten_report_value(value: Any) -> str:
    """Convert structured payload values into compact report-table text."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, list):
        flattened_items = [_flatten_report_value(item) for item in value]
        flattened_items = [item for item in flattened_items if item]
        return "; ".join(flattened_items)
    if isinstance(value, dict):
        if value.get("label") and value.get("url"):
            return f"{value['label']} ({value['url']})"
        if value.get("paper") and value.get("finding"):
            return f"{value['paper']}: {value['finding']}"
        if value.get("variant") and value.get("summary"):
            return f"{value['variant']}: {value['summary']}"
        if value.get("location_group") and value.get("label"):
            return f"{value['location_group']} - {value['label']}"
        return json.dumps(value, ensure_ascii=True)
    return str(value)


def _report_df_from_rows(
    rows: list[dict[str, Any]],
    column_map: list[tuple[str, str]],
) -> pd.DataFrame:
    """Build a report-ready DataFrame from structured rows and user-facing labels."""
    normalized_rows: list[dict[str, str]] = []
    for row in rows:
        normalized_rows.append(
            {
                label: _flatten_report_value(row.get(source_key))
                for source_key, label in column_map
            }
        )
    if not normalized_rows:
        return pd.DataFrame(columns=[label for _, label in column_map])
    return pd.DataFrame(normalized_rows)


def _render_report_paragraphs(
    title: str,
    paragraphs: list[str],
    *,
    extra_markup: str = "",
) -> str:
    """Render one report section made of paragraphs and optional nested markup."""
    rendered_paragraphs = "".join(
        f"<p>{html.escape(paragraph)}</p>" for paragraph in paragraphs if str(paragraph).strip()
    )
    if not rendered_paragraphs and not extra_markup:
        return ""
    return f"<section><h2>{html.escape(title)}</h2>{rendered_paragraphs}{extra_markup}</section>"


def _render_variant_interpretation_report(
    variant_interpretations: dict[str, Any],
    population_insights: dict[str, Any],
) -> str:
    """Render the richer variant interpretation block for the exported HTML report."""
    if not variant_interpretations:
        return ""

    nested_sections: list[str] = []
    sample_rows = variant_interpretations.get("sample_highlights", {}).get("result_table_rows", [])
    if sample_rows:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    sample_rows,
                    [
                        ("variant_label", "Variant label"),
                        ("change", "Change"),
                        ("linked_to", "Linked to"),
                    ],
                ),
                "Sample Results",
                rows=None,
            )
        )

    matched_records = variant_interpretations.get("matched_records", [])
    if matched_records:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    matched_records,
                    [
                        ("variant", "Curated variant"),
                        ("variant_label", "Variant label"),
                        ("observed_variant", "Observed in run"),
                        ("change", "Change"),
                        ("linked_to", "Linked to"),
                        ("interpretation_scope", "Scope"),
                        ("clinical_significance", "Clinical significance"),
                        ("clinical_interpretation", "Clinical interpretation"),
                        ("methylation_interpretation", "Methylation interpretation"),
                        ("associated_conditions", "Associated conditions"),
                        ("functional_effects", "Functional effects"),
                        ("research_context", "Research context"),
                        ("relevant_probe_ids", "Relevant probes"),
                        ("evidence", "Evidence"),
                    ],
                ),
                "Matched Variant Interpretations",
                rows=None,
            )
        )

    curated_markers = variant_interpretations.get("curated_named_markers", [])
    if curated_markers:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    curated_markers,
                    [
                        ("variant", "Variant"),
                        ("common_name", "Common name"),
                        ("observed_in_run", "Observed in run"),
                        ("observed_variants", "Observed labels"),
                        ("region_class", "Curated class"),
                        ("clinical_significance", "Clinical significance"),
                        ("summary", "Interpretation"),
                        ("associated_conditions", "Associated conditions"),
                        ("functional_effects", "Functional effects"),
                        ("research_context", "Research context"),
                    ],
                ),
                "Curated Marker Catalog",
                rows=None,
            )
        )

    variant_effect_overview = variant_interpretations.get("variant_effect_overview", [])
    if variant_effect_overview:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    [{"item": item} for item in variant_effect_overview],
                    [("item", "Variant interpretation overview")],
                ),
                "How This Gene's Variants Are Usually Interpreted",
                rows=None,
            )
        )

    condition_research_overview = variant_interpretations.get("condition_research_overview", [])
    if condition_research_overview:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    [{"item": item} for item in condition_research_overview],
                    [("item", "Condition or research theme")],
                ),
                "Conditions and Research Themes",
                rows=None,
            )
        )

    region_recommendations = variant_interpretations.get("region_recommendations", [])
    if region_recommendations:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    region_recommendations,
                    [
                        ("title", "View"),
                        ("region", "Region"),
                        ("purpose", "Purpose"),
                    ],
                ),
                "Region Recommendations",
                rows=None,
            )
        )

    promoter_analysis = variant_interpretations.get("promoter_analysis", {})
    if promoter_analysis.get("found_variants"):
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    promoter_analysis["found_variants"],
                    [
                        ("display", "Observed variant"),
                        ("variant_label", "Variant label"),
                        ("change", "Change"),
                        ("position", "Position"),
                        ("linked_to", "Linked to"),
                        ("clinical_significance", "Clinical significance"),
                        ("summary", "Interpretation"),
                    ],
                ),
                promoter_analysis.get("label", "Promoter Analysis") + " - Observed Variants",
                rows=None,
            )
        )
    if promoter_analysis.get("known_variants"):
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    promoter_analysis["known_variants"],
                    [
                        ("variant", "Variant"),
                        ("common_name", "Common name"),
                        ("position", "Position"),
                        ("clinical_significance", "Clinical significance"),
                        ("summary", "Interpretation"),
                        ("associated_conditions", "Associated conditions"),
                    ],
                ),
                promoter_analysis.get("label", "Promoter Analysis") + " - Curated Variants",
                rows=None,
            )
        )

    gene_analysis = variant_interpretations.get("gene_analysis", {})
    if gene_analysis.get("found_variants"):
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    gene_analysis["found_variants"],
                    [
                        ("display", "Observed variant"),
                        ("variant_label", "Variant label"),
                        ("change", "Change"),
                        ("position", "Position"),
                        ("linked_to", "Linked to"),
                        ("clinical_significance", "Clinical significance"),
                        ("summary", "Interpretation"),
                    ],
                ),
                gene_analysis.get("label", "Gene Analysis") + " - Observed Variants",
                rows=None,
            )
        )
    if gene_analysis.get("known_variants"):
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    gene_analysis["known_variants"],
                    [
                        ("variant", "Variant"),
                        ("common_name", "Common name"),
                        ("position", "Position"),
                        ("clinical_significance", "Clinical significance"),
                        ("summary", "Interpretation"),
                        ("associated_conditions", "Associated conditions"),
                    ],
                ),
                gene_analysis.get("label", "Gene Analysis") + " - Curated Variants",
                rows=None,
            )
        )

    if population_insights and population_insights.get("variant_population_records"):
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    population_insights["variant_population_records"],
                    [
                        ("display_name", "Variant"),
                        ("observed_in_run", "Observed in run"),
                        ("observed_variants", "Observed labels"),
                        ("effect_allele", "Effect allele"),
                        ("effect_summary", "Population summary"),
                        ("associated_conditions", "Associated conditions"),
                        ("functional_effects", "Functional effects"),
                        ("population_extremes", "Population extremes"),
                    ],
                ),
                "Population Reference Data",
                rows=None,
            )
        )

    if population_insights and population_insights.get("sources"):
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    population_insights["sources"],
                    [
                        ("label", "Source"),
                        ("url", "URL"),
                    ],
                ),
                "Population Sources",
                rows=None,
            )
        )

    if population_insights and population_insights.get("gene_population_patterns"):
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    population_insights["gene_population_patterns"],
                    [
                        ("variant", "Variant or locus"),
                        ("location_group", "Location group"),
                        ("summary", "Summary"),
                    ],
                ),
                "Gene-wide Population Patterns",
                rows=None,
            )
        )

    return _render_report_paragraphs(
        "Variant Interpretation",
        [
            variant_interpretations.get("sample_highlights", {}).get("summary", ""),
            variant_interpretations.get("summary", ""),
            variant_interpretations.get("gene_summary", ""),
            variant_interpretations.get("clinical_context", ""),
            variant_interpretations.get("curated_named_markers_summary", ""),
            population_insights.get("summary", "") if population_insights else "",
            population_insights.get("gene_population_patterns_intro", "") if population_insights else "",
        ],
        extra_markup="".join(nested_sections),
    )


def _render_methylation_interpretation_report(methylation_insights: dict[str, Any]) -> str:
    """Render the richer methylation interpretation block for the exported HTML report."""
    if not methylation_insights:
        return ""

    nested_sections: list[str] = []
    summary_metric_rows = methylation_insights.get("summary_metric_rows", [])
    if summary_metric_rows:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    summary_metric_rows,
                    [
                        ("metric", "Metric"),
                        ("mean_beta_display", "Mean beta"),
                        ("numeric_values", "Numeric values"),
                        ("summary", "Summary"),
                    ],
                ),
                "Methylation Summary Metrics",
                rows=None,
            )
        )

    whitelist_probe_statuses = methylation_insights.get("whitelist_probe_statuses", [])
    if whitelist_probe_statuses:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    whitelist_probe_statuses,
                    [
                        ("probe_id", "Whitelist probe"),
                        ("observed_in_run", "Observed in run"),
                    ],
                ),
                "Whitelist Probe Status",
                rows=None,
            )
        )

    methylation_effects = methylation_insights.get("methylation_effects", [])
    if methylation_effects:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    [{"item": item} for item in methylation_effects],
                    [("item", "Likely biological effect")],
                ),
                "Likely Biological Effects",
                rows=None,
            )
        )

    methylation_condition_research = methylation_insights.get("methylation_condition_research", [])
    if methylation_condition_research:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    [{"item": item} for item in methylation_condition_research],
                    [("item", "Condition or research setting")],
                ),
                "Methylation Research Context",
                rows=None,
            )
        )

    evidence = methylation_insights.get("evidence", [])
    if evidence:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    evidence,
                    [
                        ("label", "Source"),
                        ("url", "URL"),
                    ],
                ),
                "Methylation Evidence",
                rows=None,
            )
        )

    whitelist_probe_reference_rows = methylation_insights.get("whitelist_probe_reference_rows", [])
    if whitelist_probe_reference_rows:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    whitelist_probe_reference_rows,
                    [
                        ("probe_id", "Whitelist probe"),
                        ("observed_in_run", "Observed in run"),
                        ("beta", "Observed beta"),
                        ("probe_locus", "Probe locus"),
                        ("linked_variants", "Linked variants"),
                        ("nearby_manifest_variants", "Nearby manifest variants"),
                        ("papers", "Bundled papers"),
                    ],
                ),
                "Whitelist Probe Reference Map",
                rows=None,
            )
        )

    probe_preview = methylation_insights.get("probe_preview")
    if isinstance(probe_preview, pd.DataFrame) and not probe_preview.empty:
        nested_sections.append(
            _render_section_table(
                _prepare_methylation_table_for_output(probe_preview),
                "Observed Whitelist Probe Rows",
                rows=None,
            )
        )

    return _render_report_paragraphs(
        "Methylation Interpretation",
        [
            methylation_insights.get("summary", ""),
            methylation_insights.get("clinical_context", ""),
            methylation_insights.get("whitelist_explanation", ""),
            methylation_insights.get("whitelist_literature_context", ""),
            methylation_insights.get("gene_name_match_rule", ""),
            methylation_insights.get("whitelist_probe_reference_summary", ""),
        ],
        extra_markup="".join(nested_sections),
    )


def _render_predictive_theses_report(predictive_theses: dict[str, Any]) -> str:
    """Render the predictive-thesis panels in the exported HTML report."""
    if not predictive_theses:
        return ""

    nested_sections: list[str] = []
    variant_prediction_rows = predictive_theses.get("variant_prediction_rows", [])
    if variant_prediction_rows:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    variant_prediction_rows,
                    [
                        ("observed_signal", "Observed signal"),
                        ("source", "Source"),
                        ("prediction", "Prediction"),
                        ("research_focus", "Research focus"),
                    ],
                ),
                "Variant Prediction",
                rows=None,
            )
        )

    methylation_prediction_rows = predictive_theses.get("methylation_prediction_rows", [])
    if methylation_prediction_rows:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    methylation_prediction_rows,
                    [
                        ("metric_label", "Metric"),
                        ("mean_beta_display", "Mean beta"),
                        ("probe_count", "Numeric values"),
                        ("band_display", "Band"),
                        ("matched_case_label", "Matched case"),
                        ("prediction", "Prediction"),
                        ("research_focus", "Research focus"),
                    ],
                ),
                "Methylation Prediction",
                rows=None,
            )
        )

    matched_cases = predictive_theses.get("matched_cases", [])
    if matched_cases:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    matched_cases,
                    [
                        ("case_label", "Case"),
                        ("trigger", "Trigger"),
                        ("source", "Source"),
                        ("mean_beta_display", "Observed value"),
                        ("band", "Band"),
                        ("prediction", "Prediction"),
                        ("research_focus", "Research focus"),
                    ],
                ),
                "Synthesis",
                rows=None,
            )
        )

    seeded_markers = predictive_theses.get("seeded_markers", [])
    if seeded_markers:
        nested_sections.append(
            _render_section_table(
                _report_df_from_rows(
                    [{"marker": marker} for marker in seeded_markers],
                    [("marker", "Seeded marker")],
                ),
                "Seeded Predictive Markers",
                rows=None,
            )
        )

    return _render_report_paragraphs(
        "Predictive Theses",
        [
            predictive_theses.get("summary", ""),
            predictive_theses.get("variant_summary", ""),
            predictive_theses.get("matching_rule", ""),
            predictive_theses.get("disclaimer", ""),
        ],
        extra_markup="".join(nested_sections),
    )


def generate_report(
    variants: pd.DataFrame,
    methylation: pd.DataFrame,
    popstats: Any,
    output_path: str,
    *,
    gene_name: str = DEFAULT_GENE_NAME,
    region: str,
    methylation_output_path: Path | None = None,
    variant_interpretations: dict[str, Any] | None = None,
    methylation_insights: dict[str, Any] | None = None,
    population_insights: dict[str, Any] | None = None,
    predictive_theses: dict[str, Any] | None = None,
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
    predictive_theses : dict[str, Any] | None, optional
        Predictive thesis payload rendered into the report when available.

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

        variant_interpretation_section = _render_variant_interpretation_report(
            variant_interpretations or {},
            population_insights or {},
        )
        methylation_interpretation_section = _render_methylation_interpretation_report(
            methylation_insights or {},
        )
        predictive_theses_section = _render_predictive_theses_report(
            predictive_theses or {},
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
      width: min(98vw, 1800px);
      max-width: none;
      margin: 0 auto;
      padding: 38px 12px 64px;
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
      max-width: 78rem;
      overflow-wrap: anywhere;
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
      width: 100%;
      min-width: 0;
      overflow: hidden;
    }}
    h2 {{
      margin-top: 0;
      font-size: 1.3rem;
      letter-spacing: -0.02em;
    }}
    p, li, td, th, strong {{
      overflow-wrap: anywhere;
      word-break: break-word;
    }}
    .report-table-shell {{
      width: 100%;
      max-width: 100%;
      overflow-x: auto;
      border-radius: 18px;
      border: 1px solid var(--line);
      background: rgba(255, 255, 255, 0.54);
    }}
    .data-table {{
      width: 100%;
      min-width: 760px;
      border-collapse: collapse;
      font-size: 0.92rem;
      table-layout: fixed;
    }}
    .data-table th,
    .data-table td {{
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
      overflow-wrap: anywhere;
      word-break: break-word;
      white-space: normal;
    }}
    .data-table th {{
      background: rgba(15, 118, 110, 0.08);
    }}
    pre {{
      padding: 16px;
      overflow-x: auto;
      white-space: pre-wrap;
      overflow-wrap: anywhere;
      border-radius: 16px;
      background: #172023;
      color: #f4f0e8;
    }}
    @media (max-width: 760px) {{
      main {{
        width: min(100vw, 100%);
        padding: 18px 8px 42px;
      }}
      .hero,
      section {{
        padding: 18px;
      }}
      .data-table {{
        min-width: 680px;
      }}
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
    {_render_section_table(_prepare_variant_table_for_output(variants), "Genetic Variant Results", rows=None)}
    {variant_interpretation_section}
    {predictive_theses_section}
    {methylation_interpretation_section}
    {_render_section_table(_prepare_methylation_table_for_output(methylation), "Methylation Raw Results", rows=None)}
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
            "predictive_theses": predictive_theses or {},
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
                {
                    "metric": "predictive_thesis_matched_cases",
                    "value": (predictive_theses or {}).get("matched_case_count", 0),
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
    overwrite_general_database: bool = False,
    general_database_path: str | Path = GENERAL_ANALYSIS_DATABASE_PATH,
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
    overwrite_general_database : bool, optional
        When true, replace an existing central-database row for the analyzed
        gene. When false, an existing row is left unchanged.
    general_database_path : str | Path, optional
        Destination CSV for the central one-row-per-observed-variant database.

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
        variants = annotate_known_variant_ids(variants, knowledge_base)
        variant_interpretations = build_variant_interpretations(variants, knowledge_base, region=region)
        matched_variant_ids = {
            str(record.get("variant", "")).strip()
            for record in variant_interpretations.get("matched_records", [])
            if str(record.get("variant", "")).strip()
        }
        methylation_insights = build_methylation_insights(
            methylation,
            knowledge_base,
            matched_variant_ids=matched_variant_ids,
        )
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

    synthesis_database = load_gene_synthesis_database(normalized_gene_name)
    predictive_theses = build_predictive_theses(
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        knowledge_base=knowledge_base,
        synthesis_database=synthesis_database,
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

    general_database_result = update_general_analysis_database(
        gene_name=normalized_gene_name,
        variants=variants,
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        overwrite=overwrite_general_database,
        database_path=general_database_path,
    )

    final_report_path = generate_report(
        variants,
        methylation,
        popstats,
        str(report_path),
        gene_name=normalized_gene_name,
        region=region,
        methylation_output_path=methylation_output_path,
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        population_insights=population_insights,
        predictive_theses=predictive_theses,
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
        predictive_theses=predictive_theses,
        general_database_path=Path(general_database_result["path"]),
        general_database_status=str(general_database_result["message"]),
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
    print(f"{result.general_database_status} ({result.general_database_path})")
    return 0


if __name__ == "__main__":
    print(__doc__)
    raise SystemExit(main())
