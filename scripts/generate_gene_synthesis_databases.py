#!/usr/bin/env python3
"""Generate gene-level predictive synthesis databases from the bundled interpretation JSON files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENE_DATA_DIR = PROJECT_ROOT / "src" / "gene_data"
VERSION = "2026-04-18"
TARGET_GENES = [
    "HERC2",
    "DRD4",
    "IGF1R",
    "FOXO3",
    "MTOR",
    "RPS6",
    "SIK3",
    "FLCN",
    "SIRT6",
    "PRKAA1",
    "NAMPT",
    "CDKN2A",
    "TERT",
]

METHYLATION_SOURCES = [
    {
        "key": "whitelist",
        "label": "Whitelist mean beta",
        "description": "Uses the curated whitelist probes stored in the interpretation database for the current gene.",
    },
    {
        "key": "gene_name_related",
        "label": "Gene-name-related mean beta",
        "description": "Uses rows whose gene annotation explicitly names the current gene.",
    },
    {
        "key": "all_numeric",
        "label": "All numeric-row mean beta",
        "description": "Uses every numeric beta value that survived preprocessing for the current sample.",
    },
]

BAND_CONTEXT = {
    "high": "This pattern leans toward stronger regulatory restraint or reduced local accessibility around the locus.",
    "medium": "This pattern suggests an intermediate or mixed regulatory state rather than a strongly polarized epigenetic signal.",
    "low": "This pattern leans toward a more permissive or less methylated local regulatory state.",
}


def _clean_text(value: Any) -> str:
    """Normalize arbitrary values into compact single-line strings."""
    return " ".join(str(value or "").split())


def _first_nonempty(*values: Any) -> str:
    """Return the first cleaned non-empty value from a list of candidates."""
    for value in values:
        cleaned = _clean_text(value)
        if cleaned:
            return cleaned
    return ""


def _dedupe_text_items(values: list[str]) -> list[str]:
    """Return text items without duplicates while preserving order."""
    deduped: list[str] = []
    seen: set[str] = set()
    for value in values:
        cleaned = _clean_text(value)
        if not cleaned:
            continue
        key = cleaned.casefold()
        if key in seen:
            continue
        seen.add(key)
        deduped.append(cleaned)
    return deduped


def _candidate_interpretation_paths(gene_name: str) -> list[Path]:
    """Return likely interpretation database paths for one gene."""
    return [
        GENE_DATA_DIR / f"{gene_name.lower()}_interpretation_db.json",
        GENE_DATA_DIR / f"{gene_name}_interpretation_db.json",
        GENE_DATA_DIR / f"{gene_name.upper()}_interpretation_db.json",
    ]


def _load_interpretation_database(gene_name: str) -> dict[str, Any]:
    """Load the bundled interpretation database for a gene."""
    for path in _candidate_interpretation_paths(gene_name):
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
    raise FileNotFoundError(f"No interpretation database found for {gene_name}")


def _build_seeded_markers(knowledge_base: dict[str, Any]) -> list[str]:
    """Return a short display list of seeded markers for the synthesis UI."""
    markers: list[str] = []
    for record in knowledge_base.get("variant_records", []):
        display_name = _first_nonempty(record.get("display_name"), record.get("variant"))
        common_name = _clean_text(record.get("common_name"))
        if display_name and common_name:
            markers.append(f"{display_name} ({common_name})")
        elif display_name:
            markers.append(display_name)
    return _dedupe_text_items(markers)[:6]


def _collect_research_focus(knowledge_base: dict[str, Any]) -> list[str]:
    """Collect a concise research-focus list from the interpretation bundle."""
    gene_context = knowledge_base.get("gene_context", {})
    focus_items: list[str] = []
    focus_items.extend(gene_context.get("condition_research_overview", []))
    focus_items.extend(gene_context.get("methylation_condition_research", []))
    for record in knowledge_base.get("variant_records", []):
        focus_items.extend(record.get("associated_conditions", []))
    return _dedupe_text_items(focus_items)[:4]


def _build_base_case(
    *,
    gene_name: str,
    clinical_context: str,
    variant_focus: str,
    research_focus: list[str],
) -> dict[str, Any]:
    """Build the base variant-only synthesis case."""
    prediction = _clean_text(
        f"An observed {gene_name} variant places this sample in the curated {gene_name} research context. "
        f"{clinical_context} {variant_focus}"
    )
    return {
        "case_id": "gene_variant_found",
        "label": "Gene variant found",
        "requires_variant": True,
        "methylation_source": None,
        "methylation_band": None,
        "prediction": prediction,
        "rationale": (
            "This is the base thesis that activates as soon as any promoter or gene-body variant is visible in the current sample."
        ),
        "research_focus": research_focus[:3],
    }


def _build_combined_case(
    *,
    gene_name: str,
    source_key: str,
    source_label: str,
    source_description: str,
    band: str,
    clinical_context: str,
    variant_focus: str,
    methylation_context: str,
    research_focus: list[str],
) -> dict[str, Any]:
    """Build one of the nine variant-plus-methylation synthesis cases."""
    prediction = _clean_text(
        f"When a {gene_name} variant is paired with {band} methylation in the {source_label.lower()}, "
        f"the sample best fits a combined regulatory-context thesis. {BAND_CONTEXT[band]} "
        f"{methylation_context} {clinical_context} {variant_focus}"
    )
    return {
        "case_id": f"gene_variant_found__{source_key}__{band}",
        "label": f"Gene variant found + {band} {source_label.lower()}",
        "requires_variant": True,
        "methylation_source": source_key,
        "methylation_band": band,
        "prediction": prediction,
        "rationale": source_description,
        "research_focus": research_focus[:3],
    }


def build_synthesis_database(knowledge_base: dict[str, Any]) -> dict[str, Any]:
    """Create the 10-case predictive synthesis matrix for one gene."""
    gene_context = knowledge_base.get("gene_context", {})
    gene_name = _first_nonempty(gene_context.get("gene_name"), "UNKNOWN")
    clinical_context = _first_nonempty(
        gene_context.get("clinical_context"),
        gene_context.get("gene_summary"),
        f"The bundled {gene_name} database is intended for research context.",
    )
    variant_focus = _first_nonempty(
        *(gene_context.get("variant_effect_overview", []) or []),
        "Observed variants should be treated as locus-aware research context unless stronger external evidence exists.",
    )
    methylation_context = _first_nonempty(
        gene_context.get("methylation_interpretation"),
        *(gene_context.get("methylation_effects", []) or []),
        f"{gene_name} methylation is best interpreted as regulatory context rather than as a standalone biomarker.",
    )
    research_focus = _collect_research_focus(knowledge_base)

    cases = [
        _build_base_case(
            gene_name=gene_name,
            clinical_context=clinical_context,
            variant_focus=variant_focus,
            research_focus=research_focus,
        )
    ]

    for source in METHYLATION_SOURCES:
        for band in ("high", "medium", "low"):
            cases.append(
                _build_combined_case(
                    gene_name=gene_name,
                    source_key=source["key"],
                    source_label=source["label"],
                    source_description=source["description"],
                    band=band,
                    clinical_context=clinical_context,
                    variant_focus=variant_focus,
                    methylation_context=methylation_context,
                    research_focus=research_focus,
                )
            )

    return {
        "database_name": f"NophiGene {gene_name} Predictive Synthesis Database",
        "version": VERSION,
        "gene_name": gene_name,
        "source_interpretation_database": knowledge_base.get(
            "database_name",
            f"NophiGene {gene_name} Interpretation Database",
        ),
        "matching_rule": (
            "One base case matches when a promoter or gene-body variant is visible in the current sample. "
            "Three additional case families match when the whitelist mean beta, the gene-name-related mean beta, "
            "or the all-numeric mean beta resolves to low, medium, or high methylation."
        ),
        "disclaimer": (
            "Predictive theses in this database are literature-guided research summaries derived from the bundled gene interpretation bundle. "
            "They are designed for exploratory synthesis in the UI and should not be treated as diagnostic or therapeutic claims."
        ),
        "seeded_markers": _build_seeded_markers(knowledge_base),
        "case_count": len(cases),
        "methylation_sources": METHYLATION_SOURCES,
        "cases": cases,
    }


def main() -> int:
    """Generate one synthesis JSON file per target gene."""
    GENE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    for gene_name in TARGET_GENES:
        knowledge_base = _load_interpretation_database(gene_name)
        synthesis_database = build_synthesis_database(knowledge_base)
        output_path = GENE_DATA_DIR / f"{gene_name.lower()}_synthesis.json"
        output_path.write_text(
            json.dumps(synthesis_database, indent=2, ensure_ascii=True) + "\n",
            encoding="utf-8",
        )
        print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
