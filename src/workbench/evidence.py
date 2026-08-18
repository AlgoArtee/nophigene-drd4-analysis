"""Evidence normalization, literature deduplication, and medical gating."""

from __future__ import annotations

import hashlib
import re
from collections import defaultdict
from typing import Any, Iterable

CORE_SOURCE_KEYS = frozenset(
    {
        "ensembl", "ncbi_gene", "uniprot", "interpro", "quickgo", "human_protein_atlas",
        "gnomad", "clinvar", "gwas_catalog", "gtex", "opentargets", "encode", "screen",
        "jaspar", "unibind", "string", "reactome", "intact", "ensembl_interactions",
        "pubmed", "pmc", "clingen", "cpic", "pharmgkb", "ema", "fda_pgx", "civic",
        "clinicaltrials", "rcsb_pdb", "alphafold_db",
    }
)

CORE_SOURCE_METADATA: dict[str, dict[str, str]] = {
    "ensembl": {"name": "Ensembl VEP / MANE", "lane": "identity"},
    "ncbi_gene": {"name": "NCBI Gene", "lane": "identity"},
    "uniprot": {"name": "UniProt", "lane": "identity"},
    "interpro": {"name": "InterPro", "lane": "identity"},
    "quickgo": {"name": "QuickGO", "lane": "identity"},
    "human_protein_atlas": {"name": "Human Protein Atlas", "lane": "identity"},
    "gnomad": {"name": "gnomAD", "lane": "variant_population"},
    "clinvar": {"name": "ClinVar", "lane": "variant_population"},
    "gwas_catalog": {"name": "GWAS Catalog", "lane": "variant_population"},
    "gtex": {"name": "GTEx", "lane": "regulatory"},
    "encode": {"name": "ENCODE", "lane": "regulatory"},
    "screen": {"name": "ENCODE SCREEN", "lane": "regulatory"},
    "jaspar": {"name": "JASPAR", "lane": "regulatory"},
    "unibind": {"name": "UniBind", "lane": "regulatory"},
    "string": {"name": "STRING", "lane": "interactions"},
    "reactome": {"name": "Reactome", "lane": "interactions"},
    "intact": {"name": "IntAct", "lane": "interactions"},
    "ensembl_interactions": {"name": "Ensembl Interactions", "lane": "interactions"},
    "clingen": {"name": "ClinGen", "lane": "medical"},
    "cpic": {"name": "CPIC", "lane": "medical"},
    "pharmgkb": {"name": "PharmGKB", "lane": "medical"},
    "ema": {"name": "EMA labels", "lane": "medical"},
    "fda_pgx": {"name": "FDA pharmacogenomic labels", "lane": "medical"},
    "civic": {"name": "CIViC", "lane": "medical"},
    "opentargets": {"name": "Open Targets", "lane": "disease"},
    "clinicaltrials": {"name": "ClinicalTrials.gov", "lane": "literature"},
    "pubmed": {"name": "PubMed", "lane": "literature"},
    "pmc": {"name": "PubMed Central", "lane": "literature"},
    "rcsb_pdb": {"name": "RCSB PDB", "lane": "structure"},
    "alphafold_db": {"name": "AlphaFold DB", "lane": "structure"},
}


def list_core_source_metadata() -> list[dict[str, str]]:
    """Return the stable typed Version 2 source catalog."""
    return [
        {"key": key, **CORE_SOURCE_METADATA.get(key, {"name": key, "lane": "other"})}
        for key in sorted(CORE_SOURCE_KEYS)
    ]

ESTABLISHED_MEDICAL_AUTHORITIES = frozenset(
    {"clingen", "cpic", "dpwg", "ema", "fda", "fda_pgx", "clinvar_expert_panel", "professional_guideline"}
)

PREPRINT_SOURCES = frozenset({"biorxiv", "medrxiv", "research_square"})


def _clean(value: Any) -> str:
    return " ".join(str(value or "").split())


def publication_key(record: dict[str, Any]) -> str:
    for field in ("pmid", "pmcid", "doi"):
        value = _clean(record.get(field)).casefold()
        if value:
            return f"{field}:{value}"
    title = re.sub(r"[^a-z0-9]+", " ", _clean(record.get("title")).casefold()).strip()
    return "title:" + hashlib.sha256(title.encode("utf-8")).hexdigest() if title else ""


def deduplicate_literature(records: Iterable[dict[str, Any]], *, candidate_limit: int = 200) -> list[dict[str, Any]]:
    deduplicated: dict[str, dict[str, Any]] = {}
    for raw in records:
        record = dict(raw)
        key = publication_key(record)
        if not key:
            continue
        if key not in deduplicated:
            record["publication_key"] = key
            record["sources"] = sorted({_clean(record.get("source_key"))} - {""})
            record["preprint"] = _clean(record.get("source_key")).casefold() in PREPRINT_SOURCES
            deduplicated[key] = record
        else:
            merged_sources = set(deduplicated[key].get("sources", []))
            merged_sources.add(_clean(record.get("source_key")))
            deduplicated[key]["sources"] = sorted(merged_sources - {""})
            for field in ("abstract", "summary", "doi", "pmid", "pmcid", "url"):
                if not deduplicated[key].get(field) and record.get(field):
                    deduplicated[key][field] = record[field]
    ranked = sorted(
        deduplicated.values(),
        key=lambda item: (
            not bool(item.get("context_matched")),
            not bool(item.get("observed_entity_matched")),
            bool(item.get("preprint")),
            -int(item.get("publication_year") or 0),
            str(item.get("title") or ""),
        ),
    )
    return ranked[:candidate_limit]


def medical_gate(record: dict[str, Any]) -> dict[str, Any]:
    source = _clean(record.get("source_key")).casefold()
    authority = _clean(record.get("authority") or source).casefold()
    evidence_type = _clean(record.get("evidence_type")).casefold()
    review_status = _clean(record.get("review_status")).casefold()
    preprint = bool(record.get("preprint")) or source in PREPRINT_SOURCES
    blockers: list[str] = []
    if preprint:
        blockers.append("preprint_not_medical_evidence")
    if evidence_type in {"gwas_association", "trial_recruitment", "case_report", "adverse_event_signal"}:
        blockers.append(f"{evidence_type}_belongs_in_literature")
    established = authority in ESTABLISHED_MEDICAL_AUTHORITIES or (
        source == "clinvar" and any(token in review_status for token in ("expert panel", "practice guideline"))
    )
    if not established:
        blockers.append("no_authoritative_clinical_assertion")
    if not _clean(record.get("effective_date") or record.get("source_release")):
        blockers.append("missing_effective_date_or_release")
    return {
        "eligible": not blockers,
        "blockers": blockers,
        "destination": "medical" if not blockers else "literature",
    }


def source_coverage(statuses: Iterable[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for status in statuses:
        state = _clean(status.get("status") or "not_assessed").casefold()
        if state in {"ok", "imported", "assessed", "empty"}:
            grouped["assessed"].append(dict(status))
        elif state in {"failed", "error", "timeout"}:
            grouped["failed"].append(dict(status))
        else:
            grouped["not_assessed"].append(dict(status))
    return {
        "assessed": grouped["assessed"],
        "failed": grouped["failed"],
        "not_assessed": grouped["not_assessed"],
        "assessed_absence": bool(grouped["assessed"]) and not any(
            int(item.get("record_count") or 0) > 0 for item in grouped["assessed"]
        ),
    }
