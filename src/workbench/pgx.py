"""Conservative pharmacogenomic diplotype matching for the Predictions lane."""

from __future__ import annotations

from typing import Any


def resolve_pgx_diplotype(definition: dict[str, Any], calls: list[dict[str, Any]]) -> dict[str, Any]:
    gene = str(definition.get("gene") or "").upper()
    candidates = [item for item in definition.get("candidates", []) if isinstance(item, dict)]
    call_map = {str(item.get("locus") or item.get("rsid") or ""): item for item in calls if isinstance(item, dict)}
    defining_loci = sorted(
        {
            str(locus)
            for candidate in candidates
            for locus in (candidate.get("required_genotypes") or {})
            if str(locus)
        }
    )
    blockers: list[str] = []
    for locus in defining_loci:
        call = call_map.get(locus)
        if not call:
            blockers.append(f"missing_defining_locus:{locus}")
        elif not call.get("qc_pass") or not call.get("covered", True):
            blockers.append(f"defining_locus_not_qc_covered:{locus}")
    if definition.get("cnv_dependent") and not definition.get("cnv_result"):
        blockers.append("cnv_result_required")
    matches: list[dict[str, Any]] = []
    if not blockers:
        for candidate in candidates:
            requirements = candidate.get("required_genotypes") or {}
            if candidate.get("phase_required") and not all(call_map[locus].get("phased") for locus in requirements):
                continue
            if all(str(call_map[locus].get("genotype") or "") == str(genotype) for locus, genotype in requirements.items()):
                matches.append(candidate)
    if len(matches) == 1:
        match = matches[0]
        return {
            "gene": gene,
            "status": "resolved",
            "diplotype": match.get("diplotype"),
            "publisher_phenotype": match.get("phenotype"),
            "phenotype_status": "publisher_mapping_available" if match.get("phenotype") else "not_assessed",
            "defining_loci": defining_loci,
            "blockers": [],
            "candidate_count": 1,
            "policy": "all defining loci QC-covered and exactly one solution",
        }
    return {
        "gene": gene,
        "status": "unresolved",
        "diplotype": None,
        "publisher_phenotype": None,
        "phenotype_status": "not_assessed",
        "defining_loci": defining_loci,
        "blockers": blockers or (["multiple_solutions"] if len(matches) > 1 else ["no_exact_solution"]),
        "candidate_diplotypes": [item.get("diplotype") for item in matches],
        "candidate_count": len(matches),
        "policy": "ambiguous, unphased, CNV-dependent, or incompletely covered results remain unresolved",
    }
