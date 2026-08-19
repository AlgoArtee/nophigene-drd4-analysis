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

# Sources in this set can contribute condition or pathology context. Membership
# here does not make a record authoritative; that remains the job of
# ``medical_gate`` below.
MEDICAL_CONTEXT_SOURCE_KEYS = frozenset(
    {
        "clingen", "clinvar", "cpic", "dpwg", "ema", "fda", "fda_pgx",
        "pharmgkb", "pharmvar", "civic", "medgen", "panelapp", "opentargets",
        "clinicaltrials", "omim", "oncokb", "dgidb", "drugbank",
        "professional_guideline",
    }
)

PREPRINT_SOURCES = frozenset({"biorxiv", "medrxiv", "research_square"})

LITERATURE_SOURCE_KEYS = frozenset(
    {
        "pubmed", "pmc", "europe_pmc", "openalex", "crossref", "semantic_scholar",
        "biorxiv", "medrxiv", "research_square", "local_pdf_articles", "local_pdf",
    }
)

LITERATURE_PRIORITY_LABELS = {
    1: "Clinical relevance + direct experimental/functional evidence",
    2: "Clinical relevance",
    3: "Direct experimental/functional evidence",
    4: "Review, meta-analysis, or replicated evidence",
    5: "Human association or observational finding",
    6: "Other gene-relevant evidence",
    7: "Preprint",
}

_CLINICAL_PATTERN = re.compile(
    r"\b(?:clinical|patient|case[- ]control|diagnos\w*|disease|disorder|syndrome|risk|symptom|"
    r"treatment|response|therap\w*|drug|pharmaco\w*|schizophrenia|adhd|addiction|heroin|autism|"
    r"depression|bipolar|cancer|tumou?rs?|diabetes|cardiovascular|cohort)\b",
    re.IGNORECASE,
)
_EXPERIMENTAL_PATTERN = re.compile(
    r"\b(?:experiment|functional|reporter|luciferase|assay|transcription|expression|methylation|"
    r"mqtl|binding|biochemical|mechanis\w*|engineered|cell line|in vitro|in vivo|animal model|"
    r"postmortem|post-mortem|brain tissue|knockdown|knockout|crispr|receptor|promoter activity)\b",
    re.IGNORECASE,
)
_REVIEW_PATTERN = re.compile(r"\b(?:systematic review|review|meta-analysis|meta analysis|pooled|replicat\w*)\b", re.IGNORECASE)
_ASSOCIATION_PATTERN = re.compile(
    r"\b(?:association|observational|case[- ]control|cohort|gwas|odds ratio|genotype frequenc\w*|correlat\w*)\b",
    re.IGNORECASE,
)


def _clean(value: Any) -> str:
    return " ".join(str(value or "").split())


def _normalized_source(value: Any) -> str:
    source = re.sub(r"[^a-z0-9]+", "_", _clean(value).casefold()).strip("_")
    aliases = {
        "pubmed_central": "pmc",
        "europe_pmc": "europe_pmc",
        "semantic_scholar": "semantic_scholar",
        "local_pdf_article": "local_pdf_articles",
    }
    return aliases.get(source, source)


def canonical_medical_source(value: Any) -> str:
    """Normalize connector keys and human-readable authority aliases."""
    source = _normalized_source(value)
    aliases = {
        "clinical_genome_resource": "clingen",
        "clinical_genome": "clingen",
        "ncbi_clinvar": "clinvar",
        "ncbi_medgen": "medgen",
        "genomics_england_panelapp": "panelapp",
        "genomics_england_panel_app": "panelapp",
        "panel_app": "panelapp",
        "clinical_pharmacogenetics_implementation_consortium": "cpic",
        "dutch_pharmacogenetics_working_group": "dpwg",
        "european_medicines_agency": "ema",
        "fda_pharmacogenetic_associations": "fda_pgx",
        "fda_pharmacogenomic_biomarkers": "fda_pgx",
        "fda_pharmacogenomic_labels": "fda_pgx",
        "fda_pharmacogenomics": "fda_pgx",
        "food_and_drug_administration": "fda",
        "open_targets": "opentargets",
        "clinicaltrials_gov": "clinicaltrials",
    }
    return aliases.get(source, source)


def _publication_identifiers(record: dict[str, Any]) -> dict[str, str]:
    material = " ".join(
        _clean(record.get(field))
        for field in ("pmid", "pmcid", "doi", "source_id", "title", "paper", "label", "url", "citation")
    )
    identifiers: dict[str, str] = {}
    pmid = re.sub(r"\D", "", _clean(record.get("pmid")))
    pmcid = _clean(record.get("pmcid")).upper()
    doi = _clean(record.get("doi"))
    source_key = _normalized_source(record.get("source_key") or record.get("source"))
    source_id = _clean(record.get("source_id"))
    if not pmid:
        match = re.search(r"(?:PMID\s*[:#]?\s*|pubmed\.ncbi\.nlm\.nih\.gov/)(\d{5,10})", material, re.IGNORECASE)
        if match:
            pmid = match.group(1)
        elif source_key == "pubmed" and source_id.isdigit():
            pmid = source_id
    if not pmcid:
        match = re.search(r"\bPMC\s*(\d{4,10})\b", material, re.IGNORECASE)
        if match:
            pmcid = f"PMC{match.group(1)}"
        elif source_key == "pmc" and source_id:
            pmcid = source_id.upper() if source_id.upper().startswith("PMC") else f"PMC{source_id}"
    if not doi:
        match = re.search(r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+", material, re.IGNORECASE)
        if match:
            doi = match.group(0).rstrip(".,;]")
    if pmid:
        identifiers["pmid"] = pmid
    if pmcid:
        normalized_pmcid = re.sub(r"\s+", "", pmcid.upper())
        identifiers["pmcid"] = normalized_pmcid if normalized_pmcid.startswith("PMC") else f"PMC{normalized_pmcid}"
    if doi:
        normalized_doi = doi.casefold().removeprefix("https://doi.org/").removeprefix("http://doi.org/")
        identifiers["doi"] = normalized_doi
    return identifiers


def normalize_literature_record(
    raw: dict[str, Any],
    *,
    gene: str = "",
    evidence_origin: str = "",
    parent_variant: str = "",
) -> dict[str, Any]:
    """Normalize publication-like source records without synthesizing scientific claims."""
    record = dict(raw)
    identifiers = _publication_identifiers(record)
    source_key = _normalized_source(record.get("source_key") or record.get("source"))
    url = _clean(record.get("url"))
    if not source_key:
        if "pubmed.ncbi.nlm.nih.gov" in url:
            source_key = "pubmed"
        elif "pmc.ncbi.nlm.nih.gov" in url:
            source_key = "pmc"
        elif "doi.org" in url:
            source_key = "crossref"
    title = _clean(record.get("title") or record.get("label") or record.get("paper"))
    paper = _clean(record.get("paper") or record.get("citation"))
    finding = _clean(record.get("finding") or record.get("assertion") or record.get("snippet"))
    normalized = {
        **record,
        **identifiers,
        "gene": _clean(record.get("gene") or record.get("gene_name") or gene).upper(),
        "title": title,
        "paper": paper,
        "finding": finding,
        "phenotype": _clean(record.get("phenotype") or record.get("condition")),
        "genotypes": _clean(record.get("genotypes") or record.get("genotype")),
        "variant": _clean(record.get("variant") or record.get("rsid") or parent_variant),
        "url": url,
        "source_key": source_key or "unknown",
        "evidence_origin": _clean(record.get("evidence_origin") or evidence_origin or "source_record"),
    }
    year = record.get("publication_year") or record.get("year")
    if not year:
        match = re.search(r"\b(?:19|20)\d{2}\b", " ".join((title, paper)))
        year = int(match.group(0)) if match else None
    try:
        normalized["publication_year"] = int(year) if year else None
    except (TypeError, ValueError):
        normalized["publication_year"] = None
    category = _clean(record.get("category")).casefold()
    normalized["is_publication_record"] = bool(
        identifiers
        or record.get("title")
        or record.get("paper")
        or finding
        or source_key in LITERATURE_SOURCE_KEYS
        or category == "literature"
    )
    identity = "|".join(
        filter(
            None,
            (
                normalized["gene"],
                source_key,
                normalized["evidence_origin"],
                identifiers.get("pmid"),
                identifiers.get("pmcid"),
                identifiers.get("doi"),
                title.casefold(),
                finding.casefold(),
                normalized["variant"].casefold(),
            ),
        )
    )
    normalized["record_id"] = _clean(record.get("record_id")) or hashlib.sha256(identity.encode("utf-8")).hexdigest()[:32]
    normalized["entity_type"] = _clean(record.get("entity_type")) or ("variant" if normalized["variant"] else "gene")
    normalized["entity_key"] = _clean(record.get("entity_key")) or normalized["variant"] or normalized["gene"]
    normalized["evidence_type"] = _clean(record.get("evidence_type") or record.get("publication_type")) or "literature"
    return normalized


def publication_key(record: dict[str, Any]) -> str:
    for field in ("pmid", "pmcid", "doi"):
        value = _clean(record.get(field)).casefold()
        if value:
            return f"{field}:{value}"
    title = re.sub(r"[^a-z0-9]+", " ", _clean(record.get("title")).casefold()).strip()
    return "title:" + hashlib.sha256(title.encode("utf-8")).hexdigest() if title else ""


def _publication_aliases(record: dict[str, Any]) -> set[str]:
    aliases = {
        f"{field}:{_clean(record.get(field)).casefold()}"
        for field in ("pmid", "pmcid", "doi")
        if _clean(record.get(field))
    }
    title = re.sub(r"[^a-z0-9]+", " ", _clean(record.get("title")).casefold()).strip()
    if title:
        aliases.add("title:" + hashlib.sha256(title.encode("utf-8")).hexdigest())
    return aliases


def _merge_publication(target: dict[str, Any], record: dict[str, Any]) -> None:
    target["sources"] = sorted(
        set(target.get("sources", [])) | {_clean(record.get("source_key"))} - {""}
    )
    target["evidence_origins"] = sorted(
        set(target.get("evidence_origins", [])) | {_clean(record.get("evidence_origin"))} - {""}
    )
    target["preprint"] = bool(target.get("preprint")) or bool(record.get("preprint")) or (
        _normalized_source(record.get("source_key")) in PREPRINT_SOURCES
    )
    for field in (
        "abstract", "summary", "doi", "pmid", "pmcid", "url", "paper", "publication_type",
        "publication_year", "phenotype", "genotypes", "variant", "gene",
    ):
        if not target.get(field) and record.get(field):
            target[field] = record[field]
    target["context_matched"] = bool(target.get("context_matched")) or bool(record.get("context_matched"))
    target["observed_entity_matched"] = bool(target.get("observed_entity_matched")) or bool(record.get("observed_entity_matched"))


def _group_publications(records: Iterable[dict[str, Any]]) -> list[tuple[dict[str, Any], list[dict[str, Any]]]]:
    groups: list[dict[str, Any] | None] = []
    alias_to_group: dict[str, int] = {}
    for raw in records:
        record = dict(raw)
        aliases = _publication_aliases(record)
        if not aliases:
            continue
        matching = sorted({alias_to_group[alias] for alias in aliases if alias in alias_to_group and groups[alias_to_group[alias]]})
        if not matching:
            index = len(groups)
            merged = dict(record)
            merged["sources"] = sorted({_clean(record.get("source_key"))} - {""})
            merged["evidence_origins"] = sorted({_clean(record.get("evidence_origin"))} - {""})
            merged["preprint"] = bool(record.get("preprint")) or _normalized_source(record.get("source_key")) in PREPRINT_SOURCES
            groups.append({"merged": merged, "members": [record], "aliases": set(aliases)})
        else:
            index = matching[0]
            group = groups[index]
            assert group is not None
            for duplicate_index in matching[1:]:
                duplicate = groups[duplicate_index]
                if duplicate is None:
                    continue
                _merge_publication(group["merged"], duplicate["merged"])
                group["members"].extend(duplicate["members"])
                group["aliases"].update(duplicate["aliases"])
                groups[duplicate_index] = None
            _merge_publication(group["merged"], record)
            group["members"].append(record)
            group["aliases"].update(aliases)
        group = groups[index]
        assert group is not None
        for alias in group["aliases"]:
            alias_to_group[alias] = index
    output: list[tuple[dict[str, Any], list[dict[str, Any]]]] = []
    for group in groups:
        if group is None:
            continue
        merged = group["merged"]
        merged["publication_key"] = publication_key(merged) or sorted(group["aliases"])[0]
        output.append((merged, group["members"]))
    return output


def classify_literature_finding(record: dict[str, Any]) -> dict[str, Any]:
    evidence_text = " ".join(
        _clean(record.get(field))
        for field in (
            "title", "paper", "finding", "summary", "abstract", "phenotype", "publication_type",
            "evidence_type", "interpretation_scope", "clinical_significance",
        )
    )
    preprint = bool(record.get("preprint")) or _normalized_source(record.get("source_key")) in PREPRINT_SOURCES
    clinical = bool(_CLINICAL_PATTERN.search(evidence_text))
    experimental = bool(_EXPERIMENTAL_PATTERN.search(evidence_text))
    review = bool(_REVIEW_PATTERN.search(evidence_text))
    association = bool(_ASSOCIATION_PATTERN.search(evidence_text))
    finding_available = bool(_clean(record.get("finding") or record.get("summary") or record.get("abstract")))
    tags: list[str] = []
    reasons: list[str] = []
    if clinical:
        tags.append("clinical_relevance")
        reasons.append("The reported phenotype or study context has clinical relevance.")
    if experimental:
        tags.append("experimental_functional")
        reasons.append("The record describes direct experimental, functional, molecular, or tissue evidence.")
    if review:
        tags.append("review_meta_analysis")
        reasons.append("The record reports review, meta-analysis, pooled, or replicated evidence.")
    if association:
        tags.append("human_association_observational")
        reasons.append("The record reports human association or observational evidence.")
    if preprint:
        tags.append("preprint")
        reasons.append("Preprints are ranked after reviewed publications.")
        tier = 7
    elif not finding_available:
        tier = 6
    elif clinical and experimental:
        tier = 1
    elif clinical:
        tier = 2
    elif experimental:
        tier = 3
    elif review:
        tier = 4
    elif association:
        tier = 5
    else:
        tier = 6
    if not finding_available:
        tags.append("citation_only")
        reasons.append("Only citation metadata is available; no finding was synthesized.")
    return {
        "priority_tier": tier,
        "priority_label": LITERATURE_PRIORITY_LABELS[tier],
        "evidence_tags": tags,
        "rank_reasons": reasons,
        "finding_status": "finding_available" if finding_available else "citation_metadata_only",
        "preprint": preprint,
    }


def deduplicate_literature(
    records: Iterable[dict[str, Any]], *, candidate_limit: int | None = None
) -> list[dict[str, Any]]:
    ranked = sorted(
        (merged for merged, _members in _group_publications(records)),
        key=lambda item: (
            not bool(item.get("context_matched")),
            not bool(item.get("observed_entity_matched")),
            bool(item.get("preprint")),
            -int(item.get("publication_year") or 0),
            str(item.get("title") or ""),
        ),
    )
    return ranked[:candidate_limit] if candidate_limit is not None else ranked


def build_literature_catalog(records: Iterable[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return publication-level citations and distinct, completely ranked findings."""
    publications: list[dict[str, Any]] = []
    findings: list[dict[str, Any]] = []
    for publication, members in _group_publications(records):
        rich_members = [
            member for member in members
            if _clean(member.get("finding") or member.get("summary") or member.get("abstract"))
        ]
        finding_members = rich_members or [publication]
        publication_findings: list[dict[str, Any]] = []
        seen_findings: set[str] = set()
        for member in finding_members:
            finding_text = _clean(member.get("finding") or member.get("summary") or member.get("abstract"))
            entity_context = _clean(member.get("variant") or member.get("entity_key"))
            distinct_key = hashlib.sha256(
                "|".join((publication["publication_key"], finding_text.casefold(), entity_context.casefold())).encode("utf-8")
            ).hexdigest()
            if distinct_key in seen_findings:
                continue
            seen_findings.add(distinct_key)
            finding = dict(member)
            for field in ("pmid", "pmcid", "doi", "url", "title", "publication_year", "gene"):
                if not finding.get(field) and publication.get(field):
                    finding[field] = publication[field]
            finding["publication_key"] = publication["publication_key"]
            finding["sources"] = list(publication.get("sources", []))
            finding["evidence_origins"] = list(publication.get("evidence_origins", []))
            finding["preprint"] = bool(publication.get("preprint"))
            finding.update(classify_literature_finding(finding))
            finding["finding_id"] = distinct_key[:32]
            finding["record_id"] = finding["finding_id"]
            publication_findings.append(finding)
        publication["finding_count"] = len(publication_findings)
        aggregate_tags = sorted({tag for finding in publication_findings for tag in finding["evidence_tags"]})
        publication["evidence_tags"] = aggregate_tags
        publication["priority_tier"] = min(finding["priority_tier"] for finding in publication_findings)
        publication["priority_label"] = LITERATURE_PRIORITY_LABELS[publication["priority_tier"]]
        publications.append(publication)
        findings.extend(publication_findings)
    findings.sort(
        key=lambda item: (
            int(item.get("priority_tier") or 99),
            item.get("finding_status") == "citation_metadata_only",
            not bool(item.get("observed_entity_matched")),
            not bool(item.get("context_matched")),
            -int(item.get("publication_year") or 0),
            _clean(item.get("paper") or item.get("title")).casefold(),
            _clean(item.get("variant")).casefold(),
            _clean(item.get("finding")).casefold(),
        )
    )
    for rank, finding in enumerate(findings, start=1):
        finding["rank"] = rank
    publications.sort(
        key=lambda item: (
            int(item.get("priority_tier") or 99),
            -int(item.get("publication_year") or 0),
            _clean(item.get("title")).casefold(),
        )
    )
    return publications, findings


def medical_gate(record: dict[str, Any]) -> dict[str, Any]:
    source = canonical_medical_source(record.get("source_key") or record.get("source"))
    authority = canonical_medical_source(record.get("authority") or source)
    evidence_type = _clean(record.get("evidence_type")).casefold().replace(" ", "_")
    review_status = _clean(record.get("review_status")).casefold()
    preprint = bool(record.get("preprint")) or source in PREPRINT_SOURCES
    blockers: list[str] = []
    if preprint:
        blockers.append("preprint_not_medical_evidence")
    if evidence_type in {"gwas_association", "trial_recruitment", "case_report", "adverse_event_signal"}:
        blockers.append(f"{evidence_type}_belongs_in_literature")
    expert_panel_clinvar = source == "clinvar" and any(
        token in review_status for token in ("expert panel", "practice guideline")
    )
    established = authority in ESTABLISHED_MEDICAL_AUTHORITIES or expert_panel_clinvar
    if not established:
        blockers.append("no_authoritative_clinical_assertion")

    category = _clean(record.get("category") or record.get("record_type")).casefold().replace(" ", "_")
    metadata_only = (
        category in {"source_metadata", "metadata", "linkout", "linkout_only"}
        or evidence_type in {"source_metadata", "metadata", "linkout", "linkout_only"}
        or bool(record.get("metadata_only"))
    )
    if metadata_only:
        blockers.append("metadata_only_or_linkout_record")

    effective_date = _clean(
        record.get("effective_date")
        or record.get("source_release")
        or record.get("date")
        or record.get("published_date")
        or record.get("last_updated")
        or record.get("modification_date")
        or record.get("last_evaluated")
        or record.get("updatedate")
        or record.get("date_release")
    )
    if not effective_date:
        blockers.append("missing_effective_date_or_release")

    substantive_fields = (
        "assertion", "classification", "actionability", "actionability_statement",
        "clinical_significance", "recommendation", "guideline", "label_statement",
        "haploinsufficiency", "triplosensitivity", "intervention", "overall",
    )
    substantive = any(_clean(record.get(field)) for field in substantive_fields)
    substantive = substantive or evidence_type in {
        "professional_guideline", "practice_guideline", "drug_label", "expert_panel_assertion",
    }
    if not substantive:
        blockers.append("missing_substantive_medical_assertion")

    blockers = list(dict.fromkeys(blockers))
    return {
        "eligible": not blockers,
        "blockers": blockers,
        "destination": "medical" if not blockers else "literature",
        "canonical_source_key": source,
        "effective_date": effective_date or None,
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
