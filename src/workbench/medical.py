"""Layered queried-gene medical, condition, cohort, and pathology context."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Iterable

from .evidence import (
    MEDICAL_CONTEXT_SOURCE_KEYS,
    canonical_medical_source,
    medical_gate,
)


def _clean(value: Any) -> str:
    return " ".join(str(value or "").split())


def _items(value: Any) -> list[Any]:
    if value in (None, ""):
        return []
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def _text_items(value: Any) -> list[str]:
    output: list[str] = []
    for item in _items(value):
        if isinstance(item, dict):
            text = _clean(
                item.get("name")
                or item.get("label")
                or item.get("title")
                or item.get("disease")
                or item.get("condition")
                or item.get("id")
            )
        else:
            text = _clean(item)
        if text and text not in output:
            output.append(text)
    return output


def _normalized_text(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", " ", _clean(value).casefold()).strip()


def _safe_int(value: Any) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def _match_rank(item: dict[str, Any]) -> int:
    basis = set(item.get("match_basis") or [])
    if "observed_variant" in basis:
        return 0
    if "observed_methylation_probe" in basis:
        return 1
    return 2


def _tagged_context(
    *,
    evidence_class: str = "research_context",
    authority_status: str = "not_authoritative",
    observed_match: bool = False,
    match_basis: Iterable[str] = ("queried_gene",),
) -> dict[str, Any]:
    return {
        "evidence_class": evidence_class,
        "authority_status": authority_status,
        "research_context_only": authority_status != "authoritative",
        "observed_match": observed_match,
        "match_basis": list(dict.fromkeys(_clean(item) for item in match_basis if _clean(item))),
    }


def _observed_context(payload: dict[str, Any]) -> tuple[dict[str, list[str]], set[str]]:
    variant_interpretations = (
        payload.get("variant_interpretations")
        if isinstance(payload.get("variant_interpretations"), dict)
        else {}
    )
    matched_variants: dict[str, list[str]] = defaultdict(list)
    for record in variant_interpretations.get("matched_records", []):
        if not isinstance(record, dict):
            continue
        observed = _clean(
            record.get("observed_variant")
            or record.get("variant_label")
            or record.get("variant")
            or record.get("rsid")
        )
        for value in (
            record.get("variant"),
            record.get("rsid"),
            record.get("variant_label"),
        ):
            key = _normalized_text(value)
            if key and observed and observed not in matched_variants[key]:
                matched_variants[key].append(observed)
    methylation_insights = (
        payload.get("methylation_insights")
        if isinstance(payload.get("methylation_insights"), dict)
        else {}
    )
    observed_probes = {_clean(item) for item in methylation_insights.get("probe_ids", []) if _clean(item)}
    for status in methylation_insights.get("whitelist_probe_statuses", []):
        if isinstance(status, dict) and status.get("observed_in_run") and _clean(status.get("probe_id")):
            observed_probes.add(_clean(status.get("probe_id")))
    return dict(matched_variants), observed_probes


def _variant_match(
    variant: Any,
    variant_record: dict[str, Any] | None,
    matched_variants: dict[str, list[str]],
    observed_probes: set[str],
) -> dict[str, Any]:
    record = variant_record or {}
    keys = {
        _normalized_text(value)
        for value in (
            variant,
            record.get("variant"),
            record.get("display_name"),
            record.get("common_name"),
            *_items(record.get("lookup_keys")),
        )
        if _normalized_text(value)
    }
    observed_variants = sorted(
        {
            observed
            for key in keys
            for observed in matched_variants.get(key, [])
        }
    )
    linked_probes = {
        _clean(item) for item in _items(record.get("relevant_methylation_probe_ids")) if _clean(item)
    }
    matched_probes = sorted(linked_probes & observed_probes)
    basis: list[str] = []
    if observed_variants:
        basis.append("observed_variant")
    if matched_probes:
        basis.append("observed_methylation_probe")
    if not basis:
        basis.append("queried_gene")
    return {
        "observed_match": bool(observed_variants or matched_probes),
        "match_basis": basis,
        "observed_variants": observed_variants,
        "observed_methylation_probes": matched_probes,
    }


def _gene_context(payload: dict[str, Any], gene: str) -> dict[str, Any]:
    knowledge_base = payload.get("knowledge_base") if isinstance(payload.get("knowledge_base"), dict) else {}
    context = knowledge_base.get("gene_context") if isinstance(knowledge_base.get("gene_context"), dict) else {}
    variant_interpretations = (
        payload.get("variant_interpretations")
        if isinstance(payload.get("variant_interpretations"), dict)
        else {}
    )
    methylation_insights = (
        payload.get("methylation_insights")
        if isinstance(payload.get("methylation_insights"), dict)
        else {}
    )
    result = {
        "gene": gene,
        "gene_summary": _clean(context.get("gene_summary") or variant_interpretations.get("gene_summary")),
        "clinical_context": _clean(context.get("clinical_context") or variant_interpretations.get("clinical_context")),
        "condition_research_overview": _text_items(
            context.get("condition_research_overview") or variant_interpretations.get("condition_research_overview")
        ),
        "methylation_condition_research": _text_items(
            context.get("methylation_condition_research") or methylation_insights.get("methylation_condition_research")
        ),
        "variant_effect_overview": _text_items(
            context.get("variant_effect_overview") or variant_interpretations.get("variant_effect_overview")
        ),
        "methylation_effects": _text_items(context.get("methylation_effects") or methylation_insights.get("methylation_effects")),
        "methylation_interpretation": _clean(
            context.get("methylation_interpretation") or methylation_insights.get("whitelist_literature_context")
        ),
        "provenance": "curated_gene_knowledge_base" if context else "run_interpretation_context",
        **_tagged_context(match_basis=("queried_gene",)),
    }
    result["availability"] = "available" if any(
        value for key, value in result.items() if key not in {"gene", "provenance", "availability", "match_basis"}
        and key not in {"evidence_class", "authority_status", "research_context_only", "observed_match"}
    ) else "no_data"
    return result


def _variant_contexts(
    payload: dict[str, Any],
    matched_variants: dict[str, list[str]],
    observed_probes: set[str],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    knowledge_base = payload.get("knowledge_base") if isinstance(payload.get("knowledge_base"), dict) else {}
    records = knowledge_base.get("variant_records") if isinstance(knowledge_base.get("variant_records"), list) else []
    output: list[dict[str, Any]] = []
    lookup: dict[str, dict[str, Any]] = {}
    for raw in records:
        if not isinstance(raw, dict):
            continue
        variant = _clean(raw.get("variant") or raw.get("display_name"))
        match = _variant_match(variant, raw, matched_variants, observed_probes)
        item = {
            "variant": variant,
            "display_name": _clean(raw.get("display_name") or variant),
            "interpretation_scope": _clean(raw.get("interpretation_scope")),
            "clinical_significance": _clean(raw.get("clinical_significance")),
            "clinical_interpretation": _clean(raw.get("clinical_interpretation")),
            "associated_conditions": _text_items(raw.get("associated_conditions")),
            "clinical_parameter_summary": _clean(raw.get("clinical_parameter_summary")) or None,
            "functional_effects": _text_items(raw.get("functional_effects")),
            "methylation_interpretation": _clean(raw.get("methylation_interpretation")),
            "relevant_methylation_probe_ids": _text_items(raw.get("relevant_methylation_probe_ids")),
            "provenance": "curated_variant_record",
            **match,
            **_tagged_context(observed_match=match["observed_match"], match_basis=match["match_basis"]),
        }
        output.append(item)
        for key_value in (
            raw.get("variant"), raw.get("display_name"), raw.get("common_name"), *_items(raw.get("lookup_keys"))
        ):
            key = _normalized_text(key_value)
            if key:
                lookup[key] = raw
    output.sort(key=lambda item: (_match_rank(item), item["variant"].casefold()))
    return output, lookup


def _finding_match(
    finding: dict[str, Any],
    variant_lookup: dict[str, dict[str, Any]],
    matched_variants: dict[str, list[str]],
    observed_probes: set[str],
) -> dict[str, Any]:
    variant = _clean(finding.get("variant"))
    return _variant_match(variant, variant_lookup.get(_normalized_text(variant)), matched_variants, observed_probes)


_DESIGN_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("clinical_trial", re.compile(r"\b(?:clinical trial|randomi[sz]ed|placebo)\b", re.IGNORECASE)),
    ("meta_analysis", re.compile(r"\b(?:meta-analysis|meta analysis|pooled analysis)\b", re.IGNORECASE)),
    ("case_control", re.compile(r"\b(?:case-control|case control|cases?\b.+\bcontrols?)\b", re.IGNORECASE)),
    ("imaging", re.compile(r"\b(?:imaging|mri|fmri|connectivity|gray[- ]matter|grey[- ]matter)\b", re.IGNORECASE)),
    ("tissue_expression", re.compile(r"\b(?:postmortem|post-mortem|tissue|mrna|gene expression)\b", re.IGNORECASE)),
    ("functional_assay", re.compile(r"\b(?:reporter|luciferase|functional assay|cell line|in vitro|engineered)\b", re.IGNORECASE)),
    ("methylation_study", re.compile(r"\b(?:methylation|mqtl|cpg)\b", re.IGNORECASE)),
    ("cohort", re.compile(r"\bcohort\b", re.IGNORECASE)),
    ("association_study", re.compile(r"\b(?:association|correlation|risk|odds ratio)\b", re.IGNORECASE)),
)

_TISSUE_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("postmortem brain tissue", re.compile(r"\b(?:postmortem|post-mortem).{0,30}brain|brain.{0,30}(?:postmortem|post-mortem)\b", re.IGNORECASE)),
    ("superior frontal gyrus", re.compile(r"\bsuperior frontal gyrus\b", re.IGNORECASE)),
    ("lymphoblastoid cell line", re.compile(r"\blymphoblastoid cell lines?\b", re.IGNORECASE)),
    ("neural cell line", re.compile(r"\bneural cell lines?\b", re.IGNORECASE)),
    ("brain tissue", re.compile(r"\bbrain tissue\b", re.IGNORECASE)),
    ("blood", re.compile(r"\b(?:whole )?blood\b", re.IGNORECASE)),
    ("saliva", re.compile(r"\bsaliva(?:ry)?\b", re.IGNORECASE)),
    ("neuroimaging", re.compile(r"\b(?:neuroimaging|resting-state imaging|fmri|mri)\b", re.IGNORECASE)),
)

_SAMPLE_PATTERN = re.compile(
    r"\b(?:n\s*=\s*)?\d[\d,]*(?:\s+[a-z][a-z-]*){0,4}\s+"
    r"(?:cases?|controls?|patients?|participants?|subjects?|samples?|datasets?|cohorts?|adults?|children)\b",
    re.IGNORECASE,
)
_N_PATTERN = re.compile(r"\b[nN]\s*=\s*\d[\d,]*\b")
_EFFECT_PATTERNS = (
    re.compile(r"\bOR\s*[=:]?\s*\d+(?:\.\d+)?", re.IGNORECASE),
    re.compile(r"\b\d{1,3}%\s*CI\s*\d+(?:\.\d+)?\s*[-–]\s*\d+(?:\.\d+)?", re.IGNORECASE),
    re.compile(r"\b(?:about\s+)?\d+(?:\.\d+)?%", re.IGNORECASE),
)


def _study_parameters(text: str) -> dict[str, Any]:
    designs = [label for label, pattern in _DESIGN_PATTERNS if pattern.search(text)]
    if not designs:
        designs = ["other_or_unspecified"]
    sample_mentions = []
    for match in [*_SAMPLE_PATTERN.finditer(text), *_N_PATTERN.finditer(text)]:
        value = _clean(match.group(0))
        if value and value not in sample_mentions:
            sample_mentions.append(value)
    tissue_or_model = [label for label, pattern in _TISSUE_PATTERNS if pattern.search(text)]
    effect_mentions: list[str] = []
    for pattern in _EFFECT_PATTERNS:
        for match in pattern.finditer(text):
            value = _clean(match.group(0))
            if value and value not in effect_mentions:
                effect_mentions.append(value)
    return {
        "study_designs": designs,
        "sample_size_mentions": sample_mentions or None,
        "tissue_or_model": tissue_or_model or None,
        "reported_effect_measure_mentions": effect_mentions or None,
    }


def _cohort_studies(
    literature: dict[str, Any],
    variant_lookup: dict[str, dict[str, Any]],
    matched_variants: dict[str, list[str]],
    observed_probes: set[str],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for finding in literature.get("findings", []):
        if not isinstance(finding, dict) or finding.get("finding_status") == "citation_metadata_only":
            continue
        finding_text = _clean(finding.get("finding") or finding.get("summary") or finding.get("abstract"))
        phenotype = _clean(finding.get("phenotype"))
        paper = _clean(finding.get("paper") or finding.get("title"))
        genotype_groups = _clean(finding.get("genotypes"))
        combined = " ".join((phenotype, genotype_groups, finding_text, paper))
        parameters = _study_parameters(combined)
        match = _finding_match(finding, variant_lookup, matched_variants, observed_probes)
        unavailable: dict[str, str] = {}
        if not parameters["sample_size_mentions"]:
            unavailable["sample_size_mentions"] = "No explicit sample-size or group-count phrase is present in the available record."
        if not parameters["tissue_or_model"]:
            unavailable["tissue_or_model"] = "No explicit tissue or experimental-model phrase is present in the available record."
        if not parameters["reported_effect_measure_mentions"]:
            unavailable["reported_effect_measure_mentions"] = "No supported effect-measure phrase is present in the available record."
        if not phenotype:
            unavailable["cohort_or_model_context"] = "No phenotype, population, cohort, tissue, or model description is present."
        if not genotype_groups:
            unavailable["genotype_or_comparator_groups"] = "No genotype or comparator groups are present in the available record."
        item = {
            "rank": finding.get("rank"),
            "paper": paper,
            "pmid": finding.get("pmid"),
            "pmcid": finding.get("pmcid"),
            "doi": finding.get("doi"),
            "url": _clean(finding.get("url")),
            "variant": _clean(finding.get("variant")),
            "cohort_or_model_context": phenotype or None,
            "genotype_or_comparator_groups": genotype_groups or None,
            "finding": finding_text,
            "source_key": _clean(finding.get("source_key")),
            "sources": list(finding.get("sources") or []),
            "provenance": [
                value
                for value in (finding.get("evidence_origins") or [finding.get("evidence_origin")])
                if _clean(value)
            ],
            **parameters,
            "unavailable_parameters": unavailable,
            **match,
            **_tagged_context(observed_match=match["observed_match"], match_basis=match["match_basis"]),
        }
        output.append(item)
    output.sort(
        key=lambda item: (
            _match_rank(item),
            int(item.get("rank") or 10**9),
            item["paper"].casefold(),
        )
    )
    return output


def _record_conditions(record: dict[str, Any]) -> list[str]:
    conditions: list[str] = []
    for field in (
        "disease", "disease_name", "condition", "conditions", "trait", "phenotype", "phenotypes",
        "relevant_disorders", "indication", "disease_names", "panel_disease_group", "panel_disease_sub_group",
    ):
        for value in _text_items(record.get(field)):
            if value and value not in conditions:
                conditions.append(value)
    if not conditions and _clean(record.get("category")) == "clinical_condition":
        title = _clean(record.get("title") or record.get("label"))
        if title:
            conditions.append(title)
    return conditions


def _investigated_conditions(
    payload: dict[str, Any],
    records: list[dict[str, Any]],
    literature: dict[str, Any],
    variant_contexts: list[dict[str, Any]],
    variant_lookup: dict[str, dict[str, Any]],
    matched_variants: dict[str, list[str]],
    observed_probes: set[str],
) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}

    def add(
        label: Any,
        *,
        summary: Any = "",
        source_key: str = "curated_gene_knowledge_base",
        source_url: Any = "",
        variant: Any = "",
        provenance: str,
        evidence_class: str = "research_context",
        authority_status: str = "not_authoritative",
        match: dict[str, Any] | None = None,
    ) -> None:
        text = _clean(label)
        key = _normalized_text(text)
        if not key:
            return
        match = match or {"observed_match": False, "match_basis": ["queried_gene"]}
        evidence_link = {
            "source_key": source_key,
            "source_url": _clean(source_url),
            "variant": _clean(variant),
            "summary": _clean(summary),
            "provenance": provenance,
            "evidence_class": evidence_class,
            "authority_status": authority_status,
        }
        if key not in grouped:
            grouped[key] = {
                "condition_or_topic": text,
                "variants": [],
                "sources": [],
                "evidence_links": [],
                "provenance": [],
                **_tagged_context(
                    evidence_class=evidence_class,
                    authority_status=authority_status,
                    observed_match=bool(match.get("observed_match")),
                    match_basis=match.get("match_basis") or ("queried_gene",),
                ),
            }
        item = grouped[key]
        if evidence_link not in item["evidence_links"]:
            item["evidence_links"].append(evidence_link)
        for value, field in ((variant, "variants"), (source_key, "sources"), (provenance, "provenance")):
            cleaned = _clean(value)
            if cleaned and cleaned not in item[field]:
                item[field].append(cleaned)
        if match.get("observed_match"):
            item["observed_match"] = True
        item["match_basis"] = sorted(set(item["match_basis"]) | set(match.get("match_basis") or []))
        if authority_status == "authoritative":
            item["authority_status"] = "authoritative"
            item["research_context_only"] = False
            item["evidence_class"] = "authoritative_assertion"

    gene_context = _gene_context(payload, _clean(payload.get("gene") or payload.get("gene_name")).upper())
    for field in ("condition_research_overview", "methylation_condition_research"):
        for summary in gene_context.get(field, []):
            label = summary.split(":", 1)[0] if ":" in summary else summary
            add(label, summary=summary, provenance=field)
    for context in variant_contexts:
        match = {"observed_match": context["observed_match"], "match_basis": context["match_basis"]}
        for condition in context["associated_conditions"]:
            add(
                condition,
                summary=context["clinical_interpretation"],
                variant=context["variant"],
                provenance="curated_variant_associated_condition",
                match=match,
            )
    for finding in literature.get("findings", []):
        if not isinstance(finding, dict) or finding.get("finding_status") == "citation_metadata_only":
            continue
        match = _finding_match(finding, variant_lookup, matched_variants, observed_probes)
        if _clean(finding.get("phenotype")):
            add(
                finding.get("phenotype"),
                summary=finding.get("finding"),
                source_key=_clean(finding.get("source_key") or "literature"),
                source_url=finding.get("url"),
                variant=finding.get("variant"),
                provenance="study_phenotype",
                match=match,
            )
    for record in records:
        source_key = _clean(record.get("source_key") or record.get("source") or "unknown")
        canonical = canonical_medical_source(source_key)
        if canonical not in MEDICAL_CONTEXT_SOURCE_KEYS and not _record_conditions(record):
            continue
        gate = medical_gate(record)
        match = _variant_match(
            record.get("variant") or record.get("rsid"),
            variant_lookup.get(_normalized_text(record.get("variant") or record.get("rsid"))),
            matched_variants,
            observed_probes,
        )
        for condition in _record_conditions(record):
            add(
                condition,
                summary=record.get("summary") or record.get("definition") or record.get("assertion"),
                source_key=source_key,
                source_url=record.get("url"),
                variant=record.get("variant") or record.get("rsid"),
                provenance="dynamic_clinical_source",
                evidence_class="authoritative_assertion" if gate["eligible"] else "clinical_database_context",
                authority_status="authoritative" if gate["eligible"] else "not_authoritative",
                match=match,
            )
    output = list(grouped.values())
    output.sort(key=lambda item: (_match_rank(item), item["condition_or_topic"].casefold()))
    return output


def _pathology_contexts(
    gene_context: dict[str, Any],
    variant_contexts: list[dict[str, Any]],
    records: list[dict[str, Any]],
    variant_lookup: dict[str, dict[str, Any]],
    matched_variants: dict[str, list[str]],
    observed_probes: set[str],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()

    def add(
        context_type: str,
        summary: Any,
        *,
        variant: Any = "",
        source_key: str = "curated_gene_knowledge_base",
        source_url: Any = "",
        provenance: str,
        evidence_class: str = "research_context",
        authority_status: str = "not_authoritative",
        match: dict[str, Any] | None = None,
    ) -> None:
        text = _clean(summary)
        variant_text = _clean(variant)
        key = (context_type, _normalized_text(text), _normalized_text(variant_text))
        if not text or key in seen:
            return
        seen.add(key)
        match = match or {"observed_match": False, "match_basis": ["queried_gene"]}
        output.append(
            {
                "context_type": context_type,
                "summary": text,
                "variant": variant_text,
                "source_key": source_key,
                "source_url": _clean(source_url),
                "provenance": provenance,
                **_tagged_context(
                    evidence_class=evidence_class,
                    authority_status=authority_status,
                    observed_match=bool(match.get("observed_match")),
                    match_basis=match.get("match_basis") or ("queried_gene",),
                ),
            }
        )

    add("gene_function", gene_context.get("gene_summary"), provenance="gene_summary")
    add("clinical_context", gene_context.get("clinical_context"), provenance="gene_clinical_context")
    for summary in gene_context.get("variant_effect_overview", []):
        add("variant_mechanism_overview", summary, provenance="gene_variant_effect_overview")
    for summary in gene_context.get("methylation_effects", []):
        add("methylation_regulation", summary, provenance="gene_methylation_effect")
    add(
        "methylation_regulation",
        gene_context.get("methylation_interpretation"),
        provenance="gene_methylation_interpretation",
    )
    for context in variant_contexts:
        match = {"observed_match": context["observed_match"], "match_basis": context["match_basis"]}
        add(
            "variant_clinical_interpretation",
            context.get("clinical_interpretation"),
            variant=context["variant"],
            provenance="curated_variant_interpretation",
            match=match,
        )
        for summary in context.get("functional_effects", []):
            add(
                "variant_functional_effect",
                summary,
                variant=context["variant"],
                provenance="curated_variant_functional_effect",
                match=match,
            )
        add(
            "variant_methylation_context",
            context.get("methylation_interpretation"),
            variant=context["variant"],
            provenance="curated_variant_methylation_interpretation",
            match=match,
        )
    for record in records:
        if not _record_conditions(record) and canonical_medical_source(record.get("source_key") or record.get("source")) not in MEDICAL_CONTEXT_SOURCE_KEYS:
            continue
        summary = record.get("definition") or record.get("summary") or record.get("assertion") or record.get("classification")
        if not _clean(summary):
            continue
        gate = medical_gate(record)
        variant = record.get("variant") or record.get("rsid")
        match = _variant_match(
            variant,
            variant_lookup.get(_normalized_text(variant)),
            matched_variants,
            observed_probes,
        )
        add(
            "clinical_database_context",
            summary,
            variant=variant,
            source_key=_clean(record.get("source_key") or record.get("source") or "unknown"),
            source_url=record.get("url"),
            provenance="dynamic_clinical_source",
            evidence_class="authoritative_assertion" if gate["eligible"] else "clinical_database_context",
            authority_status="authoritative" if gate["eligible"] else "not_authoritative",
            match=match,
        )
    output.sort(key=lambda item: (_match_rank(item), item["context_type"], item["summary"].casefold()))
    return output


def _medical_source_states(statuses: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    assessed: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []
    not_assessed: list[dict[str, Any]] = []
    for raw in statuses:
        source_key = _clean(raw.get("source_key") or raw.get("source"))
        canonical = canonical_medical_source(source_key)
        if canonical not in MEDICAL_CONTEXT_SOURCE_KEYS:
            continue
        state = _clean(raw.get("status") or "not_assessed").casefold()
        item = {
            "source_key": source_key,
            "canonical_source_key": canonical,
            "status": state,
            "record_count": _safe_int(raw.get("record_count")),
            "snapshot_date": _clean(raw.get("retrieved_at") or raw.get("snapshot_date")),
        }
        if state in {"ok", "empty", "assessed", "imported"}:
            assessed.append(item)
        elif state in {"failed", "error", "timeout"}:
            failed.append(item)
        else:
            not_assessed.append(item)
    return assessed, failed, not_assessed


def build_medical_context(
    *,
    payload: dict[str, Any],
    gene: str,
    records: list[dict[str, Any]],
    statuses: list[dict[str, Any]],
    literature: dict[str, Any],
) -> dict[str, Any]:
    """Build layered medical context while retaining a strict authoritative record lane."""
    matched_variants, observed_probes = _observed_context(payload)
    gene_context = _gene_context(payload, gene)
    variant_contexts, variant_lookup = _variant_contexts(payload, matched_variants, observed_probes)
    cohort_studies = _cohort_studies(
        literature,
        variant_lookup,
        matched_variants,
        observed_probes,
    )
    established_records = []
    for record in records:
        gate = medical_gate(record)
        if gate["eligible"]:
            match = _variant_match(
                record.get("variant") or record.get("rsid"),
                variant_lookup.get(_normalized_text(record.get("variant") or record.get("rsid"))),
                matched_variants,
                observed_probes,
            )
            established_records.append(
                {
                    **record,
                    "medical_gate": gate,
                    "effective_date": _clean(record.get("effective_date")) or gate.get("effective_date"),
                    "provenance": _clean(record.get("evidence_origin")) or "authoritative_clinical_source",
                    **match,
                    **_tagged_context(
                        evidence_class="authoritative_assertion",
                        authority_status="authoritative",
                        observed_match=match["observed_match"],
                        match_basis=match["match_basis"],
                    ),
                }
            )
    established_records.sort(key=lambda item: (_match_rank(item), canonical_medical_source(item.get("source_key"))))
    investigated_conditions = _investigated_conditions(
        payload,
        records,
        literature,
        variant_contexts,
        variant_lookup,
        matched_variants,
        observed_probes,
    )
    pathology_contexts = _pathology_contexts(
        gene_context,
        variant_contexts,
        records,
        variant_lookup,
        matched_variants,
        observed_probes,
    )
    checked_sources, failed_sources, not_assessed_sources = _medical_source_states(statuses)
    if established_records:
        established_status = "assessed"
    elif checked_sources:
        established_status = "assessed_absence"
    elif failed_sources:
        established_status = "source_failed"
    else:
        established_status = "not_assessed"
    has_context = any(
        (
            gene_context.get("availability") == "available",
            investigated_conditions,
            cohort_studies,
            variant_contexts,
            pathology_contexts,
        )
    )
    context_status = "partial" if has_context and failed_sources else "available" if has_context else "no_data"
    observed_match_count = sum(
        int(bool(item.get("observed_match")))
        for collection in (investigated_conditions, cohort_studies, variant_contexts, pathology_contexts)
        for item in collection
    )
    counts = {
        "investigated_condition_count": len(investigated_conditions),
        "cohort_study_count": len(cohort_studies),
        "variant_context_count": len(variant_contexts),
        "pathology_context_count": len(pathology_contexts),
        "observed_match_count": observed_match_count,
        "authoritative_assertion_count": len(established_records),
    }
    inclusion_policy = "authoritative guidelines, labels, expert panels, and reviewed classifications only"
    established_evidence = {
        "status": established_status,
        "records": established_records,
        "checked_sources": checked_sources,
        "failed_sources": failed_sources,
        "not_assessed_sources": not_assessed_sources,
        "inclusion_policy": inclusion_policy,
    }
    return {
        "status": established_status,
        "context_status": context_status,
        "scope": "queried_gene",
        "selection_policy": "all_gene_context_observed_matches_first",
        "counts": counts,
        "gene_context": gene_context,
        "investigated_conditions": investigated_conditions,
        "cohort_studies": cohort_studies,
        "variant_context": variant_contexts,
        "pathology_context": pathology_contexts,
        "established_evidence": established_evidence,
        # Compatibility fields remain authoritative-only.
        "records": established_records,
        "checked_sources": checked_sources,
        "failed_sources": failed_sources,
        "not_assessed_sources": not_assessed_sources,
        "research_use_only": True,
        "disclaimer": (
            "Research-use-only gene and condition context; investigated associations and mechanisms are not "
            "diagnoses, causal findings, pathogenic classifications, or treatment recommendations."
        ),
        "inclusion_policy": inclusion_policy,
        "research_inclusion_policy": (
            "All available queried-gene condition, cohort, variant, methylation, and pathology context is shown; "
            "exact observed variant or methylation matches are ranked first."
        ),
    }
