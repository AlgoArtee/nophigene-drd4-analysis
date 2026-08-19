"""Canonical schema-3 report builder and compact evidence-first HTML renderer."""

from __future__ import annotations

import html
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .evidence import (
    LITERATURE_SOURCE_KEYS,
    build_literature_catalog,
    medical_gate,
    normalize_literature_record,
    source_coverage,
)
from .interactions import INTERACTION_SOURCE_KEYS, build_interaction_graph, bundled_interaction_records
from .medical import build_medical_context
from .model_registry import inspect_model_inputs, list_model_manifests
from .pgx import resolve_pgx_diplotype
from .predictions import (
    mark_exact_source_matches,
    observed_variant_alleles,
    select_alphagenome_variants,
    source_native_annotations,
)
from .statistics import build_single_person_statistics

REPORT_SCHEMA_VERSION = "3.0"
DEFAULT_PAGE_SIZE = 20


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, pd.DataFrame):
        return [_jsonable(row) for row in value.to_dict(orient="records")]
    if isinstance(value, pd.Series):
        return _jsonable(value.to_dict())
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if value is pd.NA:
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _column(frame: pd.DataFrame, *names: str) -> str | None:
    lookup = {str(name).casefold(): str(name) for name in frame.columns}
    for name in names:
        if name.casefold() in lookup:
            return lookup[name.casefold()]
    return None


def _text(value: Any) -> str:
    if value is None or value is pd.NA:
        return ""
    return str(value).strip()


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _variant_records(frame: pd.DataFrame) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if frame is None or frame.empty:
        return [], []
    filter_col = _column(frame, "FILTER", "filter", "filter_status")
    gq_col = _column(frame, "GQ", "genotype_quality")
    dp_col = _column(frame, "DP", "depth")
    gt_col = _column(frame, "GT", "gt_raw", "genotype")
    alt_col = _column(frame, "ALT", "alt")
    primary: list[dict[str, Any]] = []
    all_rows: list[dict[str, Any]] = []
    for row in frame.to_dict(orient="records"):
        reasons: list[str] = []
        filter_value = _text(row.get(filter_col)) if filter_col else ""
        if filter_value and filter_value not in {"PASS", "."}:
            reasons.append(f"filter:{filter_value}")
        gq = _number(row.get(gq_col)) if gq_col else None
        dp = _number(row.get(dp_col)) if dp_col else None
        if gq is not None and gq < 20:
            reasons.append("GQ<20")
        if dp is not None and dp < 10:
            reasons.append("DP<10")
        genotype = _text(row.get(gt_col)) if gt_col else ""
        non_reference = genotype not in {"", ".", "./.", ".|.", "0/0", "0|0"}
        if not genotype:
            reasons.append("genotype_missing")
        output = {**_jsonable(row), "qc_pass": not reasons, "qc_reasons": reasons, "non_reference": non_reference}
        if alt_col and "," in _text(row.get(alt_col)):
            output["normalization_status"] = "requires_multiallelic_split"
        all_rows.append(output)
        if output["qc_pass"] and non_reference:
            primary.append(output)
    primary.sort(
        key=lambda row: (
            _number(row.get(gq_col)) if gq_col and _number(row.get(gq_col)) is not None else -1.0,
            _number(row.get(dp_col)) if dp_col and _number(row.get(dp_col)) is not None else -1.0,
        ),
        reverse=True,
    )
    for rank, row in enumerate(primary, start=1):
        row["objective_rank"] = rank
    return primary, all_rows


def _methylation_records(frame: pd.DataFrame) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if frame is None or frame.empty:
        return [], []
    detection_col = _column(frame, "detection_p", "Detection Pval", "detection_p_value")
    bead_col = _column(frame, "bead_count", "NBeads", "beads")
    flag_columns = [
        column
        for column in frame.columns
        if any(token in str(column).casefold() for token in ("cross_reactive", "snp_overlap", "probe_flag"))
    ]
    primary: list[dict[str, Any]] = []
    all_rows: list[dict[str, Any]] = []
    for row in frame.to_dict(orient="records"):
        reasons: list[str] = []
        detection = _number(row.get(detection_col)) if detection_col else None
        beads = _number(row.get(bead_col)) if bead_col else None
        if detection is not None and detection > 0.01:
            reasons.append("detection_p>0.01")
        if beads is not None and beads < 3:
            reasons.append("bead_count<3")
        for column in flag_columns:
            value = row.get(column)
            if value is True or _text(value).casefold() in {"true", "yes", "1", "flagged"}:
                reasons.append(str(column))
        output = {**_jsonable(row), "qc_pass": not reasons, "qc_reasons": reasons}
        all_rows.append(output)
        if output["qc_pass"]:
            primary.append(output)
    primary.sort(
        key=lambda row: (
            _number(row.get(detection_col)) if detection_col and _number(row.get(detection_col)) is not None else float("inf"),
            -(_number(row.get(bead_col)) if bead_col and _number(row.get(bead_col)) is not None else -1.0),
        )
    )
    for rank, row in enumerate(primary, start=1):
        row["objective_rank"] = rank
    return primary, all_rows


def _source_records(payload: dict[str, Any], gene: str = "") -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    def append_source_record(
        item: dict[str, Any], *, evidence_origin: str, parent_variant: str = "", publications_only: bool = False
    ) -> None:
        if (
            item.get("partner_gene")
            or item.get("target_gene") and item.get("edge_type")
            or _text(item.get("category")).casefold() in {
                "interaction", "protein_interaction", "protein_functional_association"
            }
        ):
            if not publications_only:
                records.append({**item, "evidence_origin": item.get("evidence_origin") or evidence_origin})
            return
        normalized = normalize_literature_record(
            item,
            gene=gene,
            evidence_origin=evidence_origin,
            parent_variant=parent_variant,
        )
        if normalized["is_publication_record"]:
            records.append(normalized)
        elif not publications_only:
            records.append(dict(item))

    containers = (
        (payload, "run_source_record"),
        (payload.get("knowledge_base", {}), "curated_source_record"),
        (payload.get("dynamic_knowledge_base", {}), "dynamic_source_record"),
    )
    for container, origin in containers:
        if not isinstance(container, dict):
            continue
        for key in ("source_records", "records", "literature_records", "interaction_records"):
            value = container.get(key)
            if isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        append_source_record(item, evidence_origin=origin)
        gene_context = container.get("gene_context")
        if isinstance(gene_context, dict):
            for item in gene_context.get("evidence", []):
                if isinstance(item, dict):
                    append_source_record(
                        item,
                        evidence_origin="curated_gene_context",
                        publications_only=True,
                    )
        variant_records = container.get("variant_records")
        if isinstance(variant_records, list):
            for variant_record in variant_records:
                if not isinstance(variant_record, dict):
                    continue
                parent_variant = _text(
                    variant_record.get("variant")
                    or variant_record.get("display_name")
                    or variant_record.get("common_name")
                )
                variant_gene = _text(variant_record.get("gene_name") or gene)
                for item in variant_record.get("evidence", []):
                    if isinstance(item, dict):
                        append_source_record(
                            {**item, "gene": variant_gene},
                            evidence_origin="curated_variant_evidence",
                            parent_variant=parent_variant,
                            publications_only=True,
                        )
                for item in variant_record.get("literature_findings", []):
                    if isinstance(item, dict):
                        append_source_record(
                            {**item, "gene": variant_gene},
                            evidence_origin="curated_variant_finding",
                            parent_variant=parent_variant,
                            publications_only=True,
                        )
    dynamic = payload.get("dynamic_knowledge_base", {})
    if isinstance(dynamic, dict):
        local_articles = dynamic.get("local_article_evidence", {})
        if isinstance(local_articles, dict) and isinstance(local_articles.get("records"), list):
            for item in local_articles["records"]:
                if isinstance(item, dict):
                    append_source_record(item, evidence_origin="local_pdf")
    return records


def _provider_statuses(payload: dict[str, Any]) -> list[dict[str, Any]]:
    statuses: list[dict[str, Any]] = []
    for container in (payload, payload.get("knowledge_base", {}), payload.get("dynamic_knowledge_base", {})):
        if isinstance(container, dict) and isinstance(container.get("provider_statuses"), list):
            statuses.extend(dict(item) for item in container["provider_statuses"] if isinstance(item, dict))
    snapshot = payload.get("interpretation", {}).get("evidence_snapshot", {})
    if isinstance(snapshot, dict) and isinstance(snapshot.get("providers"), list):
        statuses.extend(dict(item) for item in snapshot["providers"] if isinstance(item, dict))
    dynamic = payload.get("dynamic_knowledge_base", {})
    if isinstance(dynamic, dict):
        for workflow in dynamic.get("workflow_runs", []):
            if isinstance(workflow, dict) and isinstance(workflow.get("provider_statuses"), list):
                statuses.extend(dict(item) for item in workflow["provider_statuses"] if isinstance(item, dict))
    return statuses


def _literature_and_medical(
    records: list[dict[str, Any]], statuses: list[dict[str, Any]]
) -> tuple[dict[str, Any], dict[str, Any]]:
    literature_candidates = [
        record
        for record in records
        if record.get("is_publication_record")
        or record.get("title")
        or record.get("pmid")
        or record.get("pmcid")
        or record.get("doi")
    ]
    publications, findings = build_literature_catalog(literature_candidates)
    literature_categories: dict[str, list[dict[str, Any]]] = {
        "primary_studies": [],
        "reviews": [],
        "preprints": [],
        "case_reports": [],
        "database_assertions": [],
    }
    for record in publications:
        publication_type = _text(record.get("publication_type") or record.get("evidence_type")).casefold()
        if record.get("preprint"):
            category = "preprints"
        elif "review" in publication_type or "meta-analysis" in publication_type:
            category = "reviews"
        elif "case report" in publication_type or publication_type == "case_report":
            category = "case_reports"
        elif "database" in publication_type or publication_type == "database_assertion":
            category = "database_assertions"
        else:
            category = "primary_studies"
        literature_categories[category].append(record)
    priority_category_counts = {
        str(tier): {
            "label": label,
            "count": sum(1 for finding in findings if int(finding.get("priority_tier") or 0) == tier),
        }
        for tier, label in (
            (1, "Clinical relevance + direct experimental/functional evidence"),
            (2, "Clinical relevance"),
            (3, "Direct experimental/functional evidence"),
            (4, "Review, meta-analysis, or replicated evidence"),
            (5, "Human association or observational finding"),
            (6, "Other gene-relevant evidence"),
            (7, "Preprint"),
        )
    }
    medical: list[dict[str, Any]] = []
    for record in records:
        gate = medical_gate(record)
        if gate["eligible"]:
            medical.append({**record, "medical_gate": gate})
    compact_states = [
        {
            "source_key": _text(item.get("source_key") or item.get("source")),
            "status": _text(item.get("status") or "not_assessed"),
            "record_count": int(item.get("record_count") or 0),
            "snapshot_date": _text(item.get("retrieved_at") or item.get("snapshot_date")),
        }
        for item in statuses
    ]
    assessed_states = [item for item in compact_states if item["status"].casefold() in {"ok", "empty", "assessed", "imported"}]
    failed_states = [item for item in compact_states if item["status"].casefold() in {"failed", "error", "timeout"}]
    not_assessed_states = [
        item
        for item in compact_states
        if item not in assessed_states and item not in failed_states
    ]
    if medical:
        medical_status = "assessed"
    elif assessed_states:
        medical_status = "assessed_absence"
    elif failed_states:
        medical_status = "source_failed"
    else:
        medical_status = "not_assessed"
    failed_source_names = sorted(
        {
            item["source_key"]
            for item in failed_states
            if item["source_key"].casefold().replace(" ", "_") in LITERATURE_SOURCE_KEYS
        }
    )
    unavailable_source_names = sorted(
        {
            item["source_key"]
            for item in not_assessed_states
            if item["source_key"].casefold().replace(" ", "_") in LITERATURE_SOURCE_KEYS
        }
    )
    limitations = [
        "Coverage includes all gene-related evidence available to this run; it is not an exhaustive internet-wide search.",
        "Priority reflects evidence type and clinical relevance, not proof of causality or medical applicability.",
    ]
    if failed_source_names:
        limitations.append("Unavailable or failed literature sources: " + ", ".join(failed_source_names) + ".")
    if unavailable_source_names:
        limitations.append("Literature sources not assessed in this run: " + ", ".join(unavailable_source_names) + ".")
    return (
        {
            "status": "assessed" if literature_candidates else "not_assessed",
            "coverage_scope": "all_available_gene_related_evidence",
            "candidate_limit": len(publications),
            "detailed_review_limit": len(findings),
            "candidate_count": len(publications),
            "publication_count": len(publications),
            "finding_count": len(findings),
            "source_count": len(
                {
                    source
                    for finding in findings
                    for source in (finding.get("sources") or [finding.get("source_key")])
                    if source
                }
            ),
            "findings": findings,
            "publications": publications,
            "detailed_records": [record for record in findings if not record.get("preprint")],
            "candidate_index": publications,
            "categories": literature_categories,
            "priority_category_counts": priority_category_counts,
            "preprints": [record for record in findings if record.get("preprint")],
            "preprints_separated": True,
            "available_records_complete": True,
            "ranking_policy": {
                "method": "deterministic_evidence_priority",
                "tiers": priority_category_counts,
                "tie_breakers": [
                    "finding before citation-only metadata",
                    "observed entity and gene-context match",
                    "publication year descending",
                    "paper, variant, and finding text",
                ],
                "disclaimer": "Direct experimental/functional evidence does not by itself establish causality or clinical validity.",
            },
            "limitations": limitations,
            "generative_synthesis": {"status": "not_requested", "claims": []},
        },
        {
            "status": medical_status,
            "records": medical,
            "checked_sources": assessed_states,
            "failed_sources": failed_states,
            "research_use_only": True,
            "disclaimer": "Research-use-only evidence review; not a diagnosis or treatment recommendation.",
            "inclusion_policy": "authoritative guidelines, labels, expert panels, and reviewed classifications only",
        },
    )


def _interaction_edges(records: Iterable[dict[str, Any]], query_gene: str) -> list[dict[str, Any]]:
    edges_by_key: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    allowed = {
        "physical_binding", "genetic_interaction", "regulatory", "pathway_co_membership",
        "coexpression", "functional_association", "predicted_functional_association",
    }
    for record in records:
        edge_type = _text(record.get("edge_type") or record.get("interaction_type")).casefold().replace(" ", "_")
        partner = _text(record.get("partner_gene") or record.get("target_gene")).upper()
        source_gene = _text(record.get("source_gene") or query_gene).upper()
        if edge_type not in allowed or not partner:
            continue
        directed = bool(record.get("directed", False))
        if not directed and partner == query_gene.upper() and source_gene != query_gene.upper():
            source_gene, partner = partner, source_gene
        score_value = record.get("native_score")
        if score_value is None:
            score_value = record.get("score")
        edge = {
            "source_gene": source_gene,
            "target_gene": partner,
            "edge_type": edge_type,
            "directed": directed,
            "source_key": _text(record.get("source_key") or "unknown"),
            "source_record_id": _text(record.get("record_id") or record.get("source_record_id")),
            "native_score": _number(score_value),
            "native_score_label": _text(record.get("score_type")),
            "tissue": _text(record.get("tissue")),
            "association_scope": _text(record.get("association_scope"))
            or "functional association; not necessarily direct physical binding",
            "summary": _text(record.get("summary")),
            "url": _text(record.get("url")),
            "evidence": _jsonable(record),
        }
        key = (source_gene, partner, edge_type, edge["source_key"].casefold())
        current = edges_by_key.get(key)
        current_origin = _text((current or {}).get("evidence", {}).get("evidence_origin"))
        candidate_origin = _text(record.get("evidence_origin"))
        prefer_candidate = current is None or (
            current_origin == "bundled_versioned_interaction_snapshot"
            and candidate_origin != "bundled_versioned_interaction_snapshot"
        )
        if not prefer_candidate and current is not None:
            prefer_candidate = (edge.get("native_score") or -1.0) > (current.get("native_score") or -1.0)
        if prefer_candidate:
            edges_by_key[key] = edge
    return list(edges_by_key.values())


def _prediction_section(
    payload: dict[str, Any], context: dict[str, Any], primary_variants: list[dict[str, Any]]
) -> dict[str, Any]:
    interpretation = payload.get("interpretation", {}) if isinstance(payload.get("interpretation"), dict) else {}
    legacy_assessments = interpretation.get("model_assessments", [])
    requested = payload.get("requested_models") or []
    requested_ids: list[str] = []
    for item in requested:
        if isinstance(item, str):
            requested_ids.append(item)
        elif isinstance(item, dict) and item.get("model_id"):
            requested_ids.append(str(item["model_id"]))
    manifests = list_model_manifests()
    manifest_ids = {item["id"] for item in manifests}
    executable_requested_ids = list(dict.fromkeys(item for item in requested_ids if item in manifest_ids))
    input_availability = dict(context.get("model_inputs") or {})
    input_availability.setdefault("genome_build", str(payload.get("genome_build") or ""))
    input_assessments = [inspect_model_inputs(model_id, input_availability) for model_id in executable_requested_ids]
    observed_alleles = observed_variant_alleles(
        primary_variants,
        genome_build=str(payload.get("genome_build") or ""),
        scope_regions=dict(payload.get("scope_regions") or {}),
    )
    observed_alleles = mark_exact_source_matches(payload, observed_alleles)
    native_annotations, unmatched_predictor_records = source_native_annotations(payload, observed_alleles)
    variant_selection = select_alphagenome_variants(observed_alleles)
    pgx = None
    if isinstance(payload.get("pgx_definition"), dict):
        pgx = resolve_pgx_diplotype(payload["pgx_definition"], list(payload.get("pgx_calls") or []))
    if native_annotations:
        status = "available"
    elif executable_requested_ids:
        status = "blocked" if input_assessments and all(not item.get("eligible") for item in input_assessments) else "not_requested"
    elif observed_alleles:
        status = "not_requested"
    else:
        status = "no_data"
    available_models = []
    for manifest in manifests:
        applicability = inspect_model_inputs(manifest["id"], input_availability)
        if manifest["id"] == "alphagenome-api":
            runtime_availability = "configuration_required"
        elif manifest.get("status") in {
            "unsupported_for_epic_5mc", "blocked_pending_verified_release", "framework_only"
        }:
            runtime_availability = "scientifically_blocked"
        elif not applicability.get("eligible"):
            runtime_availability = "unsupported_input"
        else:
            runtime_availability = "missing_adapter_assets"
        available_models.append(
            {
                **manifest,
                "runtime_availability": runtime_availability,
                "input_blockers": list(applicability.get("blockers") or []),
            }
        )
    return {
        "status": status,
        "scope": "single_person",
        "selection_policy": "exact_source_alleles_plus_independent_opt_in_models",
        "requested_models": executable_requested_ids,
        "assessments": _jsonable(legacy_assessments),
        "model_input_assessments": input_assessments,
        "model_metadata_assessments": _jsonable(legacy_assessments),
        "legacy_metadata_assessments": _jsonable(legacy_assessments),
        "source_native_annotations": native_annotations,
        "source_native_context_only_record_count": unmatched_predictor_records,
        "variant_selection": variant_selection,
        "model_runs": [],
        "predictions": [],
        "counts": {
            "observed_allele_count": len(observed_alleles),
            "source_native_annotation_count": len(native_annotations),
            "requested_model_count": len(executable_requested_ids),
            "model_run_count": 0,
            "completed_prediction_count": 0,
        },
        "blockers": variant_selection.get("blockers", []),
        "failures": [],
        "pgx_diplotype": pgx,
        "consensus": None,
        "policy": "independent model outputs; no consensus or app-defined low/medium/high labels",
        "limitations": [
            "Source-native predictor scales and directions are source-specific and are not comparable.",
            "Model outputs are molecular research predictions, not statistical significance measurements, diagnoses, disease risk, or clinical classifications.",
            "An empty model response is not evidence that a variant has no biological effect.",
        ],
        "available_models": available_models,
    }


def _statistics_section(
    payload: dict[str, Any],
    *,
    gene: str,
    variants: pd.DataFrame,
    primary_variants: list[dict[str, Any]],
    all_variants: list[dict[str, Any]],
    methylation: pd.DataFrame,
    primary_methylation: list[dict[str, Any]],
    all_methylation: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build statistics only from the current person's observed result rows."""
    knowledge_base = payload.get("knowledge_base") if isinstance(payload.get("knowledge_base"), dict) else {}
    gene_context = knowledge_base.get("gene_context") if isinstance(knowledge_base.get("gene_context"), dict) else {}
    methylation_insights = (
        payload.get("methylation_insights") if isinstance(payload.get("methylation_insights"), dict) else {}
    )
    curated_probe_ids = gene_context.get("relevant_methylation_probe_ids") or methylation_insights.get("probe_ids") or []
    result = build_single_person_statistics(
        gene=gene,
        variants=variants,
        all_variant_rows=all_variants,
        primary_variant_rows=primary_variants,
        methylation=methylation,
        all_methylation_rows=all_methylation,
        primary_methylation_rows=primary_methylation,
        scope_regions=dict(payload.get("scope_regions") or {}),
        analysis_scope=_text(payload.get("analysis_scope")),
        active_region=_text(payload.get("region")),
        curated_probe_ids=curated_probe_ids,
    )
    result["records"][0].update(
        {
            "genome_build": _text(payload.get("genome_build")),
            "analysis_scope": _text(payload.get("analysis_scope")),
            "region": _text(payload.get("region")),
            "scope_regions": _jsonable(payload.get("scope_regions") or {}),
        }
    )
    return result


def build_canonical_report(payload: dict[str, Any]) -> dict[str, Any]:
    gene = _text(payload.get("gene") or payload.get("gene_name")).upper()
    variants = payload.get("variants")
    methylation = payload.get("methylation")
    if not isinstance(variants, pd.DataFrame):
        variants = pd.DataFrame(variants or [])
    if not isinstance(methylation, pd.DataFrame):
        methylation = pd.DataFrame(methylation or [])
    primary_variants, all_variants = _variant_records(variants)
    primary_methylation, all_methylation = _methylation_records(methylation)
    records = _source_records(payload, gene)
    existing_record_ids = {
        _text(record.get("record_id")) for record in records if _text(record.get("record_id"))
    }
    for record in bundled_interaction_records(gene):
        if _text(record.get("record_id")) not in existing_record_ids:
            records.append(record)
            existing_record_ids.add(_text(record.get("record_id")))
    statuses = _provider_statuses(payload)
    coverage = source_coverage(statuses)
    literature, _legacy_medical = _literature_and_medical(records, statuses)
    medical = build_medical_context(
        payload=payload,
        gene=gene,
        records=records,
        statuses=statuses,
        literature=literature,
    )
    context = dict(payload.get("sample_context") or {})
    if not context:
        context = dict(
            payload.get("interpretation", {})
            .get("interpretation_context", {})
            .get("sample_context", {})
        )
    interaction_evidence_edges = _interaction_edges(records, gene)
    graph = build_interaction_graph(gene, interaction_evidence_edges, max_hops=1, node_cap=150)
    interaction_sources = sorted({edge["source_key"] for edge in interaction_evidence_edges})
    interaction_source_statuses = [
        dict(status)
        for status in statuses
        if _text(status.get("source_key") or status.get("source")).casefold().replace(" ", "_")
        in INTERACTION_SOURCE_KEYS
    ]
    reported_interaction_source_keys = {
        _text(status.get("source_key") or status.get("source")).casefold().replace(" ", "_")
        for status in interaction_source_statuses
    }
    for source_key in interaction_sources:
        if source_key.casefold() in reported_interaction_source_keys:
            continue
        source_edges = [edge for edge in interaction_evidence_edges if edge["source_key"] == source_key]
        evidence_rows = [edge.get("evidence", {}) for edge in source_edges]
        bundled = bool(evidence_rows) and all(
            row.get("evidence_origin") == "bundled_versioned_interaction_snapshot"
            for row in evidence_rows
        )
        interaction_source_statuses.append(
            {
                "source_key": source_key,
                "status": "bundled_snapshot" if bundled else "available_in_report",
                "record_count": len(source_edges),
                "source_release": next(
                    (_text(row.get("source_release")) for row in evidence_rows if row.get("source_release")),
                    "",
                ),
                "snapshot_date": next(
                    (_text(row.get("snapshot_date")) for row in evidence_rows if row.get("snapshot_date")),
                    "",
                ),
            }
        )
    failed_interaction_sources = [
        status
        for status in interaction_source_statuses
        if _text(status.get("status")).casefold() in {"failed", "error", "timeout"}
    ]
    interaction_assessed = bool(graph["edge_count"]) or any(
        _text(status.get("status")).casefold() not in {"failed", "error", "timeout", "not_assessed", "skipped"}
        for status in interaction_source_statuses
    )
    missing_context = [
        key for key in ("tissue", "platform", "normalization", "ancestry", "phenotype") if not context.get(key)
    ]
    snapshot = payload.get("interpretation", {}).get("evidence_snapshot", {})
    report = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "generated_at": utc_now_iso(),
        "run": {
            "id": _text(payload.get("run_id") or payload.get("job_id")),
            "gene": gene,
            "genes": _jsonable(payload.get("genes") or ([gene] if gene else [])),
            "genome_build": _text(payload.get("genome_build")),
            "region": _text(payload.get("region")),
            "analysis_scope": _text(payload.get("analysis_scope")),
            "scope_regions": _jsonable(payload.get("scope_regions") or {}),
            "sample_context": _jsonable(context),
            "evidence_snapshot_id": _text(snapshot.get("snapshot_id")) if isinstance(snapshot, dict) else "",
            "status": _text(payload.get("status") or "succeeded"),
        },
        "summary": {
            "qc": {
                "variant_pass_count": len(primary_variants),
                "variant_total_count": len(all_variants),
                "methylation_pass_count": len(primary_methylation),
                "methylation_total_count": len(all_methylation),
            },
            "measured_findings": {
                "non_reference_variants": len(primary_variants),
                "qc_passing_methylation_probes": len(primary_methylation),
            },
            "evidence_coverage": {
                "assessed_source_count": len(coverage["assessed"]),
                "failed_source_count": len(coverage["failed"]),
                "not_assessed_source_count": len(coverage["not_assessed"]),
                "assessed_absence": coverage["assessed_absence"],
            },
            "missing_context": missing_context,
            "research_use_only": True,
        },
        "sections": {
            "objective_data": {
                "status": "assessed",
                "default_page_size": DEFAULT_PAGE_SIZE,
                "variants": primary_variants[:DEFAULT_PAGE_SIZE],
                "methylation": primary_methylation[:DEFAULT_PAGE_SIZE],
                "all_variants": all_variants,
                "all_methylation": all_methylation,
                "unsupported_variant_classes": ["structural_variants", "copy_number_variants", "VNTRs"],
                "policy": "measured observations only; no phenotype inference",
            },
            "statistics": _statistics_section(
                payload,
                gene=gene,
                variants=variants,
                primary_variants=primary_variants,
                all_variants=all_variants,
                methylation=methylation,
                primary_methylation=primary_methylation,
                all_methylation=all_methylation,
            ),
            "literature": literature,
            "medical": medical,
            "interactions": {
                "status": "assessed" if interaction_assessed else "not_assessed",
                "coverage_status": (
                    "partial" if graph["edge_count"] and failed_interaction_sources
                    else "available" if graph["edge_count"]
                    else "source_failed" if failed_interaction_sources
                    else "no_data"
                ),
                "initial_graph": graph,
                "direct_partner_count": graph["hop_counts"].get(1, 0),
                "source_count": len(interaction_sources),
                "sources": interaction_sources,
                "source_statuses": interaction_source_statuses,
                "failed_sources": failed_interaction_sources,
                "max_hops": 3,
                "node_cap": 150,
                "hop_semantics": "gene_to_gene",
                "association_policy": (
                    "Functional associations are source-backed gene/protein relationships and do not necessarily "
                    "represent direct physical binding, directionality, causality, or a person-specific effect."
                ),
                "selection_policy": "direct_queried_gene_partners_ranked_by_source_tier_corroboration_and_native_score",
            },
            "predictions": _prediction_section(payload, context, primary_variants),
        },
        "run_details": {
            "warnings": _jsonable(payload.get("warnings") or []),
            "source_provenance": _jsonable(payload.get("source_provenance") or {}),
            "source_coverage": _jsonable(coverage),
            "evidence_snapshot": _jsonable(snapshot),
            "checksums": _jsonable(payload.get("checksums") or {}),
            "artifacts": _jsonable(payload.get("artifacts") or {}),
            "raw_interpretation": _jsonable(payload.get("interpretation") or {}),
            "dynamic_knowledge_base": _jsonable(payload.get("dynamic_knowledge_base") or {}),
            "normalized_evidence_records": _jsonable(records),
            "provider_statuses": _jsonable(statuses),
            "interaction_evidence_edges": _jsonable(interaction_evidence_edges),
        },
    }
    assert "predictive_theses" not in report
    return report


def _cell(value: Any) -> str:
    if isinstance(value, (dict, list)):
        value = json.dumps(value, ensure_ascii=False, sort_keys=True)
    return html.escape(_text(value))


def _table(rows: Any, *, limit: int | None = DEFAULT_PAGE_SIZE) -> str:
    if not isinstance(rows, list) or not rows:
        return '<p class="empty">No eligible records were available for this section.</p>'
    selected = rows if limit is None else rows[:limit]
    normalized = [row if isinstance(row, dict) else {"value": row} for row in selected]
    columns: list[str] = []
    for row in normalized:
        for key in row:
            if key not in columns and key not in {"raw_fields", "evidence", "abstract"}:
                columns.append(key)
        if len(columns) >= 10:
            break
    header = "".join(f"<th scope=\"col\">{html.escape(key.replace('_', ' ').title())}</th>" for key in columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{_cell(row.get(key))}</td>" for key in columns) + "</tr>"
        for row in normalized
    )
    return f'<div class="table-shell"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>'


def _metric_cards(items: Iterable[tuple[str, Any]]) -> str:
    return '<div class="metrics">' + "".join(
        f"<article><span>{html.escape(label)}</span><strong>{_cell(value)}</strong></article>"
        for label, value in items
    ) + "</div>"


def render_evidence_first_html(report: dict[str, Any]) -> str:
    run = report["run"]
    summary = report["summary"]
    sections = report["sections"]
    qc = summary["qc"]
    tab_specs = (
        ("summary", "Summary"),
        ("objective_data", "Objective Data"),
        ("statistics", "Statistics"),
        ("literature", "Scientific Literature"),
        ("medical", "Medical Information"),
        ("interactions", "Interactions"),
        ("predictions", "Predictions"),
    )
    tabs = "".join(
        f'<button class="tab" role="tab" aria-selected="{str(key == "summary").lower()}" '
        f'aria-controls="panel-{key}" id="tab-{key}" data-tab="{key}">{html.escape(label)}</button>'
        for key, label in tab_specs
    )
    missing_context = ", ".join(summary.get("missing_context", [])) or "none"
    summary_panel = f"""
      <section id="panel-summary" class="panel active" role="tabpanel" aria-labelledby="tab-summary">
        <h2>Run validity and measured findings</h2>
        <div class="metrics">
          <article><span>QC variants</span><strong>{qc['variant_pass_count']} / {qc['variant_total_count']}</strong></article>
          <article><span>QC methylation probes</span><strong>{qc['methylation_pass_count']} / {qc['methylation_total_count']}</strong></article>
          <article><span>Assessed evidence sources</span><strong>{summary['evidence_coverage']['assessed_source_count']}</strong></article>
          <article><span>Source failures</span><strong>{summary['evidence_coverage']['failed_source_count']}</strong></article>
        </div>
        <p><strong>Context unavailable:</strong> {html.escape(missing_context)}. Medical applicability or external model eligibility may remain unassessed; descriptive statistics are still calculated from the observed rows.</p>
        <p class="notice">Research-use-only workbench. Measured observations, statistics, established medical evidence, and predictions are intentionally separated.</p>
      </section>
    """
    objective = sections["objective_data"]
    objective_panel = f"""
      <section id="panel-objective_data" class="panel" role="tabpanel" aria-labelledby="tab-objective_data" hidden>
        <h2>Objective Data</h2><p>QC-passing non-reference calls and methylation measurements. Raw and failed-QC rows remain in JSON/CSV exports.</p>
        <h3>Variants</h3>{_table(objective['variants'])}
        <h3>Methylation</h3>{_table(objective['methylation'])}
      </section>
    """
    statistics = sections["statistics"]
    variant_statistics = statistics["variant_statistics"]
    methylation_statistics = statistics["methylation_statistics"]
    variant_counts = variant_statistics["counts"]
    methylation_counts = methylation_statistics["counts"]
    variant_region_counts = {row["category"]: row["count"] for row in variant_statistics["by_region"]}
    variant_cards = _metric_cards(
        (
            ("QC non-reference variants", variant_counts["qc_passing_non_reference_count"]),
            ("Promoter variants", variant_region_counts.get("promoter", 0)),
            ("Gene-body variants", variant_region_counts.get("gene_body", 0)),
            ("Named variant IDs", variant_counts["named_variant_count"]),
        )
    )
    methylation_cards = _metric_cards(
        (
            ("Valid beta values", methylation_counts["valid_beta_count"]),
            ("Gene-named probes", methylation_counts["gene_named_probe_count"]),
            ("Promoter probes", methylation_counts["promoter_probe_count"]),
            ("Gene-body probes", methylation_counts["gene_body_probe_count"]),
        )
    )
    extremes = [
        {"extreme": "lowest", **row} for row in methylation_statistics["extremes"]["lowest"]
    ] + [{"extreme": "highest", **row} for row in methylation_statistics["extremes"]["highest"]]
    stats_panel = f"""
      <section id="panel-statistics" class="panel" role="tabpanel" aria-labelledby="tab-statistics" hidden>
        <h2>Single-person descriptive statistics</h2>
        <p>Status: <strong>{html.escape(statistics['status'])}</strong>. Every value below is calculated only from the variant and methylation rows in this person's gene analysis.</p>
        <h3>Variants</h3>{variant_cards}
        <h4>Region distribution</h4>{_table(variant_statistics['by_region'])}
        <h4>Variant and genotype distributions</h4>{_table(variant_statistics['by_type'])}{_table(variant_statistics['genotypes'])}
        <h4>A/C/G/T composition</h4>{_table(variant_statistics['reference_bases'])}{_table(variant_statistics['alternate_bases'])}{_table(variant_statistics['dosage_weighted_alternate_copies'])}
        <h4>Substitutions and density</h4>{_table(variant_statistics['substitutions'])}{_table(variant_statistics['density'])}
        <h4>Variant quality summaries</h4>{_table(variant_statistics['quality'])}
        <h3>Methylation</h3>{methylation_cards}
        <h4>Beta summaries by subset</h4>{_table(methylation_statistics['subsets'])}
        <h4>Genomic and manifest annotations</h4>{_table(methylation_statistics['by_region'])}{_table(methylation_statistics['by_refgene_group'])}{_table(methylation_statistics['by_cpg_island_relation'])}
        <h4>Beta distribution and measurement quality</h4>{_table(methylation_statistics['beta_histogram'])}{_table(methylation_statistics['quality'])}
        <h4>Annotation coverage and probe extremes</h4>{_table(methylation_statistics['annotation_coverage'])}{_table(extremes)}
        <details><summary><strong>Raw field coverage</strong></summary><h4>Variant fields</h4>{_table(variant_statistics['raw_field_coverage'])}<h4>Methylation fields</h4>{_table(methylation_statistics['raw_field_coverage'])}</details>
      </section>
    """
    literature = sections["literature"]
    literature_findings = literature.get("findings") or literature.get("detailed_records") or []
    ranking_policy = literature.get("ranking_policy") or {
        "disclaimer": "Evidence priority does not establish causality or clinical validity.",
        "tie_breakers": [],
    }
    literature_rows = [
        {
            "rank": record.get("rank"),
            "evidence_priority": record.get("priority_label"),
            "finding": record.get("finding") or "Citation metadata only; no finding was synthesized.",
            "paper": record.get("paper") or record.get("title"),
            "phenotype": record.get("phenotype"),
            "genotype": record.get("genotypes"),
            "variant": record.get("variant"),
            "sources": record.get("sources"),
            "citation_identifiers": " · ".join(
                filter(
                    None,
                    (
                        f"PMID {record['pmid']}" if record.get("pmid") else "",
                        _text(record.get("pmcid")),
                        f"DOI {record['doi']}" if record.get("doi") else "",
                    ),
                )
            ),
            "url": record.get("url"),
        }
        for record in literature_findings
    ]
    literature_limitations = "".join(
        f"<li>{html.escape(_text(item))}</li>" for item in literature.get("limitations", [])
    )
    literature_panel = f"""
      <section id="panel-literature" class="panel" role="tabpanel" aria-labelledby="tab-literature" hidden>
        <h2>Scientific Literature</h2>
        <p>All gene-related evidence available to this run is shown below in ranked order. This is not an exhaustive internet-wide search.</p>
        {_metric_cards((('Distinct findings', literature.get('finding_count', len(literature_findings))), ('Publications', literature.get('publication_count', literature.get('candidate_count', 0))), ('Sources', literature.get('source_count', 'not recorded')), ('Preprint findings', len(literature.get('preprints', [])))))}
        <p class="notice">{html.escape(_text(ranking_policy.get('disclaimer')))} Clinically interesting literature remains separate from authoritative Medical Information.</p>
        {_table(literature_rows, limit=None)}
        <details><summary><strong>Coverage and ranking details</strong></summary><ul>{literature_limitations}</ul><p><strong>Ranking tie-breakers:</strong> {html.escape('; '.join(ranking_policy.get('tie_breakers', [])))}.</p></details>
      </section>
    """
    medical = sections["medical"]
    medical_counts = medical.get("counts") or {}
    gene_medical = medical.get("gene_context") or {}
    investigated_conditions = medical.get("investigated_conditions") or []
    cohort_studies = medical.get("cohort_studies") or []
    variant_context = medical.get("variant_context") or []
    pathology_context = medical.get("pathology_context") or []
    gene_overview_parts = []
    for label, field in (
        ("Gene function", "gene_summary"),
        ("Clinical research context", "clinical_context"),
        ("Methylation interpretation", "methylation_interpretation"),
    ):
        value = _text(gene_medical.get(field))
        if value:
            gene_overview_parts.append(f"<h4>{html.escape(label)}</h4><p>{html.escape(value)}</p>")
    for label, field in (
        ("Condition-research overview", "condition_research_overview"),
        ("Methylation-condition research", "methylation_condition_research"),
        ("Variant-effect overview", "variant_effect_overview"),
        ("Methylation mechanisms", "methylation_effects"),
    ):
        values = gene_medical.get(field) or []
        if values:
            items = "".join(f"<li>{html.escape(_text(value))}</li>" for value in values)
            gene_overview_parts.append(f"<h4>{html.escape(label)}</h4><ul>{items}</ul>")
    condition_rows = [
        {
            "match": "observed match" if item.get("observed_match") else "queried-gene context",
            "match_basis": item.get("match_basis"),
            "condition_or_topic": item.get("condition_or_topic"),
            "variants": item.get("variants"),
            "evidence_summaries": [
                link.get("summary") for link in item.get("evidence_links", []) if link.get("summary")
            ],
            "sources": item.get("sources"),
            "authority_status": item.get("authority_status"),
        }
        for item in investigated_conditions
    ]
    cohort_rows = []
    for item in cohort_studies:
        unavailable = item.get("unavailable_parameters") or {}
        cohort_rows.append(
            {
                "match": "observed match" if item.get("observed_match") else "queried-gene context",
                "paper": item.get("paper"),
                "identifiers": {
                    key: item.get(key) for key in ("pmid", "pmcid", "doi") if item.get(key)
                },
                "url": item.get("url"),
                "variant": item.get("variant"),
                "phenotype_cohort_or_model": item.get("cohort_or_model_context")
                or unavailable.get("cohort_or_model_context"),
                "study_designs": item.get("study_designs"),
                "sample_parameters": item.get("sample_size_mentions")
                or unavailable.get("sample_size_mentions"),
                "genotype_or_comparators": item.get("genotype_or_comparator_groups")
                or unavailable.get("genotype_or_comparator_groups"),
                "tissue_or_model": item.get("tissue_or_model") or unavailable.get("tissue_or_model"),
                "reported_effect": item.get("reported_effect_measure_mentions")
                or unavailable.get("reported_effect_measure_mentions"),
                "finding": item.get("finding"),
            }
        )
    variant_rows = [
        {
            "match": "observed match" if item.get("observed_match") else "queried-gene context",
            "match_basis": item.get("match_basis"),
            "variant": item.get("display_name") or item.get("variant"),
            "observed_variants": item.get("observed_variants"),
            "observed_methylation_probes": item.get("observed_methylation_probes"),
            "scope": item.get("interpretation_scope"),
            "associated_conditions": item.get("associated_conditions"),
            "clinical_research_interpretation": item.get("clinical_interpretation"),
            "functional_effects": item.get("functional_effects"),
            "methylation_context": item.get("methylation_interpretation"),
            "clinical_parameters": item.get("clinical_parameter_summary")
            or "Not available in the current record.",
        }
        for item in variant_context
    ]
    pathology_rows = [
        {
            "match": "observed match" if item.get("observed_match") else "queried-gene context",
            "context_type": item.get("context_type"),
            "variant": item.get("variant"),
            "summary": item.get("summary"),
            "source": item.get("source_key"),
            "authority_status": item.get("authority_status"),
        }
        for item in pathology_context
    ]
    authoritative_records = medical.get("records") or []
    authoritative_empty = "" if authoritative_records else (
        '<p class="notice">No authoritative guideline, regulatory label, expert-panel assertion, or reviewed '
        "classification was available to this run. The research context remains visible but is not promoted "
        "into medical guidance.</p>"
    )
    medical_panel = f"""
      <section id="panel-medical" class="panel" role="tabpanel" aria-labelledby="tab-medical" hidden>
        <h2>Medical Information</h2>
        <p>Queried-gene research context is ranked with this person's observed variant or methylation matches first. Research associations remain separate from authoritative assertions.</p>
        {_metric_cards((('Investigated conditions', medical_counts.get('investigated_condition_count', len(investigated_conditions))), ('Cohort studies', medical_counts.get('cohort_study_count', len(cohort_studies))), ('Variant contexts', medical_counts.get('variant_context_count', len(variant_context))), ('Observed-context matches', medical_counts.get('observed_match_count', 0)), ('Established assertions', medical_counts.get('authoritative_assertion_count', len(authoritative_records))))) }
        <p class="notice">{html.escape(_text(medical.get('disclaimer')))}</p>
        <h3>1. Gene-level medical and pathology overview</h3>{''.join(gene_overview_parts) or '<p class="empty">No gene-level context was available.</p>'}
        <h3>2. Investigated conditions and phenotypes</h3>{_table(condition_rows, limit=None)}
        <h3>3. Cohort and study parameters</h3>{_table(cohort_rows, limit=None)}
        <h3>4. Variant and methylation context</h3>{_table(variant_rows, limit=None)}
        <details><summary><strong>Pathology and mechanism context ({len(pathology_rows)})</strong></summary>{_table(pathology_rows, limit=None)}</details>
        <h3>5. Established medical evidence and source coverage</h3>
        <p><strong>Authority assessment:</strong> {html.escape(_text(medical.get('status')).replace('_', ' '))}. <strong>Inclusion policy:</strong> {html.escape(_text(medical.get('inclusion_policy')))}.</p>
        {authoritative_empty}{_table(authoritative_records, limit=None)}
        <h4>Checked medical sources</h4>{_table(medical.get('checked_sources') or [], limit=None)}
        <h4>Medical sources that could not be assessed</h4>{_table(medical.get('failed_sources') or [], limit=None)}
        <details><summary><strong>Medical sources not assessed</strong></summary>{_table(medical.get('not_assessed_sources') or [], limit=None)}</details>
      </section>
    """
    interactions = sections["interactions"]
    graph = interactions["initial_graph"]
    interaction_rows = [
        {
            "hop": edge.get("hop"),
            "source_gene": edge.get("source_gene"),
            "target_gene": edge.get("target_gene"),
            "relationship": edge.get("edge_type"),
            "source": edge.get("source_key"),
            "native_score": edge.get("native_score"),
            "score_label": edge.get("native_score_label"),
            "evidence_summary": edge.get("summary"),
            "scope": edge.get("association_scope"),
            "url": edge.get("url"),
        }
        for edge in graph["edges"]
    ]
    interactions_panel = f"""
      <section id="panel-interactions" class="panel" role="tabpanel" aria-labelledby="tab-interactions" hidden>
        <h2>Interactions</h2><p>Direct source-backed partners are shown first. Gene-to-gene expansion supports three hops and at most 150 nodes.</p>
        {_metric_cards((('Direct partners', interactions.get('direct_partner_count', graph.get('hop_counts', {}).get(1, 0))), ('Visible nodes', graph['node_count']), ('Visible edges', graph['edge_count']), ('Interaction sources', interactions.get('source_count', len(interactions.get('sources', []))))))}
        <p class="notice">{html.escape(_text(interactions.get('association_policy') or 'Functional associations do not necessarily represent direct physical binding, directionality, causality, or a person-specific effect.'))}</p>
        <p><strong>Sources:</strong> {html.escape(', '.join(interactions.get('sources', [])) or 'none assessed')}.</p>
        {_table(interaction_rows, limit=None)}
        <details><summary><strong>Interaction source coverage</strong></summary>{_table(interactions.get('source_statuses') or [], limit=None)}{_table(interactions.get('failed_sources') or [], limit=None)}</details>
      </section>
    """
    predictions = sections["predictions"]
    prediction_counts = predictions.get("counts") or {}
    source_native_rows = [
        {
            "variant": item.get("variant"),
            "rsID": item.get("rsid"),
            "predictor": item.get("predictor"),
            "native value": item.get("native_value_original"),
            "native label": item.get("native_label"),
            "source": item.get("source") or item.get("source_key"),
            "release": item.get("source_release"),
            "match basis": item.get("match_basis"),
            "URL": item.get("url"),
        }
        for item in predictions.get("source_native_annotations", [])
    ]
    prediction_rows = [
        {
            "variant": item.get("variant"),
            "output type": item.get("output_type") or item.get("output_name"),
            "raw score": item.get("raw_score"),
            "quantile score": item.get("quantile_score"),
            "scorer": item.get("scorer"),
            "track": item.get("track_name"),
            "gene": item.get("gene_name") or item.get("gene_id"),
            "ontology": item.get("ontology_curie"),
            "biosample": item.get("biosample_name"),
            "model interval": item.get("model_interval_0_based_half_open"),
        }
        for item in predictions.get("predictions", [])
    ]
    predictions_panel = f"""
      <section id="panel-predictions" class="panel" role="tabpanel" aria-labelledby="tab-predictions" hidden>
        <h2>Predictions</h2><p>Exact-allele source-native annotations and explicitly requested molecular model runs are reported independently. No consensus or clinical score is calculated.</p>
        {_metric_cards((('Observed alleles', prediction_counts.get('observed_allele_count', 0)), ('Source-native annotations', prediction_counts.get('source_native_annotation_count', 0)), ('Model runs', prediction_counts.get('model_run_count', 0)), ('Completed model scores', prediction_counts.get('completed_prediction_count', 0))))}
        <h3>Exact-allele source-native annotations</h3><p>Values retain the source predictor's native scale and direction and are not normalized into pathogenicity or phenotype conclusions.</p>
        {_table(source_native_rows, limit=None)}
        <h3>Opt-in model jobs</h3>{_table(predictions.get('model_runs') or [], limit=None)}
        <h3>Completed molecular model scores</h3>{_table(prediction_rows, limit=None)}
        <details><summary><strong>Model applicability and metadata-only assessments</strong></summary>{_table(predictions.get('model_input_assessments') or [], limit=None)}{_table(predictions.get('model_metadata_assessments') or [], limit=None)}{_table(predictions.get('available_models') or [], limit=None)}</details>
        <details><summary><strong>Blockers, failures, notices, and limitations</strong></summary>{_table(predictions.get('failures') or [], limit=None)}{_table(predictions.get('notices') or [], limit=None)}<pre>{html.escape(json.dumps({'blockers': predictions.get('blockers', []), 'limitations': predictions.get('limitations', [])}, indent=2, ensure_ascii=False, default=str))}</pre></details>
        <p class="notice">Molecular research predictions are not statistical significance measurements, diagnoses, causal proof, disease-risk estimates, pathogenicity classifications, or authoritative medical evidence. Empty provider output does not establish no biological effect.</p>
      </section>
    """
    details_json = html.escape(json.dumps(report["run_details"], indent=2, ensure_ascii=False, default=str))
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(run.get('gene') or 'Gene')} evidence-first report</title>
<style>
:root{{--bg:#f8e7e6;--bg-2:#f2d3d6;--panel:rgba(255,247,247,.90);--panel-strong:rgba(255,252,252,.96);--ink:#2a1118;--muted:#6b4b56;--accent:#a1143d;--accent-2:#d11f4f;--line:rgba(118,19,45,.14);--warn:rgba(209,31,79,.10);--shadow:0 24px 70px rgba(90,12,36,.16)}}
*{{box-sizing:border-box}}body{{margin:0;background:radial-gradient(circle at 0 0,rgba(209,31,79,.30),transparent 28rem),radial-gradient(circle at 100% 0,rgba(124,13,45,.24),transparent 26rem),linear-gradient(180deg,var(--bg),var(--bg-2));color:var(--ink);font:15px/1.5 system-ui,sans-serif}}main{{max-width:1500px;margin:auto;padding:28px}}
header,.panel,details{{background:var(--panel);border:1px solid var(--line);border-radius:16px;padding:22px;margin-bottom:18px}}h1{{margin:.1em 0}}.muted{{color:var(--muted)}}
.tabs{{display:flex;gap:8px;overflow:auto;margin:18px 0}}.tab{{border:1px solid var(--line);background:var(--panel-strong);color:var(--muted);padding:10px 14px;border-radius:999px;white-space:nowrap;cursor:pointer}}.tab[aria-selected=true]{{background:linear-gradient(135deg,var(--accent),var(--accent-2));color:#fff9fa}}
.metrics{{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px;margin:15px 0}}.metrics article{{border:1px solid var(--line);border-radius:12px;padding:14px}}.metrics span{{display:block;color:var(--muted)}}.metrics strong{{font-size:1.5rem}}
.notice{{background:var(--warn);border-radius:10px;padding:12px}}.empty{{color:var(--muted);font-style:italic}}.table-shell{{overflow:auto;border:1px solid var(--line);border-radius:12px;margin-bottom:20px}}table{{border-collapse:collapse;width:100%;font-size:.88rem}}th,td{{text-align:left;padding:9px;border-bottom:1px solid var(--line);vertical-align:top;max-width:360px;overflow-wrap:anywhere}}th{{background:rgba(209,31,79,.10);position:sticky;top:0}}pre{{white-space:pre-wrap;max-height:36rem;overflow:auto}}
</style></head><body><main>
<header><p class="muted">NophiGene Version 2 · schema {REPORT_SCHEMA_VERSION}</p><h1>{html.escape(run.get('gene') or 'Gene')} evidence-first report</h1><p>{html.escape(run.get('genome_build') or 'build not declared')} · {html.escape(run.get('region') or 'region not declared')}</p></header>
<nav class="tabs" role="tablist" aria-label="Result sections">{tabs}</nav>
{summary_panel}{objective_panel}{stats_panel}{literature_panel}{medical_panel}{interactions_panel}{predictions_panel}
<details><summary><strong>Run Details</strong></summary><p>Technical provenance, failures, URLs, checksums, and raw interpretation are kept out of the primary result.</p><pre>{details_json}</pre></details>
</main><script>
document.querySelectorAll('[data-tab]').forEach(button=>button.addEventListener('click',()=>{{
 document.querySelectorAll('[data-tab]').forEach(item=>item.setAttribute('aria-selected','false'));
 document.querySelectorAll('.panel').forEach(panel=>{{panel.hidden=true;panel.classList.remove('active')}});
 button.setAttribute('aria-selected','true');const panel=document.getElementById('panel-'+button.dataset.tab);panel.hidden=false;panel.classList.add('active');button.focus();
}}));
</script></body></html>"""
