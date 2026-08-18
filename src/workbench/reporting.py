"""Canonical schema-3 report builder and compact evidence-first HTML renderer."""

from __future__ import annotations

import html
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .evidence import deduplicate_literature, medical_gate, source_coverage
from .interactions import build_interaction_graph
from .model_registry import inspect_model_inputs, list_model_manifests
from .pgx import resolve_pgx_diplotype
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


def _source_records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for container in (
        payload,
        payload.get("knowledge_base", {}),
        payload.get("dynamic_knowledge_base", {}),
    ):
        if not isinstance(container, dict):
            continue
        for key in ("source_records", "records", "literature_records"):
            value = container.get(key)
            if isinstance(value, list):
                records.extend(dict(item) for item in value if isinstance(item, dict))
    dynamic = payload.get("dynamic_knowledge_base", {})
    if isinstance(dynamic, dict):
        local_articles = dynamic.get("local_article_evidence", {})
        if isinstance(local_articles, dict) and isinstance(local_articles.get("records"), list):
            records.extend(dict(item) for item in local_articles["records"] if isinstance(item, dict))
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
        if record.get("title") or record.get("pmid") or record.get("pmcid") or record.get("doi")
    ]
    literature = deduplicate_literature(literature_candidates, candidate_limit=200)
    literature_categories: dict[str, list[dict[str, Any]]] = {
        "primary_studies": [],
        "reviews": [],
        "preprints": [],
        "case_reports": [],
        "database_assertions": [],
    }
    for record in literature:
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
    if medical:
        medical_status = "assessed"
    elif assessed_states:
        medical_status = "assessed_absence"
    elif failed_states:
        medical_status = "source_failed"
    else:
        medical_status = "not_assessed"
    return (
        {
            "status": "assessed" if literature_candidates else "not_assessed",
            "candidate_limit": 200,
            "detailed_review_limit": 50,
            "candidate_count": len(literature),
            "detailed_records": [record for record in literature[:50] if not record.get("preprint")],
            "candidate_index": literature,
            "categories": literature_categories,
            "preprints": [record for record in literature[:50] if record.get("preprint")],
            "preprints_separated": True,
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
    edges: list[dict[str, Any]] = []
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
        edges.append(
            {
                "source_gene": source_gene,
                "target_gene": partner,
                "edge_type": edge_type,
                "directed": bool(record.get("directed", False)),
                "source_key": _text(record.get("source_key") or "unknown"),
                "source_record_id": _text(record.get("record_id")),
                "native_score": _number(record.get("native_score") or record.get("score")),
                "native_score_label": _text(record.get("score_type")),
                "tissue": _text(record.get("tissue")),
                "evidence": _jsonable(record),
            }
        )
    return edges


def _prediction_section(payload: dict[str, Any], context: dict[str, Any]) -> dict[str, Any]:
    interpretation = payload.get("interpretation", {}) if isinstance(payload.get("interpretation"), dict) else {}
    legacy_assessments = interpretation.get("model_assessments", [])
    requested = payload.get("requested_models") or []
    if not requested and isinstance(legacy_assessments, list):
        requested = [item.get("model_id") for item in legacy_assessments if isinstance(item, dict) and item.get("model_id")]
    requested_ids = []
    for item in requested:
        if isinstance(item, str):
            requested_ids.append(item)
        elif isinstance(item, dict) and item.get("model_id"):
            requested_ids.append(str(item["model_id"]))
    input_availability = dict(context.get("model_inputs") or {})
    input_availability.setdefault("genome_build", str(payload.get("genome_build") or ""))
    assessments = [inspect_model_inputs(model_id, input_availability) for model_id in requested_ids]
    pgx = None
    if isinstance(payload.get("pgx_definition"), dict):
        pgx = resolve_pgx_diplotype(payload["pgx_definition"], list(payload.get("pgx_calls") or []))
    return {
        "status": "assessed" if requested_ids or legacy_assessments else "not_requested",
        "requested_models": requested_ids,
        "assessments": assessments,
        "legacy_metadata_assessments": _jsonable(legacy_assessments),
        "predictions": [],
        "pgx_diplotype": pgx,
        "consensus": None,
        "policy": "independent model outputs; no consensus or app-defined low/medium/high labels",
        "available_models": list_model_manifests(),
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
    records = _source_records(payload)
    statuses = _provider_statuses(payload)
    coverage = source_coverage(statuses)
    literature, medical = _literature_and_medical(records, statuses)
    context = dict(payload.get("sample_context") or {})
    if not context:
        context = dict(
            payload.get("interpretation", {})
            .get("interpretation_context", {})
            .get("sample_context", {})
        )
    interaction_evidence_edges = _interaction_edges(records, gene)
    graph = build_interaction_graph(gene, interaction_evidence_edges, max_hops=1, node_cap=150)
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
                "status": "assessed" if graph["edge_count"] else "not_assessed",
                "initial_graph": graph,
                "max_hops": 3,
                "node_cap": 150,
                "hop_semantics": "gene_to_gene",
            },
            "predictions": _prediction_section(payload, context),
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


def _table(rows: Any, *, limit: int = DEFAULT_PAGE_SIZE) -> str:
    if not isinstance(rows, list) or not rows:
        return '<p class="empty">No eligible records were available for this section.</p>'
    normalized = [row if isinstance(row, dict) else {"value": row} for row in rows[:limit]]
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
    literature_panel = f"""
      <section id="panel-literature" class="panel" role="tabpanel" aria-labelledby="tab-literature" hidden>
        <h2>Scientific Literature</h2><p>{literature['candidate_count']} deduplicated candidates; detailed review is capped at {literature['detailed_review_limit']}.</p>
        {_table(literature['detailed_records'])}
      </section>
    """
    medical = sections["medical"]
    medical_panel = f"""
      <section id="panel-medical" class="panel" role="tabpanel" aria-labelledby="tab-medical" hidden>
        <h2>Medical Information</h2><p class="notice">{html.escape(medical['disclaimer'])}</p>{_table(medical['records'])}
      </section>
    """
    interactions = sections["interactions"]
    graph = interactions["initial_graph"]
    interactions_panel = f"""
      <section id="panel-interactions" class="panel" role="tabpanel" aria-labelledby="tab-interactions" hidden>
        <h2>Interactions</h2><p>Direct partners are shown first. Gene-to-gene expansion supports three hops and at most 150 nodes.</p>
        <div class="metrics"><article><span>Nodes</span><strong>{graph['node_count']}</strong></article><article><span>Edges</span><strong>{graph['edge_count']}</strong></article></div>
        {_table(graph['edges'])}
      </section>
    """
    predictions = sections["predictions"]
    predictions_panel = f"""
      <section id="panel-predictions" class="panel" role="tabpanel" aria-labelledby="tab-predictions" hidden>
        <h2>Predictions</h2><p>Each model is reported independently. Unsupported inputs produce eligibility blockers, never a substitute score.</p>
        {_table(predictions['assessments'])}
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
