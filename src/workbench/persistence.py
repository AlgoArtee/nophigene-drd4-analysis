"""Idempotent persistence of canonical schema-3 reports into normalized tables."""

from __future__ import annotations

import hashlib
import json
from typing import Any

from sqlalchemy import select
from sqlalchemy.orm import Session

from .audit import append_audit_event
from .models import (
    EvidenceRecord,
    EvidenceSnapshot,
    Gene,
    GenomicLocus,
    InteractionEdge,
    MedicalAssertion,
    MethylationMeasurement,
    Run,
    StatisticalResult,
    VariantCall,
)


def _first(row: dict[str, Any], *keys: str, default: Any = "") -> Any:
    for key in keys:
        if row.get(key) not in (None, ""):
            return row[key]
    return default


def _float(value: Any) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _int(value: Any) -> int | None:
    try:
        return int(float(value)) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    ).hexdigest()


def _gene(session: Session, symbol: str) -> Gene | None:
    if not symbol:
        return None
    gene = session.scalar(select(Gene).where(Gene.symbol == symbol))
    if gene is None:
        gene = Gene(symbol=symbol)
        session.add(gene)
        session.flush()
    return gene


def _variant_locus(
    session: Session,
    *,
    row: dict[str, Any],
    assembly: str,
    gene_id: str | None,
) -> GenomicLocus | None:
    chromosome = str(_first(row, "CHROM", "chrom", "chromosome")).removeprefix("chr")
    position = _int(_first(row, "POS", "pos", "position"))
    ref = str(_first(row, "REF", "ref"))
    alt = str(_first(row, "ALT", "alt"))
    if not chromosome or position is None:
        return None
    end = position + max(1, len(ref)) - 1
    locus = session.scalar(
        select(GenomicLocus).where(
            GenomicLocus.assembly == assembly,
            GenomicLocus.chromosome == chromosome,
            GenomicLocus.start == position,
            GenomicLocus.end == end,
            GenomicLocus.ref == ref,
            GenomicLocus.alt == alt,
        )
    )
    if locus is None:
        locus = GenomicLocus(
            gene_id=gene_id,
            assembly=assembly,
            chromosome=chromosome,
            start=position,
            end=end,
            ref=ref,
            alt=alt,
            locus_type="variant",
        )
        session.add(locus)
        session.flush()
    return locus


def persist_canonical_report(session: Session, report: dict[str, Any]) -> dict[str, int]:
    """Persist observations/evidence once while keeping report generation reproducible."""
    run_id = str(report.get("run", {}).get("id") or "")
    run = session.get(Run, run_id) if run_id else None
    if run is None:
        return {"skipped": 1}
    counts = {
        "variants": 0,
        "methylation": 0,
        "statistics": 0,
        "evidence": 0,
        "medical": 0,
        "interactions": 0,
    }
    run.schema_version = "3.0"
    report_status = str(report.get("run", {}).get("status") or "")
    if run.status not in {"partial", "failed", "cancelled"} and report_status:
        run.status = report_status
    gene = _gene(session, str(report.get("run", {}).get("gene") or "").upper())
    gene_id = gene.id if gene else None
    assembly = str(report.get("run", {}).get("genome_build") or "")
    objective = report.get("sections", {}).get("objective_data", {})
    for row in objective.get("all_variants", []):
        if not isinstance(row, dict):
            continue
        locus = _variant_locus(session, row=row, assembly=assembly, gene_id=gene_id)
        if locus is None:
            continue
        genotype = str(_first(row, "GT", "gt_raw", "genotype"))
        exists = session.scalar(
            select(VariantCall.id).where(
                VariantCall.run_id == run_id,
                VariantCall.locus_id == locus.id,
                VariantCall.genotype == genotype,
            )
        )
        if exists:
            continue
        session.add(
            VariantCall(
                run_id=run_id,
                sample_id=run.sample_id,
                gene_id=gene_id,
                locus_id=locus.id,
                rsid=str(_first(row, "rsid", "ID", "id")),
                genotype=genotype,
                zygosity=str(row.get("zygosity") or ""),
                filter_status=str(_first(row, "FILTER", "filter", "filter_status")),
                genotype_quality=_float(_first(row, "GQ", "genotype_quality", default=None)),
                depth=_int(_first(row, "DP", "depth", default=None)),
                qc_pass=bool(row.get("qc_pass")),
                qc_reasons=list(row.get("qc_reasons") or []),
                consequence=dict(row.get("consequence") or {}),
                raw_fields=dict(row),
            )
        )
        counts["variants"] += 1
    for row in objective.get("all_methylation", []):
        if not isinstance(row, dict):
            continue
        probe_id = str(_first(row, "probe_id", "IlmnID", "Name"))
        if not probe_id or session.scalar(
            select(MethylationMeasurement.id).where(
                MethylationMeasurement.run_id == run_id,
                MethylationMeasurement.probe_id == probe_id,
            )
        ):
            continue
        session.add(
            MethylationMeasurement(
                run_id=run_id,
                sample_id=run.sample_id,
                gene_id=gene_id,
                probe_id=probe_id,
                beta_value=_float(_first(row, "beta_value", "beta", "Beta_value", default=None)),
                m_value=_float(_first(row, "m_value", "M_value", default=None)),
                detection_p=_float(_first(row, "detection_p", "detection_p_value", "Detection Pval", default=None)),
                bead_count=_int(_first(row, "bead_count", "NBeads", default=None)),
                normalization=str(row.get("normalization") or ""),
                manifest_version=str(row.get("manifest_version") or ""),
                qc_pass=bool(row.get("qc_pass")),
                qc_reasons=list(row.get("qc_reasons") or []),
                raw_fields=dict(row),
            )
        )
        counts["methylation"] += 1
    for row in report.get("sections", {}).get("statistics", {}).get("records", []):
        if not isinstance(row, dict):
            continue
        entity_key = str(row.get("entity_key") or "")
        family = str(row.get("family") or "unspecified")
        if session.scalar(
            select(StatisticalResult.id).where(
                StatisticalResult.run_id == run_id,
                StatisticalResult.entity_key == entity_key,
                StatisticalResult.family == family,
            )
        ):
            continue
        session.add(
            StatisticalResult(
                run_id=run_id,
                entity_type=str(row.get("entity_type") or "measurement"),
                entity_key=entity_key,
                family=family,
                method=str(row.get("method") or ""),
                effect_size=_float(row.get("effect_size")),
                percentile=_float(row.get("percentile")),
                raw_p=_float(row.get("raw_p")),
                q_value=_float(row.get("q_value")),
                status=str(row.get("status") or "not_assessed"),
                limitations=list(row.get("limitations") or []),
                details=dict(row),
            )
        )
        counts["statistics"] += 1

    details = report.get("run_details", {})
    records = [item for item in details.get("normalized_evidence_records", []) if isinstance(item, dict)]
    statuses = [item for item in details.get("provider_statuses", []) if isinstance(item, dict)]
    edges = [item for item in details.get("interaction_evidence_edges", []) if isinstance(item, dict)]
    snapshot_material = {"records": records, "statuses": statuses, "run": run_id}
    checksum = _digest(snapshot_material)
    snapshot = session.scalar(select(EvidenceSnapshot).where(EvidenceSnapshot.checksum_sha256 == checksum))
    if snapshot is None:
        snapshot = EvidenceSnapshot(
            id=checksum,
            run_id=run_id,
            checksum_sha256=checksum,
            query_context={"gene": gene.symbol if gene else "", "build": assembly, "provider_statuses": statuses},
        )
        session.add(snapshot)
        session.flush()
    evidence_lookup: dict[tuple[str, str], EvidenceRecord] = {}
    for record in records:
        source = str(record.get("source_key") or record.get("source") or "unknown")
        record_id = str(record.get("record_id") or record.get("id") or _digest(record)[:32])
        existing = session.scalar(
            select(EvidenceRecord).where(
                EvidenceRecord.snapshot_id == snapshot.id,
                EvidenceRecord.source_key == source,
                EvidenceRecord.record_id == record_id,
            )
        )
        if existing is None:
            citations = [
                {key: record[key] for key in ("pmid", "pmcid", "doi") if record.get(key)}
            ]
            existing = EvidenceRecord(
                snapshot_id=snapshot.id,
                source_key=source,
                source_release=str(record.get("source_release") or record.get("release") or ""),
                record_id=record_id,
                evidence_type=str(record.get("evidence_type") or record.get("record_type") or "source_record"),
                entity_type=str(record.get("entity_type") or "gene"),
                entity_key=str(record.get("entity_key") or (gene.symbol if gene else "")),
                assertion=str(record.get("assertion") or record.get("summary") or record.get("title") or ""),
                effect=str(record.get("effect") or record.get("finding") or ""),
                evidence_level=str(record.get("evidence_level") or record.get("review_status") or ""),
                context=dict(record.get("context") or {}),
                citations=[item for item in citations if item],
                license_status=str(record.get("license_status") or "not_recorded"),
                raw_checksum_sha256=str(record.get("raw_checksum_sha256") or record.get("raw_response_hash") or _digest(record)),
                status=str(record.get("status") or "assessed"),
            )
            session.add(existing)
            session.flush()
            counts["evidence"] += 1
        evidence_lookup[(source, record_id)] = existing
    for record in report.get("sections", {}).get("medical", {}).get("records", []):
        if not isinstance(record, dict):
            continue
        source = str(record.get("source_key") or record.get("source") or "unknown")
        record_id = str(record.get("record_id") or record.get("id") or _digest(record)[:32])
        evidence = evidence_lookup.get((source, record_id))
        if evidence is None or session.scalar(
            select(MedicalAssertion.id).where(MedicalAssertion.evidence_record_id == evidence.id)
        ):
            continue
        session.add(
            MedicalAssertion(
                evidence_record_id=evidence.id,
                authority=str(record.get("authority") or source),
                jurisdiction=str(record.get("jurisdiction") or ""),
                assertion_type=str(record.get("assertion_type") or record.get("evidence_type") or ""),
                review_status=str(record.get("review_status") or ""),
                effective_date=str(record.get("effective_date") or ""),
                population=str(record.get("population") or ""),
                applicability=dict(record.get("applicability") or {}),
                conflict_group=str(record.get("conflict_group") or ""),
            )
        )
        counts["medical"] += 1
    for edge in edges:
        key = {
            "source_gene": str(edge.get("source_gene") or "").upper(),
            "target_gene": str(edge.get("target_gene") or "").upper(),
            "edge_type": str(edge.get("edge_type") or "functional_association"),
            "source_key": str(edge.get("source_key") or "unknown"),
            "source_record_id": str(edge.get("source_record_id") or ""),
        }
        if not key["source_gene"] or not key["target_gene"] or session.scalar(
            select(InteractionEdge.id).where(
                InteractionEdge.snapshot_id == snapshot.id,
                InteractionEdge.source_gene == key["source_gene"],
                InteractionEdge.target_gene == key["target_gene"],
                InteractionEdge.edge_type == key["edge_type"],
                InteractionEdge.source_key == key["source_key"],
                InteractionEdge.source_record_id == key["source_record_id"],
            )
        ):
            continue
        session.add(
            InteractionEdge(
                snapshot_id=snapshot.id,
                **key,
                directed=bool(edge.get("directed")),
                native_score=_float(edge.get("native_score")),
                native_score_label=str(edge.get("native_score_label") or ""),
                tissue=str(edge.get("tissue") or ""),
                evidence_count=int(edge.get("corroborating_source_count") or 1),
                details=dict(edge),
            )
        )
        counts["interactions"] += 1
    if any(counts.values()):
        append_audit_event(
            session,
            "canonical_report_persisted",
            entity_type="run",
            entity_id=run_id,
            payload={"schema_version": "3.0", "inserted": counts},
        )
    return counts
