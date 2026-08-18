"""Normalized SQLAlchemy schema for NophiGene Version 2.

The schema deliberately separates observations, statistical results, curated
evidence, medical assertions, interactions, and model predictions.  JSON
columns contain source-specific details, never the canonical identity fields
used for joins and deduplication.
"""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from typing import Any

from sqlalchemy import (
    Boolean,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    JSON,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def new_id() -> str:
    return uuid.uuid4().hex


class Base(DeclarativeBase):
    pass


class TimestampMixin:
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utc_now, nullable=False)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utc_now, onupdate=utc_now, nullable=False
    )


class Gene(Base, TimestampMixin):
    __tablename__ = "genes"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    symbol: Mapped[str] = mapped_column(String(64), unique=True, index=True)
    hgnc_id: Mapped[str | None] = mapped_column(String(32), unique=True)
    ensembl_id: Mapped[str | None] = mapped_column(String(32), index=True)
    ncbi_gene_id: Mapped[str | None] = mapped_column(String(32), index=True)
    name: Mapped[str] = mapped_column(Text, default="")
    aliases: Mapped[list[str]] = mapped_column(JSON, default=list)


class Transcript(Base, TimestampMixin):
    __tablename__ = "transcripts"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    gene_id: Mapped[str] = mapped_column(ForeignKey("genes.id", ondelete="CASCADE"), index=True)
    accession: Mapped[str] = mapped_column(String(64), index=True)
    version: Mapped[str] = mapped_column(String(16), default="")
    assembly: Mapped[str] = mapped_column(String(16), index=True)
    is_mane_select: Mapped[bool] = mapped_column(Boolean, default=False)
    protein_accession: Mapped[str] = mapped_column(String(64), default="")
    metadata_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)

    gene: Mapped[Gene] = relationship()
    __table_args__ = (UniqueConstraint("accession", "version", "assembly"),)


class GenomicLocus(Base, TimestampMixin):
    __tablename__ = "genomic_loci"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    gene_id: Mapped[str | None] = mapped_column(ForeignKey("genes.id", ondelete="SET NULL"), index=True)
    assembly: Mapped[str] = mapped_column(String(16), index=True)
    chromosome: Mapped[str] = mapped_column(String(16), index=True)
    start: Mapped[int] = mapped_column(Integer)
    end: Mapped[int] = mapped_column(Integer)
    ref: Mapped[str] = mapped_column(Text, default="")
    alt: Mapped[str] = mapped_column(Text, default="")
    locus_type: Mapped[str] = mapped_column(String(32), default="variant")

    __table_args__ = (
        Index("ix_locus_identity", "assembly", "chromosome", "start", "end"),
    )


class CrossBuildMapping(Base, TimestampMixin):
    __tablename__ = "cross_build_mappings"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    source_locus_id: Mapped[str] = mapped_column(ForeignKey("genomic_loci.id", ondelete="CASCADE"))
    target_locus_id: Mapped[str | None] = mapped_column(ForeignKey("genomic_loci.id", ondelete="SET NULL"))
    chain_name: Mapped[str] = mapped_column(String(128))
    tool: Mapped[str] = mapped_column(String(64))
    tool_version: Mapped[str] = mapped_column(String(32), default="")
    status: Mapped[str] = mapped_column(String(32), index=True)
    ref_validated: Mapped[bool] = mapped_column(Boolean, default=False)
    ambiguous: Mapped[bool] = mapped_column(Boolean, default=False)
    details: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class Sample(Base, TimestampMixin):
    __tablename__ = "samples"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    pseudonym: Mapped[str] = mapped_column(String(128), unique=True, index=True)
    user_label: Mapped[str] = mapped_column(Text, default="")
    context: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    deleted_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))


class InputArtifact(Base, TimestampMixin):
    __tablename__ = "input_artifacts"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    sample_id: Mapped[str | None] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"), index=True)
    kind: Mapped[str] = mapped_column(String(32), index=True)
    source_path: Mapped[str] = mapped_column(Text, default="")
    managed_path: Mapped[str] = mapped_column(Text, default="")
    checksum_sha256: Mapped[str] = mapped_column(String(64), index=True)
    size_bytes: Mapped[int] = mapped_column(Integer, default=0)
    assembly: Mapped[str] = mapped_column(String(16), default="")
    metadata_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class Run(Base, TimestampMixin):
    __tablename__ = "runs"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    parent_run_id: Mapped[str | None] = mapped_column(ForeignKey("runs.id", ondelete="SET NULL"))
    sample_id: Mapped[str | None] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"), index=True)
    schema_version: Mapped[str] = mapped_column(String(16), default="3.0")
    status: Mapped[str] = mapped_column(String(32), index=True)
    stage: Mapped[str] = mapped_column(String(64), default="queued")
    progress_percent: Mapped[int] = mapped_column(Integer, default=0)
    genes: Mapped[list[str]] = mapped_column(JSON, default=list)
    context: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    configuration: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    started_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    finished_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    error: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class QcEvent(Base, TimestampMixin):
    __tablename__ = "qc_events"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    run_id: Mapped[str] = mapped_column(ForeignKey("runs.id", ondelete="CASCADE"), index=True)
    entity_type: Mapped[str] = mapped_column(String(32), index=True)
    entity_key: Mapped[str] = mapped_column(String(256), default="")
    check: Mapped[str] = mapped_column(String(128))
    status: Mapped[str] = mapped_column(String(32), index=True)
    value: Mapped[str] = mapped_column(Text, default="")
    threshold: Mapped[str] = mapped_column(Text, default="")
    reason: Mapped[str] = mapped_column(Text, default="")


class VariantCall(Base, TimestampMixin):
    __tablename__ = "variant_calls"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    run_id: Mapped[str] = mapped_column(ForeignKey("runs.id", ondelete="CASCADE"), index=True)
    sample_id: Mapped[str | None] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"), index=True)
    gene_id: Mapped[str | None] = mapped_column(ForeignKey("genes.id", ondelete="SET NULL"), index=True)
    locus_id: Mapped[str] = mapped_column(ForeignKey("genomic_loci.id", ondelete="RESTRICT"), index=True)
    rsid: Mapped[str] = mapped_column(String(64), default="", index=True)
    genotype: Mapped[str] = mapped_column(String(32), default="")
    zygosity: Mapped[str] = mapped_column(String(32), default="")
    filter_status: Mapped[str] = mapped_column(String(64), default="")
    genotype_quality: Mapped[float | None] = mapped_column(Float)
    depth: Mapped[int | None] = mapped_column(Integer)
    qc_pass: Mapped[bool] = mapped_column(Boolean, default=False, index=True)
    qc_reasons: Mapped[list[str]] = mapped_column(JSON, default=list)
    consequence: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    raw_fields: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)

    __table_args__ = (UniqueConstraint("run_id", "locus_id", "genotype"),)


class MethylationMeasurement(Base, TimestampMixin):
    __tablename__ = "methylation_measurements"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    run_id: Mapped[str] = mapped_column(ForeignKey("runs.id", ondelete="CASCADE"), index=True)
    sample_id: Mapped[str | None] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"), index=True)
    gene_id: Mapped[str | None] = mapped_column(ForeignKey("genes.id", ondelete="SET NULL"), index=True)
    locus_id: Mapped[str | None] = mapped_column(ForeignKey("genomic_loci.id", ondelete="SET NULL"), index=True)
    probe_id: Mapped[str] = mapped_column(String(64), index=True)
    beta_value: Mapped[float | None] = mapped_column(Float)
    m_value: Mapped[float | None] = mapped_column(Float)
    detection_p: Mapped[float | None] = mapped_column(Float)
    bead_count: Mapped[int | None] = mapped_column(Integer)
    normalization: Mapped[str] = mapped_column(String(64), default="")
    manifest_version: Mapped[str] = mapped_column(String(64), default="")
    qc_pass: Mapped[bool] = mapped_column(Boolean, default=False, index=True)
    qc_reasons: Mapped[list[str]] = mapped_column(JSON, default=list)
    raw_fields: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)

    __table_args__ = (UniqueConstraint("run_id", "probe_id"),)


class ReferenceCohort(Base, TimestampMixin):
    __tablename__ = "reference_cohorts"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    name: Mapped[str] = mapped_column(Text)
    source_type: Mapped[str] = mapped_column(String(32))
    release: Mapped[str] = mapped_column(String(64), default="")
    sample_count: Mapped[int] = mapped_column(Integer)
    context: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    checksum_sha256: Mapped[str] = mapped_column(String(64), default="")


class StatisticalDataset(Base, TimestampMixin):
    """Immutable catalog entry for cohort-level statistical inputs."""

    __tablename__ = "statistical_datasets"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    name: Mapped[str] = mapped_column(Text)
    description: Mapped[str] = mapped_column(Text, default="")
    method: Mapped[str] = mapped_column(String(64), index=True, default="dandelion")
    phenotype: Mapped[str] = mapped_column(String(256), index=True)
    exposure_type: Mapped[str] = mapped_column(String(16), index=True)
    assembly: Mapped[str] = mapped_column(String(16), index=True)
    gene_namespace: Mapped[str] = mapped_column(String(64))
    storage_mode: Mapped[str] = mapped_column(String(32), default="registered_path")
    context: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    applicability: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    manifest_checksum_sha256: Mapped[str] = mapped_column(String(64), unique=True, index=True)


class StatisticalDatasetFile(Base, TimestampMixin):
    __tablename__ = "statistical_dataset_files"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    dataset_id: Mapped[str] = mapped_column(ForeignKey("statistical_datasets.id", ondelete="CASCADE"), index=True)
    role: Mapped[str] = mapped_column(String(64))
    relative_path: Mapped[str] = mapped_column(Text)
    managed_path: Mapped[str] = mapped_column(Text, default="")
    file_format: Mapped[str] = mapped_column(String(16))
    checksum_sha256: Mapped[str] = mapped_column(String(64), index=True)
    size_bytes: Mapped[int] = mapped_column(Integer)
    mapping: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    inspection: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)

    __table_args__ = (UniqueConstraint("dataset_id", "role"),)


class StatisticalAnalysisRun(Base, TimestampMixin):
    """Execution metadata for a cohort statistical method, separate from AI models."""

    __tablename__ = "statistical_analysis_runs"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    run_id: Mapped[str] = mapped_column(ForeignKey("runs.id", ondelete="CASCADE"), unique=True, index=True)
    dataset_id: Mapped[str] = mapped_column(ForeignKey("statistical_datasets.id", ondelete="RESTRICT"), index=True)
    method: Mapped[str] = mapped_column(String(64), index=True)
    method_version: Mapped[str] = mapped_column(String(64))
    worker_job_id: Mapped[str] = mapped_column(String(32), unique=True, index=True)
    status: Mapped[str] = mapped_column(String(32), index=True)
    parameters: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    resource_estimate: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    result_checksum_sha256: Mapped[str] = mapped_column(String(64), default="")
    error: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class StatisticalResult(Base, TimestampMixin):
    __tablename__ = "statistical_results"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    run_id: Mapped[str] = mapped_column(ForeignKey("runs.id", ondelete="CASCADE"), index=True)
    analysis_run_id: Mapped[str | None] = mapped_column(
        ForeignKey("statistical_analysis_runs.id", ondelete="CASCADE"), index=True
    )
    cohort_id: Mapped[str | None] = mapped_column(ForeignKey("reference_cohorts.id", ondelete="SET NULL"))
    entity_type: Mapped[str] = mapped_column(String(32), index=True)
    entity_key: Mapped[str] = mapped_column(String(256), index=True)
    family: Mapped[str] = mapped_column(String(64), index=True)
    method: Mapped[str] = mapped_column(String(128))
    effect_size: Mapped[float | None] = mapped_column(Float)
    percentile: Mapped[float | None] = mapped_column(Float)
    raw_p: Mapped[float | None] = mapped_column(Float)
    q_value: Mapped[float | None] = mapped_column(Float)
    status: Mapped[str] = mapped_column(String(32), index=True)
    limitations: Mapped[list[str]] = mapped_column(JSON, default=list)
    details: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class EvidenceSnapshot(Base, TimestampMixin):
    __tablename__ = "evidence_snapshots"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    run_id: Mapped[str | None] = mapped_column(ForeignKey("runs.id", ondelete="SET NULL"), index=True)
    generated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utc_now)
    checksum_sha256: Mapped[str] = mapped_column(String(64), unique=True)
    refresh_policy: Mapped[str] = mapped_column(String(64), default="explicit_user_triggered")
    query_context: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class EvidenceRecord(Base, TimestampMixin):
    __tablename__ = "evidence_records"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    snapshot_id: Mapped[str] = mapped_column(ForeignKey("evidence_snapshots.id", ondelete="CASCADE"), index=True)
    source_key: Mapped[str] = mapped_column(String(64), index=True)
    source_release: Mapped[str] = mapped_column(String(128), default="")
    record_id: Mapped[str] = mapped_column(String(256), default="")
    evidence_type: Mapped[str] = mapped_column(String(64), index=True)
    entity_type: Mapped[str] = mapped_column(String(32), index=True)
    entity_key: Mapped[str] = mapped_column(String(256), index=True)
    assertion: Mapped[str] = mapped_column(Text, default="")
    effect: Mapped[str] = mapped_column(Text, default="")
    evidence_level: Mapped[str] = mapped_column(String(64), default="")
    context: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    citations: Mapped[list[dict[str, Any]]] = mapped_column(JSON, default=list)
    license_status: Mapped[str] = mapped_column(String(64), default="")
    raw_checksum_sha256: Mapped[str] = mapped_column(String(64), default="")
    status: Mapped[str] = mapped_column(String(32), default="assessed")
    raw_artifact_id: Mapped[str | None] = mapped_column(ForeignKey("artifacts.id", ondelete="SET NULL"))

    __table_args__ = (UniqueConstraint("snapshot_id", "source_key", "record_id"),)


class MedicalAssertion(Base, TimestampMixin):
    __tablename__ = "medical_assertions"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    evidence_record_id: Mapped[str] = mapped_column(ForeignKey("evidence_records.id", ondelete="CASCADE"), index=True)
    authority: Mapped[str] = mapped_column(String(128))
    jurisdiction: Mapped[str] = mapped_column(String(64))
    assertion_type: Mapped[str] = mapped_column(String(64))
    review_status: Mapped[str] = mapped_column(String(128))
    effective_date: Mapped[str] = mapped_column(String(32), default="")
    population: Mapped[str] = mapped_column(Text, default="")
    applicability: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    conflict_group: Mapped[str] = mapped_column(String(128), default="")


class InteractionEdge(Base, TimestampMixin):
    __tablename__ = "interaction_edges"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    snapshot_id: Mapped[str] = mapped_column(ForeignKey("evidence_snapshots.id", ondelete="CASCADE"), index=True)
    source_gene: Mapped[str] = mapped_column(String(64), index=True)
    target_gene: Mapped[str] = mapped_column(String(64), index=True)
    edge_type: Mapped[str] = mapped_column(String(64), index=True)
    directed: Mapped[bool] = mapped_column(Boolean, default=False)
    source_key: Mapped[str] = mapped_column(String(64), index=True)
    source_record_id: Mapped[str] = mapped_column(String(256), default="")
    native_score: Mapped[float | None] = mapped_column(Float)
    native_score_label: Mapped[str] = mapped_column(String(64), default="")
    tissue: Mapped[str] = mapped_column(String(128), default="")
    evidence_count: Mapped[int] = mapped_column(Integer, default=1)
    details: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)

    __table_args__ = (
        UniqueConstraint("snapshot_id", "source_gene", "target_gene", "edge_type", "source_key", "source_record_id"),
    )


class ModelDefinition(Base, TimestampMixin):
    __tablename__ = "model_definitions"

    id: Mapped[str] = mapped_column(String(128), primary_key=True)
    name: Mapped[str] = mapped_column(String(128))
    version: Mapped[str] = mapped_column(String(64))
    wave: Mapped[str] = mapped_column(String(16), index=True)
    task: Mapped[str] = mapped_column(String(128))
    execution_mode: Mapped[str] = mapped_column(String(32))
    status: Mapped[str] = mapped_column(String(32), index=True)
    input_contract: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    output_contract: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    supported_builds: Mapped[list[str]] = mapped_column(JSON, default=list)
    resource_requirements: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    license: Mapped[str] = mapped_column(Text, default="")
    citation: Mapped[str] = mapped_column(Text, default="")
    limitations: Mapped[list[str]] = mapped_column(JSON, default=list)
    manifest_checksum_sha256: Mapped[str] = mapped_column(String(64), default="")


class ModelRun(Base, TimestampMixin):
    __tablename__ = "model_runs"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    run_id: Mapped[str] = mapped_column(ForeignKey("runs.id", ondelete="CASCADE"), index=True)
    model_id: Mapped[str] = mapped_column(ForeignKey("model_definitions.id", ondelete="RESTRICT"), index=True)
    status: Mapped[str] = mapped_column(String(32), index=True)
    blockers: Mapped[list[str]] = mapped_column(JSON, default=list)
    input_checksum_sha256: Mapped[str] = mapped_column(String(64), default="")
    container_digest: Mapped[str] = mapped_column(String(128), default="")
    started_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    finished_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    error: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class Prediction(Base, TimestampMixin):
    __tablename__ = "predictions"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    model_run_id: Mapped[str] = mapped_column(ForeignKey("model_runs.id", ondelete="CASCADE"), index=True)
    entity_type: Mapped[str] = mapped_column(String(32), index=True)
    entity_key: Mapped[str] = mapped_column(String(256), index=True)
    output_name: Mapped[str] = mapped_column(String(128))
    raw_score: Mapped[float | None] = mapped_column(Float)
    validated_label: Mapped[str] = mapped_column(String(128), default="")
    calibration: Mapped[str] = mapped_column(Text, default="")
    applicability: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    limitations: Mapped[list[str]] = mapped_column(JSON, default=list)
    raw_artifact_id: Mapped[str | None] = mapped_column(ForeignKey("artifacts.id", ondelete="SET NULL"))


class Artifact(Base, TimestampMixin):
    __tablename__ = "artifacts"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    run_id: Mapped[str | None] = mapped_column(ForeignKey("runs.id", ondelete="CASCADE"), index=True)
    kind: Mapped[str] = mapped_column(String(64), index=True)
    relative_path: Mapped[str] = mapped_column(Text)
    checksum_sha256: Mapped[str] = mapped_column(String(64), index=True)
    size_bytes: Mapped[int] = mapped_column(Integer)
    media_type: Mapped[str] = mapped_column(String(128), default="application/octet-stream")
    sensitive: Mapped[bool] = mapped_column(Boolean, default=True)
    metadata_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class AuditEvent(Base):
    __tablename__ = "audit_events"

    sequence: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    occurred_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utc_now, nullable=False)
    event_type: Mapped[str] = mapped_column(String(64), index=True)
    actor: Mapped[str] = mapped_column(String(128), default="local-user")
    entity_type: Mapped[str] = mapped_column(String(64), default="")
    entity_id: Mapped[str] = mapped_column(String(128), default="")
    payload: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    previous_hash: Mapped[str] = mapped_column(String(64), default="")
    event_hash: Mapped[str] = mapped_column(String(64), unique=True)


class DeletionRecord(Base):
    __tablename__ = "deletion_records"

    id: Mapped[str] = mapped_column(String(32), primary_key=True, default=new_id)
    occurred_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utc_now)
    entity_type: Mapped[str] = mapped_column(String(64))
    tombstone_hash: Mapped[str] = mapped_column(String(64), unique=True)
    deleted_counts: Mapped[dict[str, int]] = mapped_column(JSON, default=dict)
    retained_counts: Mapped[dict[str, int]] = mapped_column(JSON, default=dict)
