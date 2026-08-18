"""Add DANDELION statistical dataset and analysis entities.

Revision ID: 20260818_0002
Revises: 20260817_0001
"""

from alembic import op
import sqlalchemy as sa


revision = "20260818_0002"
down_revision = "20260817_0001"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # The original baseline revision was derived from live ORM metadata. A
    # brand-new installation can therefore already contain these objects when
    # this revision follows it, while an existing Version 2 database does not.
    # Keep the upgrade transactional and idempotent for both paths.
    inspector = sa.inspect(op.get_bind())
    tables = set(inspector.get_table_names())
    result_columns = (
        {column["name"] for column in inspector.get_columns("statistical_results")}
        if "statistical_results" in tables
        else set()
    )
    dandelion_tables = {
        "statistical_datasets",
        "statistical_dataset_files",
        "statistical_analysis_runs",
    }
    if dandelion_tables.issubset(tables):
        if "analysis_run_id" not in result_columns:
            with op.batch_alter_table("statistical_results") as batch:
                batch.add_column(sa.Column("analysis_run_id", sa.String(length=32), nullable=True))
                batch.create_foreign_key(
                    "fk_statistical_results_analysis_run_id",
                    "statistical_analysis_runs",
                    ["analysis_run_id"],
                    ["id"],
                    ondelete="CASCADE",
                )
                batch.create_index("ix_statistical_results_analysis_run_id", ["analysis_run_id"])
        return
    op.create_table(
        "statistical_datasets",
        sa.Column("id", sa.String(length=32), nullable=False),
        sa.Column("name", sa.Text(), nullable=False),
        sa.Column("description", sa.Text(), nullable=False),
        sa.Column("method", sa.String(length=64), nullable=False),
        sa.Column("phenotype", sa.String(length=256), nullable=False),
        sa.Column("exposure_type", sa.String(length=16), nullable=False),
        sa.Column("assembly", sa.String(length=16), nullable=False),
        sa.Column("gene_namespace", sa.String(length=64), nullable=False),
        sa.Column("storage_mode", sa.String(length=32), nullable=False),
        sa.Column("context", sa.JSON(), nullable=False),
        sa.Column("applicability", sa.JSON(), nullable=False),
        sa.Column("manifest_checksum_sha256", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), nullable=False),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("manifest_checksum_sha256"),
    )
    op.create_index("ix_statistical_datasets_method", "statistical_datasets", ["method"])
    op.create_index("ix_statistical_datasets_phenotype", "statistical_datasets", ["phenotype"])
    op.create_index("ix_statistical_datasets_exposure_type", "statistical_datasets", ["exposure_type"])
    op.create_index("ix_statistical_datasets_assembly", "statistical_datasets", ["assembly"])
    op.create_index(
        "ix_statistical_datasets_manifest_checksum_sha256",
        "statistical_datasets",
        ["manifest_checksum_sha256"],
    )
    op.create_table(
        "statistical_dataset_files",
        sa.Column("id", sa.String(length=32), nullable=False),
        sa.Column("dataset_id", sa.String(length=32), nullable=False),
        sa.Column("role", sa.String(length=64), nullable=False),
        sa.Column("relative_path", sa.Text(), nullable=False),
        sa.Column("managed_path", sa.Text(), nullable=False),
        sa.Column("file_format", sa.String(length=16), nullable=False),
        sa.Column("checksum_sha256", sa.String(length=64), nullable=False),
        sa.Column("size_bytes", sa.Integer(), nullable=False),
        sa.Column("mapping", sa.JSON(), nullable=False),
        sa.Column("inspection", sa.JSON(), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), nullable=False),
        sa.ForeignKeyConstraint(["dataset_id"], ["statistical_datasets.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("dataset_id", "role"),
    )
    op.create_index("ix_statistical_dataset_files_dataset_id", "statistical_dataset_files", ["dataset_id"])
    op.create_index("ix_statistical_dataset_files_checksum_sha256", "statistical_dataset_files", ["checksum_sha256"])
    op.create_table(
        "statistical_analysis_runs",
        sa.Column("id", sa.String(length=32), nullable=False),
        sa.Column("run_id", sa.String(length=32), nullable=False),
        sa.Column("dataset_id", sa.String(length=32), nullable=False),
        sa.Column("method", sa.String(length=64), nullable=False),
        sa.Column("method_version", sa.String(length=64), nullable=False),
        sa.Column("worker_job_id", sa.String(length=32), nullable=False),
        sa.Column("status", sa.String(length=32), nullable=False),
        sa.Column("parameters", sa.JSON(), nullable=False),
        sa.Column("resource_estimate", sa.JSON(), nullable=False),
        sa.Column("result_checksum_sha256", sa.String(length=64), nullable=False),
        sa.Column("error", sa.JSON(), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), nullable=False),
        sa.ForeignKeyConstraint(["dataset_id"], ["statistical_datasets.id"], ondelete="RESTRICT"),
        sa.ForeignKeyConstraint(["run_id"], ["runs.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("run_id"),
        sa.UniqueConstraint("worker_job_id"),
    )
    op.create_index("ix_statistical_analysis_runs_run_id", "statistical_analysis_runs", ["run_id"], unique=True)
    op.create_index("ix_statistical_analysis_runs_dataset_id", "statistical_analysis_runs", ["dataset_id"])
    op.create_index("ix_statistical_analysis_runs_method", "statistical_analysis_runs", ["method"])
    op.create_index("ix_statistical_analysis_runs_worker_job_id", "statistical_analysis_runs", ["worker_job_id"], unique=True)
    op.create_index("ix_statistical_analysis_runs_status", "statistical_analysis_runs", ["status"])
    with op.batch_alter_table("statistical_results") as batch:
        batch.add_column(sa.Column("analysis_run_id", sa.String(length=32), nullable=True))
        batch.create_foreign_key(
            "fk_statistical_results_analysis_run_id",
            "statistical_analysis_runs",
            ["analysis_run_id"],
            ["id"],
            ondelete="CASCADE",
        )
        batch.create_index("ix_statistical_results_analysis_run_id", ["analysis_run_id"])


def downgrade() -> None:
    with op.batch_alter_table("statistical_results") as batch:
        batch.drop_index("ix_statistical_results_analysis_run_id")
        batch.drop_column("analysis_run_id")
    op.drop_table("statistical_analysis_runs")
    op.drop_table("statistical_dataset_files")
    op.drop_table("statistical_datasets")
