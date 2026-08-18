"""Acceptance tests for the schema-3 evidence-first workbench."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from src.workbench.audit import append_audit_event, verify_audit_chain
from src.workbench.artifacts import ArtifactStore
from src.workbench.database import DatabaseSecurityError, create_database_engine, ensure_schema
from src.workbench.evidence import deduplicate_literature, medical_gate, source_coverage
from src.workbench.interactions import build_interaction_graph
from src.workbench.legacy_migration import import_legacy_reports, inventory_legacy_stores
from src.workbench.model_registry import inspect_model_inputs, list_model_manifests
from src.workbench.models import EvidenceRecord, MethylationMeasurement, Run, VariantCall
from src.workbench.persistence import persist_canonical_report
from src.workbench.pgx import resolve_pgx_diplotype
from src.workbench.reporting import build_canonical_report, render_evidence_first_html
from src.workbench.statistics import apply_family_fdr, empirical_single_sample_result, reference_compatibility


def test_report_schema_separates_six_sections_and_limits_primary_tables() -> None:
    variants = pd.DataFrame(
        [
            {"CHROM": "11", "POS": 100 + index, "REF": "A", "ALT": "G", "GT": "0/1", "FILTER": "PASS", "GQ": 30, "DP": 15}
            for index in range(25)
        ]
        + [{"CHROM": "11", "POS": 999, "REF": "A", "ALT": "T", "GT": "0/1", "FILTER": "LowQual", "GQ": 8, "DP": 2}]
    )
    methylation = pd.DataFrame(
        [{"probe_id": f"cg{index:08d}", "beta_value": 0.25, "detection_p": 0.001} for index in range(24)]
        + [{"probe_id": "cg_failed", "beta_value": 0.9, "detection_p": 0.2}]
    )
    report = build_canonical_report({"gene": "DRD4", "genome_build": "GRCh38", "variants": variants, "methylation": methylation})

    assert report["schema_version"] == "3.0"
    assert list(report["sections"]) == ["objective_data", "statistics", "literature", "medical", "interactions", "predictions"]
    assert len(report["sections"]["objective_data"]["variants"]) == 20
    assert len(report["sections"]["objective_data"]["methylation"]) == 20
    assert len(report["sections"]["objective_data"]["all_variants"]) == 26
    assert report["summary"]["qc"]["variant_pass_count"] == 25
    assert "phenotype_prediction" not in str(report)


def test_medical_gate_keeps_exploratory_records_in_literature() -> None:
    approved = medical_gate(
        {
            "source_key": "cpic",
            "authority": "CPIC",
            "evidence_type": "professional_guideline",
            "effective_date": "2025-01-01",
        }
    )
    rejected = medical_gate(
        {
            "source_key": "biorxiv",
            "authority": "biorxiv",
            "evidence_type": "gwas_association",
            "effective_date": "2025-01-01",
            "preprint": True,
        }
    )
    assert approved["eligible"] is True
    assert rejected["eligible"] is False
    assert rejected["destination"] == "literature"


def test_literature_dedup_and_provider_states_are_explicit() -> None:
    records = deduplicate_literature(
        [
            {"pmid": "123", "title": "A study", "source_key": "pubmed"},
            {"pmid": "123", "title": "A study", "source_key": "pmc", "pmcid": "PMC1"},
            {"doi": "10.1/example", "title": "A preprint", "source_key": "biorxiv"},
        ]
    )
    coverage = source_coverage(
        [
            {"source_key": "clinvar", "status": "empty", "record_count": 0},
            {"source_key": "gtex", "status": "failed"},
            {"source_key": "string", "status": "not_assessed"},
        ]
    )
    assert len(records) == 2
    assert records[1]["preprint"] is True
    assert coverage["assessed_absence"] is True
    assert coverage["failed"][0]["source_key"] == "gtex"
    assert coverage["not_assessed"][0]["source_key"] == "string"


def test_statistics_require_compatible_context_and_thirty_references() -> None:
    exact = {"tissue": "blood", "platform": "EPIC", "normalization": "SeSAMe", "genome_build": "GRCh38"}
    assert reference_compatibility(exact, dict(exact))["compatible"] is True
    mismatch = reference_compatibility(exact, {**exact, "tissue": "brain"})
    assert mismatch["compatible"] is False
    assert mismatch["mismatches"] == ["tissue"]
    assert empirical_single_sample_result(0.2, [0.1] * 29)["status"] == "descriptive_only"
    exploratory = empirical_single_sample_result(0.2, [index / 100 for index in range(30)])
    assert exploratory["status"] == "exploratory"
    rows = apply_family_fdr(
        [
            {"family": "methylation", "raw_p": 0.01},
            {"family": "methylation", "raw_p": 0.04},
            {"family": "variant", "raw_p": 0.02},
        ]
    )
    assert rows[0]["q_value"] == 0.02
    assert rows[1]["q_value"] == 0.04
    assert rows[2]["q_value"] == 0.02


def test_interaction_graph_uses_gene_hops_and_hard_cap() -> None:
    edges = [
        {
            "source_gene": "DRD4" if index == 0 else f"G{index}",
            "target_gene": f"G{index + 1}",
            "edge_type": "physical_binding",
            "source_key": "intact",
            "native_score": 0.8,
            "directed": False,
        }
        for index in range(200)
    ]
    graph = build_interaction_graph("DRD4", edges, max_hops=3, node_cap=150)
    assert graph["max_hops"] == 3
    assert graph["node_count"] <= 150
    assert all(edge["hop"] <= 3 for edge in graph["edges"])
    assert all("not a probability" in edge["combined_rank_policy"] for edge in graph["edges"])


def test_model_contracts_block_invalid_inputs_and_never_offer_consensus() -> None:
    manifests = {item["id"]: item for item in list_model_manifests()}
    assert {"alphagenome-api", "borzoi-local", "methylbert-local", "melody-local"} <= set(manifests)
    alpha = inspect_model_inputs(
        "alphagenome-api",
        {
            "genome_build": "GRCh37",
            "grch38_variant": True,
            "reference_validated": True,
            "requested_modalities": True,
        },
    )
    assert alpha["eligible"] is False
    assert "external_transfer_not_approved" in alpha["blockers"]
    assert "unsupported_build:GRCh37" in alpha["blockers"]
    methylbert = inspect_model_inputs("methylbert-local", {"epic_beta_table": True})
    assert methylbert["eligible"] is False
    assert "bismark_bam_with_xm" in methylbert["blockers"]


def test_database_schema_and_audit_chain_detect_tampering() -> None:
    engine = create_database_engine(test_url="sqlite+pysqlite:///:memory:")
    ensure_schema(engine)
    with Session(engine) as session:
        first = append_audit_event(session, "run_created", entity_type="run", entity_id="r1", payload={"gene": "DRD4"})
        append_audit_event(session, "export_created", entity_type="run", entity_id="r1")
        session.commit()
        assert verify_audit_chain(session)["valid"] is True
        first.payload = {"gene": "HERC2"}
        session.commit()
        assert verify_audit_chain(session)["valid"] is False


def test_production_database_fails_closed_without_key(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("NOPHIGENE_REQUIRE_ENCRYPTION", "1")
    monkeypatch.delenv("NOPHIGENE_DATABASE_KEY_FILE", raising=False)
    try:
        create_database_engine(database_path=tmp_path / "workbench.db")
    except DatabaseSecurityError:
        pass
    else:
        raise AssertionError("Production database accepted a plaintext configuration")


def test_html_has_accessible_seven_result_tabs_and_run_details() -> None:
    report = build_canonical_report({"gene": "DRD4", "variants": [], "methylation": []})
    rendered = render_evidence_first_html(report)
    assert rendered.count('role="tab"') == 7
    assert "Objective Data" in rendered
    assert "Scientific Literature" in rendered
    assert "Medical Information" in rendered
    assert "Run Details" in rendered
    assert "Predictive Theses" not in rendered


def test_canonical_report_runs_only_compatible_raw_reference_comparisons() -> None:
    context = {"tissue": "blood", "platform": "EPIC", "normalization": "SeSAMe", "genome_build": "GRCh38"}
    report = build_canonical_report(
        {
            "gene": "DRD4",
            "sample_context": context,
            "statistical_comparisons": [
                {
                    "entity_type": "methylation_probe",
                    "entity_key": "cg1",
                    "family": "methylation",
                    "test_value": 0.4,
                    "reference_values": [index / 100 for index in range(30)],
                    "cohort_id": "public-a",
                    "cohort_source_type": "public",
                    "cohort_context": context,
                },
                {
                    "entity_type": "methylation_probe",
                    "entity_key": "cg2",
                    "family": "methylation",
                    "test_value": 0.4,
                    "reference_values": [index / 100 for index in range(30)],
                    "cohort_id": "user-b",
                    "cohort_source_type": "user",
                    "cohort_context": {**context, "tissue": "brain"},
                },
            ],
        }
    )
    rows = report["sections"]["statistics"]["records"]
    assert rows[0]["status"] == "exploratory"
    assert rows[0]["q_value"] is not None
    assert rows[1]["status"] == "not_assessed"
    assert rows[1]["q_value"] is None


def test_pgx_requires_complete_qc_coverage_and_one_solution() -> None:
    definition = {
        "gene": "TEST",
        "candidates": [
            {"diplotype": "*1/*1", "required_genotypes": {"rs1": "A/A"}, "phenotype": "normal"},
            {"diplotype": "*1/*2", "required_genotypes": {"rs1": "A/G"}, "phenotype": "intermediate"},
        ],
    }
    unresolved = resolve_pgx_diplotype(definition, [{"locus": "rs1", "genotype": "A/G", "qc_pass": False}])
    resolved = resolve_pgx_diplotype(definition, [{"locus": "rs1", "genotype": "A/G", "qc_pass": True, "covered": True}])
    assert unresolved["status"] == "unresolved"
    assert unresolved["phenotype_status"] == "not_assessed"
    assert resolved["status"] == "resolved"
    assert resolved["diplotype"] == "*1/*2"


def test_canonical_report_persists_normalized_rows_idempotently() -> None:
    engine = create_database_engine(test_url="sqlite+pysqlite:///:memory:")
    ensure_schema(engine)
    with Session(engine) as session:
        session.add(Run(id="r1", status="succeeded", stage="complete", genes=["DRD4"]))
        session.commit()
        report = build_canonical_report(
            {
                "run_id": "r1",
                "gene": "DRD4",
                "genome_build": "GRCh38",
                "variants": [{"CHROM": "11", "POS": 100, "REF": "A", "ALT": "G", "GT": "0/1", "FILTER": "PASS", "GQ": 30, "DP": 20}],
                "methylation": [{"probe_id": "cg1", "beta_value": 0.2, "m_value": -2.0, "detection_p": 0.001}],
                "source_records": [{"source_key": "pubmed", "record_id": "pmid:1", "title": "A study", "pmid": "1", "source_release": "2026"}],
            }
        )
        first = persist_canonical_report(session, report)
        session.commit()
        second = persist_canonical_report(session, report)
        session.commit()
        assert first["variants"] == first["methylation"] == first["evidence"] == 1
        assert second["variants"] == second["methylation"] == second["evidence"] == 0
        assert session.scalar(select(func.count()).select_from(VariantCall)) == 1
        assert session.scalar(select(func.count()).select_from(MethylationMeasurement)) == 1
        assert session.scalar(select(func.count()).select_from(EvidenceRecord)) == 1


def test_legacy_import_archives_synthesis_and_marks_linkouts_not_assessed(tmp_path: Path) -> None:
    gene_data = tmp_path / "src" / "gene_data"
    gene_data.mkdir(parents=True)
    archived_gene_data = tmp_path / "version1" / "gene_data" / "loose"
    archived_gene_data.mkdir(parents=True)
    (archived_gene_data / "test_synthesis.json").write_text(
        json.dumps({"cases": [{"prediction": "unsupported"}]}), encoding="utf-8"
    )
    (gene_data / "test_interpretation_db.json").write_text(
        json.dumps(
            {
                "version": "legacy-1",
                "gene_context": {
                    "gene_name": "TEST",
                    "concrete_variant_prediction": "must not be imported",
                    "evidence": [{"label": "PubMed record", "url": "https://pubmed.ncbi.nlm.nih.gov/1/"}],
                },
            }
        ),
        encoding="utf-8",
    )
    inventory = inventory_legacy_stores(tmp_path)
    synthetic = next(item for item in inventory.candidates if item.kind == "synthetic_synthesis")
    assert synthetic.status == "archive_only"
    assert synthetic.excluded_fields == ["entire_file"]

    engine = create_database_engine(test_url="sqlite+pysqlite:///:memory:")
    ensure_schema(engine)
    with Session(engine) as session:
        result = import_legacy_reports(session, inventory, artifact_store=ArtifactStore(tmp_path / "artifacts"))
        session.commit()
        record = session.scalar(select(EvidenceRecord))
        assert result["evidence_references_imported"] == 1
        assert record is not None
        assert record.status == "not_assessed"
        assert record.evidence_type == "bibliographic_linkout"
        assert record.assertion == "PubMed record"
        assert "prediction" not in record.assertion.casefold()
