"""Acceptance tests for the schema-3 evidence-first workbench."""

from __future__ import annotations

import json
import re
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
from src.workbench.statistics import describe_numeric


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


def test_numeric_description_reports_quantiles_population_std_and_invalid_values() -> None:
    summary = describe_numeric([0.0, 0.2, 0.8, 1.0, 1.2, None], valid_range=(0.0, 1.0))

    assert summary == {
        "row_count": 6,
        "valid_count": 4,
        "missing_count": 1,
        "invalid_count": 1,
        "minimum": 0.0,
        "p10": 0.06,
        "q1": 0.15,
        "median": 0.5,
        "mean": 0.5,
        "q3": 0.85,
        "p90": 0.94,
        "maximum": 1.0,
        "iqr": 0.7,
        "population_std": 0.412311,
        "availability": "available",
        "explanation": None,
    }


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
    for forbidden in ("reference sample", "cohort comparison", "p-value", "FDR", "≥30"):
        assert forbidden.casefold() not in rendered.casefold()


def test_statistics_handle_empty_data_single_observations_and_unavailable_promoter_coordinates() -> None:
    empty = build_canonical_report({"gene": "EMPTY", "variants": [], "methylation": []})
    assert empty["sections"]["statistics"]["status"] == "no_data"

    report = build_canonical_report(
        {
            "gene": "REVERSE",
            "region": "19:100-250",
            "analysis_scope": "gene_only",
            "scope_regions": {"promoter_only": "", "gene_only": "19:100-250", "promoter_plus_gene": "19:100-250"},
            "variants": [
                {"CHROM": "chr19", "POS": 250, "REF": "T", "ALT": "C", "GT": "0/1", "FILTER": "PASS"}
            ],
            "methylation": [{"probe_id": "cg1", "beta_value": 0.4}],
        }
    )
    statistics = report["sections"]["statistics"]
    region_counts = {row["category"]: row["count"] for row in statistics["variant_statistics"]["by_region"]}
    promoter_density = next(row for row in statistics["variant_statistics"]["density"] if row["region"] == "promoter")
    beta = next(row for row in statistics["methylation_statistics"]["subsets"] if row["subset"] == "all_rows")

    assert region_counts["gene_body"] == 1
    assert promoter_density["variants_per_kb"] is None
    assert promoter_density["explanation"] == "Region coordinates are unavailable."
    assert beta["population_std"] == 0.0
    assert next(row for row in statistics["variant_statistics"]["quality"] if row["metric"] == "genotype_quality")["explanation"] == "The field has no valid numeric values."


def test_version_two_uses_the_version_one_color_palette() -> None:
    stylesheet = (Path(__file__).parents[1] / "src" / "static" / "v2.css").read_text(encoding="utf-8")
    report = render_evidence_first_html(
        build_canonical_report({"gene": "DRD4", "variants": [], "methylation": []})
    )
    version_one_colors = {"#f8e7e6", "#f2d3d6", "#2a1118", "#6b4b56", "#a1143d", "#d11f4f", "#7c0d2d"}

    assert version_one_colors <= set(color.lower() for color in re.findall(r"#[0-9a-fA-F]{6}", stylesheet))
    assert {"#f8e7e6", "#f2d3d6", "#2a1118", "#6b4b56", "#a1143d", "#d11f4f"} <= set(
        color.lower() for color in re.findall(r"#[0-9a-fA-F]{6}", report)
    )
    assert "#126a55" not in stylesheet
    assert "#176b58" not in report


def test_canonical_report_ignores_reference_comparisons_and_builds_single_person_statistics() -> None:
    variants = pd.DataFrame(
        [
            {"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "G", "GT": "0/1", "FILTER": "PASS", "GQ": 30, "DP": 20, "QUAL": 10, "ID": "rs1"},
            {"CHROM": "1", "POS": 150, "REF": "G", "ALT": "C", "GT": "0/1", "FILTER": "PASS", "GQ": 30, "DP": 20, "QUAL": 20, "ID": "rs2"},
            {"CHROM": "1", "POS": 200, "REF": "C", "ALT": "T", "GT": "1/1", "FILTER": "PASS", "GQ": 30, "DP": 20, "QUAL": 30, "ID": "."},
            {"CHROM": "1", "POS": 250, "REF": "T", "ALT": "TA", "GT": "0/1", "FILTER": "PASS", "GQ": 30, "DP": 20, "QUAL": 40, "ID": "custom"},
            {"CHROM": "1", "POS": 300, "REF": "A", "ALT": "G,T", "GT": "1/2", "FILTER": "PASS", "GQ": 30, "DP": 20, "QUAL": 50, "ID": "rs5"},
            {"CHROM": "1", "POS": 310, "REF": "A", "ALT": "G", "GT": "0/0", "FILTER": "PASS", "GQ": 30, "DP": 20, "QUAL": 60, "ID": "rs6"},
            {"CHROM": "1", "POS": 320, "REF": "A", "ALT": "G", "GT": "./.", "FILTER": "PASS", "GQ": 30, "DP": 20, "QUAL": 70, "ID": "rs7"},
            {"CHROM": "1", "POS": 330, "REF": "A", "ALT": "G", "GT": "0/1", "FILTER": "LowQual", "GQ": 30, "DP": 20, "QUAL": 80, "ID": "rs8"},
        ]
    )
    methylation = pd.DataFrame(
        [
            {"probe_id": "cg1", "CHR": "chr1", "MAPINFO": 100, "beta_value": 0.0, "m_value": -3, "detection_p": 0.001, "bead_count": 5, "UCSC_RefGene_Name": "TEST", "UCSC_RefGene_Group": "TSS200", "Relation_to_UCSC_CpG_Island": "Island"},
            {"probe_id": "cg2", "CHR": "1", "MAPINFO": 150, "beta_value": 0.2, "m_value": -2, "detection_p": 0.001, "bead_count": 5, "UCSC_RefGene_Name": "TEST;OTHER", "UCSC_RefGene_Group": "TSS1500", "Relation_to_UCSC_CpG_Island": "Shore"},
            {"probe_id": "cg3", "CHR": "1", "MAPINFO": 200, "beta_value": 0.8, "m_value": 2, "detection_p": 0.001, "bead_count": 5, "UCSC_RefGene_Name": "OTHER", "UCSC_RefGene_Group": "Body", "Relation_to_UCSC_CpG_Island": "Island"},
            {"probe_id": "cg4", "CHR": "1", "MAPINFO": 300, "beta_value": 1.0, "m_value": 3, "detection_p": 0.001, "bead_count": 5, "UCSC_RefGene_Name": "TEST", "UCSC_RefGene_Group": "Body", "Relation_to_UCSC_CpG_Island": "OpenSea"},
            {"probe_id": "cg5", "CHR": "1", "MAPINFO": 250, "beta_value": 1.2, "m_value": 4, "detection_p": 0.001, "bead_count": 5, "UCSC_RefGene_Name": "TEST", "UCSC_RefGene_Group": "Body", "Relation_to_UCSC_CpG_Island": "OpenSea"},
            {"probe_id": "cg6", "CHR": "1", "MAPINFO": 260, "beta_value": None, "m_value": None, "detection_p": 0.2, "bead_count": 2, "UCSC_RefGene_Name": "", "UCSC_RefGene_Group": "", "Relation_to_UCSC_CpG_Island": ""},
        ]
    )
    report = build_canonical_report(
        {
            "gene": "TEST",
            "region": "1:50-350",
            "analysis_scope": "promoter_plus_gene",
            "scope_regions": {"promoter_only": "1:50-150", "gene_only": "1:150-350", "promoter_plus_gene": "1:50-350"},
            "variants": variants,
            "methylation": methylation,
            "knowledge_base": {"gene_context": {"relevant_methylation_probe_ids": ["cg2", "cg4"]}},
            "statistical_comparisons": [
                {"entity_key": "must-be-ignored", "reference_values": [0.1] * 100, "raw_p": 0.001}
            ],
        }
    )
    statistics = report["sections"]["statistics"]
    variant = statistics["variant_statistics"]
    methylation_stats = statistics["methylation_statistics"]
    substitutions = {row["category"]: row["count"] for row in variant["substitutions"]}
    regions = {row["category"]: row["count"] for row in variant["by_region"]}
    types = {row["category"]: row["count"] for row in variant["by_type"]}
    dosage = {row["base"]: row["allele_copies"] for row in variant["dosage_weighted_alternate_copies"]}
    all_beta = next(row for row in methylation_stats["subsets"] if row["subset"] == "all_rows")

    assert statistics["scope"] == "single_person"
    assert statistics["status"] == "descriptive"
    assert variant["counts"] == {
        "total_count": 8,
        "raw_row_count": 8,
        "filter_pass_count": 7,
        "filter_non_pass_count": 1,
        "qc_passing_non_reference_count": 5,
        "non_reference_count": 6,
        "reference_count": 1,
        "missing_genotype_count": 1,
        "named_rsid_count": 6,
        "named_variant_count": 6,
        "unlabeled_count": 1,
        "unlabeled_variant_count": 1,
    }
    assert regions == {"promoter": 1, "gene_body": 3, "promoter_and_gene": 1, "other_analyzed_region": 0, "unclassified": 0}
    assert types["snv"] == 3 and types["insertion"] == 1 and types["multiallelic"] == 1
    assert dosage == {"A": 0.0, "C": 1.0, "G": 2.0, "T": 3.0}
    assert substitutions["A>G"] == 2 and substitutions["A>T"] == 1 and substitutions["G>C"] == 1 and substitutions["C>T"] == 1
    assert variant["transition_transversion"] == {"transition_count": 3, "transversion_count": 2, "ratio": 1.5}
    assert next(row for row in variant["quality"] if row["metric"] == "qual")["median"] == 30.0
    assert methylation_stats["counts"]["valid_beta_count"] == 4
    assert methylation_stats["counts"]["invalid_beta_count"] == 1
    assert methylation_stats["counts"]["missing_beta_count"] == 1
    assert methylation_stats["counts"]["gene_named_probe_count"] == 4
    assert methylation_stats["counts"]["curated_whitelist_probe_count"] == 2
    assert all_beta["mean"] == all_beta["median"] == 0.5
    assert sum(row["count"] for row in methylation_stats["beta_histogram"]) == 4
    assert methylation_stats["extremes"]["highest"][0]["probe_id"] == "cg4"
    assert methylation_stats["extremes"]["lowest"][0]["probe_id"] == "cg1"
    assert "must-be-ignored" not in json.dumps(statistics)
    assert not {"raw_p", "q_value", "percentile"} & statistics["records"][0].keys()


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
