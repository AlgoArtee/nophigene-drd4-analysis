"""Coverage for the local asynchronous gene workflow API."""

from __future__ import annotations

import threading
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from flask import Flask
from sqlalchemy.orm import Session

from src.api.jobs import JobManager
from src.api.profiles import ProfileStore
from src.api.routes import api_v1
from src.api.v2_routes import api_v2
from src.api.serialization import write_json_atomic
from src.api.workflow_runner import WorkflowRunner, normalize_job_request
from src.workbench.database import create_database_engine, ensure_schema
from src.workbench.models import Run
from src.workflow import select_profile_variant_source


def _profile_files(tmp_path: Path) -> dict[str, Path]:
    sample = tmp_path / "sample"
    for suffix in ("_Grn.idat", "_Red.idat"):
        sample.with_name(sample.name + suffix).write_bytes(b"idat")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("IlmnID,CHR,MAPINFO,UCSC_RefGene_Name\n", encoding="utf-8")
    vcf19 = tmp_path / "sample_hg19.vcf"
    vcf19.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    vcf38 = tmp_path / "sample_hg38.vcf"
    vcf38.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    bam = tmp_path / "sample_hg38.bam"
    bam.write_bytes(b"bam")
    return {
        "sample": sample,
        "manifest": manifest,
        "vcf19": vcf19,
        "vcf38": vcf38,
        "bam": bam,
    }


def _profile_payload(files: dict[str, Path], *, profile_id: str = "sample-one") -> dict[str, object]:
    return {
        "id": profile_id,
        "display_name": "Sample One",
        "default_genome_build": "hg19",
        "idat_prefix": str(files["sample"]),
        "manifest_path": str(files["manifest"]),
        "population_statistics_path": "",
        "vcf_sources": [
            {"path": str(files["vcf19"]), "genome_build": "hg19"},
            {"path": str(files["vcf38"]), "genome_build": "hg38"},
        ],
        "bam_sources": [{"path": str(files["bam"]), "genome_build": "hg38"}],
    }


def _test_app(profile_store: ProfileStore, manager: JobManager) -> Flask:
    app = Flask(__name__)
    app.config.update(
        TESTING=True,
        NOPHIGENE_PROFILE_STORE=profile_store,
        NOPHIGENE_JOB_MANAGER=manager,
    )
    app.register_blueprint(api_v1)
    return app


def test_v1_keeps_read_endpoints_and_freezes_mutations(tmp_path: Path) -> None:
    store = ProfileStore(tmp_path / "profiles.json")
    manager = JobManager(jobs_root=tmp_path / "jobs", profile_store=store)
    client = _test_app(store, manager).test_client()

    listing = client.get("/api/v1/knowledge-sources")
    assert listing.status_code == 200
    assert listing.get_json()["count"] >= 1
    assert client.get("/api/v1/knowledge-workflows").status_code == 200

    for endpoint, payload in (
        ("/api/v1/knowledge-sources/test", {"source_key": "clinvar"}),
        ("/api/v1/profiles", {"id": "sample-one"}),
        ("/api/v1/jobs", {"operation": "resolve_regions", "genes": ["DRD4"]}),
    ):
        response = client.post(endpoint, json=payload)
        assert response.status_code == 410
        body = response.get_json()
        assert body["error"]["code"] == "api_v1_read_only"
        assert body["error"]["details"]["replacement"] == "/api/v2"


def test_v2_exposes_schema_three_and_enforces_batch_limit(tmp_path: Path) -> None:
    store = ProfileStore(tmp_path / "profiles.json")
    manager = JobManager(jobs_root=tmp_path / "jobs", profile_store=store)
    engine = create_database_engine(test_url="sqlite+pysqlite:///:memory:")
    ensure_schema(engine)
    app = Flask(__name__)
    app.config.update(
        TESTING=True,
        NOPHIGENE_PROFILE_STORE=store,
        NOPHIGENE_JOB_MANAGER=manager,
        NOPHIGENE_DATABASE_ENGINE=engine,
    )
    app.register_blueprint(api_v2)
    client = app.test_client()

    index = client.get("/api/v2/")
    assert index.status_code == 200
    assert index.get_json()["report_schema_version"] == "3.0"
    assert client.get("/api/v2/evidence/sources").get_json()["count"] >= 20
    assert client.get("/api/v2/models").get_json()["count"] >= 10

    oversized = client.post(
        "/api/v2/runs",
        json={"operation": "resolve_regions", "genes": [f"GENE{i}" for i in range(101)]},
    )
    assert oversized.status_code == 422
    assert oversized.get_json()["error"]["code"] == "too_many_genes"

    submitted = client.post(
        "/api/v2/runs",
        json={"operation": "resolve_regions", "gene": "DRD4", "sample_context": {"tissue": "whole blood"}},
    )
    assert submitted.status_code == 202
    run_id = submitted.get_json()["id"]
    assert submitted.get_json()["schema_version"] == "3.0"
    with Session(engine) as session:
        persisted = session.get(Run, run_id)
        assert persisted is not None
        assert persisted.genes == ["DRD4"]


def test_job_validation_rejects_empty_and_oversized_gene_lists() -> None:
    with pytest.raises(Exception) as empty:
        normalize_job_request({"operation": "resolve_regions", "genes": []})
    assert getattr(empty.value, "code", "") == "invalid_genes"

    with pytest.raises(Exception) as oversized:
        normalize_job_request(
            {
                "operation": "resolve_regions",
                "genes": [f"GENE{i}" for i in range(101)],
            }
        )
    assert getattr(oversized.value, "code", "") == "too_many_genes"

    with pytest.raises(Exception) as traversal_gene:
        normalize_job_request({"operation": "resolve_regions", "genes": [".."]})
    assert getattr(traversal_gene.value, "code", "") == "invalid_gene"

    normalized = normalize_job_request(
        {
            "operation": "build_knowledge_bases",
            "genes": ["DRD4"],
            "profile_id": "sample-one",
            "options": {
                "knowledge_workflows": "clinical_variant_triage,population_frequency_association",
                "knowledge_sources": ["clinvar"],
                "use_local_article_evidence": True,
                "article_pdf_folder": "data/articles",
                "article_pdf_recursive": False,
                "max_article_pdfs": "25",
                "interpretation_mode": "dual",
                "evidence_snapshot_id": "evidence-123",
                "requested_models": [{"model_id": "PGS000001"}],
            },
        }
    )
    assert normalized["options"]["knowledge_workflows"] == [
        "clinical_variant_triage",
        "population_frequency_association",
    ]
    assert normalized["options"]["use_local_article_evidence"] is True
    assert normalized["options"]["article_pdf_folder"] == "data/articles"
    assert normalized["options"]["article_pdf_recursive"] is False
    assert normalized["options"]["max_article_pdfs"] == 25
    assert normalized["options"]["interpretation_mode"] == "dual"
    assert normalized["options"]["evidence_snapshot_id"] == "evidence-123"
    assert normalized["options"]["requested_models"] == [{"model_id": "PGS000001"}]

    with pytest.raises(Exception) as bad_mode:
        normalize_job_request(
            {
                "operation": "build_knowledge_bases",
                "genes": ["DRD4"],
                "profile_id": "sample-one",
                "options": {"interpretation_mode": "diagnostic"},
            }
        )
    assert getattr(bad_mode.value, "code", "") == "invalid_interpretation_mode"

    with pytest.raises(Exception) as bad_article_limit:
        normalize_job_request(
            {
                "operation": "build_knowledge_bases",
                "genes": ["DRD4"],
                "profile_id": "sample-one",
                "options": {"max_article_pdfs": "many"},
            }
        )
    assert getattr(bad_article_limit.value, "code", "") == "invalid_max_article_pdfs"


def test_running_jobs_are_marked_interrupted_and_queued_jobs_can_be_cancelled(tmp_path: Path) -> None:
    jobs_root = tmp_path / "jobs"
    running_id = "a" * 32
    queued_id = "b" * 32
    base_job = {
        "operation": "resolve_regions",
        "genes": ["DRD4"],
        "stage": "starting",
        "progress": {"completed": 0, "total": 1, "percent": 0},
        "outcomes": [],
        "error": None,
        "created_at": "2026-01-01T00:00:00Z",
        "updated_at": "2026-01-01T00:00:00Z",
        "started_at": None,
        "finished_at": None,
    }
    write_json_atomic(
        jobs_root / running_id / "job.json",
        {"id": running_id, "status": "running", **base_job},
    )
    write_json_atomic(
        jobs_root / queued_id / "job.json",
        {"id": queued_id, "status": "queued", **base_job},
    )
    write_json_atomic(
        jobs_root / queued_id / "request.json",
        normalize_job_request({"operation": "resolve_regions", "genes": ["DRD4"]}),
    )

    blocker = threading.Event()

    class BlockingRunner:
        def execute(self, *_args, **_kwargs):
            blocker.wait(timeout=2)
            return {
                "status": "succeeded",
                "genes": [],
            }

    manager = JobManager(
        jobs_root=jobs_root,
        profile_store=ProfileStore(tmp_path / "profiles.json"),
        runner=BlockingRunner(),
    )
    manager.start()
    assert manager.get(running_id)["error"]["code"] == "interrupted"

    # The recovered queued job may already be running, so create a second queued
    # record directly to test the cancellation contract deterministically.
    cancellable_id = "c" * 32
    write_json_atomic(
        jobs_root / cancellable_id / "job.json",
        {"id": cancellable_id, "status": "queued", **base_job},
    )
    cancelled = manager.cancel(cancellable_id)
    assert cancelled["status"] == "cancelled"
    blocker.set()


def test_variant_source_prefers_matching_vcf_then_bam(tmp_path: Path) -> None:
    profile = {
        "id": "sample",
        "vcf_sources": [{"path": str(tmp_path / "hg19.vcf"), "genome_build": "hg19"}],
        "bam_sources": [{"path": str(tmp_path / "hg38.bam"), "genome_build": "hg38"}],
    }
    assert select_profile_variant_source(profile, "GRCh37")["type"] == "vcf"
    assert select_profile_variant_source(profile, "hg38")["type"] == "bam"


def test_all_workflow_operations_and_multi_gene_methylprep_reuse(
    monkeypatch,
    tmp_path: Path,
) -> None:
    files = _profile_files(tmp_path)
    profile = _profile_payload(files)
    profile["id"] = "sample-one"

    class StaticStore:
        def get(self, profile_id: str):
            assert profile_id == "sample-one"
            return profile

    jobs_root = tmp_path / "jobs"
    runner = WorkflowRunner(StaticStore(), jobs_root)
    calls = {"manifest": 0, "methylprep": 0, "analysis": 0}

    full_manifest = pd.DataFrame(
        [
            {
                "IlmnID": "cg1",
                "CHR": "1",
                "MAPINFO": 150,
                "CHR_hg38": "chr1",
                "Start_hg38": 150,
                "UCSC_RefGene_Name": "GENE1",
            },
            {
                "IlmnID": "cg2",
                "CHR": "2",
                "MAPINFO": 350,
                "CHR_hg38": "chr2",
                "Start_hg38": 350,
                "UCSC_RefGene_Name": "GENE2",
            },
        ]
    )

    def fake_manifest(_path: str) -> pd.DataFrame:
        calls["manifest"] += 1
        return full_manifest

    def fake_beta(_prefix: str, manifest_filepath: str | None = None) -> pd.DataFrame:
        calls["methylprep"] += 1
        assert manifest_filepath
        return pd.DataFrame([{"probe_id": "cg1", "beta": 0.4}, {"probe_id": "cg2", "beta": 0.8}])

    def fake_resolve(gene: str, **_kwargs):
        chrom = "1" if gene == "GENE1" else "2"
        build = "hg19" if gene == "GENE1" else "hg38"
        region = f"{'chr' if build == 'hg38' else ''}{chrom}:{100 if gene == 'GENE1' else 300}-{200 if gene == 'GENE1' else 400}"
        return {
            "gene": gene,
            "genome_build": build,
            "region": region,
            "scope": "promoter_plus_gene",
            "scope_regions": {"promoter_plus_gene": region, "promoter_only": "", "gene_only": region},
            "selected_gene_region": region,
            "selected_sources": ["test"],
            "candidate_regions": [],
            "curated_coordinates": True,
        }

    def fake_variants(path: str, region: str) -> pd.DataFrame:
        return pd.DataFrame(
            [{"sample": "sample", "chrom": region.split(":")[0], "id": "rs1", "pos": 150, "ref": "A", "alt": "G"}]
        )

    def fake_analysis(**kwargs):
        calls["analysis"] += 1
        variants = kwargs["variants"]
        methylation = kwargs["methylation"]
        return SimpleNamespace(
            variants=variants,
            methylation=methylation,
            popstats=kwargs.get("popstats"),
            analysis_scope="promoter_plus_gene",
            analysis_scope_label="Promoter + gene",
            variant_interpretations={"gene_name": kwargs["gene_name"], "matched_records": []},
            methylation_insights={
                "gene_name": kwargs["gene_name"],
                "probe_preview": methylation.copy(),
            },
            knowledge_base={"gene_context": {"gene_name": kwargs["gene_name"]}},
            population_insights={},
            population_database={},
            general_database_path=tmp_path / "general.csv",
            general_database_status=(
                "updated" if kwargs["update_general_database_enabled"] else "not requested"
            ),
        )

    def fake_report(_variants, _methylation, _popstats, output_path: str, **_kwargs):
        path = Path(output_path)
        path.write_text("<html>report</html>", encoding="utf-8")
        return path

    monkeypatch.setattr("src.api.workflow_runner.load_full_methylation_manifest", fake_manifest)
    monkeypatch.setattr("src.api.workflow_runner.load_methylation_beta_values", fake_beta)
    monkeypatch.setattr("src.api.workflow_runner.resolve_gene_region", fake_resolve)
    monkeypatch.setattr("src.api.workflow_runner.load_variants", fake_variants)
    monkeypatch.setattr("src.api.workflow_runner.analyze_prepared_data", fake_analysis)
    monkeypatch.setattr("src.api.workflow_runner.generate_report", fake_report)

    def run(operation: str, job_id: str, *, source_job_id: str = ""):
        request_payload = normalize_job_request(
            {
                "operation": operation,
                "genes": ["GENE1", "GENE2"],
                "profile_id": "sample-one" if operation != "render_reports" else "",
                "source_job_id": source_job_id,
            }
        )
        job_dir = jobs_root / job_id
        return runner.execute(job_id, request_payload, job_dir, lambda *_args: None)

    assert run("resolve_regions", "1" * 32)["status"] == "succeeded"
    assert run("prepare_manifests", "2" * 32)["status"] == "succeeded"
    assert run("extract_variants", "3" * 32)["status"] == "succeeded"
    analysis = run("analyze", "4" * 32)
    assert analysis["status"] == "succeeded"
    rendered = run("render_reports", "5" * 32, source_job_id="4" * 32)
    assert rendered["status"] == "succeeded"
    full = run("full_workflow", "6" * 32)
    assert full["status"] == "succeeded"

    assert calls["manifest"] == 3  # prepare, analyze, full
    assert calls["methylprep"] == 2  # once for analyze and once for full
    assert calls["analysis"] == 4  # two genes in analyze and two in full
    for gene in ("GENE1", "GENE2"):
        gene_dir = jobs_root / ("6" * 32) / "genes" / gene
        assert (gene_dir / "report.html").is_file()
        assert (gene_dir / "report.json").is_file()
        assert (gene_dir / "report_summary.csv").is_file()
        assert (gene_dir / "variants.csv").is_file()
        assert (gene_dir / "methylation.csv").is_file()
        report_payload = (gene_dir / "report.json").read_text(encoding="utf-8")
        assert '"schema_version": "3.0"' in report_payload
        assert '"source_provenance"' in report_payload
    assert (jobs_root / ("6" * 32) / "artifacts.zip").is_file()


def test_batch_continues_after_per_gene_failure(monkeypatch, tmp_path: Path) -> None:
    runner = WorkflowRunner(ProfileStore(tmp_path / "profiles.json"), tmp_path / "jobs")

    def resolve(gene: str, **_kwargs):
        if gene == "BAD":
            raise ValueError("No region")
        return {
            "gene": gene,
            "genome_build": "hg19",
            "region": "1:1-10",
            "scope": "promoter_plus_gene",
            "scope_regions": {"promoter_plus_gene": "1:1-10"},
            "selected_gene_region": "1:1-10",
            "selected_sources": ["test"],
            "candidate_regions": [],
            "curated_coordinates": False,
        }

    monkeypatch.setattr("src.api.workflow_runner.resolve_gene_region", resolve)
    request_payload = normalize_job_request(
        {"operation": "resolve_regions", "genes": ["GOOD", "BAD"]}
    )
    result = runner.execute(
        "d" * 32,
        request_payload,
        tmp_path / "jobs" / ("d" * 32),
        lambda *_args: None,
    )
    assert result["status"] == "partial"
    assert result["counts"] == {"requested": 2, "succeeded": 1, "failed": 1}
    assert result["genes"][1]["error"]["code"] == "gene_workflow_failed"
