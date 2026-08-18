"""Scientific, persistence, queue, and API tests for the DANDELION lane."""

from __future__ import annotations

import json
from pathlib import Path

import pyzipper
import pytest
from flask import Flask
from sqlalchemy.orm import Session

from src.api.v2_routes import api_v2
from src.workbench.dandelion import (
    DandelionValidationError,
    build_report,
    normalize_analysis_parameters,
    normalize_dataset_manifest,
)
from src.workbench.database import create_database_engine, ensure_schema
from src.workbench.models import StatisticalAnalysisRun, StatisticalDataset
from src.workbench.statistical_jobs import StatisticalJobError, StatisticalJobManager


def _cohort_files(root: Path) -> dict[str, object]:
    root.mkdir(parents=True)
    (root / "trans.tsv").write_text("gene\tE1\tE2\nG1\t0.01\t0.2\nG2\t0.02\t0.3\n", encoding="utf-8")
    (root / "burden.tsv").write_text("gene\tp_value\nG1\t0.001\nG2\t0.2\n", encoding="utf-8")
    (root / "genes.tsv").write_text(
        "gene_name\ttype\tChromosome\tstart\tend\n"
        "G1\tprotein_coding\tchr1\t100\t200\n"
        "G2\tprotein_coding\tchr2\t100\t200\n"
        "E1\tprotein_coding\tchr3\t100\t200\n"
        "E2\tprotein_coding\tchr4\t100\t200\n",
        encoding="utf-8",
    )
    return {
        "name": "Test cohort",
        "phenotype": "EFO:test",
        "exposure_type": "Gene",
        "assembly": "GRCh38",
        "gene_namespace": "HGNC symbol",
        "context": {"tissue": "blood", "ancestry": "test ancestry"},
        "files": {
            "trans_matrix": {"path": "trans.tsv", "mapping": {"row_id": "gene"}},
            "gene_association": {"path": "burden.tsv"},
            "gene_annotation": {"path": "genes.tsv"},
        },
    }


def test_dataset_manifest_is_immutable_context_aware_and_path_bounded(tmp_path: Path) -> None:
    payload = _cohort_files(tmp_path / "imports")
    manifest = normalize_dataset_manifest(payload, tmp_path / "imports")
    assert manifest["method"] == "dandelion"
    assert manifest["applicability"]["status"] == "complete"
    assert manifest["files"]["trans_matrix"]["inspection"]["row_count"] == 2
    assert len(manifest["manifest_checksum_sha256"]) == 64

    payload["context"] = {}
    limited = normalize_dataset_manifest(payload, tmp_path / "imports")
    assert limited["applicability"] == {"status": "limited", "missing_context": ["tissue", "ancestry"]}

    payload["files"]["trans_matrix"]["path"] = "../outside.tsv"
    with pytest.raises(DandelionValidationError, match="relative"):
        normalize_dataset_manifest(payload, tmp_path / "imports")


def test_snp_mode_requires_reference_and_parameters_are_bounded(tmp_path: Path) -> None:
    payload = _cohort_files(tmp_path / "imports")
    payload["exposure_type"] = "SNP"
    with pytest.raises(DandelionValidationError, match="snp_reference"):
        normalize_dataset_manifest(payload, tmp_path / "imports")
    assert normalize_analysis_parameters({})["multiple_testing"].startswith("Benjamini-Hochberg")
    with pytest.raises(DandelionValidationError, match="target_fdr"):
        normalize_analysis_parameters({"target_fdr": 0})


def test_signed_queue_rejects_tampering_and_supports_cancellation(tmp_path: Path) -> None:
    key_file = tmp_path / "runner.key"
    key_file.write_bytes(b"k" * 48)
    manager = StatisticalJobManager(tmp_path / "work", key_file=key_file)
    state = manager.submit(analysis_id="a" * 32, dataset={"id": "dataset"}, parameters={})
    queue_path = manager.queue / f"{state['id']}.json"
    envelope = json.loads(queue_path.read_text(encoding="utf-8"))
    assert manager.verify_envelope(envelope)["analysis_id"] == "a" * 32
    envelope["payload"]["analysis_id"] = "b" * 32
    with pytest.raises(StatisticalJobError, match="authentication"):
        manager.verify_envelope(envelope)
    cancelled = manager.cancel(state["id"])
    assert cancelled and cancelled["cancel_requested"] is True
    assert (manager.jobs / state["id"] / "cancel.requested").is_file()


def test_report_keeps_dandelion_out_of_medical_and_predictions() -> None:
    report = build_report(
        analysis_id="a" * 32,
        dataset={
            "id": "dataset",
            "name": "Cohort",
            "phenotype": "EFO:test",
            "exposure_type": "Gene",
            "assembly": "GRCh38",
            "gene_namespace": "HGNC symbol",
            "context": {"tissue": "blood"},
            "applicability": {"status": "limited", "missing_context": ["ancestry"]},
            "manifest_checksum_sha256": "c" * 64,
        },
        parameters=normalize_analysis_parameters({}),
        normalized_result={
            "tested_pair_count": 20,
            "records": [
                {
                    "exposure": "E1",
                    "source_gene": "E1",
                    "candidate_gene": "G1",
                    "trans_p_value": 0.01,
                    "gene_association_p_value": 0.001,
                    "p_value": 0.002,
                    "q_value": 0.01,
                    "significant": True,
                }
            ],
        },
        job={"status": "completed"},
    )
    assert list(report["sections"]) == ["objective_data", "statistics", "literature", "medical", "interactions", "predictions"]
    assert report["sections"]["objective_data"]["status"] == "not_applicable"
    assert report["sections"]["medical"]["status"] == "not_assessed"
    assert report["sections"]["predictions"]["status"] == "not_requested"
    assert report["sections"]["interactions"]["edges"][0]["native_score_label"].startswith("BH q-value")
    assert report["summary"]["tested_count"] == 20


def _api_app(tmp_path: Path):
    import_root = tmp_path / "imports"
    manifest = _cohort_files(import_root)
    engine = create_database_engine(test_url="sqlite+pysqlite:///:memory:")
    ensure_schema(engine)
    runner_key = tmp_path / "runner.key"
    runner_key.write_bytes(b"r" * 48)
    artifact_key = tmp_path / "artifact.key"
    artifact_key.write_bytes(b"a" * 48)
    manager = StatisticalJobManager(tmp_path / "work", key_file=runner_key)
    app = Flask(__name__)
    app.config.update(
        TESTING=True,
        NOPHIGENE_DATABASE_ENGINE=engine,
        NOPHIGENE_STATISTICAL_JOB_MANAGER=manager,
        NOPHIGENE_DANDELION_IMPORT_ROOT=import_root,
        NOPHIGENE_DANDELION_ARTIFACT_KEY_FILE=artifact_key,
    )
    app.register_blueprint(api_v2)
    return app, engine, manager, manifest


def test_dataset_analysis_and_result_api_contract(tmp_path: Path) -> None:
    app, engine, manager, manifest = _api_app(tmp_path)
    client = app.test_client()
    registered = client.post("/api/v2/statistical-datasets", json=manifest)
    assert registered.status_code == 201
    dataset_id = registered.get_json()["id"]
    assert client.get("/api/v2/statistical-datasets").get_json()["count"] == 1

    submitted = client.post(
        "/api/v2/statistical-analyses",
        json={"method": "dandelion", "dataset_id": dataset_id, "parameters": {"target_fdr": 0.1}},
    )
    assert submitted.status_code == 202
    analysis_id = submitted.get_json()["id"]
    with Session(engine) as session:
        analysis = session.get(StatisticalAnalysisRun, analysis_id)
        assert analysis is not None
        job_id = analysis.worker_job_id

    result = {
        "contract_version": "1.0",
        "tested_pair_count": 2,
        "records": [
            {
                "exposure": "E1",
                "source_gene": "E1",
                "source_node_type": "gene",
                "candidate_gene": "G1",
                "trans_p_value": 0.01,
                "gene_association_p_value": 0.001,
                "p_value": 0.002,
                "q_value": 0.004,
                "significant": True,
            }
        ],
    }
    job_dir = manager.jobs / job_id
    (job_dir / "result.json").write_text(json.dumps(result), encoding="utf-8")
    state = manager.get(job_id)
    assert state is not None
    state.update(status="completed", stage="completed", progress_percent=100)
    (job_dir / "status.json").write_text(json.dumps(state), encoding="utf-8")

    response = client.get(f"/api/v2/statistical-analyses/{analysis_id}/result")
    assert response.status_code == 200
    report = response.get_json()
    assert report["schema_version"] == "3.0"
    assert report["sections"]["statistics"]["status"] == "assessed"
    assert report["sections"]["medical"]["status"] == "not_assessed"
    graph = client.get(f"/api/v2/statistical-analyses/{analysis_id}/interactions?node_cap=2").get_json()
    assert graph["node_count"] == 2
    assert graph["score_direction"] == "lower BH q-value is stronger"
    retry = client.post(f"/api/v2/statistical-analyses/{analysis_id}/retry", json={})
    assert retry.status_code == 202
    assert retry.get_json()["parent_analysis_id"] == analysis_id


def test_managed_copy_is_new_immutable_aes_dataset(tmp_path: Path) -> None:
    app, engine, _manager, manifest = _api_app(tmp_path)
    client = app.test_client()
    source = client.post("/api/v2/statistical-datasets", json=manifest).get_json()
    copied = client.post(f"/api/v2/statistical-datasets/{source['id']}/managed-copy", json={})
    assert copied.status_code == 201
    value = copied.get_json()
    assert value["id"] != source["id"]
    assert value["storage_mode"] == "managed_encrypted_copy"
    assert value["manifest_checksum_sha256"] != source["manifest_checksum_sha256"]
    archive = tmp_path / "work" / value["files"]["trans_matrix"]["managed_path"]
    with pyzipper.AESZipFile(archive) as bundle:
        bundle.setpassword(b"a" * 48)
        assert bundle.read(bundle.namelist()[0]).startswith(b"gene\tE1")
    with Session(engine) as session:
        assert session.get(StatisticalDataset, source["id"]).storage_mode == "registered_path"
