"""Prediction-contract, isolated-worker, and consent-gate coverage."""

from __future__ import annotations

import importlib.util
import json
import sys
import types
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from flask import Flask
from sqlalchemy.orm import Session

from src.api.jobs import JobManager
from src.api.profiles import ProfileStore
from src.api.serialization import utc_now
from src.api.v2_routes import api_v2, persist_completed_model_job
from src.workbench.database import create_database_engine, ensure_schema
from src.workbench.model_credentials import ModelCredentialStore
from src.workbench.model_jobs import ModelJobManager, canonical_json
from src.workbench.models import Artifact, ModelRun, Prediction, Run
from src.workbench.predictions import (
    observed_variant_alleles,
    prepare_alphagenome_disclosure,
    select_alphagenome_variants,
    validate_reference_alleles,
)
from src.workbench.reporting import build_canonical_report, render_evidence_first_html


def _write_fasta(path: Path, chromosome: str, sequence: str) -> None:
    header = f">{chromosome}\n"
    path.write_bytes(f"{header}{sequence}\n".encode("ascii"))
    Path(f"{path}.fai").write_bytes(
        f"{chromosome}\t{len(sequence)}\t{len(header)}\t{len(sequence)}\t{len(sequence) + 1}\n".encode("ascii")
    )


def _drd4_prediction_report() -> dict:
    return build_canonical_report(
        {
            "run_id": "1" * 32,
            "gene": "DRD4",
            "genome_build": "GRCh38",
            "scope_regions": {"promoter_only": "11:636000-637268", "gene_only": "11:637269-640706"},
            "variants": [
                {
                    "CHROM": "chr11",
                    "POS": 636433,
                    "REF": "G",
                    "ALT": "A",
                    "GT": "0/1",
                    "FILTER": "PASS",
                    "GQ": 99,
                    "DP": 40,
                }
            ],
            "methylation": [],
            "dynamic_knowledge_base": {
                "variant_records": [
                    {
                        "gene": "DRD4",
                        "reference_genome": "GRCh38",
                        "chromosome": "11",
                        "position": 636433,
                        "reference": "G",
                        "alternate": "A",
                        "rsid": "rs4987059",
                        "source_key": "gnomad",
                        "source": "gnomAD",
                        "source_id": "11-636433-G-A",
                        "release": "v4",
                        "in_silico_predictors": [
                            {"id": "CADD", "value": "6.46"},
                            {"id": "phyloP", "value": "-0.473"},
                        ],
                    },
                    {
                        "reference_genome": "GRCh38",
                        "chromosome": "11",
                        "position": 636433,
                        "reference": "G",
                        "alternate": "T",
                        "rsid": "rs4987059",
                        "source": "context-only",
                        "in_silico_predictors": [{"id": "CADD", "value": "99"}],
                    },
                ]
            },
        }
    )


def test_drd4_exact_allele_exposes_native_cadd_and_phylop_without_rsid_only_match() -> None:
    report = _drd4_prediction_report()
    section = report["sections"]["predictions"]
    values = {(item["predictor"], item["native_value"]) for item in section["source_native_annotations"]}
    assert values == {("CADD", 6.46), ("phyloP", -0.473)}
    assert section["counts"]["source_native_annotation_count"] == 2
    assert section["source_native_context_only_record_count"] == 1
    assert section["variant_selection"]["selected"][0]["rsid"] == "rs4987059"
    assert section["variant_selection"]["selected"][0]["exact_source_match"] is True
    assert section["consensus"] is None
    rendered = render_evidence_first_html(report)
    assert "6.46" in rendered and "-0.473" in rendered
    assert "pathogenicity classification" in rendered


def test_variant_selection_caps_at_twenty_and_reference_validation_is_exact(tmp_path: Path) -> None:
    multiallelic = observed_variant_alleles(
        [{"CHROM": "1", "POS": 2, "REF": "A", "ALT": "C,G", "GT": "0/2", "qc_pass": True, "non_reference": True}],
        genome_build="GRCh38",
    )
    assert [item["alternate"] for item in multiallelic] == ["G"]
    alleles = [
        {
            "entity_key": f"GRCh38|1|{index + 1}|A|G",
            "assembly": "GRCh38",
            "chromosome": "chr1",
            "position": index + 1,
            "reference": "A",
            "alternate": "G",
            "variant": f"chr1:{index + 1}:A>G",
            "curated_or_named": index == 24,
            "region": "gene_body",
            "gq": 20 + index,
            "dp": 10,
        }
        for index in range(25)
    ]
    selection = select_alphagenome_variants(alleles)
    assert len(selection["selected"]) == 20
    assert selection["selected"][0]["position"] == 25
    assert len(selection["omitted"]) == 5
    assert {item["omission_reason"] for item in selection["omitted"]} == {"maximum_variant_limit"}
    unsupported = select_alphagenome_variants(
        [{**alleles[0], "alternate": "<DEL>", "variant": "chr1:1:A><DEL>"}]
    )
    assert unsupported["status"] == "no_data"
    assert unsupported["omitted"][0]["omission_reason"] == "unsupported_variant_representation"

    fasta = tmp_path / "hg38.fa"
    _write_fasta(fasta, "chr1", "A" * 30)
    validated = validate_reference_alleles(selection, fasta)
    assert validated["status"] == "available"
    assert all(item["reference_allele_verified"] for item in validated["selected"])
    validated["selected"][0]["reference"] = "C"
    mismatch = validate_reference_alleles(validated, fasta)
    assert mismatch["status"] == "blocked"
    assert any(item.startswith("reference_mismatch:") for item in mismatch["blockers"])


def test_disclosure_validates_settings_and_has_stable_confirmation_digest() -> None:
    selection = {
        "selected": [
            {
                "variant": "chr1:10:A>G",
                "chromosome": "chr1",
                "position": 10,
                "reference": "A",
                "alternate": "G",
                "reference_allele_verified": True,
            }
        ],
        "omitted": [],
        "blockers": [],
    }
    first = prepare_alphagenome_disclosure(
        selection, ontology_terms=["UBERON:0000178"], modalities=["RNA_SEQ", "ATAC"], sequence_length=1_048_576
    )
    second = prepare_alphagenome_disclosure(
        selection, ontology_terms=["UBERON:0000178"], modalities=["RNA_SEQ", "ATAC"], sequence_length=1_048_576
    )
    assert first["status"] == "ready"
    assert first["payload_sha256"] == second["payload_sha256"]
    blocked = prepare_alphagenome_disclosure(
        selection, ontology_terms=[], modalities=["UNSUPPORTED"], sequence_length=1000
    )
    assert blocked["status"] == "blocked"
    assert set(blocked["blockers"]) == {
        "ontology_terms_must_contain_1_to_5_items",
        "unsupported_or_missing_modalities",
        "unsupported_sequence_length",
    }


def test_model_credential_is_encrypted_and_job_manifest_is_signed_and_idempotent(tmp_path: Path) -> None:
    key_file = tmp_path / "credential.key"
    key_file.write_bytes(b"credential-encryption-key-contains-more-than-thirty-two-bytes")
    store = ModelCredentialStore(tmp_path / "credentials", key_file=key_file)
    state = store.save("alphagenome-api", "top-secret-provider-key")
    assert state["status"] == "configured"
    assert b"top-secret-provider-key" not in (tmp_path / "credentials" / "alphagenome-api.aes.zip").read_bytes()
    assert store.read("alphagenome-api") == "top-secret-provider-key"
    assert "top-secret-provider-key" not in json.dumps(store.status("alphagenome-api"))

    signing_key = tmp_path / "runner.key"
    signing_key.write_bytes(b"model-runner-signing-key-contains-more-than-thirty-two-bytes")
    manager = ModelJobManager(tmp_path / "model-jobs", key_file=signing_key)
    inputs = {"operation": "metadata", "credential_revision": state["updated_at"]}
    first = manager.submit(run_id="model-settings", model_id="alphagenome-api", inputs=inputs)
    second = manager.submit(run_id="model-settings", model_id="alphagenome-api", inputs=inputs)
    assert first["id"] == second["id"]
    envelope = json.loads((manager.queue / f"{first['id']}.json").read_text(encoding="utf-8"))
    assert manager.verify_envelope(envelope)["inputs"] == inputs
    assert "top-secret-provider-key" not in json.dumps(envelope)
    assert manager.cancel(first["id"])["status"] == "queued"
    replacement = manager.submit(run_id="model-settings", model_id="alphagenome-api", inputs=inputs)
    assert replacement["id"] != first["id"]


def _load_worker(monkeypatch, tmp_path: Path):
    worker_dir = Path(__file__).parents[1] / "docker" / "alphagenome"
    monkeypatch.syspath_prepend(str(worker_dir))
    spec = importlib.util.spec_from_file_location("test_alphagenome_worker", worker_dir / "worker.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    module.WORK_ROOT = tmp_path / "work"
    module.QUEUE = module.WORK_ROOT / "queue"
    module.JOBS = module.WORK_ROOT / "jobs"
    module.QUEUE.mkdir(parents=True)
    module.JOBS.mkdir(parents=True)
    module.version = lambda _name: "0.8.0"
    module.api_key = lambda: "fake-key"
    return module


def test_alphagenome_worker_fake_sdk_preserves_raw_quantile_and_partial_failures(tmp_path: Path, monkeypatch) -> None:
    worker = _load_worker(monkeypatch, tmp_path)
    class SpliceJunctionScorer:
        def to_proto(self):
            return "splice_junction {}"

    assert worker.scorer_matches(SpliceJunctionScorer(), {"SPLICE_JUNCTIONS"}) is True
    assert worker.tidy_rows(None, requested_modalities={"RNA_SEQ"}, ontology_terms=set()) == []
    fasta = tmp_path / "hg38.fa"
    _write_fasta(fasta, "chr1", "AGAAAAAAAAAAAAAAAAAA")
    worker.REFERENCE_FASTA = fasta
    job_id = "a" * 32
    job_dir = worker.JOBS / job_id
    job_dir.mkdir()
    (job_dir / "status.json").write_bytes(canonical_json({"id": job_id, "status": "running"}))

    class FakeScorer:
        def to_proto(self):
            return "RNA_SEQ"

    class FakeClient:
        def score_variant(self, **_kwargs):
            return object()

    dna_client = types.ModuleType("alphagenome.models.dna_client")
    dna_client.Organism = SimpleNamespace(HOMO_SAPIENS=SimpleNamespace(to_proto=lambda: "human"))
    dna_client.create = lambda *_args, **_kwargs: FakeClient()
    variant_scorers = types.ModuleType("alphagenome.models.variant_scorers")
    variant_scorers.get_recommended_scorers = lambda _organism: [FakeScorer()]
    variant_scorers.tidy_scores = lambda _scores: pd.DataFrame(
        [
            {
                "variant_id": "chr1:2:G>A",
                "output_type": "RNA_SEQ",
                "variant_scorer": "RNA scorer",
                "track_name": "track-1",
                "gene_id": "ENSG1",
                "gene_name": "GENE1",
                "ontology_curie": "UBERON:0000178",
                "biosample_name": "blood",
                "raw_score": -0.25,
                "quantile_score": 0.91,
            }
        ]
    )
    genome = types.ModuleType("alphagenome.data.genome")
    genome.Interval = lambda *args: args
    genome.Variant = lambda **kwargs: kwargs
    packages = {
        "alphagenome": types.ModuleType("alphagenome"),
        "alphagenome.data": types.ModuleType("alphagenome.data"),
        "alphagenome.data.genome": genome,
        "alphagenome.models": types.ModuleType("alphagenome.models"),
        "alphagenome.models.dna_client": dna_client,
        "alphagenome.models.variant_scorers": variant_scorers,
    }
    packages["alphagenome.data"].genome = genome
    packages["alphagenome.models"].dna_client = dna_client
    packages["alphagenome.models"].variant_scorers = variant_scorers
    for name, value in packages.items():
        monkeypatch.setitem(sys.modules, name, value)

    def variant(position: int, reference: str) -> dict:
        length = 131_072
        start = max(0, position - 1 - length // 2)
        return {
            "variant": f"chr1:{position}:{reference}>A",
            "assembly": "GRCh38",
            "chromosome": "chr1",
            "position_1_based": position,
            "reference": reference,
            "alternate": "A",
            "reference_allele_verified": True,
            "model_interval_0_based_half_open": {"start": start, "end": start + length, "length": length},
        }

    result = worker.execute(
        {
            "contract_version": "1.0",
            "job_id": job_id,
            "model_id": "alphagenome-api",
            "inputs": {
                "operation": "predictions",
                "requested_modalities": ["RNA_SEQ"],
                "ontology_terms": ["UBERON:0000178"],
                "sequence_length": 131_072,
                "explicit_transfer_consent": True,
                "variants": [variant(2, "G"), variant(3, "G")],
            },
        }
    )
    assert len(result["predictions"]) == 1
    assert result["predictions"][0]["raw_score"] == -0.25
    assert result["predictions"][0]["quantile_score"] == 0.91
    assert result["predictions"][0]["ontology_curie"] == "UBERON:0000178"
    assert result["failures"][0]["code"] == "worker_reference_mismatch"
    runner_key = tmp_path / "runner.key"
    runner_key.write_bytes(b"model-runner-signing-key-contains-more-than-thirty-two-bytes")
    monkeypatch.setenv("NOPHIGENE_MODEL_RUNNER_KEY_FILE", str(runner_key))
    recovery_result = {**result, "predictions": result["predictions"], "failures": []}
    import hashlib

    recovery_result["result_checksum_sha256"] = hashlib.sha256(canonical_json(recovery_result)).hexdigest()
    (job_dir / "result.json").write_bytes(canonical_json(worker.sign(recovery_result)))
    assert worker.recover_completed_result(job_id) is True
    recovered_status = json.loads((job_dir / "status.json").read_text(encoding="utf-8"))
    assert recovered_status["status"] == "succeeded"
    assert recovered_status["stage"] == "completed_after_worker_restart"


def test_prediction_api_requires_preview_hash_and_consent_and_never_serializes_credential(tmp_path: Path) -> None:
    profile_store = ProfileStore(tmp_path / "profiles.json")
    workflow_jobs = JobManager(jobs_root=tmp_path / "jobs", profile_store=profile_store)
    model_key = tmp_path / "model.key"
    model_key.write_bytes(b"model-runner-signing-key-contains-more-than-thirty-two-bytes")
    model_jobs = ModelJobManager(tmp_path / "model-jobs", key_file=model_key)
    credential_key = tmp_path / "credential.key"
    credential_key.write_bytes(b"credential-encryption-key-contains-more-than-thirty-two-bytes")
    credentials = ModelCredentialStore(tmp_path / "credentials", key_file=credential_key)
    credentials.save("alphagenome-api", "top-secret-provider-key")
    credentials.mark_status("alphagenome-api", "verified")
    (credentials.root / "alphagenome-api.metadata.json").write_text(
        json.dumps({"ontology_terms": [{"ontology_curie": "UBERON:0000178", "biosample_name": "blood"}]}),
        encoding="utf-8",
    )
    (model_jobs.root / "worker-heartbeat.json").write_text(
        json.dumps({"status": "ready", "updated_at": utc_now(), "sdk_version": "0.8.0"}), encoding="utf-8"
    )
    fasta = tmp_path / "hg38.fa"
    _write_fasta(fasta, "chr1", "AGA")

    run_id = "b" * 32
    report = build_canonical_report(
        {
            "run_id": run_id,
            "gene": "GENE1",
            "genome_build": "GRCh38",
            "region": "1:1-3",
            "variants": [{"CHROM": "1", "POS": 2, "REF": "G", "ALT": "A", "GT": "0/1", "FILTER": "PASS", "GQ": 50, "DP": 20}],
            "methylation": [],
        }
    )
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    engine = create_database_engine(test_url=f"sqlite+pysqlite:///{(tmp_path / 'db.sqlite').as_posix()}")
    ensure_schema(engine)
    with Session(engine) as database:
        database.add(
            Run(
                id=run_id,
                status="succeeded",
                stage="analysis",
                progress_percent=100,
                genes=["GENE1"],
                configuration={"kind": "single_person_gene_analysis", "canonical_report_path": str(report_path)},
                started_at=datetime.now(timezone.utc),
                finished_at=datetime.now(timezone.utc),
            )
        )
        database.commit()

    app = Flask(__name__)
    app.config.update(
        TESTING=True,
        NOPHIGENE_PROFILE_STORE=profile_store,
        NOPHIGENE_JOB_MANAGER=workflow_jobs,
        NOPHIGENE_MODEL_JOB_MANAGER=model_jobs,
        NOPHIGENE_MODEL_CREDENTIAL_STORE=credentials,
        NOPHIGENE_MODEL_ARTIFACT_ROOT=str(tmp_path / "artifacts"),
        NOPHIGENE_DATABASE_ENGINE=engine,
        NOPHIGENE_HG38_REFERENCE_FASTA=str(fasta),
    )
    app.register_blueprint(api_v2)
    client = app.test_client()
    settings = {"ontology_terms": ["UBERON:0000178"], "modalities": ["RNA_SEQ"], "sequence_length": 131_072}
    preview_response = client.post(f"/api/v2/runs/{run_id}/predictions/preview", json=settings)
    assert preview_response.status_code == 200, preview_response.get_json()
    preview = preview_response.get_json()
    assert preview["status"] == "ready", preview
    assert preview["payload"]["variants"][0]["reference_allele_verified"] is True

    missing_consent = client.post(
        f"/api/v2/runs/{run_id}/predictions", json={**settings, "payload_sha256": preview["payload_sha256"]}
    )
    assert missing_consent.status_code == 422
    stale = client.post(
        f"/api/v2/runs/{run_id}/predictions",
        json={**settings, "payload_sha256": "0" * 64, "external_transfer_consent": True},
    )
    assert stale.status_code == 409
    submitted = client.post(
        f"/api/v2/runs/{run_id}/predictions",
        json={**settings, "payload_sha256": preview["payload_sha256"], "external_transfer_consent": True},
    )
    assert submitted.status_code == 202
    duplicate = client.post(
        f"/api/v2/runs/{run_id}/predictions",
        json={**settings, "payload_sha256": preview["payload_sha256"], "external_transfer_consent": True},
    )
    assert duplicate.get_json()["id"] == submitted.get_json()["id"]
    manifest_text = (model_jobs.queue / f"{submitted.get_json()['id']}.json").read_text(encoding="utf-8")
    assert "top-secret-provider-key" not in manifest_text
    legacy = client.post(f"/api/v2/runs/{run_id}/models/alphagenome-api", json={})
    assert legacy.status_code == 409

    job_id = submitted.get_json()["id"]
    result = {
        "kind": "predictions",
        "model_id": "alphagenome-api",
        "sdk_version": "0.8.0",
        "sdk_source_commit": "71a6beb8c30832f121309a81c2530efa5af7986a",
        "provider_model_revision": None,
        "retrieved_at": utc_now(),
        "external_transfer_occurred": True,
        "predictions": [
            {
                "entity_type": "observed_variant",
                "entity_key": "chr1:2:G>A",
                "output_name": "RNA_SEQ:stable",
                "variant": "chr1:2:G>A",
                "output_type": "RNA_SEQ",
                "raw_score": -0.2,
                "quantile_score": 0.8,
                "scorer": "RNA scorer",
                "track_name": "track",
                "ontology_curie": "UBERON:0000178",
                "limitations": ["molecular research output"],
            }
        ],
        "failures": [],
    }
    import hashlib

    result["result_checksum_sha256"] = hashlib.sha256(canonical_json(result)).hexdigest()
    (model_jobs.jobs / job_id / "result.json").write_bytes(canonical_json(model_jobs._signed_envelope(result)))
    status_path = model_jobs.jobs / job_id / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status.update(status="succeeded", stage="completed", progress_percent=100, finished_at=utc_now())
    status_path.write_bytes(canonical_json(status))
    persist_completed_model_job(app, model_jobs, job_id)
    persist_completed_model_job(app, model_jobs, job_id)
    with Session(engine) as database:
        model_run = database.get(ModelRun, job_id)
        assert model_run is not None and model_run.status == "succeeded"
        assert database.query(Prediction).filter_by(model_run_id=job_id).count() == 1
        artifact = database.query(Artifact).filter_by(run_id=run_id, kind="model_prediction_raw").one()
        assert artifact.sensitive is True
    assert all(
        "top-secret-provider-key" not in path.read_text(encoding="utf-8")
        for path in (tmp_path / "artifacts").rglob("*")
        if path.is_file()
    )
