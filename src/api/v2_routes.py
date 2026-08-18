"""NophiGene schema-v3 REST API exposed under /api/v2."""

from __future__ import annotations

import json
import hashlib
import shutil
import uuid
import os
from pathlib import Path
from typing import Any

import pyzipper
from flask import Blueprint, current_app, jsonify, request, send_file
from sqlalchemy import select
from werkzeug.exceptions import HTTPException

from .errors import APIError
from .serialization import read_json, utc_now
from .workflow_runner import MAX_GENES_PER_JOB, normalize_job_request

try:
    from ..variant_knowledge.registry import list_source_specs
    from ..workbench.audit import append_audit_event
    from ..workbench.database import get_default_engine, session_scope
    from ..workbench.dandelion import (
        DandelionValidationError,
        PACKAGE_VERSION as DANDELION_VERSION,
        build_report as build_dandelion_report,
        estimate_resources as estimate_dandelion_resources,
        normalize_analysis_parameters,
        normalize_dataset_manifest,
        sha256_file,
        validate_runner_result,
    )
    from ..workbench.evidence import CORE_SOURCE_KEYS, list_core_source_metadata
    from ..workbench.exports import create_export_bundle
    from ..workbench.interactions import build_interaction_graph
    from ..workbench.model_registry import (
        estimate_model_resources,
        inspect_model_inputs,
        list_model_manifests,
        model_installation_plan,
    )
    from ..workbench.models import (
        DeletionRecord,
        ModelRun,
        Prediction,
        Run,
        Sample,
        StatisticalAnalysisRun,
        StatisticalDataset,
        StatisticalDatasetFile,
        StatisticalResult,
    )
    from ..workbench.persistence import persist_canonical_report
except ImportError:
    from variant_knowledge.registry import list_source_specs
    from workbench.audit import append_audit_event
    from workbench.database import get_default_engine, session_scope
    from workbench.dandelion import (
        DandelionValidationError,
        PACKAGE_VERSION as DANDELION_VERSION,
        build_report as build_dandelion_report,
        estimate_resources as estimate_dandelion_resources,
        normalize_analysis_parameters,
        normalize_dataset_manifest,
        sha256_file,
        validate_runner_result,
    )
    from workbench.evidence import CORE_SOURCE_KEYS, list_core_source_metadata
    from workbench.exports import create_export_bundle
    from workbench.interactions import build_interaction_graph
    from workbench.model_registry import (
        estimate_model_resources,
        inspect_model_inputs,
        list_model_manifests,
        model_installation_plan,
    )
    from workbench.models import (
        DeletionRecord,
        ModelRun,
        Prediction,
        Run,
        Sample,
        StatisticalAnalysisRun,
        StatisticalDataset,
        StatisticalDatasetFile,
        StatisticalResult,
    )
    from workbench.persistence import persist_canonical_report

api_v2 = Blueprint("api_v2", __name__, url_prefix="/api/v2")


def _jobs():
    return current_app.config["NOPHIGENE_JOB_MANAGER"]


def _profiles():
    return current_app.config["NOPHIGENE_PROFILE_STORE"]


def _model_jobs():
    return current_app.config["NOPHIGENE_MODEL_JOB_MANAGER"]


def _statistical_jobs():
    return current_app.config["NOPHIGENE_STATISTICAL_JOB_MANAGER"]


def _dandelion_import_root() -> Path:
    configured = current_app.config.get("NOPHIGENE_DANDELION_IMPORT_ROOT") or os.environ.get(
        "NOPHIGENE_DANDELION_IMPORT_ROOT", "data/dandelion"
    )
    return Path(configured)


def _dandelion_artifact_key() -> bytes:
    configured = current_app.config.get("NOPHIGENE_DANDELION_ARTIFACT_KEY_FILE") or os.environ.get(
        "NOPHIGENE_DANDELION_ARTIFACT_KEY_FILE"
    )
    if not configured:
        raise APIError("managed_copy_unavailable", "The managed-copy encryption key is not configured.", 503)
    key = Path(configured).read_bytes().strip()
    if len(key) < 32:
        raise APIError("managed_copy_unavailable", "The managed-copy encryption key is invalid.", 503)
    return key


def _engine():
    return current_app.config.get("NOPHIGENE_DATABASE_ENGINE") or get_default_engine()


def _json_body() -> dict[str, Any]:
    payload = request.get_json(silent=True)
    if not isinstance(payload, dict):
        raise APIError("invalid_json", "Request body must contain a JSON object.", 400)
    return payload


@api_v2.errorhandler(APIError)
def handle_api_error(error: APIError):
    return error.to_response()


@api_v2.errorhandler(Exception)
def handle_unexpected_error(error: Exception):
    if isinstance(error, HTTPException):
        return APIError("http_error", error.description, error.code or 500).to_response()
    current_app.logger.exception("Unhandled API v2 error", exc_info=error)
    return APIError("internal_error", "The API could not complete the request.", 500).to_response()


@api_v2.get("")
@api_v2.get("/")
def index():
    return jsonify(
        {
            "name": "NophiGene Evidence-First API",
            "version": "2.0",
            "report_schema_version": "3.0",
            "runs": "/api/v2/runs",
            "evidence_sources": "/api/v2/evidence/sources",
            "models": "/api/v2/models",
            "statistical_datasets": "/api/v2/statistical-datasets",
            "statistical_analyses": "/api/v2/statistical-analyses",
            "health": "/api/v2/health",
            "openapi": "/api/v2/openapi.json",
        }
    )


@api_v2.get("/openapi.json")
def openapi_document():
    operations = {
        "/api/v2/runs": {"post": "Submit a one-to-100-gene run", "get": "List runs"},
        "/api/v2/runs/{id}": {"get": "Read run state, progress, blockers, and partial status"},
        "/api/v2/runs/{id}/result": {"get": "Read canonical report schema 3.0"},
        "/api/v2/runs/{id}/retry": {"post": "Create an immutable targeted retry"},
        "/api/v2/runs/{id}/evidence/payload-preview": {"get": "Preview exact external disclosure payload"},
        "/api/v2/runs/{id}/evidence/refresh": {"post": "Explicitly refresh selected evidence sources"},
        "/api/v2/runs/{id}/interactions": {"get": "Expand one to three gene hops"},
        "/api/v2/runs/{id}/exports": {"post": "Create a sensitive full-fidelity export"},
        "/api/v2/evidence/sources": {"get": "List the typed vetted source catalog"},
        "/api/v2/models": {"get": "List immutable model manifests"},
        "/api/v2/statistical-datasets": {
            "post": "Register an immutable local cohort dataset",
            "get": "List cohort datasets",
        },
        "/api/v2/statistical-datasets/{id}": {"get": "Read an immutable cohort dataset manifest"},
        "/api/v2/statistical-datasets/{id}/managed-copy": {
            "post": "Create a new immutable dataset backed by AES-encrypted managed copies"
        },
        "/api/v2/statistical-analyses": {
            "post": "Queue an offline cohort statistical analysis",
            "get": "List cohort statistical analyses",
        },
        "/api/v2/statistical-analyses/{id}": {"get": "Read cohort analysis state and resource estimate"},
        "/api/v2/statistical-analyses/{id}/cancel": {"post": "Request cancellation"},
        "/api/v2/statistical-analyses/{id}/retry": {"post": "Create a new immutable retry job"},
        "/api/v2/statistical-analyses/{id}/result": {"get": "Read canonical report schema 3.0"},
        "/api/v2/statistical-analyses/{id}/interactions": {
            "get": "Read typed DANDELION prioritisation edges without treating them as physical interactions"
        },
        "/api/v2/models/{id}/eligibility": {"post": "Inspect input eligibility"},
        "/api/v2/models/{id}/estimate": {"post": "Estimate execution resources"},
        "/api/v2/models/{id}/installation": {"get": "Read installation gate", "post": "Confirm installation manifest"},
        "/api/v2/runs/{run_id}/models/{model_id}": {"post": "Queue an independent model job"},
        "/api/v2/model-jobs/{id}": {"get": "Read model job state"},
        "/api/v2/model-jobs/{id}/cancel": {"post": "Cancel a queued model job"},
        "/api/v2/model-jobs/{id}/result": {"get": "Read independent normalized predictions"},
        "/api/v2/samples/{id}": {"delete": "Delete a sample and sample-linked derivatives"},
    }
    return jsonify(
        {
            "openapi": "3.1.0",
            "info": {"title": "NophiGene Evidence-First API", "version": "2.0.0"},
            "servers": [{"url": "http://127.0.0.1:8766"}],
            "paths": {
                path: {
                    method: {
                        "summary": summary,
                        "responses": {"200": {"description": "Success"}},
                    }
                    for method, summary in methods.items()
                }
                for path, methods in operations.items()
            },
            "components": {
                "schemas": {
                    "CanonicalReport": {
                        "type": "object",
                        "required": ["schema_version", "run", "summary", "sections", "run_details"],
                        "properties": {
                            "schema_version": {"const": "3.0"},
                            "run": {"type": "object"},
                            "summary": {"type": "object"},
                            "sections": {
                                "type": "object",
                                "required": ["objective_data", "statistics", "literature", "medical", "interactions", "predictions"],
                            },
                            "run_details": {"type": "object"},
                        },
                    }
                }
            },
        }
    )


@api_v2.get("/profiles")
def list_profiles():
    profiles = _profiles().list()
    return jsonify({"profiles": profiles, "count": len(profiles)})


@api_v2.post("/profiles")
def create_profile():
    profile = _profiles().create(_json_body())
    return jsonify(profile), 201


def _dataset_payload(session, dataset: StatisticalDataset, *, include_files: bool = True) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "id": dataset.id,
        "method": dataset.method,
        "name": dataset.name,
        "description": dataset.description,
        "phenotype": dataset.phenotype,
        "exposure_type": dataset.exposure_type,
        "assembly": dataset.assembly,
        "gene_namespace": dataset.gene_namespace,
        "storage_mode": dataset.storage_mode,
        "context": dataset.context or {},
        "applicability": dataset.applicability or {},
        "manifest_checksum_sha256": dataset.manifest_checksum_sha256,
        "created_at": dataset.created_at.isoformat() if dataset.created_at else None,
    }
    if include_files:
        files = session.scalars(
            select(StatisticalDatasetFile)
            .where(StatisticalDatasetFile.dataset_id == dataset.id)
            .order_by(StatisticalDatasetFile.role)
        ).all()
        payload["files"] = {
            item.role: {
                "role": item.role,
                "relative_path": item.relative_path,
                "managed_path": item.managed_path,
                "format": item.file_format,
                "checksum_sha256": item.checksum_sha256,
                "size_bytes": item.size_bytes,
                "mapping": item.mapping or {},
                "inspection": item.inspection or {},
            }
            for item in files
        }
    return payload


@api_v2.post("/statistical-datasets")
def register_statistical_dataset():
    payload = _json_body()
    try:
        manifest = normalize_dataset_manifest(payload, _dandelion_import_root())
    except (DandelionValidationError, OSError) as exc:
        raise APIError("invalid_statistical_dataset", str(exc), 422) from exc
    if manifest["storage_mode"] != "registered_path":
        raise APIError(
            "managed_copy_requires_import",
            "Register the dataset by path first; an encrypted managed copy is created through the dataset import action.",
            422,
        )
    with session_scope(_engine()) as session:
        existing = session.scalar(
            select(StatisticalDataset).where(
                StatisticalDataset.manifest_checksum_sha256 == manifest["manifest_checksum_sha256"]
            )
        )
        if existing is not None:
            result = _dataset_payload(session, existing)
            result["already_registered"] = True
            return jsonify(result), 200
        dataset = StatisticalDataset(
            name=manifest["name"],
            description=manifest["description"],
            method=manifest["method"],
            phenotype=manifest["phenotype"],
            exposure_type=manifest["exposure_type"],
            assembly=manifest["assembly"],
            gene_namespace=manifest["gene_namespace"],
            storage_mode=manifest["storage_mode"],
            context=manifest["context"],
            applicability=manifest["applicability"],
            manifest_checksum_sha256=manifest["manifest_checksum_sha256"],
        )
        session.add(dataset)
        session.flush()
        for role, item in manifest["files"].items():
            session.add(
                StatisticalDatasetFile(
                    dataset_id=dataset.id,
                    role=role,
                    relative_path=item["relative_path"],
                    file_format=item["format"],
                    checksum_sha256=item["checksum_sha256"],
                    size_bytes=item["size_bytes"],
                    mapping=item["mapping"],
                    inspection=item["inspection"],
                )
            )
        append_audit_event(
            session,
            "statistical_dataset_registered",
            entity_type="statistical_dataset",
            entity_id=dataset.id,
            payload={
                "method": "dandelion",
                "manifest_checksum_sha256": manifest["manifest_checksum_sha256"],
                "file_checksums": {role: item["checksum_sha256"] for role, item in manifest["files"].items()},
            },
        )
        session.flush()
        result = _dataset_payload(session, dataset)
    return jsonify(result), 201


@api_v2.get("/statistical-datasets")
def list_statistical_datasets():
    with session_scope(_engine()) as session:
        datasets = session.scalars(select(StatisticalDataset).order_by(StatisticalDataset.created_at.desc())).all()
        values = [_dataset_payload(session, item, include_files=False) for item in datasets]
    return jsonify({"datasets": values, "count": len(values)})


@api_v2.get("/statistical-datasets/<dataset_id>")
def get_statistical_dataset(dataset_id: str):
    with session_scope(_engine()) as session:
        dataset = session.get(StatisticalDataset, dataset_id)
        if dataset is None:
            raise APIError("statistical_dataset_not_found", f"Dataset '{dataset_id}' was not found.", 404)
        result = _dataset_payload(session, dataset)
    return jsonify(result)


@api_v2.post("/statistical-datasets/<dataset_id>/managed-copy")
def create_statistical_dataset_managed_copy(dataset_id: str):
    key = _dandelion_artifact_key()
    new_id = uuid.uuid4().hex
    with session_scope(_engine()) as session:
        source = session.get(StatisticalDataset, dataset_id)
        if source is None:
            raise APIError("statistical_dataset_not_found", f"Dataset '{dataset_id}' was not found.", 404)
        if source.storage_mode != "registered_path":
            raise APIError("managed_copy_exists", "This dataset is already backed by managed encrypted copies.", 409)
        source_payload = _dataset_payload(session, source)

    managed_root = (_statistical_jobs().root / "managed" / new_id).resolve()
    work_root = _statistical_jobs().root.resolve()
    managed_root.relative_to(work_root)
    managed_root.mkdir(parents=True, exist_ok=False)
    managed_files: dict[str, dict[str, Any]] = {}
    try:
        import_root = _dandelion_import_root().resolve(strict=True)
        for role, item in source_payload["files"].items():
            relative = Path(item["relative_path"])
            current = import_root
            for part in relative.parts:
                current = current / part
                if current.is_symlink():
                    raise APIError("managed_copy_source_invalid", f"Source path for '{role}' contains a symbolic link.", 422)
            source_path = (import_root / relative).resolve(strict=True)
            source_path.relative_to(import_root)
            if not source_path.is_file():
                raise APIError("managed_copy_source_invalid", f"Source file for '{role}' is unavailable.", 422)
            if sha256_file(source_path) != item["checksum_sha256"]:
                raise APIError("managed_copy_source_changed", f"Source file for '{role}' changed after registration.", 409)
            archive = managed_root / f"{role}.zip"
            temporary = managed_root / f".{role}.tmp"
            with pyzipper.AESZipFile(
                temporary,
                "w",
                compression=pyzipper.ZIP_DEFLATED,
                encryption=pyzipper.WZ_AES,
            ) as bundle:
                bundle.setpassword(key)
                bundle.setencryption(pyzipper.WZ_AES, nbits=256)
                bundle.write(source_path, arcname=f"payload.{item['format']}")
            temporary.replace(archive)
            inspection = dict(item.get("inspection") or {})
            inspection["managed_archive_sha256"] = sha256_file(archive)
            managed_files[role] = {
                **item,
                "managed_path": archive.relative_to(work_root).as_posix(),
                "inspection": inspection,
            }
    except Exception:
        shutil.rmtree(managed_root, ignore_errors=True)
        raise

    managed_checksum = hashlib.sha256(
        json.dumps(
            {
                "source_manifest_checksum_sha256": source_payload["manifest_checksum_sha256"],
                "storage_mode": "managed_encrypted_copy",
                "archives": {
                    role: item["inspection"]["managed_archive_sha256"] for role, item in managed_files.items()
                },
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    with session_scope(_engine()) as session:
        copied = StatisticalDataset(
            id=new_id,
            name=f"{source_payload['name']} (managed copy)",
            description=source_payload["description"],
            method="dandelion",
            phenotype=source_payload["phenotype"],
            exposure_type=source_payload["exposure_type"],
            assembly=source_payload["assembly"],
            gene_namespace=source_payload["gene_namespace"],
            storage_mode="managed_encrypted_copy",
            context=source_payload["context"],
            applicability=source_payload["applicability"],
            manifest_checksum_sha256=managed_checksum,
        )
        session.add(copied)
        for role, item in managed_files.items():
            session.add(
                StatisticalDatasetFile(
                    dataset_id=new_id,
                    role=role,
                    relative_path=item["relative_path"],
                    managed_path=item["managed_path"],
                    file_format=item["format"],
                    checksum_sha256=item["checksum_sha256"],
                    size_bytes=item["size_bytes"],
                    mapping=item["mapping"],
                    inspection=item["inspection"],
                )
            )
        append_audit_event(
            session,
            "statistical_dataset_managed_copy_created",
            entity_type="statistical_dataset",
            entity_id=new_id,
            payload={"source_dataset_id": dataset_id, "manifest_checksum_sha256": managed_checksum, "encryption": "AES-256"},
        )
        session.flush()
        result = _dataset_payload(session, copied)
    return jsonify(result), 201


def _analysis_payload(session, analysis: StatisticalAnalysisRun, run: Run, state: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": analysis.id,
        "run_id": run.id,
        "dataset_id": analysis.dataset_id,
        "method": analysis.method,
        "method_version": analysis.method_version,
        "status": state.get("status", analysis.status),
        "stage": state.get("stage", run.stage),
        "progress_percent": int(state.get("progress_percent") or 0),
        "parameters": analysis.parameters or {},
        "resource_estimate": analysis.resource_estimate or {},
        "cancel_requested": bool(state.get("cancel_requested")),
        "submitted_at": state.get("submitted_at"),
        "updated_at": state.get("updated_at"),
        "result_url": f"/api/v2/statistical-analyses/{analysis.id}/result",
    }


def _load_analysis(analysis_id: str) -> tuple[dict[str, Any], dict[str, Any]]:
    with session_scope(_engine()) as session:
        analysis = session.get(StatisticalAnalysisRun, analysis_id)
        if analysis is None:
            raise APIError("statistical_analysis_not_found", f"Analysis '{analysis_id}' was not found.", 404)
        run = session.get(Run, analysis.run_id)
        dataset = session.get(StatisticalDataset, analysis.dataset_id)
        if run is None or dataset is None:
            raise APIError("statistical_analysis_invalid", "Analysis persistence is incomplete.", 500)
        state = _statistical_jobs().get(analysis.worker_job_id)
        if state is None:
            raise APIError("statistical_job_not_found", "The offline worker job is unavailable.", 500)
        analysis.status = str(state.get("status") or analysis.status)
        run.status = analysis.status
        run.stage = str(state.get("stage") or run.stage)
        run.progress_percent = int(state.get("progress_percent") or 0)
        if analysis.status in {"failed", "cancelled"}:
            analysis.error = dict(state.get("error") or {})
            run.error = analysis.error
        analysis_payload = _analysis_payload(session, analysis, run, state)
        dataset_payload = _dataset_payload(session, dataset)
    return analysis_payload, dataset_payload


@api_v2.post("/statistical-analyses")
def submit_statistical_analysis():
    payload = _json_body()
    if str(payload.get("method") or "dandelion").lower() != "dandelion":
        raise APIError("unsupported_statistical_method", "Only the allowlisted DANDELION method is available.", 422)
    dataset_id = str(payload.get("dataset_id") or "")
    try:
        parameters = normalize_analysis_parameters(payload.get("parameters"))
    except DandelionValidationError as exc:
        raise APIError("invalid_analysis_parameters", str(exc), 422) from exc
    analysis_id = uuid.uuid4().hex
    with session_scope(_engine()) as session:
        dataset = session.get(StatisticalDataset, dataset_id)
        if dataset is None:
            raise APIError("statistical_dataset_not_found", f"Dataset '{dataset_id}' was not found.", 404)
        dataset_payload = _dataset_payload(session, dataset)
    resources = estimate_dandelion_resources(dataset_payload, parameters)
    state = _statistical_jobs().submit(
        analysis_id=analysis_id,
        dataset=dataset_payload,
        parameters=parameters,
    )
    with session_scope(_engine()) as session:
        session.add(
            Run(
                id=analysis_id,
                sample_id=None,
                status=state["status"],
                stage=state["stage"],
                progress_percent=0,
                genes=[],
                context=dict(dataset_payload.get("context") or {}),
                configuration={"kind": "cohort_statistical_analysis", "method": "dandelion", "parameters": parameters},
            )
        )
        session.add(
            StatisticalAnalysisRun(
                id=analysis_id,
                run_id=analysis_id,
                dataset_id=dataset_id,
                method="dandelion",
                method_version=DANDELION_VERSION,
                worker_job_id=state["id"],
                status=state["status"],
                parameters=parameters,
                resource_estimate=resources,
            )
        )
        append_audit_event(
            session,
            "statistical_analysis_submitted",
            entity_type="statistical_analysis",
            entity_id=analysis_id,
            payload={"dataset_id": dataset_id, "method": "dandelion", "parameters": parameters},
        )
    response = jsonify(
        {
            "id": analysis_id,
            "dataset_id": dataset_id,
            "method": "dandelion",
            "method_version": DANDELION_VERSION,
            "status": state["status"],
            "stage": state["stage"],
            "resource_estimate": resources,
            "result_url": f"/api/v2/statistical-analyses/{analysis_id}/result",
        }
    )
    response.status_code = 202
    response.headers["Location"] = f"/api/v2/statistical-analyses/{analysis_id}"
    return response


@api_v2.get("/statistical-analyses")
def list_statistical_analyses():
    values = []
    with session_scope(_engine()) as session:
        analyses = session.scalars(
            select(StatisticalAnalysisRun).order_by(StatisticalAnalysisRun.created_at.desc())
        ).all()
        for analysis in analyses:
            run = session.get(Run, analysis.run_id)
            state = _statistical_jobs().get(analysis.worker_job_id) or {"status": "unavailable"}
            if run is not None:
                values.append(_analysis_payload(session, analysis, run, state))
    return jsonify({"analyses": values, "count": len(values)})


@api_v2.get("/statistical-analyses/<analysis_id>")
def get_statistical_analysis(analysis_id: str):
    analysis, _dataset = _load_analysis(analysis_id)
    return jsonify(analysis)


@api_v2.post("/statistical-analyses/<analysis_id>/cancel")
def cancel_statistical_analysis(analysis_id: str):
    with session_scope(_engine()) as session:
        analysis = session.get(StatisticalAnalysisRun, analysis_id)
        if analysis is None:
            raise APIError("statistical_analysis_not_found", f"Analysis '{analysis_id}' was not found.", 404)
        state = _statistical_jobs().cancel(analysis.worker_job_id)
        if state is None:
            raise APIError("statistical_job_not_found", "The offline worker job is unavailable.", 500)
        append_audit_event(
            session,
            "statistical_analysis_cancellation_requested",
            entity_type="statistical_analysis",
            entity_id=analysis_id,
        )
    return jsonify({"id": analysis_id, **state}), 202


@api_v2.post("/statistical-analyses/<analysis_id>/retry")
def retry_statistical_analysis(analysis_id: str):
    with session_scope(_engine()) as session:
        source = session.get(StatisticalAnalysisRun, analysis_id)
        if source is None:
            raise APIError("statistical_analysis_not_found", f"Analysis '{analysis_id}' was not found.", 404)
        source_state = _statistical_jobs().get(source.worker_job_id)
        if source_state is None or source_state.get("status") not in {"completed", "failed", "cancelled"}:
            raise APIError("statistical_analysis_not_retryable", "Only terminal cohort analyses can be retried.", 409)
        dataset = session.get(StatisticalDataset, source.dataset_id)
        if dataset is None:
            raise APIError("statistical_dataset_not_found", "The source dataset is unavailable.", 404)
        dataset_payload = _dataset_payload(session, dataset)
        parameters = dict(source.parameters or {})
        resources = dict(source.resource_estimate or estimate_dandelion_resources(dataset_payload, parameters))
    retry_id = uuid.uuid4().hex
    state = _statistical_jobs().submit(analysis_id=retry_id, dataset=dataset_payload, parameters=parameters)
    with session_scope(_engine()) as session:
        session.add(
            Run(
                id=retry_id,
                parent_run_id=analysis_id,
                sample_id=None,
                status=state["status"],
                stage=state["stage"],
                progress_percent=0,
                genes=[],
                context=dict(dataset_payload.get("context") or {}),
                configuration={"kind": "cohort_statistical_analysis", "method": "dandelion", "parameters": parameters},
            )
        )
        session.add(
            StatisticalAnalysisRun(
                id=retry_id,
                run_id=retry_id,
                dataset_id=dataset_payload["id"],
                method="dandelion",
                method_version=DANDELION_VERSION,
                worker_job_id=state["id"],
                status=state["status"],
                parameters=parameters,
                resource_estimate=resources,
            )
        )
        append_audit_event(
            session,
            "statistical_analysis_retried",
            entity_type="statistical_analysis",
            entity_id=retry_id,
            payload={"source_analysis_id": analysis_id, "dataset_id": dataset_payload["id"]},
        )
    response = jsonify(
        {
            "id": retry_id,
            "parent_analysis_id": analysis_id,
            "dataset_id": dataset_payload["id"],
            "method": "dandelion",
            "status": state["status"],
            "stage": state["stage"],
            "resource_estimate": resources,
            "result_url": f"/api/v2/statistical-analyses/{retry_id}/result",
        }
    )
    response.status_code = 202
    response.headers["Location"] = f"/api/v2/statistical-analyses/{retry_id}"
    return response


@api_v2.get("/statistical-analyses/<analysis_id>/result")
def get_statistical_analysis_result(analysis_id: str):
    analysis_payload, dataset = _load_analysis(analysis_id)
    if analysis_payload["status"] != "completed":
        raise APIError("statistical_analysis_not_finished", "The cohort analysis has not completed.", 409)
    with session_scope(_engine()) as session:
        analysis = session.get(StatisticalAnalysisRun, analysis_id)
        if analysis is None:
            raise APIError("statistical_analysis_not_found", f"Analysis '{analysis_id}' was not found.", 404)
        raw_result = _statistical_jobs().result(analysis.worker_job_id)
        if raw_result is None:
            raise APIError("statistical_result_missing", "The worker completed without a readable result.", 500)
        try:
            normalized = validate_runner_result(raw_result)
        except DandelionValidationError as exc:
            raise APIError("invalid_statistical_result", str(exc), 500) from exc
        report = build_dandelion_report(
            analysis_id=analysis_id,
            dataset=dataset,
            parameters=analysis.parameters or {},
            normalized_result=normalized,
            job=analysis_payload,
        )
        existing = session.scalar(
            select(StatisticalResult.id).where(StatisticalResult.analysis_run_id == analysis_id).limit(1)
        )
        if existing is None:
            for row in normalized.get("records", []):
                session.add(
                    StatisticalResult(
                        run_id=analysis.run_id,
                        analysis_run_id=analysis.id,
                        entity_type="gene_pair",
                        entity_key=f"{row.get('exposure', '')}>{row.get('candidate_gene', '')}",
                        family=f"dandelion:{row.get('exposure', '')}",
                        method=f"DANDELION {DANDELION_VERSION}",
                        effect_size=None,
                        raw_p=row.get("p_value") if isinstance(row.get("p_value"), (int, float)) else None,
                        q_value=row.get("q_value") if isinstance(row.get("q_value"), (int, float)) else None,
                        status="significant" if row.get("significant") else "not_significant",
                        limitations=["Hypothesis-generating cohort association; not causal or medical evidence."],
                        details=dict(row),
                    )
                )
            analysis.result_checksum_sha256 = hashlib.sha256(
                json.dumps(normalized, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
            ).hexdigest()
            append_audit_event(
                session,
                "statistical_result_normalized",
                entity_type="statistical_analysis",
                entity_id=analysis_id,
                payload={"record_count": len(normalized.get("records", [])), "checksum_sha256": analysis.result_checksum_sha256},
            )
    return jsonify(report)


@api_v2.get("/statistical-analyses/<analysis_id>/interactions")
def get_statistical_analysis_interactions(analysis_id: str):
    report_response = get_statistical_analysis_result(analysis_id)
    report = report_response.get_json()
    try:
        node_cap = max(2, min(150, int(request.args.get("node_cap", 150))))
    except ValueError as exc:
        raise APIError("invalid_graph_options", "node_cap must be an integer.", 422) from exc
    requested = str(request.args.get("source") or "").strip()
    edges = list(report.get("sections", {}).get("interactions", {}).get("edges", []))
    if requested:
        edges = [edge for edge in edges if edge.get("source_gene") == requested or edge.get("target_gene") == requested]
    ranked = sorted(
        edges,
        key=lambda edge: (
            edge.get("native_score") if isinstance(edge.get("native_score"), (int, float)) else 2.0,
            str(edge.get("source_gene") or ""),
            str(edge.get("target_gene") or ""),
        ),
    )
    nodes: dict[str, dict[str, Any]] = {}
    accepted = []
    for edge in ranked:
        source = str(edge.get("source_gene") or "")
        target = str(edge.get("target_gene") or "")
        needed = [node for node in (source, target) if node and node not in nodes]
        if len(nodes) + len(needed) > node_cap:
            continue
        nodes.setdefault(
            source,
            {"id": source, "node_type": edge.get("details", {}).get("source_node_type", "gene")},
        )
        nodes.setdefault(target, {"id": target, "node_type": "gene"})
        accepted.append(edge)
    return jsonify(
        {
            "analysis_id": analysis_id,
            "semantics": "trans-regulatory disease-prioritisation hypotheses; not a physical interaction network",
            "score_direction": "lower BH q-value is stronger",
            "nodes": list(nodes.values()),
            "edges": accepted,
            "node_count": len(nodes),
            "edge_count": len(accepted),
            "truncated": len(accepted) < len(edges),
            "node_cap": node_cap,
        }
    )


def _v2_job_request(payload: dict[str, Any]) -> dict[str, Any]:
    genes = payload.get("genes")
    if genes is None and payload.get("gene"):
        genes = [payload["gene"]]
    options = dict(payload.get("options") or {})
    if payload.get("sample_context") is not None:
        options["sample_context"] = payload["sample_context"]
    if payload.get("source_set") is not None:
        options["knowledge_sources"] = payload["source_set"]
        options["use_dynamic_knowledge_base"] = bool(payload["source_set"])
    if payload.get("external_consents") is not None:
        options["external_consents"] = payload["external_consents"]
    if payload.get("models") is not None:
        options["requested_models"] = [
            {"model_id": item} if isinstance(item, str) else dict(item)
            for item in payload["models"]
        ]
    return normalize_job_request(
        {
            "operation": payload.get("operation") or "full_workflow",
            "genes": genes,
            "profile_id": payload.get("profile_id") or payload.get("sample_id"),
            "source_job_id": payload.get("source_run_id") or payload.get("source_job_id"),
            "analysis_scope": payload.get("analysis_scope") or "promoter_plus_gene",
            "genome_build": payload.get("genome_build") or "auto",
            "region_overrides": payload.get("region_overrides") or {},
            "options": options,
        }
    )


def _persist_submitted_run(job: dict[str, Any], normalized: dict[str, Any], *, parent_run_id: str = "") -> None:
    with session_scope(_engine()) as session:
        sample_id = normalized.get("profile_id") or None
        if sample_id and session.get(Sample, sample_id) is None:
            session.add(Sample(id=sample_id, pseudonym=sample_id, user_label=sample_id, context=normalized["options"].get("sample_context", {})))
        session.add(
            Run(
                id=job["id"],
                parent_run_id=parent_run_id or None,
                sample_id=sample_id,
                status=job["status"],
                stage=job["stage"],
                progress_percent=0,
                genes=list(normalized["genes"]),
                context=dict(normalized["options"].get("sample_context", {})),
                configuration=normalized,
            )
        )
        append_audit_event(
            session,
            "run_submitted",
            entity_type="run",
            entity_id=job["id"],
            payload={"genes": normalized["genes"], "external_sources": normalized["options"].get("knowledge_sources", [])},
        )


@api_v2.post("/runs")
def submit_run():
    normalized = _v2_job_request(_json_body())
    job = _jobs().submit(normalized)
    _persist_submitted_run(job, normalized)
    response = jsonify(
        {
            **job,
            "schema_version": "3.0",
            "result_url": f"/api/v2/runs/{job['id']}/result",
            "run_url": f"/api/v2/runs/{job['id']}",
        }
    )
    response.status_code = 202
    response.headers["Location"] = f"/api/v2/runs/{job['id']}"
    return response


def _sync_run(job: dict[str, Any]) -> None:
    with session_scope(_engine()) as session:
        run = session.get(Run, job["id"])
        if run is None:
            return
        run.status = job["status"]
        run.stage = job["stage"]
        run.progress_percent = int(job.get("progress", {}).get("percent") or 0)
        run.error = dict(job.get("error") or {})


@api_v2.get("/runs")
def list_runs():
    jobs = _jobs().list()
    return jsonify({"runs": [{**job, "result_url": f"/api/v2/runs/{job['id']}/result"} for job in jobs], "count": len(jobs)})


@api_v2.get("/runs/<run_id>")
def get_run(run_id: str):
    job = _jobs().get(run_id)
    _sync_run(job)
    return jsonify({**job, "schema_version": "3.0", "result_url": f"/api/v2/runs/{run_id}/result"})


def _canonical_results(run_id: str) -> list[dict[str, Any]]:
    manager = _jobs()
    result_manifest = manager.result(run_id)
    reports: list[dict[str, Any]] = []
    for outcome in result_manifest.get("genes", []):
        if outcome.get("status") != "succeeded":
            continue
        report = read_json(manager.jobs_root / run_id / "genes" / str(outcome.get("gene")) / "report.json")
        if isinstance(report, dict) and report.get("schema_version") == "3.0":
            reports.append(report)
    if reports:
        with session_scope(_engine()) as session:
            for report in reports:
                persist_canonical_report(session, report)
    return reports


@api_v2.get("/runs/<run_id>/result")
def get_result(run_id: str):
    job = _jobs().get(run_id)
    if job["status"] not in {"succeeded", "partial"}:
        raise APIError("run_not_finished", f"Run '{run_id}' has not produced a result.", 409)
    reports = _canonical_results(run_id)
    if len(reports) == 1:
        return jsonify(reports[0])
    return jsonify(
        {
            "schema_version": "3.0",
            "run": {"id": run_id, "genes": job["genes"], "status": job["status"]},
            "results": reports,
            "partial_failures": [outcome for outcome in job.get("outcomes", []) if outcome.get("status") != "succeeded"],
        }
    )


@api_v2.post("/runs/<run_id>/retry")
def retry_run(run_id: str):
    manager = _jobs()
    source = manager.get(run_id)
    if source["status"] not in {"succeeded", "partial", "failed"}:
        raise APIError("run_not_retryable", "Only terminal runs can be retried.", 409)
    original = read_json(manager.jobs_root / run_id / "request.json")
    if not isinstance(original, dict):
        raise APIError("run_request_missing", "The original immutable run request is unavailable.", 500)
    retry_options = _json_body() if request.data else {}
    cloned = json.loads(json.dumps(original))
    cloned["source_job_id"] = ""
    cloned["options"]["retry_targets"] = retry_options.get("targets", [])
    job = manager.submit(cloned)
    _persist_submitted_run(job, cloned, parent_run_id=run_id)
    return jsonify({**job, "parent_run_id": run_id, "result_url": f"/api/v2/runs/{job['id']}/result"}), 202


def _refresh_payload_preview(run_id: str, source_keys: list[str]) -> dict[str, Any]:
    reports = _canonical_results(run_id)
    if not reports:
        raise APIError("result_not_found", "A completed schema-v3 result is required before evidence refresh.", 404)
    disclosures: list[dict[str, Any]] = []
    for report in reports:
        objective = report.get("sections", {}).get("objective_data", {})
        variants = []
        for row in objective.get("all_variants", []):
            variants.append(
                {
                    key: row.get(key)
                    for key in ("CHROM", "chrom", "POS", "pos", "REF", "ref", "ALT", "alt", "ID", "id", "rsid")
                    if row.get(key) not in (None, "")
                }
            )
        disclosures.append(
            {
                "gene": report.get("run", {}).get("gene"),
                "genome_build": report.get("run", {}).get("genome_build"),
                "region": report.get("run", {}).get("region"),
                "variants": variants,
                "phenotype_terms": report.get("run", {}).get("sample_context", {}).get("phenotype_terms", []),
            }
        )
    preview = {"run_id": run_id, "sources": sorted(set(source_keys)), "payloads": disclosures}
    preview["payload_sha256"] = hashlib.sha256(
        json.dumps(preview, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    ).hexdigest()
    return preview


@api_v2.get("/runs/<run_id>/evidence/payload-preview")
def evidence_payload_preview(run_id: str):
    sources = [item.strip() for item in request.args.get("sources", "").split(",") if item.strip()]
    invalid = sorted(set(sources) - set(CORE_SOURCE_KEYS))
    if not sources or invalid:
        raise APIError("invalid_sources", "Supply one or more vetted source keys.", 422)
    return jsonify(_refresh_payload_preview(run_id, sources))


@api_v2.post("/runs/<run_id>/evidence/refresh")
def refresh_evidence(run_id: str):
    payload = _json_body()
    sources = [str(item).strip() for item in payload.get("sources", []) if str(item).strip()]
    invalid = sorted(set(sources) - set(CORE_SOURCE_KEYS))
    if not sources or invalid:
        raise APIError("invalid_sources", "Evidence refresh requires vetted source keys only.", 422)
    consents = payload.get("external_consents") or {}
    imports = payload.get("source_imports") or {}
    unapproved = [source for source in sources if not consents.get(source) and source not in imports]
    if unapproved:
        raise APIError("external_consent_required", "Consent is missing for: " + ", ".join(unapproved), 422)
    preview = _refresh_payload_preview(run_id, sources)
    if payload.get("payload_sha256") != preview["payload_sha256"]:
        raise APIError(
            "payload_preview_confirmation_required",
            "Retrieve and confirm the current exact payload preview before external transfer.",
            409,
        )
    manager = _jobs()
    original = read_json(manager.jobs_root / run_id / "request.json")
    if not isinstance(original, dict):
        raise APIError("run_request_missing", "The immutable source request is unavailable.", 500)
    cloned = json.loads(json.dumps(original))
    cloned["source_job_id"] = ""
    cloned["options"].update(
        use_dynamic_knowledge_base=True,
        knowledge_sources=sources,
        external_consents={key: bool(value) for key, value in consents.items()},
        knowledge_source_imports={key: str(value) for key, value in imports.items()},
        evidence_refresh_parent=run_id,
        disclosed_payload_sha256=preview["payload_sha256"],
    )
    job = manager.submit(cloned)
    _persist_submitted_run(job, cloned, parent_run_id=run_id)
    with session_scope(_engine()) as session:
        append_audit_event(
            session,
            "external_evidence_disclosure_approved",
            entity_type="run",
            entity_id=job["id"],
            payload={"sources": sources, "payload_sha256": preview["payload_sha256"], "payloads": preview["payloads"]},
        )
    return jsonify({**job, "parent_run_id": run_id, "payload_sha256": preview["payload_sha256"]}), 202


@api_v2.get("/evidence/sources")
def evidence_sources():
    registered = {spec.key: spec for spec in list_source_specs()}
    sources = []
    for metadata in list_core_source_metadata():
        spec = registered.get(metadata["key"])
        live = bool(spec and spec.connector_kind not in {"metadata", "auth_metadata", "licensed_metadata"})
        ingestion_modes = list(spec.ingestion_modes) if spec else []
        sources.append(
            {
                **metadata,
                "live_connector": live,
                "adapter_status": "available" if live else ("import_only" if "user_export" in ingestion_modes else "not_installed"),
                "ingestion_modes": ingestion_modes,
                "requires_explicit_consent": live,
                "counts_as_assessed": live,
                "license_note": spec.license_note if spec else "Adapter installation and upstream terms must be verified before use.",
            }
        )
    return jsonify({"sources": sources, "count": len(sources), "policy": "core vetted allowlist; linkouts never count as assessed"})


@api_v2.get("/models")
def models():
    manifests = list_model_manifests()
    return jsonify({"models": manifests, "count": len(manifests), "consensus_enabled": False})


@api_v2.post("/models/<model_id>/eligibility")
def model_eligibility(model_id: str):
    return jsonify(inspect_model_inputs(model_id, _json_body()))


@api_v2.post("/models/<model_id>/estimate")
def model_estimate(model_id: str):
    estimate = estimate_model_resources(model_id, _json_body())
    return jsonify(estimate), (404 if estimate.get("status") == "unknown_model" else 200)


@api_v2.get("/models/<model_id>/installation")
def model_installation(model_id: str):
    plan = model_installation_plan(model_id)
    return jsonify(plan), (404 if plan.get("status") == "unknown_model" else 200)


@api_v2.post("/models/<model_id>/installation")
def confirm_model_installation(model_id: str):
    plan = model_installation_plan(model_id)
    if plan.get("status") == "unknown_model":
        raise APIError("model_not_found", f"Model '{model_id}' is not allowlisted.", 404)
    payload = _json_body()
    confirmed = (
        payload.get("manifest_checksum_sha256") == plan["manifest_checksum_sha256"]
        and payload.get("license_acknowledged") is True
        and payload.get("asset_checksums_acknowledged") is True
    )
    if not plan["installable"]:
        raise APIError("model_installation_blocked", "This model is intentionally blocked by its scientific gate.", 409)
    if not confirmed:
        raise APIError("model_installation_confirmation_required", "Confirm the manifest, license, and asset checksums.", 422)
    return jsonify(
        {
            **plan,
            "status": "runner_required",
            "message": "The web service cannot install models. Submit this signed manifest to the privilege-separated local runner.",
        }
    ), 202


def _persist_model_job(job: dict[str, Any], *, inputs: dict[str, Any] | None = None, predictions: list[Any] | None = None) -> None:
    with session_scope(_engine()) as session:
        model_run = session.get(ModelRun, job["id"])
        if model_run is None:
            input_checksum = hashlib.sha256(
                json.dumps(inputs or {}, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
            ).hexdigest()
            model_run = ModelRun(
                id=job["id"],
                run_id=job["run_id"],
                model_id=job["model_id"],
                status=job["status"],
                blockers=list(job.get("eligibility", {}).get("blockers") or []),
                input_checksum_sha256=input_checksum,
                error=dict(job.get("error") or {}),
            )
            session.add(model_run)
        else:
            model_run.status = job["status"]
            model_run.blockers = list(job.get("eligibility", {}).get("blockers") or model_run.blockers or [])
            model_run.error = dict(job.get("error") or {})
        if predictions:
            for index, item in enumerate(predictions):
                row = item if isinstance(item, dict) else {"raw_output": item}
                entity_type = str(row.get("entity_type") or "model_output")
                entity_key = str(row.get("entity_key") or index)
                output_name = str(row.get("output_name") or row.get("name") or "score")
                exists = session.scalar(
                    select(Prediction.id).where(
                        Prediction.model_run_id == job["id"],
                        Prediction.entity_type == entity_type,
                        Prediction.entity_key == entity_key,
                        Prediction.output_name == output_name,
                    )
                )
                if exists:
                    continue
                session.add(
                    Prediction(
                        model_run_id=job["id"],
                        entity_type=entity_type,
                        entity_key=entity_key,
                        output_name=output_name,
                        raw_score=row.get("raw_score") if isinstance(row.get("raw_score"), (int, float)) else None,
                        validated_label=str(row.get("validated_label") or ""),
                        calibration=str(row.get("calibration") or ""),
                        applicability=dict(row.get("applicability") or {}),
                        limitations=list(row.get("limitations") or []),
                    )
                )


@api_v2.post("/runs/<run_id>/models/<model_id>")
def submit_model_job(run_id: str, model_id: str):
    _jobs().get(run_id)
    inputs = _json_body()
    eligibility = inspect_model_inputs(model_id, inputs)
    if eligibility.get("status") == "unknown_model":
        raise APIError("model_not_found", f"Model '{model_id}' is not allowlisted.", 404)
    job = _model_jobs().submit(run_id=run_id, model_id=model_id, inputs=inputs)
    _persist_model_job(job, inputs=inputs)
    status_code = 202 if job["status"] == "queued" else 422
    return jsonify(job), status_code


@api_v2.get("/model-jobs/<job_id>")
def get_model_job(job_id: str):
    try:
        job = _model_jobs().get(job_id)
        _persist_model_job(job)
        return jsonify(job)
    except (ValueError, FileNotFoundError) as exc:
        raise APIError("model_job_not_found", f"Model job '{job_id}' was not found.", 404) from exc


@api_v2.post("/model-jobs/<job_id>/cancel")
def cancel_model_job(job_id: str):
    try:
        job = _model_jobs().cancel(job_id)
        _persist_model_job(job)
        return jsonify(job)
    except FileNotFoundError as exc:
        raise APIError("model_job_not_found", f"Model job '{job_id}' was not found.", 404) from exc
    except ValueError as exc:
        raise APIError("model_job_not_cancellable", str(exc), 409) from exc


@api_v2.get("/model-jobs/<job_id>/result")
def get_model_job_result(job_id: str):
    try:
        job = _model_jobs().get(job_id)
    except (ValueError, FileNotFoundError) as exc:
        raise APIError("model_job_not_found", f"Model job '{job_id}' was not found.", 404) from exc
    if job["status"] != "succeeded":
        raise APIError("model_job_not_finished", "The model job has no validated result.", 409)
    payload = read_json(_model_jobs().root / job_id / "predictions.json", default=[])
    normalized = payload if isinstance(payload, list) else [payload]
    _persist_model_job(job, predictions=normalized)
    return jsonify({"model_job": job, "predictions": normalized, "consensus": None})


@api_v2.get("/runs/<run_id>/interactions")
def interactions(run_id: str):
    reports = _canonical_results(run_id)
    if not reports:
        raise APIError("result_not_found", "No schema-v3 report is available for this run.", 404)
    gene = str(request.args.get("gene") or reports[0].get("run", {}).get("gene") or "").upper()
    try:
        max_hops = max(1, min(3, int(request.args.get("hops", 1))))
        node_cap = max(2, min(150, int(request.args.get("node_cap", 150))))
    except ValueError as exc:
        raise APIError("invalid_graph_options", "hops and node_cap must be integers.", 422) from exc
    edges = reports[0].get("run_details", {}).get("interaction_evidence_edges", [])
    return jsonify(build_interaction_graph(gene, edges, max_hops=max_hops, node_cap=node_cap))


@api_v2.post("/runs/<run_id>/exports")
def export_run(run_id: str):
    reports = _canonical_results(run_id)
    if not reports:
        raise APIError("result_not_found", "No schema-v3 report is available for export.", 404)
    payload = _json_body()
    acknowledged = bool(payload.get("sensitive_data_acknowledged"))
    if not acknowledged:
        raise APIError(
            "sensitive_export_confirmation_required",
            "Full-fidelity exports contain sensitive genetic and sample-context data.",
            422,
        )
    password = str(payload.get("password") or "")
    export_dir = _jobs().jobs_root / run_id / "exports"
    export_path = export_dir / f"nophigene-{run_id}-{uuid.uuid4().hex[:8]}.zip"
    manifest = create_export_bundle(reports[0], export_path, password=password)
    with session_scope(_engine()) as session:
        append_audit_event(
            session,
            "full_fidelity_export_created",
            entity_type="run",
            entity_id=run_id,
            payload={"encrypted": bool(password), "path": str(export_path)},
        )
    return jsonify({**manifest, "download_url": f"/api/v2/runs/{run_id}/exports/{export_path.name}"}), 201


@api_v2.get("/runs/<run_id>/exports/<filename>")
def download_export(run_id: str, filename: str):
    if Path(filename).name != filename or not filename.endswith(".zip"):
        raise APIError("invalid_export_path", "Invalid export filename.", 400)
    path = (_jobs().jobs_root / run_id / "exports" / filename).resolve()
    export_root = (_jobs().jobs_root / run_id / "exports").resolve()
    path.relative_to(export_root)
    if not path.is_file():
        raise APIError("export_not_found", "Export bundle was not found.", 404)
    return send_file(path, as_attachment=True)


@api_v2.delete("/samples/<sample_id>")
def delete_sample(sample_id: str):
    payload = _json_body()
    if payload.get("confirmation") != sample_id:
        raise APIError("deletion_confirmation_required", "confirmation must exactly match the sample ID.", 422)
    run_ids: list[str] = []
    model_job_ids: list[str] = []
    with session_scope(_engine()) as session:
        sample = session.get(Sample, sample_id)
        if sample is None:
            raise APIError("sample_not_found", f"Sample '{sample_id}' was not found.", 404)
        run_ids = [row.id for row in session.scalars(select(Run).where(Run.sample_id == sample_id)).all()]
        run_count = len(run_ids)
        if run_ids:
            model_job_ids = [row.id for row in session.scalars(select(ModelRun).where(ModelRun.run_id.in_(run_ids))).all()]
        tombstone_hash = hashlib.sha256(f"{sample_id}:{uuid.uuid4().hex}".encode("utf-8")).hexdigest()
        session.add(
            DeletionRecord(
                entity_type="sample",
                tombstone_hash=tombstone_hash,
                deleted_counts={"samples": 1, "runs": run_count, "model_jobs": len(model_job_ids)},
                retained_counts={"public_evidence_snapshots": 1, "audit_tombstones": 1},
            )
        )
        append_audit_event(
            session,
            "sample_deleted",
            entity_type="sample_tombstone",
            entity_id="",
            payload={"deleted_run_count": run_count, "public_evidence_retained": True},
        )
        session.delete(sample)
    removed_paths: list[str] = []
    jobs_root = _jobs().jobs_root.resolve()
    for run_id in run_ids:
        target = (jobs_root / run_id).resolve()
        try:
            target.relative_to(jobs_root)
        except ValueError:
            continue
        if target.is_dir():
            shutil.rmtree(target)
            removed_paths.append(str(target))
    model_root = _model_jobs().root.resolve()
    for job_id in model_job_ids:
        target = (model_root / job_id).resolve()
        try:
            target.relative_to(model_root)
        except ValueError:
            continue
        if target.is_dir():
            shutil.rmtree(target)
            removed_paths.append(str(target))
    return jsonify(
        {
            "deleted": {
                "sample": 1,
                "runs": run_count,
                "model_jobs": len(model_job_ids),
                "sample_derivatives": "database cascade",
                "artifact_directories": removed_paths,
            },
            "retained": {"public_evidence_snapshots": True, "non_identifying_audit_tombstone": True},
        }
    )


@api_v2.get("/health")
def health():
    manager = _jobs()
    dandelion = _statistical_jobs().health()
    return jsonify(
        {
            "status": "ok" if manager.worker_alive else "degraded",
            "worker": {"alive": manager.worker_alive, "queue_depth": manager.queue_depth},
            "dandelion": dandelion,
            "database": {"schema": "3.0", "encryption_required_in_production": True},
            "gpu": {"status": "preflight_required", "message": "GPU models remain unavailable until WSL2/NVIDIA validation passes."},
        }
    )
