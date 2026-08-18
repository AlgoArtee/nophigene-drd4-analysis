"""Flask routes for the versioned local API."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from flask import Blueprint, current_app, jsonify, request, send_from_directory
from werkzeug.exceptions import HTTPException, NotFound

from .errors import APIError
from .openapi import build_openapi_document
from .workflow_runner import normalize_job_request

try:
    from ..bam_extraction import get_extraction_tool_status, get_hg38_reference_status
    from ..variant_knowledge.credentials import credential_status_for_specs
    from ..variant_knowledge.registry import get_source_spec, list_source_cards, list_source_specs
    from ..variant_knowledge.workflows import default_workflow_source_keys, list_workflow_cards
except ImportError:
    from bam_extraction import get_extraction_tool_status, get_hg38_reference_status
    from variant_knowledge.credentials import credential_status_for_specs
    from variant_knowledge.registry import get_source_spec, list_source_cards, list_source_specs
    from variant_knowledge.workflows import default_workflow_source_keys, list_workflow_cards

api_v1 = Blueprint("api_v1", __name__, url_prefix="/api/v1")


@api_v1.before_request
def freeze_legacy_mutations():
    """Schema v1 is retained only for reading legacy profiles/jobs/artifacts."""
    if request.method not in {"GET", "HEAD", "OPTIONS"}:
        return APIError(
            "api_v1_read_only",
            "API v1 is frozen for legacy reads. Submit new work through /api/v2.",
            410,
            details={"replacement": "/api/v2"},
        ).to_response()


def _profiles():
    return current_app.config["NOPHIGENE_PROFILE_STORE"]


def _jobs():
    return current_app.config["NOPHIGENE_JOB_MANAGER"]


def _json_body() -> Any:
    payload = request.get_json(silent=True)
    if payload is None:
        raise APIError("invalid_json", "Request body must contain valid JSON.", 400)
    return payload


@api_v1.errorhandler(APIError)
def handle_api_error(error: APIError):
    return error.to_response()


@api_v1.errorhandler(Exception)
def handle_unexpected_error(error: Exception):
    if isinstance(error, HTTPException):
        return APIError(
            "http_error",
            error.description,
            error.code or 500,
        ).to_response()
    current_app.logger.exception("Unhandled API error", exc_info=error)
    return APIError(
        "internal_error",
        "The API could not complete the request.",
        500,
    ).to_response()


@api_v1.get("")
@api_v1.get("/")
def api_index():
    return jsonify(
        {
            "name": "NophiGene Local Gene Workflow API",
            "version": "1.0",
            "health": "/api/v1/health",
            "openapi": "/api/v1/openapi.json",
            "profiles": "/api/v1/profiles",
            "jobs": "/api/v1/jobs",
            "knowledge_sources": "/api/v1/knowledge-sources",
            "knowledge_workflows": "/api/v1/knowledge-workflows",
        }
    )


@api_v1.get("/profiles")
def list_profiles():
    profiles = _profiles().list()
    return jsonify({"profiles": profiles, "count": len(profiles)})


@api_v1.post("/profiles")
def create_profile():
    return jsonify(_profiles().create(_json_body())), 201


@api_v1.get("/profiles/<profile_id>")
def get_profile(profile_id: str):
    return jsonify(_profiles().get(profile_id))


@api_v1.put("/profiles/<profile_id>")
def update_profile(profile_id: str):
    return jsonify(_profiles().update(profile_id, _json_body()))


@api_v1.delete("/profiles/<profile_id>")
def delete_profile(profile_id: str):
    _profiles().delete(profile_id)
    return "", 204


@api_v1.get("/jobs")
def list_jobs():
    jobs = _jobs().list()
    return jsonify({"jobs": jobs, "count": len(jobs)})


@api_v1.get("/knowledge-sources")
def list_knowledge_sources():
    specs = list_source_specs()
    credential_statuses = credential_status_for_specs(specs)
    cards = list_source_cards(
        selected_keys=default_workflow_source_keys(),
        credential_statuses=credential_statuses,
    )
    return jsonify({"sources": cards, "count": len(cards)})


@api_v1.get("/knowledge-workflows")
def list_knowledge_workflows():
    cards = list_workflow_cards()
    return jsonify({"workflows": cards, "count": len(cards)})


@api_v1.post("/knowledge-sources/test")
def test_knowledge_source():
    payload = _json_body()
    source_imports = payload.get("source_imports") or {}
    if source_imports and not isinstance(source_imports, dict):
        raise APIError("invalid_source_imports", "'source_imports' must be an object mapping source keys to paths.", 422)
    raw_sources = payload.get("sources")
    if raw_sources is not None:
        if isinstance(raw_sources, str):
            source_keys = [part.strip() for part in raw_sources.split(",") if part.strip()]
        elif isinstance(raw_sources, list):
            source_keys = [str(item).strip() for item in raw_sources if str(item).strip()]
        else:
            raise APIError("invalid_knowledge_sources", "'sources' must be a list or comma-separated string.", 422)
        if not source_keys:
            raise APIError("invalid_knowledge_sources", "'sources' must include at least one source key.", 422)
        specs = []
        for source_key in source_keys:
            spec = get_source_spec(source_key)
            if spec is None:
                raise APIError("unknown_knowledge_source", f"Unknown knowledge source '{source_key}'.", 404)
            specs.append(spec)
        credential_statuses = credential_status_for_specs(specs)
        return jsonify(
            {
                "results": [
                    _knowledge_source_test_result(
                        spec,
                        credential_statuses[spec.key],
                        import_path=str(source_imports.get(spec.key, "")),
                    )
                    for spec in specs
                ],
                "count": len(specs),
            }
        )

    source_key = str(payload.get("source_key") or payload.get("key") or "").strip()
    spec = get_source_spec(source_key)
    if spec is None:
        raise APIError("unknown_knowledge_source", f"Unknown knowledge source '{source_key}'.", 404)
    status = credential_status_for_specs([spec])[spec.key]
    import_path = str(payload.get("source_import") or source_imports.get(spec.key, "") or "")
    return jsonify(_knowledge_source_test_result(spec, status, import_path=import_path))


def _source_import_ready(spec, import_path: str) -> bool:
    return "user_export" in spec.ingestion_modes and bool(import_path) and Path(import_path).exists()


def _source_readiness(spec, credential_status: str, import_path: str) -> dict[str, str]:
    official_supported = "official_api" in spec.ingestion_modes
    official_ready = official_supported and (not spec.env_var or credential_status != "missing")
    user_export_supported = "user_export" in spec.ingestion_modes
    return {
        "official_api": "ready" if official_ready else "needs_credentials" if official_supported else "not_supported",
        "user_export": (
            "ready"
            if _source_import_ready(spec, import_path)
            else "needs_export"
            if user_export_supported
            else "not_supported"
        ),
        "linkout_only": "available" if "linkout_only" in spec.ingestion_modes else "not_supported",
    }


def _knowledge_source_test_result(spec, status: str, *, import_path: str = "") -> dict[str, Any]:
    readiness = _source_readiness(spec, status, import_path)
    if readiness["official_api"] == "ready" and spec.connector_kind not in {"metadata", "auth_metadata", "licensed_metadata"}:
        query_status = "queryable"
        message = f"{spec.name} has a configured official connector."
    elif readiness["user_export"] == "ready":
        query_status = "import_ready"
        message = f"{spec.name} has a user-provided export ready for dynamic KB ingestion."
    elif spec.access_type in {"auth_api", "licensed"} and status == "missing" and "official_api" in spec.ingestion_modes:
        query_status = "needs_credentials"
        message = f"Set {spec.env_var} before live querying {spec.name}."
    elif spec.requires_export or "user_export" in spec.ingestion_modes:
        query_status = "needs_export"
        message = f"Upload a permitted CSV/JSON export for {spec.name}; scraping is not performed."
    elif spec.connector_kind in {"metadata", "auth_metadata", "licensed_metadata"}:
        query_status = "metadata_only"
        message = spec.license_note
    else:
        query_status = "queryable"
        message = f"{spec.name} has a configured connector."
    return {
        "source": spec.to_card(
            credential_status=status,
            import_status=readiness["user_export"],
            import_path=import_path,
        ),
        "key": spec.key,
        "status": query_status,
        "message": message,
        "readiness": readiness,
    }


@api_v1.post("/jobs")
def submit_job():
    normalized = normalize_job_request(_json_body())
    job = _jobs().submit(normalized)
    response = jsonify(job)
    response.status_code = 202
    response.headers["Location"] = f"/api/v1/jobs/{job['id']}"
    return response


@api_v1.get("/jobs/<job_id>")
def get_job(job_id: str):
    return jsonify(_jobs().get(job_id))


@api_v1.post("/jobs/<job_id>/cancel")
def cancel_job(job_id: str):
    return jsonify(_jobs().cancel(job_id))


@api_v1.get("/jobs/<job_id>/result")
def get_job_result(job_id: str):
    return jsonify(_jobs().result(job_id))


@api_v1.get("/jobs/<job_id>/artifacts/<path:artifact_path>")
def get_job_artifact(job_id: str, artifact_path: str):
    manager = _jobs()
    manager.get(job_id)
    job_dir = (manager.jobs_root / job_id).resolve()
    requested = (job_dir / artifact_path).resolve()
    try:
        requested.relative_to(job_dir)
    except ValueError as exc:
        raise APIError("invalid_artifact_path", "Artifact path escapes the job directory.", 400) from exc
    if not requested.is_file():
        raise APIError("artifact_not_found", "Artifact was not found.", 404)
    try:
        return send_from_directory(
            str(job_dir),
            requested.relative_to(job_dir).as_posix(),
            as_attachment=requested.suffix.lower() == ".zip",
        )
    except NotFound as exc:
        raise APIError("artifact_not_found", "Artifact was not found.", 404) from exc


@api_v1.get("/health")
def health():
    manager = _jobs()
    extraction = get_extraction_tool_status()
    reference = get_hg38_reference_status()
    return jsonify(
        {
            "status": "ok" if manager.worker_alive else "degraded",
            "worker": {
                "alive": manager.worker_alive,
                "queue_depth": manager.queue_depth,
            },
            "extraction": {
                "available": extraction["available"],
                "docker_runtime": extraction["docker_runtime"],
                "missing_tools": extraction["missing_tools"],
                "message": extraction["message"],
            },
            "hg38_reference": {
                "ready": reference["ready"],
                "message": reference["message"],
            },
        }
    )


@api_v1.get("/openapi.json")
def openapi():
    return jsonify(build_openapi_document())
