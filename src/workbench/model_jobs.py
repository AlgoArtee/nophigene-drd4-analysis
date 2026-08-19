"""Signed filesystem queue for privilege-separated biological model workers."""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from .model_registry import inspect_model_inputs


MODEL_RUNNER_CONTRACT_VERSION = "1.0"
TERMINAL = {"succeeded", "partial", "failed", "blocked", "cancelled"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, default=str).encode("utf-8")


class ModelJobError(RuntimeError):
    pass


class ModelJobManager:
    def __init__(self, root: Path, *, key_file: Path | None = None):
        self.root = Path(root)
        self.queue = self.root / "queue"
        self.jobs = self.root / "jobs"
        self.key_file = Path(key_file) if key_file else None
        self._lock = threading.RLock()
        self._hooks: list[Callable[[str, dict[str, Any]], None]] = []
        self._thread: threading.Thread | None = None
        self._delivered: set[tuple[str, str]] = set()
        self.queue.mkdir(parents=True, exist_ok=True)
        self.jobs.mkdir(parents=True, exist_ok=True)

    def _key(self) -> bytes:
        configured = self.key_file or (
            Path(os.environ["NOPHIGENE_MODEL_RUNNER_KEY_FILE"])
            if os.environ.get("NOPHIGENE_MODEL_RUNNER_KEY_FILE")
            else None
        )
        if configured:
            key = configured.read_bytes().strip()
        elif os.environ.get("NOPHIGENE_REQUIRE_ENCRYPTION") == "1":
            raise ModelJobError("The model runner signing key is not configured.")
        else:
            development_key = self.root / ".development-signing-key"
            if not development_key.exists():
                development_key.write_bytes(secrets.token_bytes(48))
            key = development_key.read_bytes().strip()
        if len(key) < 32:
            raise ModelJobError("The model runner signing key must contain at least 32 bytes.")
        return key

    def _signed_envelope(self, payload: dict[str, Any]) -> dict[str, Any]:
        signature = hmac.new(self._key(), canonical_json(payload), hashlib.sha256).hexdigest()
        return {"payload": payload, "signature": signature, "algorithm": "HMAC-SHA256"}

    def verify_envelope(self, envelope: dict[str, Any]) -> dict[str, Any]:
        payload = envelope.get("payload")
        signature = str(envelope.get("signature") or "")
        if not isinstance(payload, dict) or envelope.get("algorithm") != "HMAC-SHA256":
            raise ModelJobError("Signed model envelope is invalid.")
        expected = hmac.new(self._key(), canonical_json(payload), hashlib.sha256).hexdigest()
        if not hmac.compare_digest(signature, expected):
            raise ModelJobError("Signed model envelope failed authentication.")
        return payload

    def add_completion_hook(self, hook: Callable[[str, dict[str, Any]], None]) -> None:
        self._hooks.append(hook)

    def start(self) -> None:
        with self._lock:
            if self._thread and self._thread.is_alive():
                return
            self._thread = threading.Thread(target=self._monitor, name="nophigene-model-reconciler", daemon=True)
            self._thread.start()

    def submit(self, *, run_id: str, model_id: str, inputs: dict[str, Any]) -> dict[str, Any]:
        eligibility = (
            {
                "model_id": model_id,
                "name": "AlphaGenome",
                "status": "eligible",
                "eligible": True,
                "blockers": [],
                "input_contract": ["configured_provider_credential"],
                "execution_mode": "external_api_metadata_only",
                "network_policy": "provider_only",
                "limitations": ["No person-specific data are included in metadata verification."],
            }
            if model_id == "alphagenome-api" and inputs.get("operation") == "metadata"
            else inspect_model_inputs(model_id, inputs)
        )
        input_checksum = hashlib.sha256(canonical_json(inputs)).hexdigest()
        with self._lock:
            for existing in self.list():
                if (
                    existing.get("run_id") == run_id
                    and existing.get("model_id") == model_id
                    and existing.get("input_checksum_sha256") == input_checksum
                    and existing.get("status") not in {"failed", "cancelled", "blocked"}
                    and not existing.get("cancel_requested")
                ):
                    return existing
            job_id = secrets.token_hex(16)
            status = "queued" if eligibility.get("eligible") else "blocked"
            created = utc_now()
            job = {
                "id": job_id,
                "run_id": run_id,
                "model_id": model_id,
                "status": status,
                "stage": "awaiting_isolated_worker" if status == "queued" else "eligibility_blocked",
                "eligibility": eligibility,
                "input_checksum_sha256": input_checksum,
                "progress_percent": 0,
                "created_at": created,
                "started_at": None,
                "finished_at": created if status == "blocked" else None,
                "cancel_requested": False,
                "error": None,
            }
            directory = self.jobs / job_id
            directory.mkdir(parents=False)
            (directory / "status.json").write_bytes(canonical_json(job))
            if status == "queued":
                worker_payload = {
                    "contract_version": MODEL_RUNNER_CONTRACT_VERSION,
                    "job_id": job_id,
                    "run_id": run_id,
                    "model_id": model_id,
                    "submitted_at": created,
                    "input_checksum_sha256": input_checksum,
                    "inputs": inputs,
                }
                temporary = self.queue / f".{job_id}.tmp"
                temporary.write_bytes(canonical_json(self._signed_envelope(worker_payload)))
                temporary.replace(self.queue / f"{job_id}.json")
        self.start()
        return job

    def _validate_job_id(self, job_id: str) -> None:
        if len(job_id) != 32 or any(character not in "0123456789abcdef" for character in job_id):
            raise ValueError("Invalid model job ID.")

    def get(self, job_id: str) -> dict[str, Any]:
        self._validate_job_id(job_id)
        path = self.jobs / job_id / "status.json"
        if not path.is_file():
            raise FileNotFoundError(job_id)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ModelJobError("Model job state is invalid.")
        return payload

    def list(self) -> list[dict[str, Any]]:
        jobs = []
        for path in self.jobs.glob("*/status.json"):
            try:
                item = json.loads(path.read_text(encoding="utf-8"))
                if isinstance(item, dict):
                    jobs.append(item)
            except (OSError, ValueError, json.JSONDecodeError):
                continue
        return sorted(jobs, key=lambda item: str(item.get("created_at") or ""), reverse=True)

    def cancel(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self.get(job_id)
            if job["status"] not in {"queued", "running"}:
                raise ValueError("Only queued or running model jobs can be cancelled.")
            (self.jobs / job_id / "cancel.requested").write_text(utc_now(), encoding="utf-8")
            job.update(cancel_requested=True, stage="cancellation_requested")
            (self.jobs / job_id / "status.json").write_bytes(canonical_json(job))
            return job

    def result(self, job_id: str) -> dict[str, Any] | None:
        job = self.get(job_id)
        if job.get("status") not in {"succeeded", "partial", "failed"}:
            return None
        path = self.jobs / job_id / "result.json"
        if not path.is_file():
            return None
        envelope = json.loads(path.read_text(encoding="utf-8"))
        payload = self.verify_envelope(envelope)
        claimed_checksum = str(payload.get("result_checksum_sha256") or "")
        checksum_payload = {key: value for key, value in payload.items() if key != "result_checksum_sha256"}
        actual_checksum = hashlib.sha256(canonical_json(checksum_payload)).hexdigest()
        if not claimed_checksum or not hmac.compare_digest(claimed_checksum, actual_checksum):
            raise ModelJobError("Signed model result checksum validation failed.")
        if payload.get("kind") not in {"metadata", "predictions"}:
            raise ModelJobError("Signed model result kind is invalid.")
        if payload.get("model_id") and payload.get("model_id") != job.get("model_id"):
            raise ModelJobError("Signed model result does not match the submitted model.")
        if not isinstance(payload.get("predictions"), list) or not isinstance(payload.get("failures"), list):
            raise ModelJobError("Signed model result collections are invalid.")
        return payload

    def health(self) -> dict[str, Any]:
        heartbeat_path = self.root / "worker-heartbeat.json"
        worker = {"status": "unavailable", "reason": "The AlphaGenome worker has not written a heartbeat."}
        if heartbeat_path.is_file():
            try:
                worker = json.loads(heartbeat_path.read_text(encoding="utf-8"))
                updated = str(worker.get("updated_at") or "").replace("Z", "+00:00")
                if updated and (datetime.now(timezone.utc) - datetime.fromisoformat(updated)).total_seconds() > 20:
                    worker = {**worker, "status": "unavailable", "reason": "The worker heartbeat is stale."}
            except (OSError, ValueError, json.JSONDecodeError):
                worker = {"status": "degraded", "reason": "The worker heartbeat is unreadable."}
        return {
            "worker": worker,
            "queue_depth": len(list(self.queue.glob("*.json"))),
            "network_policy": "provider_allowlist",
            "contract_version": MODEL_RUNNER_CONTRACT_VERSION,
        }

    def _monitor(self) -> None:
        while True:
            for job in self.list():
                status = str(job.get("status") or "")
                delivery_key = (str(job.get("id") or ""), status)
                if status not in TERMINAL or delivery_key in self._delivered:
                    continue
                delivered = True
                for hook in self._hooks:
                    try:
                        hook(str(job["id"]), job)
                    except Exception:
                        delivered = False
                if delivered:
                    self._delivered.add(delivery_key)
            time.sleep(1.0)


_DEFAULT_MANAGER: ModelJobManager | None = None


def get_default_model_job_manager() -> ModelJobManager:
    global _DEFAULT_MANAGER
    if _DEFAULT_MANAGER is None:
        root = Path(os.environ.get("NOPHIGENE_MODEL_WORK_ROOT", "results/model-jobs"))
        key = os.environ.get("NOPHIGENE_MODEL_RUNNER_KEY_FILE")
        _DEFAULT_MANAGER = ModelJobManager(root, key_file=Path(key) if key else None)
    return _DEFAULT_MANAGER
