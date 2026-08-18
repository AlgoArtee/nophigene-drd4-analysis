"""Signed filesystem queue for the privilege-separated DANDELION worker."""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import threading
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

from .dandelion import RUNNER_CONTRACT_VERSION, canonical_json


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


class StatisticalJobError(RuntimeError):
    pass


class StatisticalJobManager:
    """Submit immutable signed jobs and read worker-owned status artifacts."""

    def __init__(self, root: Path, *, key_file: Path | None = None):
        self.root = Path(root)
        self.queue = self.root / "queue"
        self.jobs = self.root / "jobs"
        self._key_file = key_file
        self._lock = threading.Lock()
        self.queue.mkdir(parents=True, exist_ok=True)
        self.jobs.mkdir(parents=True, exist_ok=True)

    def _key(self) -> bytes:
        if self._key_file:
            key = self._key_file.read_bytes().strip()
        else:
            configured = os.environ.get("NOPHIGENE_DANDELION_RUNNER_KEY_FILE")
            if configured:
                key = Path(configured).read_bytes().strip()
            elif os.environ.get("NOPHIGENE_REQUIRE_ENCRYPTION") == "1":
                raise StatisticalJobError("The DANDELION runner signing key is not configured.")
            else:
                # Tests and explicit unencrypted development use a stable local
                # key file; production always supplies a Compose secret.
                local = self.root / ".development-signing-key"
                if not local.exists():
                    local.write_bytes(secrets.token_bytes(48))
                key = local.read_bytes()
        if len(key) < 32:
            raise StatisticalJobError("The DANDELION runner signing key must contain at least 32 bytes.")
        return key

    def _signed_envelope(self, payload: dict[str, Any]) -> dict[str, Any]:
        signature = hmac.new(self._key(), canonical_json(payload), hashlib.sha256).hexdigest()
        return {"payload": payload, "signature": signature, "algorithm": "HMAC-SHA256"}

    def verify_envelope(self, envelope: dict[str, Any]) -> dict[str, Any]:
        payload = envelope.get("payload")
        signature = str(envelope.get("signature") or "")
        if not isinstance(payload, dict):
            raise StatisticalJobError("Signed job payload is invalid.")
        expected = hmac.new(self._key(), canonical_json(payload), hashlib.sha256).hexdigest()
        if not hmac.compare_digest(signature, expected):
            raise StatisticalJobError("Signed job payload failed authentication.")
        return payload

    def submit(self, *, analysis_id: str, dataset: dict[str, Any], parameters: dict[str, Any]) -> dict[str, Any]:
        job_id = secrets.token_hex(16)
        job_dir = self.jobs / job_id
        with self._lock:
            if job_dir.exists():
                raise StatisticalJobError("Generated duplicate job identifier.")
            job_dir.mkdir(parents=False)
            payload = {
                "contract_version": RUNNER_CONTRACT_VERSION,
                "job_id": job_id,
                "analysis_id": analysis_id,
                "submitted_at": _now(),
                "dataset": dataset,
                "parameters": parameters,
            }
            state = {
                "id": job_id,
                "analysis_id": analysis_id,
                "status": "queued",
                "stage": "awaiting_offline_runner",
                "progress_percent": 0,
                "submitted_at": payload["submitted_at"],
                "updated_at": payload["submitted_at"],
                "cancel_requested": False,
            }
            (job_dir / "status.json").write_bytes(canonical_json(state))
            temporary = self.queue / f".{job_id}.tmp"
            temporary.write_bytes(canonical_json(self._signed_envelope(payload)))
            temporary.replace(self.queue / f"{job_id}.json")
        return state

    def _state_path(self, job_id: str) -> Path:
        if len(job_id) != 32 or any(character not in "0123456789abcdef" for character in job_id):
            raise StatisticalJobError("Invalid statistical job identifier.")
        return self.jobs / job_id / "status.json"

    def get(self, job_id: str) -> dict[str, Any] | None:
        path = self._state_path(job_id)
        if not path.exists():
            return None
        return json.loads(path.read_text(encoding="utf-8"))

    def cancel(self, job_id: str) -> dict[str, Any] | None:
        state = self.get(job_id)
        if state is None:
            return None
        if state.get("status") in {"completed", "failed", "cancelled"}:
            return state
        marker = self.jobs / job_id / "cancel.requested"
        marker.write_text(_now(), encoding="utf-8")
        state.update({"cancel_requested": True, "updated_at": _now()})
        self._state_path(job_id).write_bytes(canonical_json(state))
        return state

    def result(self, job_id: str) -> dict[str, Any] | None:
        state = self.get(job_id)
        if state is None or state.get("status") != "completed":
            return None
        path = self.jobs / job_id / "result.json"
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None

    def list(self) -> list[dict[str, Any]]:
        states = []
        for path in self.jobs.glob("*/status.json"):
            try:
                states.append(json.loads(path.read_text(encoding="utf-8")))
            except (OSError, ValueError, json.JSONDecodeError):
                continue
        return sorted(states, key=lambda item: str(item.get("submitted_at") or ""), reverse=True)

    def health(self) -> dict[str, Any]:
        heartbeat = self.root / "worker-heartbeat.json"
        worker = {"status": "unavailable", "reason": "DANDELION worker has not written a heartbeat."}
        if heartbeat.exists():
            try:
                worker = json.loads(heartbeat.read_text(encoding="utf-8"))
                updated = str(worker.get("updated_at") or "").replace("Z", "+00:00")
                if updated:
                    age = (datetime.now(timezone.utc) - datetime.fromisoformat(updated)).total_seconds()
                    if age > 15:
                        worker = {**worker, "status": "unavailable", "reason": "Worker heartbeat is stale."}
            except (OSError, json.JSONDecodeError):
                worker = {"status": "degraded", "reason": "Worker heartbeat is unreadable."}
        return {
            "method": "dandelion",
            "package_version": "0.1.0",
            "queue_depth": len(list(self.queue.glob("*.json"))),
            "worker": worker,
            "network_policy": "offline",
        }


@lru_cache(maxsize=1)
def get_default_statistical_job_manager() -> StatisticalJobManager:
    root = Path(os.environ.get("NOPHIGENE_DANDELION_WORK_ROOT", "results/dandelion"))
    key = os.environ.get("NOPHIGENE_DANDELION_RUNNER_KEY_FILE")
    return StatisticalJobManager(root, key_file=Path(key) if key else None)
