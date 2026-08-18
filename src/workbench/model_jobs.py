"""Persistent one-at-a-time biological model job queue.

The queue records eligibility and lifecycle independently of analysis jobs.
Execution adapters are explicitly registered; an allowlisted model without an
installed adapter fails closed instead of fabricating output.
"""

from __future__ import annotations

import queue
import threading
import uuid
from pathlib import Path
from typing import Any

try:
    from ..api.serialization import read_json, utc_now, write_json_atomic
except ImportError:
    from api.serialization import read_json, utc_now, write_json_atomic
from .model_registry import ModelAdapter, inspect_model_inputs


TERMINAL = {"succeeded", "failed", "blocked", "cancelled"}


class ModelJobManager:
    def __init__(self, root: Path):
        self.root = Path(root)
        self._queue: queue.Queue[str] = queue.Queue()
        self._adapters: dict[str, ModelAdapter] = {}
        self._lock = threading.RLock()
        self._thread: threading.Thread | None = None

    def register_adapter(self, adapter: ModelAdapter) -> None:
        self._adapters[adapter.manifest_id] = adapter

    def start(self) -> None:
        with self._lock:
            if self._thread and self._thread.is_alive():
                return
            self.root.mkdir(parents=True, exist_ok=True)
            for path in self.root.glob("*/job.json"):
                job = read_json(path)
                if isinstance(job, dict) and job.get("status") == "running":
                    job.update(status="failed", error={"code": "interrupted", "message": "The app restarted during model execution."}, finished_at=utc_now())
                    write_json_atomic(path, job)
                elif isinstance(job, dict) and job.get("status") == "queued":
                    self._queue.put(job["id"])
            self._thread = threading.Thread(target=self._worker, name="nophigene-model-worker", daemon=True)
            self._thread.start()

    def submit(self, *, run_id: str, model_id: str, inputs: dict[str, Any]) -> dict[str, Any]:
        self.start()
        eligibility = inspect_model_inputs(model_id, inputs)
        job_id = uuid.uuid4().hex
        status = "queued" if eligibility.get("eligible") else "blocked"
        job = {
            "id": job_id,
            "run_id": run_id,
            "model_id": model_id,
            "status": status,
            "eligibility": eligibility,
            "progress_percent": 0,
            "created_at": utc_now(),
            "started_at": None,
            "finished_at": utc_now() if status == "blocked" else None,
            "error": None,
        }
        directory = self.root / job_id
        write_json_atomic(directory / "inputs.json", inputs)
        write_json_atomic(directory / "job.json", job)
        if status == "queued":
            self._queue.put(job_id)
        return job

    def get(self, job_id: str) -> dict[str, Any]:
        if not job_id or any(character not in "0123456789abcdef" for character in job_id) or len(job_id) != 32:
            raise ValueError("Invalid model job ID.")
        job = read_json(self.root / job_id / "job.json")
        if not isinstance(job, dict):
            raise FileNotFoundError(job_id)
        return job

    def cancel(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self.get(job_id)
            if job["status"] != "queued":
                raise ValueError("Only queued model jobs can be cancelled.")
            job.update(status="cancelled", finished_at=utc_now())
            write_json_atomic(self.root / job_id / "job.json", job)
            return job

    def _worker(self) -> None:
        while True:
            job_id = self._queue.get()
            try:
                self._execute(job_id)
            finally:
                self._queue.task_done()

    def _execute(self, job_id: str) -> None:
        with self._lock:
            job = self.get(job_id)
            if job["status"] != "queued":
                return
            job.update(status="running", started_at=utc_now(), progress_percent=5)
            write_json_atomic(self.root / job_id / "job.json", job)
        adapter = self._adapters.get(job["model_id"])
        if adapter is None:
            job.update(
                status="failed",
                progress_percent=100,
                finished_at=utc_now(),
                error={
                    "code": "adapter_not_installed",
                    "message": "The model is allowlisted but its versioned execution adapter/assets are not installed.",
                },
            )
            write_json_atomic(self.root / job_id / "job.json", job)
            return
        try:
            inputs = read_json(self.root / job_id / "inputs.json", default={})
            prepared = adapter.prepare(inputs)
            raw = adapter.execute(prepared)
            normalized = adapter.normalize(raw)
            validation = adapter.validate(normalized)
            write_json_atomic(self.root / job_id / "predictions.json", normalized)
            write_json_atomic(self.root / job_id / "validation.json", validation)
            job.update(status="succeeded", progress_percent=100, finished_at=utc_now(), result_url=f"/api/v2/model-jobs/{job_id}/result")
        except Exception as exc:
            job.update(status="failed", progress_percent=100, finished_at=utc_now(), error={"code": "model_execution_failed", "message": str(exc)})
        write_json_atomic(self.root / job_id / "job.json", job)


_DEFAULT_MANAGER: ModelJobManager | None = None


def get_default_model_job_manager() -> ModelJobManager:
    global _DEFAULT_MANAGER
    if _DEFAULT_MANAGER is None:
        _DEFAULT_MANAGER = ModelJobManager(Path(__file__).resolve().parents[2] / "results" / "model-jobs")
        _DEFAULT_MANAGER.start()
    return _DEFAULT_MANAGER
