"""Persistent single-worker job queue for local API workflows."""

from __future__ import annotations

import queue
import re
import threading
import time
import uuid
import logging
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable

from .errors import APIError
from .profiles import ProfileStore, get_default_profile_store
from .serialization import read_json, utc_now, write_json_atomic
from .workflow_runner import WorkflowRunner

try:
    from ..analysis import PROJECT_ROOT
except ImportError:
    from analysis import PROJECT_ROOT

DEFAULT_JOBS_ROOT = PROJECT_ROOT / "results" / "api" / "jobs"
JOB_ID_PATTERN = re.compile(r"^[a-f0-9]{32}$")
TERMINAL_STATUSES = {"succeeded", "partial", "failed", "cancelled"}
logger = logging.getLogger(__name__)


class JobManager:
    """Manage persisted jobs and execute them through one background worker."""

    def __init__(
        self,
        *,
        jobs_root: Path = DEFAULT_JOBS_ROOT,
        profile_store: ProfileStore | None = None,
        runner: WorkflowRunner | None = None,
    ):
        self.jobs_root = Path(jobs_root)
        self.profile_store = profile_store or get_default_profile_store()
        self.runner = runner or WorkflowRunner(self.profile_store, self.jobs_root)
        self._queue: queue.Queue[str] = queue.Queue()
        self._lock = threading.RLock()
        self._thread: threading.Thread | None = None
        self._started = False
        self._completion_hooks: list[Callable[[str, dict[str, Any]], None]] = []

    def add_completion_hook(self, hook: Callable[[str, dict[str, Any]], None]) -> None:
        """Register an idempotent post-completion hook for successful/partial jobs."""
        with self._lock:
            if hook not in self._completion_hooks:
                self._completion_hooks.append(hook)

    @property
    def worker_alive(self) -> bool:
        return bool(self._thread and self._thread.is_alive())

    @property
    def queue_depth(self) -> int:
        return self._queue.qsize()

    def start(self) -> None:
        with self._lock:
            if self._started and self.worker_alive:
                return
            self.jobs_root.mkdir(parents=True, exist_ok=True)
            self._recover_jobs()
            self._thread = threading.Thread(
                target=self._worker_loop,
                name="nophigene-api-worker",
                daemon=True,
            )
            self._thread.start()
            self._started = True

    def submit(self, request_payload: dict[str, Any]) -> dict[str, Any]:
        self.start()
        if request_payload.get("profile_id"):
            self.profile_store.get(request_payload["profile_id"])
        if request_payload.get("source_job_id"):
            source = self.get(request_payload["source_job_id"])
            if source["status"] not in {"succeeded", "partial"}:
                raise APIError(
                    "source_job_not_ready",
                    f"Source job '{source['id']}' is not complete.",
                    409,
                )

        job_id = uuid.uuid4().hex
        now = utc_now()
        job = {
            "id": job_id,
            "operation": request_payload["operation"],
            "genes": request_payload["genes"],
            "status": "queued",
            "stage": "queued",
            "progress": {
                "completed": 0,
                "total": len(request_payload["genes"]),
                "percent": 0,
            },
            "outcomes": [],
            "error": None,
            "created_at": now,
            "updated_at": now,
            "started_at": None,
            "finished_at": None,
            "result_url": f"/api/v1/jobs/{job_id}/result",
        }
        job_dir = self.jobs_root / job_id
        write_json_atomic(job_dir / "request.json", request_payload)
        write_json_atomic(job_dir / "job.json", job)
        self._queue.put(job_id)
        return deepcopy(job)

    def list(self) -> list[dict[str, Any]]:
        self.jobs_root.mkdir(parents=True, exist_ok=True)
        jobs: list[dict[str, Any]] = []
        for path in self.jobs_root.glob("*/job.json"):
            payload = read_json(path)
            if isinstance(payload, dict):
                jobs.append(payload)
        return sorted(jobs, key=lambda item: item.get("created_at", ""), reverse=True)

    def get(self, job_id: str) -> dict[str, Any]:
        self._validate_job_id(job_id)
        job = read_json(self.jobs_root / job_id / "job.json")
        if not isinstance(job, dict):
            raise APIError("job_not_found", f"Job '{job_id}' was not found.", 404)
        return job

    def result(self, job_id: str) -> dict[str, Any]:
        job = self.get(job_id)
        if job["status"] not in TERMINAL_STATUSES:
            raise APIError("job_not_finished", f"Job '{job_id}' has not finished.", 409)
        result = read_json(self.jobs_root / job_id / "result.json")
        if not isinstance(result, dict):
            if job["status"] in {"failed", "cancelled"}:
                return {
                    "schema_version": "1.0",
                    "job_id": job_id,
                    "operation": job["operation"],
                    "status": job["status"],
                    "counts": {
                        "requested": len(job["genes"]),
                        "succeeded": 0,
                        "failed": len(job["genes"]),
                    },
                    "genes": job.get("outcomes", []),
                    "error": job.get("error"),
                }
            raise APIError("job_result_missing", f"Job '{job_id}' has no result manifest.", 500)
        return result

    def cancel(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self.get(job_id)
            if job["status"] != "queued":
                raise APIError(
                    "job_not_cancellable",
                    "Only queued jobs can be cancelled.",
                    409,
                )
            now = utc_now()
            job.update(
                {
                    "status": "cancelled",
                    "stage": "cancelled",
                    "updated_at": now,
                    "finished_at": now,
                }
            )
            self._write_job(job)
            return job

    def wait_for_terminal(self, job_id: str, timeout: float = 10.0) -> dict[str, Any]:
        """Wait for a job in tests and local integrations."""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            job = self.get(job_id)
            if job["status"] in TERMINAL_STATUSES:
                return job
            time.sleep(0.02)
        raise TimeoutError(f"Job '{job_id}' did not finish within {timeout} seconds.")

    def _recover_jobs(self) -> None:
        for job in self.list():
            if job["status"] == "running":
                now = utc_now()
                job.update(
                    {
                        "status": "failed",
                        "stage": "interrupted",
                        "updated_at": now,
                        "finished_at": now,
                        "error": {
                            "code": "interrupted",
                            "message": "The application restarted while this job was running.",
                        },
                    }
                )
                self._write_job(job)
            elif job["status"] == "queued":
                self._queue.put(job["id"])

    def _worker_loop(self) -> None:
        while True:
            job_id = self._queue.get()
            try:
                self._execute_job(job_id)
            finally:
                self._queue.task_done()

    def _execute_job(self, job_id: str) -> None:
        with self._lock:
            job = self.get(job_id)
            if job["status"] != "queued":
                return
            now = utc_now()
            job.update(
                {
                    "status": "running",
                    "stage": "starting",
                    "started_at": now,
                    "updated_at": now,
                }
            )
            self._write_job(job)

        job_dir = self.jobs_root / job_id
        request_payload = read_json(job_dir / "request.json")

        def progress(
            stage: str,
            completed: int,
            total: int,
            outcomes: list[dict[str, Any]],
        ) -> None:
            with self._lock:
                current = self.get(job_id)
                current["stage"] = stage
                current["progress"] = {
                    "completed": completed,
                    "total": total,
                    "percent": round((completed / total) * 100) if total else 100,
                }
                current["outcomes"] = outcomes
                current["updated_at"] = utc_now()
                self._write_job(current)

        try:
            result = self.runner.execute(
                job_id,
                request_payload,
                job_dir,
                progress,
            )
            with self._lock:
                job = self.get(job_id)
                now = utc_now()
                job.update(
                    {
                        "status": result["status"],
                        "stage": "complete",
                        "progress": {
                            "completed": len(job["genes"]),
                            "total": len(job["genes"]),
                            "percent": 100,
                        },
                        "outcomes": result["genes"],
                        "updated_at": now,
                        "finished_at": now,
                    }
                )
                self._write_job(job)
            for hook in list(self._completion_hooks):
                try:
                    hook(job_id, deepcopy(result))
                except Exception:
                    logger.exception("Job completion hook failed for %s", job_id)
        except Exception as exc:
            with self._lock:
                job = self.get(job_id)
                now = utc_now()
                job.update(
                    {
                        "status": "failed",
                        "stage": "failed",
                        "updated_at": now,
                        "finished_at": now,
                        "error": {
                            "code": "job_execution_failed",
                            "message": str(exc),
                        },
                    }
                )
                self._write_job(job)

    def _write_job(self, job: dict[str, Any]) -> None:
        write_json_atomic(self.jobs_root / job["id"] / "job.json", job)

    @staticmethod
    def _validate_job_id(job_id: str) -> None:
        if not JOB_ID_PATTERN.fullmatch(str(job_id or "")):
            raise APIError("job_not_found", f"Job '{job_id}' was not found.", 404)


_DEFAULT_JOB_MANAGER: JobManager | None = None


def get_default_job_manager() -> JobManager:
    global _DEFAULT_JOB_MANAGER
    if _DEFAULT_JOB_MANAGER is None:
        _DEFAULT_JOB_MANAGER = JobManager()
    return _DEFAULT_JOB_MANAGER
