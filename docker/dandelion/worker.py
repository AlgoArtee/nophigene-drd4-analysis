"""Minimal privilege-separated worker for signed DANDELION jobs."""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyzipper


WORK_ROOT = Path(os.environ.get("NOPHIGENE_DANDELION_WORK_ROOT", "/work"))
IMPORT_ROOT = Path(os.environ.get("NOPHIGENE_DANDELION_IMPORT_ROOT", "/imports"))
KEY_FILE = Path(os.environ.get("NOPHIGENE_DANDELION_RUNNER_KEY_FILE", "/run/secrets/dandelion_runner_key"))
RUNNER = Path(os.environ.get("NOPHIGENE_DANDELION_R_RUNNER", "/opt/nophigene/runner.R"))
ARTIFACT_KEY_FILE = Path(
    os.environ.get("NOPHIGENE_DANDELION_ARTIFACT_KEY_FILE", "/run/secrets/dandelion_artifact_key")
)
TIMEOUT_SECONDS = int(os.environ.get("NOPHIGENE_DANDELION_TIMEOUT_SECONDS", "86400"))


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_bytes(canonical_json(value))
    temporary.replace(path)


def signing_key() -> bytes:
    key = KEY_FILE.read_bytes().strip()
    if len(key) < 32:
        raise RuntimeError("Runner signing key is missing or too short.")
    return key


def verify(envelope: dict[str, Any]) -> dict[str, Any]:
    payload = envelope.get("payload")
    signature = str(envelope.get("signature") or "")
    if not isinstance(payload, dict) or envelope.get("algorithm") != "HMAC-SHA256":
        raise RuntimeError("Invalid signed job envelope.")
    expected = hmac.new(signing_key(), canonical_json(payload), hashlib.sha256).hexdigest()
    if not hmac.compare_digest(signature, expected):
        raise RuntimeError("Job signature is invalid.")
    if payload.get("contract_version") != "1.0":
        raise RuntimeError("Unsupported runner contract version.")
    return payload


def safe_import(relative_path: str) -> Path:
    relative = Path(str(relative_path).replace("\\", "/"))
    if relative.is_absolute() or ".." in relative.parts:
        raise RuntimeError("Unsafe import path in signed job.")
    root = IMPORT_ROOT.resolve(strict=True)
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise RuntimeError("Signed import path contains a symbolic link.")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if not path.is_file():
        raise RuntimeError("Signed import path is not a regular file.")
    return path


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_work_file(relative_path: str) -> Path:
    relative = Path(str(relative_path).replace("\\", "/"))
    if relative.is_absolute() or ".." in relative.parts:
        raise RuntimeError("Unsafe managed-artifact path in signed job.")
    root = WORK_ROOT.resolve(strict=True)
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise RuntimeError("Managed-artifact path contains a symbolic link.")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if not path.is_file():
        raise RuntimeError("Managed artifact is not a regular file.")
    return path


def extract_managed(item: dict[str, Any], job_id: str) -> Path:
    archive = safe_work_file(item["managed_path"])
    expected_archive = str((item.get("inspection") or {}).get("managed_archive_sha256") or "")
    if not expected_archive or checksum(archive) != expected_archive:
        raise RuntimeError(f"Managed archive checksum failed: {item['role']}")
    key = ARTIFACT_KEY_FILE.read_bytes().strip()
    if len(key) < 32:
        raise RuntimeError("Managed-artifact encryption key is missing or too short.")
    target_root = Path("/tmp") / f"dandelion-{job_id}"
    target_root.mkdir(mode=0o700, parents=True, exist_ok=True)
    target = target_root / f"{item['role']}.{item['format']}"
    with pyzipper.AESZipFile(archive) as bundle:
        bundle.setpassword(key)
        members = [name for name in bundle.namelist() if not name.endswith("/")]
        if len(members) != 1 or Path(members[0]).name != members[0]:
            raise RuntimeError("Managed archive must contain exactly one flat file.")
        with bundle.open(members[0]) as source, target.open("wb") as output:
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                output.write(chunk)
    if checksum(target) != item["checksum_sha256"]:
        target.unlink(missing_ok=True)
        raise RuntimeError(f"Decrypted input checksum failed: {item['role']}")
    return target


def prepare_runtime_request(payload: dict[str, Any], job_dir: Path) -> Path:
    request = json.loads(json.dumps(payload))
    for item in request["dataset"]["files"].values():
        if request["dataset"].get("storage_mode") == "managed_encrypted_copy":
            path = extract_managed(item, payload["job_id"])
        else:
            path = safe_import(item["relative_path"])
        if checksum(path) != item["checksum_sha256"]:
            raise RuntimeError(f"Input checksum changed after registration: {item['role']}")
        item["runtime_path"] = str(path)
    request["job_dir"] = str(job_dir)
    path = job_dir / "runtime-request.json"
    write_json(path, request)
    return path


def state_for(payload: dict[str, Any], **updates: Any) -> dict[str, Any]:
    state = {
        "id": payload["job_id"],
        "analysis_id": payload["analysis_id"],
        "status": "running",
        "stage": "validating_inputs",
        "progress_percent": 1,
        "submitted_at": payload["submitted_at"],
        "updated_at": now(),
        "cancel_requested": False,
        "runner": {"package": "DANDELION", "version": "0.1.0", "network": "disabled"},
    }
    state.update(updates)
    return state


def terminate(process: subprocess.Popen[bytes]) -> None:
    process.send_signal(signal.SIGTERM)
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=5)


def run_job(queue_path: Path) -> None:
    envelope = json.loads(queue_path.read_text(encoding="utf-8"))
    payload = verify(envelope)
    job_dir = WORK_ROOT / "jobs" / payload["job_id"]
    job_dir.mkdir(parents=True, exist_ok=True)
    status_path = job_dir / "status.json"
    cancel_path = job_dir / "cancel.requested"
    state = state_for(payload)
    write_json(status_path, state)
    started = time.monotonic()
    try:
        if cancel_path.exists():
            write_json(status_path, state_for(payload, status="cancelled", stage="cancelled", progress_percent=0, cancel_requested=True))
            return
        request_path = prepare_runtime_request(payload, job_dir)
        state.update(stage="executing_dandelion", progress_percent=5, updated_at=now())
        write_json(status_path, state)
        with (job_dir / "stdout.log").open("wb") as stdout, (job_dir / "stderr.log").open("wb") as stderr:
            process = subprocess.Popen(
                ["Rscript", "--vanilla", str(RUNNER), str(request_path), str(job_dir / "result.json")],
                stdout=stdout,
                stderr=stderr,
                cwd=job_dir,
                env={**os.environ, "R_DEFAULT_PACKAGES": "datasets,utils,grDevices,graphics,stats,methods"},
            )
            while process.poll() is None:
                if cancel_path.exists():
                    terminate(process)
                    write_json(
                        status_path,
                        state_for(payload, status="cancelled", stage="cancelled", progress_percent=state["progress_percent"], cancel_requested=True),
                    )
                    return
                if time.monotonic() - started > TIMEOUT_SECONDS:
                    terminate(process)
                    raise TimeoutError(f"DANDELION exceeded its {TIMEOUT_SECONDS}-second timeout.")
                progress_path = job_dir / "progress.json"
                if progress_path.exists():
                    try:
                        progress = json.loads(progress_path.read_text(encoding="utf-8"))
                        state.update(
                            stage=str(progress.get("stage") or state["stage"]),
                            progress_percent=max(5, min(99, int(progress.get("progress_percent") or 5))),
                            updated_at=now(),
                        )
                        write_json(status_path, state)
                    except (OSError, ValueError, json.JSONDecodeError):
                        pass
                time.sleep(1)
        if process.returncode != 0:
            raise RuntimeError(f"R runner exited with code {process.returncode}; see Run Details logs.")
        result_path = job_dir / "result.json"
        if not result_path.is_file():
            raise RuntimeError("R runner did not create result.json.")
        write_json(
            status_path,
            state_for(
                payload,
                status="completed",
                stage="completed",
                progress_percent=100,
                finished_at=now(),
                result_checksum_sha256=checksum(result_path),
            ),
        )
    except Exception as exc:
        write_json(
            status_path,
            state_for(
                payload,
                status="failed",
                stage="failed",
                progress_percent=state.get("progress_percent", 0),
                finished_at=now(),
                error={"type": type(exc).__name__, "message": str(exc)},
            ),
        )


def heartbeat() -> None:
    write_json(
        WORK_ROOT / "worker-heartbeat.json",
        {
            "status": "ready",
            "updated_at": now(),
            "package": "DANDELION",
            "package_version": "0.1.0",
            "r_version": "4.6.1",
            "network": "disabled",
        },
    )


def main() -> None:
    (WORK_ROOT / "queue").mkdir(parents=True, exist_ok=True)
    (WORK_ROOT / "jobs").mkdir(parents=True, exist_ok=True)
    signing_key()
    while True:
        heartbeat()
        candidates = sorted((WORK_ROOT / "queue").glob("*.running")) + sorted((WORK_ROOT / "queue").glob("*.json"))
        if not candidates:
            time.sleep(2)
            continue
        source = candidates[0]
        running = source if source.suffix == ".running" else source.with_suffix(".running")
        if source != running:
            try:
                source.replace(running)
            except FileNotFoundError:
                continue
        try:
            run_job(running)
            running.unlink(missing_ok=True)
        except Exception as exc:
            # An unauthenticated or malformed queue item is never executed and
            # must not crash-loop the worker on every restart.
            rejected = running.with_suffix(".rejected")
            running.replace(rejected)
            write_json(
                WORK_ROOT / "worker-heartbeat.json",
                {
                    "status": "degraded",
                    "updated_at": now(),
                    "reason": f"Rejected a queue item: {type(exc).__name__}",
                    "network": "disabled",
                },
            )


if __name__ == "__main__":
    main()
