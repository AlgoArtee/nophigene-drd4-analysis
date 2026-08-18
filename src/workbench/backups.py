"""Verified encrypted database backup rotation."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def create_encrypted_backup(database_path: Path, backup_dir: Path, *, password: str, label: str = "daily") -> dict[str, Any]:
    if not password:
        raise ValueError("A backup password is required.")
    try:
        import pyzipper
    except ImportError as exc:
        raise RuntimeError("Encrypted backups require pyzipper.") from exc
    database_path = Path(database_path)
    if not database_path.is_file():
        raise FileNotFoundError(database_path)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup_dir = Path(backup_dir)
    backup_dir.mkdir(parents=True, exist_ok=True)
    output = backup_dir / f"nophigene-{label}-{timestamp}.zip"
    members = [database_path]
    for suffix in ("-wal", "-shm"):
        companion = Path(str(database_path) + suffix)
        if companion.is_file():
            members.append(companion)
    manifest = {
        "created_at": timestamp,
        "label": label,
        "database": str(database_path),
        "members": [{"name": path.name, "sha256": _sha256(path), "size_bytes": path.stat().st_size} for path in members],
    }
    with pyzipper.AESZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, encryption=pyzipper.WZ_AES) as archive:
        archive.setpassword(password.encode("utf-8"))
        archive.writestr("backup-manifest.json", json.dumps(manifest, indent=2).encode("utf-8"))
        for path in members:
            archive.write(path, f"database/{path.name}")
    return {**manifest, "path": str(output), "archive_sha256": _sha256(output)}


def rotate_backups(backup_dir: Path, *, keep_daily: int = 7, keep_weekly: int = 4) -> dict[str, list[str]]:
    backup_dir = Path(backup_dir)
    removed: list[str] = []
    retained: list[str] = []
    for label, keep in (("daily", keep_daily), ("weekly", keep_weekly)):
        paths = sorted(backup_dir.glob(f"nophigene-{label}-*.zip"), key=lambda item: item.stat().st_mtime, reverse=True)
        retained.extend(str(path) for path in paths[:keep])
        for path in paths[keep:]:
            path.unlink()
            removed.append(str(path))
    return {"retained": retained, "removed": removed}


def verify_encrypted_backup(archive_path: Path, *, password: str) -> dict[str, Any]:
    """Verify the encrypted archive and every member checksum without restoring."""
    if not password:
        raise ValueError("A backup password is required.")
    import pyzipper

    with pyzipper.AESZipFile(archive_path, "r") as archive:
        archive.setpassword(password.encode("utf-8"))
        manifest = json.loads(archive.read("backup-manifest.json"))
        for member in manifest.get("members", []):
            payload = archive.read(f"database/{member['name']}")
            if hashlib.sha256(payload).hexdigest() != member["sha256"]:
                return {"valid": False, "failed_member": member["name"]}
    return {"valid": True, "member_count": len(manifest.get("members", [])), "manifest": manifest}


def restore_encrypted_backup(
    archive_path: Path,
    destination: Path,
    *,
    password: str,
    confirm_overwrite: bool = False,
) -> dict[str, Any]:
    """Verify, stage, and atomically restore a database backup."""
    verification = verify_encrypted_backup(archive_path, password=password)
    if not verification["valid"]:
        raise ValueError(f"Backup verification failed for {verification.get('failed_member')}.")
    destination = Path(destination)
    if destination.exists() and not confirm_overwrite:
        raise FileExistsError("Destination exists; explicit overwrite confirmation is required.")
    import pyzipper

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
        staged_root = Path(temporary)
        with pyzipper.AESZipFile(archive_path, "r") as archive:
            archive.setpassword(password.encode("utf-8"))
            for member in verification["manifest"].get("members", []):
                staged = staged_root / member["name"]
                staged.write_bytes(archive.read(f"database/{member['name']}"))
        primary_name = verification["manifest"]["members"][0]["name"]
        staged_primary = staged_root / primary_name
        os.replace(staged_primary, destination)
        for suffix in ("-wal", "-shm"):
            staged = staged_root / f"{primary_name}{suffix}"
            if staged.exists():
                os.replace(staged, Path(str(destination) + suffix))
    return {"restored": True, "destination": str(destination), "verification": verification}
