"""Encrypted, non-exportable credentials for isolated model workers."""

from __future__ import annotations

import json
import os
import secrets
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyzipper


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


class ModelCredentialError(RuntimeError):
    pass


class ModelCredentialStore:
    def __init__(self, root: Path, *, key_file: Path | None = None):
        self.root = Path(root)
        self.key_file = Path(key_file) if key_file else None
        self.root.mkdir(parents=True, exist_ok=True)

    def _key(self) -> bytes:
        configured = self.key_file or (
            Path(os.environ["NOPHIGENE_MODEL_CREDENTIAL_KEY_FILE"])
            if os.environ.get("NOPHIGENE_MODEL_CREDENTIAL_KEY_FILE")
            else None
        )
        if configured:
            key = configured.read_bytes().strip()
        elif os.environ.get("NOPHIGENE_REQUIRE_ENCRYPTION") == "1":
            raise ModelCredentialError("The model credential-encryption key is not configured.")
        else:
            development_key = self.root / ".development-credential-key"
            if not development_key.exists():
                development_key.write_bytes(secrets.token_bytes(48))
            key = development_key.read_bytes().strip()
        if len(key) < 32:
            raise ModelCredentialError("The model credential-encryption key must contain at least 32 bytes.")
        return key

    def _archive(self, model_id: str) -> Path:
        if model_id != "alphagenome-api":
            raise ModelCredentialError("Unknown model credential identifier.")
        return self.root / f"{model_id}.aes.zip"

    def _status_path(self, model_id: str) -> Path:
        return self.root / f"{model_id}.status.json"

    def save(self, model_id: str, secret: str) -> dict[str, Any]:
        value = str(secret or "").strip()
        if len(value) < 8 or len(value) > 4096:
            raise ModelCredentialError("The provider credential has an invalid length.")
        target = self._archive(model_id)
        handle, temporary_name = tempfile.mkstemp(prefix=f".{model_id}-", suffix=".tmp", dir=self.root)
        os.close(handle)
        temporary = Path(temporary_name)
        try:
            with pyzipper.AESZipFile(
                temporary,
                "w",
                compression=pyzipper.ZIP_DEFLATED,
                encryption=pyzipper.WZ_AES,
            ) as archive:
                archive.setpassword(self._key())
                archive.setencryption(pyzipper.WZ_AES, nbits=256)
                archive.writestr("credential.txt", value.encode("utf-8"))
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
        state = {"model_id": model_id, "status": "configured", "updated_at": _now()}
        self._status_path(model_id).write_text(json.dumps(state, sort_keys=True), encoding="utf-8")
        return state

    def read(self, model_id: str) -> str:
        target = self._archive(model_id)
        if not target.is_file():
            raise ModelCredentialError("The model credential is not configured.")
        try:
            with pyzipper.AESZipFile(target, "r") as archive:
                archive.setpassword(self._key())
                value = archive.read("credential.txt").decode("utf-8").strip()
        except Exception as exc:
            raise ModelCredentialError("The encrypted model credential could not be opened.") from exc
        if not value:
            raise ModelCredentialError("The encrypted model credential is empty.")
        return value

    def status(self, model_id: str) -> dict[str, Any]:
        if not self._archive(model_id).is_file():
            return {"model_id": model_id, "status": "missing", "updated_at": None}
        path = self._status_path(model_id)
        if path.is_file():
            try:
                state = json.loads(path.read_text(encoding="utf-8"))
                if isinstance(state, dict):
                    return {
                        "model_id": model_id,
                        "status": str(state.get("status") or "configured"),
                        "updated_at": state.get("updated_at"),
                    }
            except (OSError, ValueError, json.JSONDecodeError):
                pass
        return {"model_id": model_id, "status": "configured", "updated_at": None}

    def mark_status(self, model_id: str, status: str) -> dict[str, Any]:
        if status not in {"configured", "verified", "invalid"}:
            raise ModelCredentialError("Invalid credential status.")
        state = {"model_id": model_id, "status": status, "updated_at": _now()}
        self._status_path(model_id).write_text(json.dumps(state, sort_keys=True), encoding="utf-8")
        return state

    def delete(self, model_id: str) -> dict[str, Any]:
        self._archive(model_id).unlink(missing_ok=True)
        self._status_path(model_id).unlink(missing_ok=True)
        return {"model_id": model_id, "status": "missing", "updated_at": _now()}


_DEFAULT_STORE: ModelCredentialStore | None = None


def get_default_model_credential_store() -> ModelCredentialStore:
    global _DEFAULT_STORE
    if _DEFAULT_STORE is None:
        root = Path(os.environ.get("NOPHIGENE_MODEL_CREDENTIAL_ROOT", "results/model-credentials"))
        key = os.environ.get("NOPHIGENE_MODEL_CREDENTIAL_KEY_FILE")
        _DEFAULT_STORE = ModelCredentialStore(root, key_file=Path(key) if key else None)
    return _DEFAULT_STORE
