"""Content-addressed artifact storage with immutable checksums."""

from __future__ import annotations

import hashlib
import mimetypes
from pathlib import Path
from typing import BinaryIO

from sqlalchemy.orm import Session

from .models import Artifact


class ArtifactStore:
    def __init__(self, root: Path):
        self.root = Path(root)

    def _target(self, digest: str, suffix: str = "") -> Path:
        safe_suffix = suffix if suffix.startswith(".") and len(suffix) <= 16 else ""
        return self.root / digest[:2] / f"{digest}{safe_suffix}"

    def put_bytes(
        self,
        session: Session,
        content: bytes,
        *,
        kind: str,
        run_id: str | None = None,
        filename: str = "",
        media_type: str = "",
        sensitive: bool = True,
    ) -> Artifact:
        digest = hashlib.sha256(content).hexdigest()
        suffix = Path(filename).suffix.lower() if filename else ""
        target = self._target(digest, suffix)
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            temporary = target.with_suffix(target.suffix + ".tmp")
            temporary.write_bytes(content)
            temporary.replace(target)
        artifact = Artifact(
            run_id=run_id,
            kind=kind,
            relative_path=target.relative_to(self.root).as_posix(),
            checksum_sha256=digest,
            size_bytes=len(content),
            media_type=media_type or mimetypes.guess_type(filename)[0] or "application/octet-stream",
            sensitive=sensitive,
            metadata_json={"original_filename": filename} if filename else {},
        )
        session.add(artifact)
        session.flush()
        return artifact

    def put_stream(self, session: Session, stream: BinaryIO, **kwargs) -> Artifact:
        return self.put_bytes(session, stream.read(), **kwargs)

    def resolve(self, artifact: Artifact) -> Path:
        target = (self.root / artifact.relative_path).resolve()
        target.relative_to(self.root.resolve())
        return target
