"""Database bootstrap and transaction helpers.

Production containers set ``NOPHIGENE_DATABASE_KEY_FILE`` and therefore use
SQLCipher.  Plain SQLite is accepted only for explicit test/development URLs;
``NOPHIGENE_REQUIRE_ENCRYPTION=1`` makes a missing cipher a startup failure.
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
from typing import Iterator

from sqlalchemy import Engine, create_engine, event, text
from sqlalchemy.orm import Session, sessionmaker

from .models import Base

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATABASE_PATH = PROJECT_ROOT / "results" / "nophigene-v2.db"


class DatabaseSecurityError(RuntimeError):
    """Raised when an encrypted production database cannot be established."""


def _read_key_file(path: Path) -> bytes:
    key = path.read_bytes().strip()
    if len(key) < 32:
        raise DatabaseSecurityError("The SQLCipher key must contain at least 32 bytes.")
    return key


def _sqlcipher_creator(database_path: Path, key_file: Path):
    key_hex = _read_key_file(key_file).hex()

    def connect():
        try:
            from sqlcipher3 import dbapi2 as sqlcipher
        except ImportError as exc:
            raise DatabaseSecurityError(
                "SQLCipher is required but the maintained sqlcipher3 binding is unavailable. Use the supported container image."
            ) from exc
        database_path.parent.mkdir(parents=True, exist_ok=True)
        connection = sqlcipher.connect(str(database_path))
        connection.execute(f"PRAGMA key = \"x'{key_hex}'\"")
        connection.execute("PRAGMA cipher_memory_security = ON")
        connection.execute("PRAGMA foreign_keys = ON")
        connection.execute("SELECT count(*) FROM sqlite_master").fetchone()
        return connection

    return connect


def create_database_engine(
    *,
    database_path: Path | None = None,
    key_file: Path | None = None,
    test_url: str | None = None,
) -> Engine:
    if test_url:
        engine = create_engine(test_url, future=True)
    else:
        path = Path(database_path or os.environ.get("NOPHIGENE_DATABASE_PATH") or DEFAULT_DATABASE_PATH)
        configured_key = key_file or (
            Path(os.environ["NOPHIGENE_DATABASE_KEY_FILE"])
            if os.environ.get("NOPHIGENE_DATABASE_KEY_FILE")
            else None
        )
        require_encryption = os.environ.get("NOPHIGENE_REQUIRE_ENCRYPTION", "0") == "1"
        if configured_key:
            engine = create_engine(
                "sqlite://",
                creator=_sqlcipher_creator(path, configured_key),
                future=True,
                pool_pre_ping=True,
            )
        elif require_encryption:
            raise DatabaseSecurityError(
                "Encrypted persistence is required; NOPHIGENE_DATABASE_KEY_FILE is not configured."
            )
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            engine = create_engine(f"sqlite+pysqlite:///{path.as_posix()}", future=True)

    @event.listens_for(engine, "connect")
    def _configure_sqlite(connection, _record):
        cursor = connection.cursor()
        cursor.execute("PRAGMA foreign_keys = ON")
        cursor.close()

    return engine


def ensure_schema(engine: Engine) -> None:
    """Create an empty schema; deployed upgrades are performed by Alembic."""
    Base.metadata.create_all(engine)
    from .model_registry import list_model_manifests
    from .models import ModelDefinition

    with Session(engine) as session:
        for manifest in list_model_manifests():
            definition = session.get(ModelDefinition, manifest["id"])
            values = {
                "name": manifest["name"],
                "version": manifest["checksum_sha256"][:12],
                "wave": manifest["wave"],
                "task": manifest["task"],
                "execution_mode": manifest["execution_mode"],
                "status": manifest["status"],
                "input_contract": {"required_inputs": list(manifest["required_inputs"])},
                "output_contract": {"semantics": manifest["output_semantics"], "consensus": False},
                "supported_builds": list(manifest["supported_builds"]),
                "resource_requirements": {"gpu": manifest["gpu"], "network_policy": manifest["network_policy"]},
                "license": manifest["license"],
                "citation": manifest["citation"],
                "limitations": list(manifest["limitations"]),
                "manifest_checksum_sha256": manifest["checksum_sha256"],
            }
            if definition is None:
                session.add(ModelDefinition(id=manifest["id"], **values))
            else:
                for key, value in values.items():
                    setattr(definition, key, value)
        session.commit()
    with engine.connect() as connection:
        connection.execute(text("SELECT 1"))


@lru_cache(maxsize=1)
def get_default_engine() -> Engine:
    engine = create_database_engine()
    ensure_schema(engine)
    return engine


def session_factory(engine: Engine | None = None) -> sessionmaker[Session]:
    return sessionmaker(bind=engine or get_default_engine(), expire_on_commit=False, future=True)


@contextmanager
def session_scope(engine: Engine | None = None) -> Iterator[Session]:
    session = session_factory(engine)()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()
