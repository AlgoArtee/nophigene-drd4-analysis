"""Container entrypoint that launches either the CLI workflow or the web UI."""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_WEB_PORT = 8766

try:
    from .env import load_dotenv
except ImportError:
    from env import load_dotenv

load_dotenv()

def build_parser() -> argparse.ArgumentParser:
    """Build the top-level launcher parser."""
    parser = argparse.ArgumentParser(
        description="Launch the NophiGene gene-analysis app in CLI or web mode."
    )
    subparsers = parser.add_subparsers(dest="mode")

    web_parser = subparsers.add_parser("web", help="Run the browser-based interface.")
    web_parser.add_argument("--host", default="127.0.0.1", help="Host interface to bind.")
    web_parser.add_argument("--port", default=DEFAULT_WEB_PORT, type=int, help="Port to serve the UI on.")
    web_parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode.")

    cli_parser = subparsers.add_parser("cli", help="Run the original command-line pipeline.")
    cli_parser.add_argument("analysis_args", nargs=argparse.REMAINDER, help="Arguments forwarded to src/analysis.py")

    db_parser = subparsers.add_parser("db", help="Initialize, inspect, migrate, or back up schema-v3 persistence.")
    db_parser.add_argument(
        "action",
        choices=("init", "inventory", "migrate-legacy", "verify-migration", "archive-legacy", "verify-audit", "backup"),
    )
    db_parser.add_argument("--apply", action="store_true", help="Apply a legacy import; without this flag it is a dry run.")
    db_parser.add_argument("--password-file", default="", help="File containing the archive/backup password.")
    db_parser.add_argument("--label", choices=("daily", "weekly"), default="daily", help="Backup rotation class.")

    return parser


def main(argv: list[str] | None = None) -> int:
    """Dispatch to web mode or CLI mode based on the selected subcommand."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.mode in (None, "web"):
        try:
            from .webapp import run_server
        except ImportError:
            from webapp import run_server
        host = getattr(args, "host", "127.0.0.1")
        port = getattr(args, "port", DEFAULT_WEB_PORT)
        debug = getattr(args, "debug", False)
        run_server(host=host, port=port, debug=debug)
        return 0

    if args.mode == "cli":
        try:
            from . import analysis
        except ImportError:
            import analysis
        analysis_args = args.analysis_args
        if analysis_args and analysis_args[0] == "--":
            analysis_args = analysis_args[1:]
        return analysis.main(analysis_args)

    if args.mode == "db":
        try:
            from .workbench.artifacts import ArtifactStore
            from .workbench.audit import verify_audit_chain
            from .workbench.backups import create_encrypted_backup, rotate_backups
            from .workbench.database import DEFAULT_DATABASE_PATH, ensure_schema, get_default_engine, session_scope
            from .workbench.legacy_migration import create_encrypted_legacy_archive, import_legacy_reports, inventory_legacy_stores
        except ImportError:
            from workbench.artifacts import ArtifactStore
            from workbench.audit import verify_audit_chain
            from workbench.backups import create_encrypted_backup, rotate_backups
            from workbench.database import DEFAULT_DATABASE_PATH, ensure_schema, get_default_engine, session_scope
            from workbench.legacy_migration import create_encrypted_legacy_archive, import_legacy_reports, inventory_legacy_stores

        project_root = Path(__file__).resolve().parents[1]
        if args.action == "backup":
            # A pre-migration backup must not initialize or mutate the schema.
            # Read the encrypted database files directly before opening the ORM.
            if not args.password_file:
                parser.error("db backup requires --password-file")
            password = Path(args.password_file).read_text(encoding="utf-8").strip()
            database_path = Path(os.environ.get("NOPHIGENE_DATABASE_PATH") or DEFAULT_DATABASE_PATH)
            result = create_encrypted_backup(
                database_path,
                project_root / "results" / "backups",
                password=password,
                label=args.label,
            )
            result["rotation"] = rotate_backups(project_root / "results" / "backups")
            print(json.dumps(result, indent=2))
            return 0
        engine = get_default_engine()
        if args.action == "init":
            ensure_schema(engine)
            print(json.dumps({"status": "ready", "schema_version": "3.0"}, indent=2))
            return 0
        if args.action in {"inventory", "migrate-legacy"}:
            inventory = inventory_legacy_stores(project_root)
            result: dict[str, object] = {"dry_run": not args.apply, **inventory.to_dict()}
            if args.action == "migrate-legacy" and args.apply:
                with session_scope(engine) as db_session:
                    result["import"] = import_legacy_reports(
                        db_session,
                        inventory,
                        artifact_store=ArtifactStore(project_root / "results" / "artifacts"),
                    )
            print(json.dumps(result, indent=2))
            return 0
        if args.action == "verify-migration":
            with engine.connect() as connection:
                foreign_key_errors = [list(row) for row in connection.exec_driver_sql("PRAGMA foreign_key_check")]
                table_counts = {
                    table: connection.exec_driver_sql(f'SELECT COUNT(*) FROM "{table}"').scalar_one()
                    for table in ("runs", "variant_calls", "methylation_measurements", "evidence_records", "artifacts")
                }
            result = {"valid": not foreign_key_errors, "foreign_key_errors": foreign_key_errors, "table_counts": table_counts}
            print(json.dumps(result, indent=2))
            return 0 if result["valid"] else 1
        if args.action == "archive-legacy":
            if not args.apply:
                parser.error("archive-legacy requires --apply; originals are retained")
            if not args.password_file:
                parser.error("archive-legacy requires --password-file")
            password = Path(args.password_file).read_text(encoding="utf-8").strip()
            inventory = inventory_legacy_stores(project_root)
            timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            output = project_root / "results" / "legacy-archive" / f"legacy-stores-{timestamp}.aes.zip"
            result = create_encrypted_legacy_archive(inventory, output, password=password)
            result["originals_deleted"] = False
            print(json.dumps(result, indent=2))
            return 0
        if args.action == "verify-audit":
            with session_scope(engine) as db_session:
                result = verify_audit_chain(db_session)
            print(json.dumps(result, indent=2))
            return 0 if result["valid"] else 1
    parser.error(f"Unsupported mode: {args.mode}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
