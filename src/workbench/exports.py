"""Full-fidelity schema-3 export bundles with optional AES encryption."""

from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path
from typing import Any

from .reporting import render_evidence_first_html


def _section_csv(section: Any) -> bytes:
    rows: list[dict[str, Any]] = []
    if isinstance(section, list):
        rows = [item if isinstance(item, dict) else {"value": item} for item in section]
    elif isinstance(section, dict):
        for key, value in section.items():
            if isinstance(value, list) and value and all(isinstance(item, dict) for item in value):
                rows.extend({"record_group": key, **item} for item in value)
        if not rows:
            rows = [{"field": key, "value": json.dumps(value, ensure_ascii=False, default=str) if isinstance(value, (dict, list)) else value} for key, value in section.items()]
    if not rows:
        return b""
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(
            {
                key: json.dumps(value, ensure_ascii=False, default=str) if isinstance(value, (dict, list)) else value
                for key, value in row.items()
            }
        )
    return buffer.getvalue().encode("utf-8-sig")


def create_export_bundle(
    report: dict[str, Any],
    output_path: Path,
    *,
    password: str = "",
    additional_files: dict[str, bytes] | None = None,
) -> dict[str, Any]:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    files: dict[str, bytes] = {
        "report.html": render_evidence_first_html(report).encode("utf-8"),
        "report.json": json.dumps(report, indent=2, ensure_ascii=False, default=str).encode("utf-8"),
        "SENSITIVE_DATA_WARNING.txt": (
            "This full-fidelity NophiGene export can contain genetic, methylation, sample-context, "
            "path, filename, and identifier data. Handle it as sensitive research data.\n"
        ).encode("utf-8"),
    }
    for key, section in report.get("sections", {}).items():
        files[f"sections/{key}.csv"] = _section_csv(section)
    for name, content in (additional_files or {}).items():
        safe = Path(name)
        if safe.is_absolute() or ".." in safe.parts:
            raise ValueError("Additional export paths must stay inside the bundle.")
        files[safe.as_posix()] = bytes(content)
    encrypted = bool(password)
    if encrypted:
        try:
            import pyzipper
        except ImportError as exc:
            raise RuntimeError("Encrypted exports require pyzipper in the supported app container.") from exc
        with pyzipper.AESZipFile(
            output_path,
            "w",
            compression=zipfile.ZIP_DEFLATED,
            encryption=pyzipper.WZ_AES,
        ) as archive:
            archive.setpassword(password.encode("utf-8"))
            for name, content in files.items():
                archive.writestr(name, content)
    else:
        with zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, content in files.items():
                archive.writestr(name, content)
    return {
        "path": str(output_path),
        "encrypted": encrypted,
        "full_fidelity": True,
        "sensitive": True,
        "files": sorted(files),
    }
