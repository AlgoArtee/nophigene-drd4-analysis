#!/usr/bin/env python3
"""Build the deterministic compressed bundle for generated gene-data artifacts."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENE_DATA_DIR = PROJECT_ROOT / "src" / "gene_data"
BUNDLE_PATH = GENE_DATA_DIR / "gene_data_bundle.zip"
GENERATED_GENE_DATA_SUFFIXES = (
    "_interpretation_db.json",
    "_population_db.json",
    "_epigenetics_hg19.csv",
)
ZIP_TIMESTAMP = (2026, 1, 1, 0, 0, 0)


def _is_generated_gene_data_file(path: Path) -> bool:
    return path.is_file() and path.name.endswith(GENERATED_GENE_DATA_SUFFIXES)


def _normalize_payload(filename: str, payload: bytes) -> bytes:
    """Return stable archive bytes, minifying JSON while preserving CSV content."""
    if filename.endswith(".json"):
        parsed = json.loads(payload.decode("utf-8"))
        return (json.dumps(parsed, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode(
            "utf-8"
        )
    return payload


def _bundle_payload(path: Path) -> bytes:
    """Return stable archive bytes for one loose generated artifact."""
    return _normalize_payload(path.name, path.read_bytes())


def _collect_bundle_payloads(bundle_path: Path) -> dict[str, bytes]:
    """Collect generated payloads from the current bundle plus loose overrides."""
    payloads: dict[str, bytes] = {}
    if bundle_path.exists():
        with zipfile.ZipFile(bundle_path) as bundle:
            for info in sorted(bundle.infolist(), key=lambda item: item.filename):
                if info.is_dir() or not info.filename.endswith(GENERATED_GENE_DATA_SUFFIXES):
                    continue
                payloads[info.filename] = _normalize_payload(info.filename, bundle.read(info.filename))

    for path in sorted(path for path in GENE_DATA_DIR.iterdir() if _is_generated_gene_data_file(path)):
        payloads[path.name] = _bundle_payload(path)

    return payloads


def build_bundle(bundle_path: Path = BUNDLE_PATH) -> dict[str, int]:
    """Write the compressed bundle and return basic size/count metrics."""
    payloads = _collect_bundle_payloads(bundle_path)
    if not payloads:
        raise FileNotFoundError(f"No generated gene-data artifacts found in {GENE_DATA_DIR}")

    tmp_path = bundle_path.with_suffix(bundle_path.suffix + ".tmp")
    raw_size = 0
    with zipfile.ZipFile(tmp_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for filename, payload in sorted(payloads.items()):
            raw_size += len(payload)
            info = zipfile.ZipInfo(filename, ZIP_TIMESTAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            bundle.writestr(info, payload)

    tmp_path.replace(bundle_path)
    return {
        "file_count": len(payloads),
        "raw_size": raw_size,
        "bundle_size": bundle_path.stat().st_size,
    }


def main() -> None:
    metrics = build_bundle()
    print(
        "Built {path} with {file_count} files; raw={raw_size:,} bytes, bundle={bundle_size:,} bytes.".format(
            path=BUNDLE_PATH,
            **metrics,
        )
    )


if __name__ == "__main__":
    main()
