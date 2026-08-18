"""DANDELION cohort-analysis contracts and report normalization.

DANDELION is deliberately treated as an offline statistical method.  Its
outputs describe disease-association prioritisation hypotheses; they are not
clinical assertions, individual-level predictions, or proof of causality.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable


METHOD_ID = "dandelion"
PACKAGE_NAME = "DANDELION"
PACKAGE_VERSION = "0.1.0"
RUNNER_CONTRACT_VERSION = "1.0"
SUPPORTED_EXPOSURES = {"Gene", "SNP"}
SUPPORTED_FORMATS = {"csv", "tsv", "rds"}
REQUIRED_FILE_ROLES = {"trans_matrix", "gene_association", "gene_annotation"}


class DandelionValidationError(ValueError):
    """Raised when a cohort dataset or analysis contract is invalid."""


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative_file(import_root: Path, value: object) -> tuple[str, Path]:
    raw = str(value or "").strip().replace("\\", "/")
    candidate = Path(raw)
    if not raw or candidate.is_absolute() or ".." in candidate.parts:
        raise DandelionValidationError("Dataset paths must be non-empty paths relative to the configured import root.")
    root = import_root.resolve(strict=True)
    unresolved = root / candidate
    current = root
    for part in candidate.parts:
        current = current / part
        if current.is_symlink():
            raise DandelionValidationError("Dataset paths must not contain symbolic links.")
    resolved = unresolved.resolve(strict=True)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise DandelionValidationError("Dataset path escapes the configured import root.") from exc
    if not resolved.is_file():
        raise DandelionValidationError("Dataset paths must identify regular, non-symlink files.")
    return resolved.relative_to(root).as_posix(), resolved


def _format_for(path: Path, configured: object) -> str:
    value = str(configured or path.suffix.lstrip(".")).lower()
    if value not in SUPPORTED_FORMATS:
        raise DandelionValidationError(f"Unsupported dataset format '{value}'; use CSV, TSV, or RDS.")
    return value


def _inspect_tabular(path: Path, file_format: str) -> dict[str, Any]:
    delimiter = "," if file_format == "csv" else "\t"
    rows = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise DandelionValidationError(f"Dataset file is empty: {path.name}") from exc
        header = [str(item).strip() for item in header]
        if not header or any(not item for item in header):
            raise DandelionValidationError(f"Dataset file has an invalid header: {path.name}")
        if len(set(header)) != len(header):
            raise DandelionValidationError(f"Dataset file has duplicate column names: {path.name}")
        for _ in reader:
            rows += 1
    return {"columns": header, "column_count": len(header), "row_count": rows}


def normalize_dataset_manifest(payload: dict[str, Any], import_root: Path) -> dict[str, Any]:
    """Validate a local cohort manifest and return immutable file metadata."""
    if not isinstance(payload, dict):
        raise DandelionValidationError("Dataset manifest must be a JSON object.")
    name = str(payload.get("name") or "").strip()
    if not name:
        raise DandelionValidationError("Dataset name is required.")
    phenotype = str(payload.get("phenotype") or "").strip()
    if not phenotype:
        raise DandelionValidationError("A structured phenotype identifier or label is required.")
    exposure_type = str(payload.get("exposure_type") or "").strip()
    if exposure_type.lower() == "gene":
        exposure_type = "Gene"
    elif exposure_type.lower() == "snp":
        exposure_type = "SNP"
    if exposure_type not in SUPPORTED_EXPOSURES:
        raise DandelionValidationError("exposure_type must be 'Gene' or 'SNP'.")
    build = str(payload.get("assembly") or payload.get("build") or "").strip()
    if build not in {"GRCh37", "GRCh38"}:
        raise DandelionValidationError("assembly must be GRCh37 or GRCh38.")
    namespace = str(payload.get("gene_namespace") or "").strip()
    if not namespace:
        raise DandelionValidationError("gene_namespace is required (for example HGNC symbol or Ensembl Gene ID).")

    raw_files = payload.get("files")
    if not isinstance(raw_files, dict):
        raise DandelionValidationError("files must be an object keyed by dataset role.")
    missing = sorted(REQUIRED_FILE_ROLES - set(raw_files))
    if exposure_type == "SNP" and "snp_reference" not in raw_files:
        missing.append("snp_reference")
    if missing:
        raise DandelionValidationError("Missing required dataset files: " + ", ".join(missing))

    normalized_files: dict[str, dict[str, Any]] = {}
    manifest_hash = hashlib.sha256()
    for role, spec in sorted(raw_files.items()):
        if not isinstance(spec, dict):
            spec = {"path": spec}
        relative_path, resolved = _relative_file(import_root, spec.get("path"))
        file_format = _format_for(resolved, spec.get("format"))
        checksum = sha256_file(resolved)
        inspection = {} if file_format == "rds" else _inspect_tabular(resolved, file_format)
        item = {
            "role": str(role),
            "relative_path": relative_path,
            "format": file_format,
            "checksum_sha256": checksum,
            "size_bytes": resolved.stat().st_size,
            "mapping": dict(spec.get("mapping") or {}),
            "inspection": inspection,
        }
        normalized_files[str(role)] = item
        manifest_hash.update(canonical_json(item))

    context = dict(payload.get("context") or {})
    missing_context = [field for field in ("tissue", "ancestry") if not context.get(field)]
    applicability = "complete" if not missing_context else "limited"
    normalized = {
        "contract_version": RUNNER_CONTRACT_VERSION,
        "method": METHOD_ID,
        "name": name,
        "description": str(payload.get("description") or "").strip(),
        "phenotype": phenotype,
        "exposure_type": exposure_type,
        "assembly": build,
        "gene_namespace": namespace,
        "storage_mode": str(payload.get("storage_mode") or "registered_path"),
        "context": context,
        "applicability": {"status": applicability, "missing_context": missing_context},
        "files": normalized_files,
    }
    manifest_hash.update(canonical_json({key: value for key, value in normalized.items() if key != "files"}))
    normalized["manifest_checksum_sha256"] = manifest_hash.hexdigest()
    return normalized


def normalize_analysis_parameters(payload: dict[str, Any] | None) -> dict[str, Any]:
    values = dict(payload or {})
    try:
        target_fdr = float(values.get("target_fdr", 0.1))
        cis_window_bp = int(values.get("cis_window_bp", 5_000_000))
        disease_p_threshold = float(values.get("gene_association_threshold", 1.0))
        chunk_size = int(values.get("chunk_size", 250))
    except (TypeError, ValueError) as exc:
        raise DandelionValidationError("DANDELION analysis parameters must be numeric.") from exc
    if not 0 < target_fdr <= 1:
        raise DandelionValidationError("target_fdr must be greater than zero and at most one.")
    if not 0 <= cis_window_bp <= 100_000_000:
        raise DandelionValidationError("cis_window_bp must be between 0 and 100,000,000.")
    if not 0 < disease_p_threshold <= 1:
        raise DandelionValidationError("gene_association_threshold must be greater than zero and at most one.")
    if not 10 <= chunk_size <= 10_000:
        raise DandelionValidationError("chunk_size must be between 10 and 10,000 exposures.")
    return {
        "target_fdr": target_fdr,
        "cis_window_bp": cis_window_bp,
        "gene_association_threshold": disease_p_threshold,
        "chunk_size": chunk_size,
        "multiple_testing": "Benjamini-Hochberg within each exposure",
    }


def estimate_resources(dataset: dict[str, Any], parameters: dict[str, Any]) -> dict[str, Any]:
    trans = dataset["files"]["trans_matrix"]
    rows = int(trans.get("inspection", {}).get("row_count") or 0)
    columns = int(trans.get("inspection", {}).get("column_count") or 0)
    size = int(trans.get("size_bytes") or 0)
    rds = trans.get("format") == "rds"
    return {
        "execution": "offline_cpu",
        "network_required": False,
        "input_bytes": sum(int(item.get("size_bytes") or 0) for item in dataset["files"].values()),
        "trans_rows": rows or None,
        "trans_columns": columns or None,
        "chunk_size": parameters["chunk_size"],
        "memory_note": (
            "RDS input is loaded as one object; confirm available memory before execution."
            if rds
            else "Delimited trans statistics are read in exposure chunks."
        ),
        "estimated_peak_bytes": max(size * (4 if rds else 2), 256 * 1024 * 1024),
    }


def _finite_number(value: object) -> float | None:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def build_report(
    *,
    analysis_id: str,
    dataset: dict[str, Any],
    parameters: dict[str, Any],
    normalized_result: dict[str, Any],
    job: dict[str, Any],
) -> dict[str, Any]:
    """Map an offline runner result into canonical report schema 3.0."""
    records = list(normalized_result.get("records") or [])
    records.sort(
        key=lambda row: (
            _finite_number(row.get("q_value")) if _finite_number(row.get("q_value")) is not None else 2.0,
            str(row.get("candidate_gene") or ""),
        )
    )
    significant = [row for row in records if bool(row.get("significant"))]
    edges = []
    for row in significant:
        exposure = str(row.get("source_gene") or row.get("exposure") or "")
        candidate = str(row.get("candidate_gene") or "")
        if not exposure or not candidate:
            continue
        edges.append(
            {
                "source_gene": exposure,
                "target_gene": candidate,
                "edge_type": "trans_regulatory_disease_prioritization",
                "directed": True,
                "source_key": "dandelion",
                "native_score": _finite_number(row.get("q_value")),
                "native_score_label": "BH q-value (lower is stronger)",
                "evidence_count": 1,
                "tissue": str(dataset.get("context", {}).get("tissue") or ""),
                "details": {
                    "hypothesis_generating": True,
                    "original_exposure": str(row.get("exposure") or ""),
                    "source_node_type": str(row.get("source_node_type") or "gene"),
                    "trans_p_value": _finite_number(row.get("trans_p_value")),
                    "gene_association_p_value": _finite_number(row.get("gene_association_p_value")),
                    "dandelion_p_value": _finite_number(row.get("p_value")),
                    "q_value": _finite_number(row.get("q_value")),
                },
            }
        )
    limitations = [
        "DANDELION prioritisation is hypothesis-generating and does not establish causality or medical relevance.",
        "Results depend on the compatibility, power, ancestry, tissue, and phenotype definitions of the supplied cohorts.",
    ]
    limitations.extend(str(item) for item in normalized_result.get("limitations") or [])
    if dataset.get("applicability", {}).get("status") != "complete":
        limitations.append("Tissue and/or ancestry context is missing; applicability is limited.")
    return {
        "schema_version": "3.0",
        "run": {
            "id": analysis_id,
            "kind": "cohort_statistical_analysis",
            "method": METHOD_ID,
            "status": job.get("status", "completed"),
            "dataset_id": dataset.get("id"),
            "phenotype": dataset.get("phenotype"),
            "assembly": dataset.get("assembly"),
        },
        "summary": {
            "title": f"DANDELION prioritisation: {dataset.get('name', '')}",
            "candidate_count": len(significant),
            "tested_count": int(normalized_result.get("tested_pair_count") or len(records)),
            "returned_record_count": len(records),
            "interpretation": "Cohort-level disease-gene prioritisation hypotheses",
            "research_use_only": True,
        },
        "sections": {
            "objective_data": {
                "status": "not_applicable",
                "reason": "This is a cohort-level association analysis, not a sample-observation run.",
                "dataset": {key: dataset.get(key) for key in ("name", "phenotype", "exposure_type", "assembly", "gene_namespace")},
            },
            "statistics": {
                "status": "assessed",
                "method": {"name": PACKAGE_NAME, "version": PACKAGE_VERSION, "parameters": parameters},
                "applicability": dataset.get("applicability", {}),
                "results": records,
                "significant_results": significant,
                "test_family": "one Benjamini-Hochberg family per exposure",
                "limitations": limitations,
            },
            "literature": {"status": "not_assessed", "reason": "Literature retrieval was not part of this statistical job."},
            "medical": {"status": "not_assessed", "reason": "Statistical prioritisation cannot populate established medical evidence."},
            "interactions": {
                "status": "assessed" if edges else "assessed_absence",
                "semantics": "directed trans-regulatory disease-prioritisation hypotheses",
                "score_direction": "Lower q-values indicate stronger statistical support; values are not probabilities.",
                "edges": edges,
            },
            "predictions": {"status": "not_requested", "reason": "DANDELION is a statistical method, not an AI prediction model."},
        },
        "run_details": {
            "runner_contract_version": RUNNER_CONTRACT_VERSION,
            "package": {"name": PACKAGE_NAME, "version": PACKAGE_VERSION, "license": "GPL-3.0"},
            "dataset_manifest_checksum_sha256": dataset.get("manifest_checksum_sha256"),
            "job": job,
            "artifacts": list(normalized_result.get("artifacts") or []),
            "warnings": list(normalized_result.get("warnings") or []),
        },
    }


def validate_runner_result(payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("contract_version") != RUNNER_CONTRACT_VERSION:
        raise DandelionValidationError("Runner result contract version is unsupported.")
    records = payload.get("records")
    if not isinstance(records, list):
        raise DandelionValidationError("Runner result records must be a list.")
    for index, row in enumerate(records):
        if not isinstance(row, dict):
            raise DandelionValidationError(f"Runner record {index} is not an object.")
        for field in ("trans_p_value", "gene_association_p_value", "p_value", "q_value"):
            value = _finite_number(row.get(field))
            if value is not None and not 0 <= value <= 1:
                raise DandelionValidationError(f"Runner record {index} has invalid {field}.")
    return payload
