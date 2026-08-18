"""Contracts for versioned variant and methylation preprocessing workers."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class MethylationPreprocessingRequest:
    idat_prefix: str
    pipeline: str = "sesame"
    platform: str = "Illumina EPIC"
    declared_tissue: str = ""
    manifest_version: str = ""

    def validate(self) -> list[str]:
        blockers: list[str] = []
        if self.pipeline not in {"sesame", "minfi_noob"}:
            blockers.append("unsupported_pipeline")
        for suffix in ("_Grn.idat", "_Red.idat"):
            if not Path(f"{self.idat_prefix}{suffix}").is_file():
                blockers.append(f"missing:{suffix}")
        return blockers


def load_methylation_worker_result(output_dir: Path) -> dict[str, Any]:
    output_dir = Path(output_dir)
    qc_path = output_dir / "qc.json"
    measurements_path = output_dir / "measurements.csv"
    if not qc_path.is_file():
        raise ValueError("The preprocessing worker did not produce qc.json.")
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    if not qc.get("ok"):
        raise ValueError(str(qc.get("error") or "Methylation preprocessing failed."))
    frame = pd.read_csv(measurements_path)
    required = {"probe_id", "beta_value", "m_value", "detection_p", "qc_pass", "qc_reason", "normalization"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError("Preprocessing output is missing columns: " + ", ".join(sorted(missing)))
    checksum = hashlib.sha256(measurements_path.read_bytes()).hexdigest()
    return {"qc": qc, "measurements": frame, "checksum_sha256": checksum}


def variant_normalization_contract(
    *, input_vcf: Path, reference_fasta: Path, genome_build: str, output_vcf: Path
) -> dict[str, Any]:
    build = str(genome_build or "").upper()
    if build not in {"GRCH37", "HG19", "GRCH38", "HG38"}:
        raise ValueError("A native GRCh37/hg19 or GRCh38/hg38 build is required.")
    if not Path(input_vcf).is_file():
        raise FileNotFoundError(input_vcf)
    if not Path(reference_fasta).is_file():
        raise FileNotFoundError(reference_fasta)
    return {
        "tool": "bcftools norm",
        "steps": ["split_multiallelic", "left_align", "trim_alleles", "validate_reference"],
        "command": [
            "bcftools", "norm", "-f", str(reference_fasta), "-m", "-any", "-c", "e",
            "-Oz", "-o", str(output_vcf), str(input_vcf),
        ],
        "genome_build": build,
        "output_vcf": str(output_vcf),
        "unsupported_primary_classes": ["SV", "CNV", "VNTR"],
    }


def vep_annotation_contract(*, genome_build: str, external_api: bool, explicit_consent: bool) -> dict[str, Any]:
    if external_api and not explicit_consent:
        return {"eligible": False, "blockers": ["external_variant_transfer_not_approved"]}
    return {
        "eligible": True,
        "mode": "ensembl_rest" if external_api else "local_or_imported_vep",
        "genome_build": genome_build,
        "primary_transcript_policy": "MANE Select",
        "retained_consequences": "all transcripts",
        "requirements": ["normalized_alleles", "declared_build", "reference_validated"],
    }


def cross_build_mapping_gate(mapping: dict[str, Any]) -> dict[str, Any]:
    """Block adapters when liftover is ambiguous, invalid, or REF mismatched."""
    blockers: list[str] = []
    if mapping.get("status") != "mapped":
        blockers.append("mapping_not_valid")
    if mapping.get("ambiguous"):
        blockers.append("ambiguous_mapping")
    if not mapping.get("ref_validated"):
        blockers.append("target_reference_not_validated")
    if not mapping.get("chain_name") or not mapping.get("tool_version"):
        blockers.append("mapping_provenance_incomplete")
    return {
        "eligible": not blockers,
        "blockers": blockers,
        "original_coordinate": mapping.get("original_coordinate"),
        "mapped_coordinate": mapping.get("mapped_coordinate"),
        "chain_name": mapping.get("chain_name"),
        "tool": mapping.get("tool"),
        "tool_version": mapping.get("tool_version"),
    }
