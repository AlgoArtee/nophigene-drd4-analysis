"""Allowlisted biological-model manifests and eligibility contracts.

This module never treats user-submitted model metadata as an execution.  A
model is runnable only when its immutable manifest is registered, installed,
and all declared input/privacy/build requirements pass.
"""

from __future__ import annotations

import hashlib
import json
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable


@dataclass(frozen=True)
class ModelManifest:
    id: str
    name: str
    wave: str
    task: str
    execution_mode: str
    status: str
    required_inputs: tuple[str, ...]
    supported_builds: tuple[str, ...] = ()
    external_transfer: bool = False
    gpu: bool = False
    network_policy: str = "offline"
    output_semantics: str = "model-specific score"
    license: str = "Review upstream terms before installation."
    citation: str = ""
    limitations: tuple[str, ...] = ()
    notes: str = ""

    @property
    def checksum_sha256(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return {**asdict(self), "checksum_sha256": self.checksum_sha256}


MODEL_MANIFESTS: tuple[ModelManifest, ...] = (
    ModelManifest(
        id="alphagenome-api",
        name="AlphaGenome",
        wave="1",
        task="molecular variant-effect prediction",
        execution_mode="external_api",
        status="available_after_configuration",
        required_inputs=("grch38_variant", "reference_validated", "requested_modalities", "explicit_transfer_consent"),
        supported_builds=("GRCh38",),
        external_transfer=True,
        network_policy="provider_only",
        output_semantics="molecular track and variant-effect scores",
        citation="https://github.com/google-deepmind/alphagenome",
        limitations=("not_pathogenicity", "not_diagnosis", "scores_are_not_p_values"),
    ),
    ModelManifest(
        id="alphamissense-precomputed",
        name="AlphaMissense",
        wave="1",
        task="missense variant-effect lookup",
        execution_mode="precomputed",
        status="available_after_assets",
        required_inputs=("missense_variant", "exact_transcript", "protein_accession", "residue_ref_alt"),
        output_semantics="source-native missense pathogenicity score and publisher label",
        citation="https://deepmind.google/research/publications/21083/",
        limitations=("computational_evidence_only", "exact_protein_mapping_required"),
    ),
    ModelManifest(
        id="eve-precomputed",
        name="EVE",
        wave="1",
        task="evolutionary missense effect lookup",
        execution_mode="precomputed",
        status="available_after_assets",
        required_inputs=("missense_variant", "exact_transcript", "protein_accession", "residue_ref_alt"),
        output_semantics="gene-specific evolutionary score",
        citation="https://github.com/OATML-Markslab/EVE",
        limitations=("not_available_for_every_protein", "do_not_remap_isoforms_silently"),
    ),
    ModelManifest(
        id="popeve-precomputed",
        name="popEVE",
        wave="1",
        task="human missense-effect lookup",
        execution_mode="precomputed",
        status="available_after_assets",
        required_inputs=("missense_variant", "exact_transcript", "protein_accession", "residue_ref_alt"),
        output_semantics="proteome-wide human missense score",
        citation="https://github.com/debbiemarkslab/popEVE",
        limitations=("computational_evidence_only", "exact_protein_mapping_required"),
    ),
    ModelManifest(
        id="alphafold-db",
        name="AlphaFold DB",
        wave="1",
        task="predicted protein-structure context",
        execution_mode="external_api",
        status="available_after_configuration",
        required_inputs=("uniprot_accession", "explicit_transfer_consent"),
        external_transfer=True,
        network_policy="provider_only",
        output_semantics="structure metadata, residue confidence, and coordinates",
        citation="https://alphafold.ebi.ac.uk/api/",
        limitations=("structure_is_not_functional_proof", "structure_is_not_pathogenicity"),
    ),
    ModelManifest(
        id="borzoi-local",
        name="Borzoi",
        wave="2A",
        task="long-range regulatory sequence prediction",
        execution_mode="local_container",
        status="available_after_assets",
        required_inputs=("reference_sequence", "variant", "reference_validated", "model_window"),
        supported_builds=("GRCh38",),
        gpu=True,
        citation="https://github.com/calico/borzoi",
        limitations=("tissue_track_applicability_required",),
    ),
    ModelManifest(
        id="enformer-local",
        name="Enformer",
        wave="2A",
        task="long-range regulatory sequence prediction",
        execution_mode="local_container",
        status="available_after_assets",
        required_inputs=("reference_sequence", "variant", "reference_validated", "model_window"),
        supported_builds=("GRCh38",),
        gpu=True,
        citation="https://www.nature.com/articles/s41592-021-01252-x",
        limitations=("model_track_not_clinical_evidence",),
    ),
    ModelManifest(
        id="sei-local",
        name="Sei",
        wave="2A",
        task="chromatin-state sequence prediction",
        execution_mode="local_container",
        status="available_after_assets",
        required_inputs=("reference_sequence", "variant", "reference_validated", "model_window"),
        supported_builds=("GRCh38",),
        gpu=True,
        citation="https://www.nature.com/articles/s41588-022-01102-2",
        limitations=("sequence_class_is_not_observed_occupancy",),
    ),
    ModelManifest(
        id="deepsea-legacy-local",
        name="DeepSEA",
        wave="2A",
        task="legacy regulatory sequence prediction",
        execution_mode="local_container",
        status="legacy_optional",
        required_inputs=("reference_sequence", "variant", "reference_validated", "model_window"),
        supported_builds=("GRCh37",),
        gpu=True,
        citation="https://deepsea.princeton.edu/help/",
        limitations=("legacy_build_and_model", "never_silently_compare_across_builds"),
    ),
    ModelManifest(
        id="esm-variant-local",
        name="ESM variant effects",
        wave="2B",
        task="protein language-model missense scoring and embeddings",
        execution_mode="local_container",
        status="available_after_assets",
        required_inputs=("protein_sequence", "protein_accession", "residue_ref_alt"),
        gpu=True,
        output_semantics="zero-shot substitution score and optional embedding artifact",
        citation="https://github.com/facebookresearch/esm",
        limitations=("selected_checkpoint_must_be_recorded", "not_a_clinical_classification"),
    ),
    ModelManifest(
        id="methylbert-local",
        name="MethylBERT",
        wave="3",
        task="read-level methylation pattern modeling",
        execution_mode="local_container",
        status="input_gated_experimental",
        required_inputs=("bismark_bam_with_xm", "reference_fasta", "dmr_regions"),
        gpu=True,
        citation="https://github.com/CompEpigen/methylbert",
        limitations=("epic_idat_is_not_a_valid_input", "assay_and_disease_validation_required"),
    ),
    ModelManifest(
        id="methyl-gp-local",
        name="Methyl-GP",
        wave="3",
        task="sequence-site methylation-mark prediction",
        execution_mode="local_container",
        status="unsupported_for_epic_5mc",
        required_inputs=("declared_species", "declared_methylation_mark", "reference_sequence"),
        gpu=True,
        citation="https://pmc.ncbi.nlm.nih.gov/articles/PMC11952970/",
        limitations=("not_ordinary_human_epic_5mc_beta_prediction",),
    ),
    ModelManifest(
        id="methylprophet-local",
        name="MethylProphet",
        wave="3",
        task="methylome reconstruction",
        execution_mode="local_container",
        status="input_gated_experimental",
        required_inputs=("matched_gene_expression", "reference_sequence", "declared_tissue"),
        gpu=True,
        citation="https://pmc.ncbi.nlm.nih.gov/articles/PMC11839017/",
        limitations=("epic_idat_alone_is_not_a_valid_input",),
    ),
    ModelManifest(
        id="melody-local",
        name="Melody",
        wave="3",
        task="locus-specific methylation prediction from sequence",
        execution_mode="local_container",
        status="blocked_pending_verified_release",
        required_inputs=("ten_kb_reference_sequence", "declared_tissue"),
        gpu=True,
        citation="https://www.biorxiv.org/content/10.1101/2025.11.23.689975v2",
        limitations=("official_code_weights_license_and_benchmark_required", "preprint_model"),
    ),
    ModelManifest(
        id="disease-methylation-adapter",
        name="Disease-specific methylation transformer adapter",
        wave="3",
        task="model-defined",
        execution_mode="local_container",
        status="framework_only",
        required_inputs=("model_card", "disease", "training_cohort", "external_validation", "license", "input_contract"),
        gpu=True,
        limitations=("no_model_is_installable_until_every_contract_field_is_verified",),
    ),
)

MANIFEST_BY_ID = {manifest.id: manifest for manifest in MODEL_MANIFESTS}


def list_model_manifests() -> list[dict[str, Any]]:
    return [manifest.to_dict() for manifest in MODEL_MANIFESTS]


def inspect_model_inputs(model_id: str, available: dict[str, Any]) -> dict[str, Any]:
    manifest = MANIFEST_BY_ID.get(model_id)
    if manifest is None:
        return {"model_id": model_id, "status": "unknown_model", "eligible": False, "blockers": ["not_allowlisted"]}
    blockers = [name for name in manifest.required_inputs if not available.get(name)]
    raw_build = str(available.get("genome_build") or "")
    build = {"hg19": "GRCh37", "grch37": "GRCh37", "hg38": "GRCh38", "grch38": "GRCh38"}.get(
        raw_build.casefold(), raw_build
    )
    if manifest.supported_builds and build and build not in manifest.supported_builds:
        blockers.append(f"unsupported_build:{build}")
    if manifest.status in {"unsupported_for_epic_5mc", "blocked_pending_verified_release", "framework_only"}:
        blockers.append(manifest.status)
    if manifest.external_transfer and not available.get("explicit_transfer_consent"):
        blockers.append("external_transfer_not_approved")
    return {
        "model_id": model_id,
        "name": manifest.name,
        "status": "eligible" if not blockers else "model_not_eligible",
        "eligible": not blockers,
        "blockers": sorted(set(blockers)),
        "input_contract": list(manifest.required_inputs),
        "supported_builds": list(manifest.supported_builds),
        "execution_mode": manifest.execution_mode,
        "network_policy": manifest.network_policy,
        "output_semantics": manifest.output_semantics,
        "limitations": list(manifest.limitations),
        "manifest_checksum_sha256": manifest.checksum_sha256,
    }


def estimate_model_resources(model_id: str, inputs: dict[str, Any]) -> dict[str, Any]:
    manifest = MANIFEST_BY_ID.get(model_id)
    if manifest is None:
        return {"model_id": model_id, "status": "unknown_model", "eligible": False}
    sequence_bases = int(inputs.get("sequence_bases") or 0)
    records = int(inputs.get("record_count") or 1)
    return {
        "model_id": model_id,
        "status": "estimate_only",
        "execution_mode": manifest.execution_mode,
        "gpu_required": manifest.gpu,
        "network_policy": manifest.network_policy,
        "input_record_count": max(1, records),
        "sequence_bases": max(0, sequence_bases),
        "vram": "model/version-specific; verify during signed asset installation" if manifest.gpu else "not required",
        "disk": "asset manifest required before installation",
        "preflight": "WSL2/NVIDIA validation required" if manifest.gpu else "CPU/API preflight required",
        "limitations": list(manifest.limitations),
    }


def model_installation_plan(model_id: str) -> dict[str, Any]:
    manifest = MANIFEST_BY_ID.get(model_id)
    if manifest is None:
        return {"model_id": model_id, "status": "unknown_model"}
    return {
        "model_id": model_id,
        "status": "explicit_confirmation_required",
        "installable": manifest.status not in {
            "unsupported_for_epic_5mc", "blocked_pending_verified_release", "framework_only"
        },
        "license": manifest.license,
        "citation": manifest.citation,
        "manifest_checksum_sha256": manifest.checksum_sha256,
        "required_confirmation_fields": ["manifest_checksum_sha256", "license_acknowledged", "asset_checksums_acknowledged"],
        "runner_policy": "signed allowlisted manifest via privilege-separated local runner; never the web container",
    }


class ModelAdapter(ABC):
    """Contract implemented by every allowlisted external or local adapter."""

    manifest_id: str

    def inspect_inputs(self, available: dict[str, Any]) -> dict[str, Any]:
        return inspect_model_inputs(self.manifest_id, available)

    @abstractmethod
    def estimate_resources(self, inputs: dict[str, Any]) -> dict[str, Any]: ...

    @abstractmethod
    def prepare(self, inputs: dict[str, Any]) -> dict[str, Any]: ...

    @abstractmethod
    def execute(self, prepared: dict[str, Any]) -> Any: ...

    @abstractmethod
    def normalize(self, raw_output: Any) -> list[dict[str, Any]]: ...

    @abstractmethod
    def validate(self, normalized: Iterable[dict[str, Any]]) -> dict[str, Any]: ...
