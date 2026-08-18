"""Execution engine for asynchronous API workflow jobs."""

from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from .errors import APIError
from .profiles import ProfileStore
from .serialization import read_json, to_jsonable, utc_now, write_json_atomic

try:
    from ..analysis import (
        ANALYSIS_SCOPE_OPTIONS,
        AnalysisError,
        PROJECT_ROOT,
        analyze_prepared_data,
        build_gene_manifest_subset,
        build_gene_methylation_table,
        fetch_population_stats,
        generate_report,
        load_full_methylation_manifest,
        load_methylation_beta_values,
        load_variants,
        normalize_analysis_scope,
    )
    from ..bam_extraction import HG38_FASTA, extract_region_vcf
    from ..variant_knowledge.orchestrator import build_dynamic_knowledge_base
    from ..workbench.evidence import CORE_SOURCE_KEYS
    from ..workbench.reporting import build_canonical_report
    from ..workflow import (
        normalize_gene_symbol,
        normalize_genome_build,
        resolve_gene_region,
        select_profile_variant_source,
    )
except ImportError:
    from analysis import (
        ANALYSIS_SCOPE_OPTIONS,
        AnalysisError,
        PROJECT_ROOT,
        analyze_prepared_data,
        build_gene_manifest_subset,
        build_gene_methylation_table,
        fetch_population_stats,
        generate_report,
        load_full_methylation_manifest,
        load_methylation_beta_values,
        load_variants,
        normalize_analysis_scope,
    )
    from bam_extraction import HG38_FASTA, extract_region_vcf
    from variant_knowledge.orchestrator import build_dynamic_knowledge_base
    from workbench.evidence import CORE_SOURCE_KEYS
    from workbench.reporting import build_canonical_report
    from workflow import (
        normalize_gene_symbol,
        normalize_genome_build,
        resolve_gene_region,
        select_profile_variant_source,
    )

JOB_OPERATIONS = {
    "resolve_regions",
    "prepare_manifests",
    "extract_variants",
    "build_knowledge_bases",
    "analyze",
    "render_reports",
    "full_workflow",
}
MAX_GENES_PER_JOB = 100
RESULT_SCHEMA_VERSION = "2.0"
REPORT_SCHEMA_VERSION = "3.0"


def normalize_job_request(payload: Any) -> dict[str, Any]:
    """Validate and normalize a submitted job request."""
    if not isinstance(payload, dict):
        raise APIError("invalid_json", "Request body must be a JSON object.", 400)
    operation = str(payload.get("operation") or "").strip()
    if operation not in JOB_OPERATIONS:
        raise APIError(
            "invalid_operation",
            f"'operation' must be one of: {', '.join(sorted(JOB_OPERATIONS))}.",
            422,
        )
    raw_genes = payload.get("genes")
    if isinstance(raw_genes, str):
        raw_genes = [raw_genes]
    if not isinstance(raw_genes, list):
        raise APIError("invalid_genes", "'genes' must be a gene name or list of gene names.", 422)
    genes: list[str] = []
    seen: set[str] = set()
    for raw_gene in raw_genes:
        try:
            gene = normalize_gene_symbol(raw_gene)
        except ValueError as exc:
            raise APIError("invalid_gene", str(exc), 422) from exc
        if gene not in seen:
            genes.append(gene)
            seen.add(gene)
    if not genes:
        raise APIError("invalid_genes", "Submit at least one gene.", 422)
    if len(genes) > MAX_GENES_PER_JOB:
        raise APIError(
            "too_many_genes",
            f"A job may contain at most {MAX_GENES_PER_JOB} unique genes.",
            422,
        )

    analysis_scope = str(payload.get("analysis_scope") or "promoter_plus_gene").strip().lower()
    analysis_scope = analysis_scope.replace("-", "_").replace(" ", "_")
    if analysis_scope not in ANALYSIS_SCOPE_OPTIONS:
        raise APIError(
            "invalid_analysis_scope",
            f"'analysis_scope' must be one of: {', '.join(sorted(ANALYSIS_SCOPE_OPTIONS))}.",
            422,
        )
    try:
        genome_build = normalize_genome_build(
            payload.get("genome_build") or "auto",
            allow_auto=True,
        )
    except ValueError as exc:
        raise APIError("invalid_genome_build", str(exc), 422) from exc

    profile_id = str(payload.get("profile_id") or "").strip()
    source_job_id = str(payload.get("source_job_id") or "").strip()
    if operation in {"prepare_manifests", "extract_variants", "build_knowledge_bases", "analyze", "full_workflow"} and not profile_id:
        raise APIError("profile_required", f"'profile_id' is required for {operation}.", 422)
    if operation == "render_reports" and not source_job_id:
        raise APIError("source_job_required", "'source_job_id' is required for render_reports.", 422)

    raw_overrides = payload.get("region_overrides") or {}
    if not isinstance(raw_overrides, dict):
        raise APIError("invalid_region_overrides", "'region_overrides' must be an object.", 422)
    region_overrides: dict[str, str] = {}
    region_pattern = re.compile(r"^(?:chr)?[^:\s]+:[\d,]+-[\d,]+$")
    for raw_gene, raw_region in raw_overrides.items():
        try:
            gene = normalize_gene_symbol(raw_gene)
        except ValueError as exc:
            raise APIError("invalid_region_override", str(exc), 422) from exc
        region = str(raw_region or "").strip()
        if not region_pattern.fullmatch(region):
            raise APIError(
                "invalid_region_override",
                f"Invalid region override for {gene}: '{region}'.",
                422,
            )
        region_overrides[gene] = region

    options = payload.get("options") or {}
    if not isinstance(options, dict):
        raise APIError("invalid_options", "'options' must be an object.", 422)
    normalized_options = {
        "update_general_database": bool(options.get("update_general_database", False)),
        "overwrite_general_database": bool(options.get("overwrite_general_database", False)),
        "use_dynamic_knowledge_base": bool(options.get("use_dynamic_knowledge_base", False)),
        "use_local_article_evidence": bool(options.get("use_local_article_evidence", False)),
        "article_pdf_folder": str(options.get("article_pdf_folder") or "").strip(),
        "article_pdf_recursive": bool(options.get("article_pdf_recursive", True)),
    }
    interpretation_mode = str(options.get("interpretation_mode") or "research").strip().lower()
    interpretation_mode = interpretation_mode.replace("-", "_").replace(" ", "_")
    if interpretation_mode not in {"research", "clinical_support", "dual"}:
        raise APIError(
            "invalid_interpretation_mode",
            "'options.interpretation_mode' must be research, clinical_support, or dual.",
            422,
        )
    normalized_options["interpretation_mode"] = interpretation_mode
    normalized_options["evidence_snapshot_id"] = str(options.get("evidence_snapshot_id") or "").strip()
    requested_models = options.get("requested_models") or []
    if not isinstance(requested_models, list) or not all(isinstance(model, dict) for model in requested_models):
        raise APIError("invalid_requested_models", "'options.requested_models' must be a list of model objects.", 422)
    normalized_options["requested_models"] = [dict(model) for model in requested_models]
    sample_context = options.get("sample_context") or {}
    if not isinstance(sample_context, dict):
        raise APIError("invalid_sample_context", "'options.sample_context' must be an object.", 422)
    normalized_options["sample_context"] = dict(sample_context)
    external_consents = options.get("external_consents") or {}
    if not isinstance(external_consents, dict):
        raise APIError("invalid_external_consents", "'options.external_consents' must map source keys to booleans.", 422)
    normalized_options["external_consents"] = {
        str(source_key).strip(): bool(approved)
        for source_key, approved in external_consents.items()
        if str(source_key).strip()
    }
    try:
        normalized_options["max_article_pdfs"] = min(
            1000,
            max(0, int(options.get("max_article_pdfs", 100) or 0)),
        )
    except (TypeError, ValueError) as exc:
        raise APIError("invalid_max_article_pdfs", "'options.max_article_pdfs' must be an integer.", 422) from exc
    raw_sources = options.get("knowledge_sources")
    if raw_sources in (None, "", []):
        normalized_options["knowledge_sources"] = []
    elif isinstance(raw_sources, str):
        normalized_options["knowledge_sources"] = [
            item.strip() for item in raw_sources.split(",") if item.strip()
        ]
    elif isinstance(raw_sources, list):
        normalized_options["knowledge_sources"] = [str(item).strip() for item in raw_sources if str(item).strip()]
    else:
        raise APIError("invalid_knowledge_sources", "'options.knowledge_sources' must be a list or comma-separated string.", 422)
    unsupported_sources = sorted(set(normalized_options["knowledge_sources"]) - set(CORE_SOURCE_KEYS))
    if unsupported_sources:
        raise APIError(
            "unsupported_knowledge_sources",
            "Version 2 accepts only the vetted core source set: " + ", ".join(unsupported_sources),
            422,
        )
    try:
        normalized_options["external_timeout_seconds"] = min(
            1200,
            max(1, int(options.get("external_timeout_seconds", 1200))),
        )
    except (TypeError, ValueError) as exc:
        raise APIError(
            "invalid_external_timeout",
            "'options.external_timeout_seconds' must be an integer between 1 and 1200.",
            422,
        ) from exc
    raw_workflows = options.get("knowledge_workflows")
    if raw_workflows in (None, "", []):
        normalized_options["knowledge_workflows"] = []
    elif isinstance(raw_workflows, str):
        normalized_options["knowledge_workflows"] = [
            item.strip() for item in raw_workflows.split(",") if item.strip()
        ]
    elif isinstance(raw_workflows, list):
        normalized_options["knowledge_workflows"] = [
            str(item).strip() for item in raw_workflows if str(item).strip()
        ]
    else:
        raise APIError(
            "invalid_knowledge_workflows",
            "'options.knowledge_workflows' must be a list or comma-separated string.",
            422,
        )
    raw_imports = options.get("knowledge_source_imports") or {}
    if not isinstance(raw_imports, dict):
        raise APIError(
            "invalid_knowledge_source_imports",
            "'options.knowledge_source_imports' must be an object mapping source keys to CSV/JSON paths.",
            422,
        )
    normalized_options["knowledge_source_imports"] = {
        str(source_key).strip(): str(import_path).strip()
        for source_key, import_path in raw_imports.items()
        if str(source_key).strip() and str(import_path).strip()
    }
    return {
        "operation": operation,
        "genes": genes,
        "profile_id": profile_id,
        "source_job_id": source_job_id,
        "analysis_scope": analysis_scope,
        "genome_build": genome_build,
        "region_overrides": region_overrides,
        "options": normalized_options,
    }


def _artifact_url(job_id: str, relative_path: str) -> str:
    return f"/api/v1/jobs/{job_id}/artifacts/{relative_path.replace(chr(92), '/')}"


def _read_source_result(jobs_root: Path, source_job_id: str) -> dict[str, Any]:
    source_dir = jobs_root / source_job_id
    result = read_json(source_dir / "result.json")
    if not isinstance(result, dict):
        raise AnalysisError(f"Source job '{source_job_id}' has no completed result manifest.")
    return result


def _source_outcome(result: dict[str, Any], gene: str) -> dict[str, Any]:
    for outcome in result.get("genes", []):
        if outcome.get("gene") == gene:
            if outcome.get("status") != "succeeded":
                raise AnalysisError(f"Source job did not complete {gene} successfully.")
            return outcome
    raise AnalysisError(f"Source job does not contain a result for {gene}.")


def _rehydrate_methylation_insights(payload: dict[str, Any]) -> dict[str, Any]:
    insights = dict(payload or {})
    if isinstance(insights.get("probe_preview"), list):
        insights["probe_preview"] = pd.DataFrame(insights["probe_preview"])
    return insights


class WorkflowRunner:
    """Run API operations against profiles and persisted source jobs."""

    def __init__(self, profile_store: ProfileStore, jobs_root: Path):
        self.profile_store = profile_store
        self.jobs_root = Path(jobs_root)

    def execute(
        self,
        job_id: str,
        request_payload: dict[str, Any],
        job_dir: Path,
        progress: Callable[[str, int, int, list[dict[str, Any]]], None],
    ) -> dict[str, Any]:
        operation = request_payload["operation"]
        genes = request_payload["genes"]
        profile = (
            self.profile_store.get(request_payload["profile_id"])
            if request_payload.get("profile_id")
            else None
        )
        source_result = (
            _read_source_result(self.jobs_root, request_payload["source_job_id"])
            if request_payload.get("source_job_id")
            else None
        )

        manifest: pd.DataFrame | None = None
        beta_values: pd.DataFrame | None = None
        if operation in {"prepare_manifests", "build_knowledge_bases", "analyze", "full_workflow"}:
            progress("loading_manifest", 0, len(genes), [])
            manifest = load_full_methylation_manifest(profile["manifest_path"])
        if operation in {"analyze", "full_workflow"}:
            progress("processing_idat", 0, len(genes), [])
            beta_values = load_methylation_beta_values(
                profile["idat_prefix"],
                manifest_filepath=profile["manifest_path"],
            )

        outcomes: list[dict[str, Any]] = []
        for index, gene in enumerate(genes):
            progress(f"{operation}:{gene}", index, len(genes), outcomes)
            gene_dir = job_dir / "genes" / gene
            gene_dir.mkdir(parents=True, exist_ok=True)
            try:
                outcome = self._execute_gene(
                    job_id=job_id,
                    operation=operation,
                    gene=gene,
                    request_payload=request_payload,
                    profile=profile,
                    source_result=source_result,
                    manifest=manifest,
                    beta_values=beta_values,
                    gene_dir=gene_dir,
                )
            except Exception as exc:
                outcome = {
                    "gene": gene,
                    "status": "failed",
                    "stage": operation,
                    "error": {
                        "code": self._error_code(exc),
                        "message": str(exc),
                    },
                    "warnings": [],
                    "artifacts": {},
                }
            outcomes.append(outcome)
            progress(f"{operation}:{gene}:complete", index + 1, len(genes), outcomes)

        succeeded = sum(outcome["status"] == "succeeded" for outcome in outcomes)
        failed = len(outcomes) - succeeded
        status = "succeeded" if failed == 0 else "failed" if succeeded == 0 else "partial"
        result = {
            "schema_version": RESULT_SCHEMA_VERSION,
            "job_id": job_id,
            "operation": operation,
            "status": status,
            "created_at": utc_now(),
            "counts": {
                "requested": len(genes),
                "succeeded": succeeded,
                "failed": failed,
            },
            "genes": outcomes,
            "artifacts": {
                "result": _artifact_url(job_id, "result.json"),
                "zip": _artifact_url(job_id, "artifacts.zip"),
            },
        }
        write_json_atomic(job_dir / "result.json", result)
        self._write_zip(job_dir)
        return result

    def _execute_gene(
        self,
        *,
        job_id: str,
        operation: str,
        gene: str,
        request_payload: dict[str, Any],
        profile: dict[str, Any] | None,
        source_result: dict[str, Any] | None,
        manifest: pd.DataFrame | None,
        beta_values: pd.DataFrame | None,
        gene_dir: Path,
    ) -> dict[str, Any]:
        if operation == "render_reports":
            source = _source_outcome(source_result or {}, gene)
            return self._render_from_source(
                job_id=job_id,
                gene=gene,
                gene_dir=gene_dir,
                source_job_id=request_payload["source_job_id"],
                source_outcome=source,
            )

        source = _source_outcome(source_result, gene) if source_result else None
        if source:
            resolved = {
                "gene": gene,
                "genome_build": source["genome_build"],
                "region": source["region"],
                "scope": source.get("analysis_scope", request_payload["analysis_scope"]),
                "scope_regions": source.get("scope_regions", {}),
                "selected_gene_region": source.get("selected_gene_region", source["region"]),
                "selected_sources": source.get("selected_sources", ["Source job"]),
                "candidate_regions": source.get("candidate_regions", []),
                "curated_coordinates": source.get("curated_coordinates", False),
            }
        else:
            resolved = resolve_gene_region(
                gene,
                genome_build=request_payload["genome_build"],
                default_genome_build=(profile or {}).get("default_genome_build", "hg19"),
                analysis_scope=request_payload["analysis_scope"],
                region_override=request_payload["region_overrides"].get(gene),
            )

        outcome: dict[str, Any] = {
            "gene": gene,
            "status": "succeeded",
            "stage": operation,
            "genome_build": resolved["genome_build"],
            "region": resolved["region"],
            "analysis_scope": resolved["scope"],
            "scope_regions": resolved["scope_regions"],
            "selected_gene_region": resolved["selected_gene_region"],
            "selected_sources": resolved["selected_sources"],
            "candidate_regions": resolved["candidate_regions"],
            "curated_coordinates": resolved["curated_coordinates"],
            "warnings": [],
            "artifacts": {},
        }
        region_path = gene_dir / "region.json"
        write_json_atomic(region_path, resolved)
        outcome["artifacts"]["region"] = _artifact_url(job_id, f"genes/{gene}/region.json")
        if operation == "resolve_regions":
            return outcome

        manifest_subset: pd.DataFrame | None = None
        if operation in {"prepare_manifests", "build_knowledge_bases", "analyze", "full_workflow"}:
            manifest_subset = build_gene_manifest_subset(
                manifest,
                region=resolved["region"],
                genome_build=resolved["genome_build"],
            )
            manifest_path = gene_dir / "manifest.csv"
            manifest_subset.to_csv(manifest_path, index=False)
            outcome["probe_count"] = len(manifest_subset)
            outcome["artifacts"]["manifest"] = _artifact_url(job_id, f"genes/{gene}/manifest.csv")
            if manifest_subset.empty:
                outcome["warnings"].append("No methylation manifest probes were found in the active region.")
        if operation == "prepare_manifests":
            return outcome

        variant_source = self._variant_source(
            profile=profile,
            resolved=resolved,
            gene=gene,
            gene_dir=gene_dir,
            source_outcome=source,
        )
        outcome["variant_source"] = variant_source
        if variant_source.get("extracted"):
            source_path = Path(variant_source["path"]).resolve()
            try:
                source_path.relative_to(gene_dir.resolve())
                artifact = _artifact_url(job_id, f"genes/{gene}/variants.vcf.gz")
            except ValueError:
                artifact = str(variant_source.get("artifact") or "")
            if artifact:
                variant_source["artifact"] = artifact
                outcome["artifacts"]["variant_vcf"] = artifact
        if operation == "extract_variants":
            return outcome

        variants = load_variants(variant_source["path"], resolved["region"])
        dynamic_knowledge_base_path: Path | None = None
        requested_snapshot_id = request_payload["options"].get("evidence_snapshot_id", "")
        source_snapshot_path = Path(str((source or {}).get("dynamic_knowledge_base_path") or ""))
        if source_snapshot_path.is_file():
            source_snapshot_payload = read_json(source_snapshot_path, default={})
            source_snapshot_id = str(
                (source_snapshot_payload or {}).get("evidence_snapshot", {}).get("snapshot_id", "")
            )
            if requested_snapshot_id and requested_snapshot_id != source_snapshot_id:
                raise AnalysisError(
                    f"Requested evidence snapshot '{requested_snapshot_id}' does not match the source job snapshot."
                )
            dynamic_knowledge_base_path = source_snapshot_path
            outcome["dynamic_knowledge_base_path"] = str(dynamic_knowledge_base_path)
            outcome["evidence_snapshot_id"] = source_snapshot_id
            outcome["evidence_snapshot_reused"] = True
            if operation == "build_knowledge_bases":
                return outcome
        elif requested_snapshot_id:
            raise AnalysisError(
                "evidence_snapshot_id requires a source job that contains the requested dynamic knowledge-base artifact."
            )
        elif operation == "build_knowledge_bases" or request_payload["options"]["use_dynamic_knowledge_base"]:
            dynamic_payload = self._build_dynamic_knowledge_base_artifact(
                gene=gene,
                resolved=resolved,
                variants=variants,
                manifest_subset=manifest_subset,
                gene_dir=gene_dir,
                request_payload=request_payload,
            )
            dynamic_knowledge_base_path = Path(dynamic_payload["artifact_path"])
            outcome["dynamic_knowledge_base_path"] = str(dynamic_knowledge_base_path)
            outcome["artifacts"]["dynamic_knowledge_base"] = _artifact_url(
                job_id,
                f"genes/{gene}/dynamic_knowledge_base/variant_kb.json",
            )
            outcome["dynamic_provider_count"] = len(dynamic_payload.get("provider_statuses", []))
            outcome["dynamic_workflow_count"] = len(dynamic_payload.get("workflow_runs", []))
            outcome["evidence_snapshot_id"] = str(
                dynamic_payload.get("evidence_snapshot", {}).get("snapshot_id", "")
            )
            if operation == "build_knowledge_bases":
                return outcome
        methylation = build_gene_methylation_table(
            beta_values,
            manifest_subset,
            region=resolved["region"],
            genome_build=resolved["genome_build"],
        )
        popstats = None
        if profile.get("population_statistics_path"):
            popstats = fetch_population_stats(profile["population_statistics_path"], variants)
        prepared = analyze_prepared_data(
            variants=variants,
            methylation=methylation,
            gene_name=gene,
            region=resolved["region"],
            analysis_scope=resolved["scope"],
            popstats=popstats,
            update_general_database_enabled=request_payload["options"]["update_general_database"],
            overwrite_general_database=request_payload["options"]["overwrite_general_database"],
            dynamic_knowledge_base_path=dynamic_knowledge_base_path,
            genome_build=resolved["genome_build"],
            interpretation_mode=request_payload["options"]["interpretation_mode"],
            sample_context={
                **dict(profile.get("sample_context", {})),
                **dict(request_payload["options"].get("sample_context", {})),
            },
            requested_models=request_payload["options"].get("requested_models", []),
        )
        variants_path = gene_dir / "variants.csv"
        methylation_path = gene_dir / "methylation.csv"
        prepared.variants.to_csv(variants_path, index=False)
        prepared.methylation.to_csv(methylation_path, index=False)
        analysis_payload = self._analysis_payload(
            gene=gene,
            resolved=resolved,
            variant_source=variant_source,
            prepared=prepared,
            warnings=outcome["warnings"],
        )
        analysis_path = gene_dir / "analysis.json"
        write_json_atomic(analysis_path, analysis_payload)
        outcome["variant_count"] = len(prepared.variants)
        outcome["methylation_count"] = len(prepared.methylation)
        outcome["general_database_status"] = prepared.general_database_status
        outcome["artifacts"].update(
            {
                "analysis": _artifact_url(job_id, f"genes/{gene}/analysis.json"),
                "variants": _artifact_url(job_id, f"genes/{gene}/variants.csv"),
                "methylation": _artifact_url(job_id, f"genes/{gene}/methylation.csv"),
            }
        )
        if operation == "full_workflow":
            outcome["artifacts"].update(
                self._render_reports(
                    job_id=job_id,
                    gene=gene,
                    gene_dir=gene_dir,
                    resolved=resolved,
                    variant_source=variant_source,
                    variants=prepared.variants,
                    methylation=prepared.methylation,
                    analysis_payload=analysis_payload,
                    popstats=prepared.popstats,
                    warnings=outcome["warnings"],
                )
            )
        return outcome

    def _build_dynamic_knowledge_base_artifact(
        self,
        *,
        gene: str,
        resolved: dict[str, Any],
        variants: pd.DataFrame,
        manifest_subset: pd.DataFrame | None,
        gene_dir: Path,
        request_payload: dict[str, Any],
    ) -> dict[str, Any]:
        output_dir = gene_dir / "dynamic_knowledge_base"
        selected_sources = request_payload["options"].get("knowledge_sources") or []
        consents = request_payload["options"].get("external_consents") or {}
        imports = request_payload["options"].get("knowledge_source_imports") or {}
        unapproved = [
            source_key
            for source_key in selected_sources
            if not consents.get(source_key) and source_key not in imports
        ]
        if unapproved:
            raise AnalysisError(
                "External evidence transfer was not approved for: " + ", ".join(sorted(unapproved))
            )
        payload = build_dynamic_knowledge_base(
            gene=gene,
            region=resolved["region"],
            genome_build=resolved["genome_build"],
            variants=variants,
            manifest_subset=manifest_subset,
            selected_sources=selected_sources,
            selected_workflows=request_payload["options"].get("knowledge_workflows") or None,
            source_imports=request_payload["options"].get("knowledge_source_imports") or None,
            use_local_article_evidence=bool(request_payload["options"].get("use_local_article_evidence", False)),
            article_pdf_folder=request_payload["options"].get("article_pdf_folder") or None,
            article_pdf_recursive=bool(request_payload["options"].get("article_pdf_recursive", True)),
            max_article_pdfs=int(request_payload["options"].get("max_article_pdfs", 100) or 0),
            output_dir=output_dir,
            cache_dir=gene_dir / ".knowledge_cache",
        )
        return payload

    def _variant_source(
        self,
        *,
        profile: dict[str, Any],
        resolved: dict[str, Any],
        gene: str,
        gene_dir: Path,
        source_outcome: dict[str, Any] | None,
    ) -> dict[str, Any]:
        if source_outcome and source_outcome.get("variant_source"):
            source = dict(source_outcome["variant_source"])
            if Path(source.get("path", "")).exists():
                source["reused_from_source_job"] = True
                return source
        source = select_profile_variant_source(profile, resolved["genome_build"])
        if source["type"] == "vcf":
            return {
                "type": "vcf",
                "path": source["path"],
                "genome_build": resolved["genome_build"],
                "extracted": False,
            }
        reference = source.get("reference_fasta")
        if not reference and resolved["genome_build"] == "hg38":
            reference = str(HG38_FASTA)
        if not reference:
            raise AnalysisError(
                f"BAM source for {resolved['genome_build']} requires a reference_fasta in the profile."
            )
        output_vcf = gene_dir / "variants.vcf.gz"
        extraction = extract_region_vcf(
            bam_path=source["path"],
            region=resolved["region"],
            output_vcf=output_vcf,
            reference_fasta=reference,
        )
        return {
            "type": "bam_extraction",
            "path": str(extraction["output_vcf"]),
            "source_bam": source["path"],
            "reference_fasta": reference,
            "genome_build": resolved["genome_build"],
            "resolved_region": extraction["resolved_region"],
            "extracted": True,
        }

    def _analysis_payload(
        self,
        *,
        gene,
        resolved,
        variant_source,
        prepared,
        warnings: list[str],
    ) -> dict[str, Any]:
        dynamic_path = getattr(prepared, "dynamic_knowledge_base_path", None)
        dynamic_payload = read_json(Path(dynamic_path), default={}) if dynamic_path else {}
        if not isinstance(dynamic_payload, dict):
            dynamic_payload = {}
        return {
            "schema_version": REPORT_SCHEMA_VERSION,
            "gene": gene,
            "genome_build": resolved["genome_build"],
            "region": resolved["region"],
            "analysis_scope": prepared.analysis_scope,
            "analysis_scope_label": prepared.analysis_scope_label,
            "source_provenance": {
                "region_sources": resolved["selected_sources"],
                "region_candidates": resolved["candidate_regions"],
                "variant_source": variant_source,
            },
            "counts": {
                "variants": len(prepared.variants),
                "methylation_probes": len(prepared.methylation),
            },
            "warnings": list(warnings),
            "variants": prepared.variants,
            "methylation": prepared.methylation,
            "variant_interpretations": prepared.variant_interpretations,
            "methylation_insights": prepared.methylation_insights,
            "knowledge_base": prepared.knowledge_base,
            "population_statistics": prepared.popstats,
            "population_insights": prepared.population_insights,
            "population_database": prepared.population_database,
            "interpretation": getattr(prepared, "interpretation", {}),
            "general_database": {
                "path": prepared.general_database_path,
                "status": prepared.general_database_status,
            },
            "dynamic_knowledge_base": {
                **dynamic_payload,
                "path": dynamic_path,
                "status": getattr(prepared, "dynamic_knowledge_base_status", ""),
            },
        }

    def _render_from_source(
        self,
        *,
        job_id: str,
        gene: str,
        gene_dir: Path,
        source_job_id: str,
        source_outcome: dict[str, Any],
    ) -> dict[str, Any]:
        source_dir = self.jobs_root / source_job_id / "genes" / gene
        analysis_payload = read_json(source_dir / "analysis.json")
        if not isinstance(analysis_payload, dict):
            raise AnalysisError(f"Source job has no reusable analysis payload for {gene}.")
        variants = pd.read_csv(source_dir / "variants.csv")
        methylation = pd.read_csv(source_dir / "methylation.csv")
        resolved = {
            "genome_build": analysis_payload["genome_build"],
            "region": analysis_payload["region"],
            "scope": analysis_payload["analysis_scope"],
        }
        artifacts = self._render_reports(
            job_id=job_id,
            gene=gene,
            gene_dir=gene_dir,
            resolved=resolved,
            variant_source=analysis_payload.get("source_provenance", {}).get("variant_source", {}),
            variants=variants,
            methylation=methylation,
            analysis_payload=analysis_payload,
            popstats=analysis_payload.get("population_statistics"),
            warnings=list(source_outcome.get("warnings", [])),
        )
        return {
            "gene": gene,
            "status": "succeeded",
            "stage": "render_reports",
            "genome_build": resolved["genome_build"],
            "region": resolved["region"],
            "analysis_scope": resolved["scope"],
            "variant_count": len(variants),
            "methylation_count": len(methylation),
            "source_job_id": source_job_id,
            "warnings": [],
            "artifacts": artifacts,
        }

    def _render_reports(
        self,
        *,
        job_id: str,
        gene: str,
        gene_dir: Path,
        resolved: dict[str, Any],
        variant_source: dict[str, Any],
        variants: pd.DataFrame,
        methylation: pd.DataFrame,
        analysis_payload: dict[str, Any],
        popstats: Any,
        warnings: list[str],
    ) -> dict[str, str]:
        html_path = gene_dir / "report.html"
        report_json_path = gene_dir / "report.json"
        summary_path = gene_dir / "report_summary.csv"
        methylation_insights = _rehydrate_methylation_insights(
            analysis_payload.get("methylation_insights", {})
        )
        generate_report(
            variants,
            methylation,
            popstats,
            str(html_path),
            gene_name=gene,
            region=resolved["region"],
            methylation_output_path=gene_dir / "methylation.csv",
            variant_interpretations=analysis_payload.get("variant_interpretations", {}),
            methylation_insights=methylation_insights,
            population_insights=analysis_payload.get("population_insights", {}),
            analysis_scope=resolved["scope"],
            interpretation=analysis_payload.get("interpretation", {}),
        )
        artifacts = {
            "report_html": _artifact_url(job_id, f"genes/{gene}/report.html"),
            "report_json": _artifact_url(job_id, f"genes/{gene}/report.json"),
            "report_summary": _artifact_url(job_id, f"genes/{gene}/report_summary.csv"),
        }
        companion_files = {
            "analysis": "analysis.json",
            "variants": "variants.csv",
            "methylation": "methylation.csv",
            "manifest": "manifest.csv",
            "region": "region.json",
            "variant_vcf": "variants.vcf.gz",
        }
        for key, filename in companion_files.items():
            if (gene_dir / filename).is_file():
                artifacts[key] = _artifact_url(job_id, f"genes/{gene}/{filename}")
        canonical = build_canonical_report(
            {
                **analysis_payload,
                "job_id": job_id,
                "variants": variants,
                "methylation": methylation,
                "artifacts": artifacts,
                "warnings": list(warnings),
                "source_provenance": {
                    **dict(analysis_payload.get("source_provenance", {})),
                    "variant_source": variant_source,
                },
            }
        )
        write_json_atomic(report_json_path, canonical)
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
            writer.writeheader()
            writer.writerows(
                [
                    {"metric": "gene", "value": gene},
                    {"metric": "genome_build", "value": resolved["genome_build"]},
                    {"metric": "region", "value": resolved["region"]},
                    {"metric": "analysis_scope", "value": resolved["scope"]},
                    {"metric": "variant_count", "value": len(variants)},
                    {"metric": "methylation_probe_count", "value": len(methylation)},
                    {"metric": "has_population_statistics", "value": popstats is not None},
                    {
                        "metric": "interpretation_mode",
                        "value": analysis_payload.get("interpretation", {})
                        .get("interpretation_context", {})
                        .get("mode", "research"),
                    },
                    {
                        "metric": "evidence_snapshot_id",
                        "value": analysis_payload.get("interpretation", {})
                        .get("evidence_snapshot", {})
                        .get("snapshot_id", ""),
                    },
                ]
            )
        return artifacts

    @staticmethod
    def _error_code(exc: Exception) -> str:
        if isinstance(exc, APIError):
            return exc.code
        if isinstance(exc, (AnalysisError, ValueError)):
            return "gene_workflow_failed"
        return "internal_error"

    @staticmethod
    def _write_zip(job_dir: Path) -> None:
        zip_path = job_dir / "artifacts.zip"
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(job_dir.rglob("*")):
                if not path.is_file() or path == zip_path or path.name in {"job.json", "request.json"}:
                    continue
                archive.write(path, path.relative_to(job_dir).as_posix())
