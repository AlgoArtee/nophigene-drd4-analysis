"""Flask-based web UI for the gene-focused analysis workflow."""

from __future__ import annotations

import io
import os
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

import pandas as pd
from flask import Flask, jsonify, render_template, request, session

try:
    from .analysis import (
        AnalysisError,
        DEFAULT_GENE_NAME,
        DEFAULT_REGION,
        DEFAULT_REPORT_NAME,
        run_analysis,
    )
    from .gene_region_extraction import find_gene_region
    from .helper_functions.filter_manifest_region import save_filtered_manifest
    from .human_protein_catalog import FEATURED_HUMAN_PROTEIN_QUERIES, get_human_protein_catalog
except ImportError:
    from analysis import AnalysisError, DEFAULT_GENE_NAME, DEFAULT_REGION, DEFAULT_REPORT_NAME, run_analysis
    from gene_region_extraction import find_gene_region
    from helper_functions.filter_manifest_region import save_filtered_manifest
    from human_protein_catalog import FEATURED_HUMAN_PROTEIN_QUERIES, get_human_protein_catalog

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
PREPROCESSED_MANIFEST_DIR = Path(__file__).resolve().parent / "gene_data"
SESSION_PREPROCESS_KEY = "preprocess_state"

app = Flask(__name__, template_folder=str(Path(__file__).resolve().parent / "templates"))
app.secret_key = os.environ.get("NOPHIGENE_SECRET_KEY", "nophigene-local-dev")


def _as_relative_display(path: Path) -> str:
    """Render a project-relative path for form fields and UI display."""
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def _resolve_user_path(raw_path: str) -> Path:
    """Resolve a user-entered path against the project root when needed."""
    candidate = Path(raw_path.strip())
    if candidate.is_absolute():
        return candidate
    return PROJECT_ROOT / candidate


def discover_vcf_files() -> list[str]:
    """Find VCF candidates under the mounted data directory."""
    if not DATA_DIR.exists():
        return []
    matches = sorted(DATA_DIR.rglob("*.vcf.gz")) + sorted(DATA_DIR.rglob("*.vcf"))
    return [_as_relative_display(path) for path in matches]


def discover_population_stats_files() -> list[str]:
    """Find optional CSV and JSON population statistics files in ``data/``."""
    if not DATA_DIR.exists():
        return []
    matches = sorted(DATA_DIR.rglob("*.csv")) + sorted(DATA_DIR.rglob("*.json"))
    return [_as_relative_display(path) for path in matches]


def discover_manifest_files() -> list[str]:
    """Find EPIC or manifest-like CSV inputs that can be prefiltered for one gene."""
    if not DATA_DIR.exists():
        return []

    matches: list[Path] = []
    for path in DATA_DIR.rglob("*"):
        lower_name = path.name.lower()
        if not path.is_file():
            continue
        if not (
            lower_name.endswith(".csv")
            or lower_name.endswith(".csv.gz")
            or lower_name.endswith(".txt")
        ):
            continue
        if "manifest" not in lower_name and "epic" not in lower_name:
            continue
        matches.append(path)

    return [_as_relative_display(path) for path in sorted(matches)]


def discover_idat_prefixes() -> list[str]:
    """Discover IDAT sample prefixes that have both green and red channels."""
    if not DATA_DIR.exists():
        return []

    prefixes: list[str] = []
    seen: set[str] = set()
    for green_file in sorted(DATA_DIR.rglob("*_Grn.idat")):
        red_file = green_file.with_name(green_file.name.replace("_Grn.idat", "_Red.idat"))
        if not red_file.exists():
            continue

        prefix = green_file.with_name(green_file.name[: -len("_Grn.idat")])
        display = _as_relative_display(prefix)
        if display not in seen:
            seen.add(display)
            prefixes.append(display)
    return prefixes


def _render_table(df: pd.DataFrame, rows: int = 12) -> str:
    """Render a compact preview table for the result cards."""
    return df.head(rows).to_html(index=False, classes="preview-table", border=0)


def _empty_form_state() -> dict[str, str]:
    """Return the initial analysis form state used on first page load."""
    vcf_files = discover_vcf_files()
    idat_prefixes = discover_idat_prefixes()
    popstats_files = discover_population_stats_files()

    return {
        "vcf": vcf_files[0] if vcf_files else "",
        "idat": idat_prefixes[0] if idat_prefixes else "",
        "out": f"results/{DEFAULT_REPORT_NAME}",
        "region": DEFAULT_REGION,
        "popstats": "",
        "manifest_file": "",
        "suggested_popstats": popstats_files[0] if popstats_files else "",
    }


def _empty_preprocess_state(manifest_files: list[str]) -> dict[str, Any]:
    """Return the default preprocessing state."""
    return {
        "gene_name": DEFAULT_GENE_NAME,
        "region": DEFAULT_REGION,
        "manifest_source": manifest_files[0] if manifest_files else "",
        "filtered_manifest": "",
        "region_candidates": [],
        "selected_sources": [],
        "region_ready": False,
        "manifest_ready": False,
        "analysis_ready": False,
        "probe_count": 0,
        "build": "hg19",
        "logs": [],
        "region_recently_updated": False,
    }


def _load_preprocess_state(manifest_files: list[str]) -> dict[str, Any]:
    """Load the persisted preprocessing state from the Flask session."""
    state = _empty_preprocess_state(manifest_files)
    saved_state = session.get(SESSION_PREPROCESS_KEY)
    if isinstance(saved_state, dict):
        state.update(saved_state)

    if not state.get("manifest_source") and manifest_files:
        state["manifest_source"] = manifest_files[0]
    return state


def _store_preprocess_state(state: dict[str, Any]) -> None:
    """Persist preprocessing state back to the Flask session."""
    session[SESSION_PREPROCESS_KEY] = {
        "gene_name": str(state.get("gene_name", DEFAULT_GENE_NAME)),
        "region": str(state.get("region", DEFAULT_REGION)),
        "manifest_source": str(state.get("manifest_source", "")),
        "filtered_manifest": str(state.get("filtered_manifest", "")),
        "region_candidates": list(state.get("region_candidates", [])),
        "selected_sources": list(state.get("selected_sources", [])),
        "region_ready": bool(state.get("region_ready", False)),
        "manifest_ready": bool(state.get("manifest_ready", False)),
        "analysis_ready": bool(state.get("analysis_ready", False)),
        "probe_count": int(state.get("probe_count", 0)),
        "build": str(state.get("build", "hg19")),
        "logs": list(state.get("logs", []))[-160:],
        "region_recently_updated": bool(state.get("region_recently_updated", False)),
    }
    session.modified = True


def _append_preprocess_log(state: dict[str, Any], message: str, *, stream: str = "stdout") -> None:
    """Append a timestamp-free developer log line to the preprocessing state."""
    normalized = str(message).rstrip()
    if not normalized:
        return

    logs = list(state.get("logs", []))
    for line in normalized.splitlines():
        logs.append(f"[{stream}] {line}")
    state["logs"] = logs[-160:]


def _capture_preprocess_call(operation, *args, **kwargs):
    """Run a preprocessing helper while capturing any stdout/stderr it emits."""
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
        result = operation(*args, **kwargs)
    return result, stdout_buffer.getvalue(), stderr_buffer.getvalue()


def _apply_preprocessing_defaults(form: dict[str, str], preprocess_state: dict[str, Any]) -> None:
    """Propagate preprocessing results into the analysis form defaults."""
    if preprocess_state.get("region"):
        form["region"] = str(preprocess_state["region"])
    if preprocess_state.get("manifest_source"):
        form["manifest_file"] = str(preprocess_state["manifest_source"])
    elif preprocess_state.get("filtered_manifest"):
        form["manifest_file"] = str(preprocess_state["filtered_manifest"])

    default_output = f"results/{DEFAULT_REPORT_NAME}"
    if form["out"] == default_output and preprocess_state.get("gene_name"):
        form["out"] = f"results/{str(preprocess_state['gene_name']).lower()}_report.html"


def _build_field_info(
    form: dict[str, str],
    *,
    preprocess_state: dict[str, Any],
    vcf_files: list[str],
    idat_prefixes: list[str],
    popstats_files: list[str],
) -> dict[str, dict[str, str]]:
    """Build UI metadata for field examples and foldable biological explanations."""
    return {
        "vcf": {
            "example": form["vcf"] or (vcf_files[0] if vcf_files else "data/gene.vcf.gz"),
            "details": (
                "A VCF, or Variant Call Format file, lists sequence differences observed "
                "between your sample and the reference genome. In this app it provides the "
                "gene-region SNPs and other variant calls that will be summarized in the report."
            ),
        },
        "idat": {
            "example": form["idat"] or (idat_prefixes[0] if idat_prefixes else "data/202277800037_R01C01"),
            "details": (
                "An IDAT base path points to the paired red and green Illumina methylation "
                "array intensity files. Those raw signal files are what methylprep uses to "
                "calculate probe-level methylation beta values for the preprocessed gene subset."
            ),
        },
        "out": {
            "example": form["out"] or f"results/{DEFAULT_REPORT_NAME}",
            "details": (
                "This controls where the generated report artifact is written so you can review "
                "or share the analysis output after the run."
            ),
        },
        "region": {
            "example": form["region"] or str(preprocess_state.get("region", DEFAULT_REGION)),
            "details": (
                "The genomic region limits the analysis to a specific coordinate interval. "
                "In this workflow it should usually come from preprocessing, where the gene name "
                "is resolved to the active genome interval."
            ),
        },
        "popstats": {
            "example": form["popstats"]
            or form["suggested_popstats"]
            or (popstats_files[0] if popstats_files else "data/reference_population.json"),
            "details": (
                "Population statistics files add cohort-level context such as allele frequency "
                "or prevalence summaries. They help interpret whether a gene-region variant looks rare, "
                "common, or enriched in reference populations."
            ),
        },
        "manifest_file": {
            "example": form["manifest_file"]
            or str(preprocess_state.get("manifest_source", "data/infinium-methylationepic-manifest-file.csv")),
            "details": (
                "This should point to the full EPIC manifest source file. During analysis the app saves a "
                "gene-specific subset like src/gene_data/GENE_epigenetics_hg19.csv before methylprep runs, "
                "and that saved subset is then reused for the probe annotation join."
            ),
        },
    }


def _build_preprocess_field_info(
    preprocess_state: dict[str, Any],
    *,
    manifest_files: list[str],
) -> dict[str, dict[str, str]]:
    """Build field hints for the preprocessing tab."""
    return {
        "gene_name": {
            "example": str(preprocess_state.get("gene_name") or DEFAULT_GENE_NAME),
            "details": (
                "Enter an HGNC-style gene symbol. The preprocessing workflow will look up the genomic "
                "coordinates for that gene and use the resulting interval in the downstream analysis."
            ),
        },
        "region": {
            "example": str(preprocess_state.get("region") or DEFAULT_REGION),
            "details": (
                "This field is filled by the gene-region lookup step and can also be adjusted manually if "
                "you want to use a custom span such as promoter-plus-gene coverage."
            ),
        },
        "manifest_source": {
            "example": str(preprocess_state.get("manifest_source") or (manifest_files[0] if manifest_files else "data/infinium-methylationepic-manifest.csv")),
            "details": (
                "Use the full EPIC manifest here. The app will filter it down to just the selected gene interval "
                "and save a much smaller CSV subset into src/gene_data for downstream methylation processing."
            ),
        },
    }


def _build_preprocess_result(preprocess_state: dict[str, Any]) -> dict[str, Any] | None:
    """Build the preprocessing result payload used by the lower status panel."""
    if not any(
        [
            preprocess_state.get("region_ready"),
            preprocess_state.get("manifest_ready"),
            preprocess_state.get("filtered_manifest"),
            preprocess_state.get("logs"),
        ]
    ):
        return None

    filtered_manifest = str(preprocess_state.get("filtered_manifest", "")).strip()
    preview_html = None
    if filtered_manifest:
        preview_path = _resolve_user_path(filtered_manifest)
        if preview_path.exists():
            try:
                preview_df = pd.read_csv(preview_path).head(12)
                preview_html = _render_table(preview_df, rows=12)
            except Exception:
                preview_html = None

    return {
        "gene_name": str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)),
        "region": str(preprocess_state.get("region", "")),
        "manifest_source": str(preprocess_state.get("manifest_source", "")),
        "filtered_manifest": filtered_manifest,
        "region_ready": bool(preprocess_state.get("region_ready", False)),
        "manifest_ready": bool(preprocess_state.get("manifest_ready", False)),
        "analysis_ready": bool(preprocess_state.get("analysis_ready", False)),
        "probe_count": int(preprocess_state.get("probe_count", 0)),
        "selected_sources": list(preprocess_state.get("selected_sources", [])),
        "region_candidates": list(preprocess_state.get("region_candidates", [])),
        "build": str(preprocess_state.get("build", "hg19")),
        "preview_html": preview_html,
        "logs": list(preprocess_state.get("logs", [])),
        "region_recently_updated": bool(preprocess_state.get("region_recently_updated", False)),
        "progress_percent": (
            100
            if preprocess_state.get("manifest_ready")
            else 50
            if preprocess_state.get("region_ready")
            else 0
        ),
        "steps": [
            {
                "title": "Find Region from Gene Name",
                "status": "complete" if preprocess_state.get("region_ready") else "pending",
                "summary": (
                    f"Resolved to {preprocess_state.get('region', DEFAULT_REGION)}"
                    if preprocess_state.get("region_ready")
                    else "Waiting for a gene-symbol lookup."
                ),
            },
            {
                "title": "Select Methylation Data",
                "status": "complete" if preprocess_state.get("manifest_ready") else "pending",
                "summary": (
                    f"{preprocess_state.get('probe_count', 0)} probes saved to {filtered_manifest}"
                    if preprocess_state.get("manifest_ready")
                    else "Waiting for the filtered EPIC manifest subset."
                ),
            },
        ],
    }


@app.get("/api/human-proteins")
def human_proteins_api() -> Any:
    """Return a page of the live human protein catalog for the UI tab."""
    query = request.args.get("q", "").strip()
    cursor = request.args.get("cursor", "").strip() or None
    longevity_page_raw = request.args.get("longevity_page", "1").strip()
    reviewed_only_raw = request.args.get("reviewed_only", "1").strip().lower()
    longevity_only_raw = request.args.get("longevity_only", "0").strip().lower()
    reviewed_only = reviewed_only_raw not in {"0", "false", "no"}
    longevity_only = longevity_only_raw in {"1", "true", "yes"}
    try:
        longevity_page = max(1, int(longevity_page_raw))
    except ValueError:
        longevity_page = 1

    payload = get_human_protein_catalog(
        query=query,
        reviewed_only=reviewed_only,
        cursor=cursor,
        longevity_only=longevity_only,
        longevity_page=longevity_page,
    )
    return jsonify(payload)


@app.route("/", methods=["GET", "POST"])
def index() -> str:
    """Render the landing page and handle preprocessing and analysis submissions."""
    result = None
    analysis_error = None
    preprocess_error = None
    preprocess_notice = None
    initial_tab = "overview"

    vcf_files = discover_vcf_files()
    idat_prefixes = discover_idat_prefixes()
    popstats_files = discover_population_stats_files()
    manifest_files = discover_manifest_files()

    preprocess_state = _load_preprocess_state(manifest_files)
    analysis_unlocked = bool(preprocess_state.get("analysis_ready", False))

    form = _empty_form_state()
    _apply_preprocessing_defaults(form, preprocess_state)

    if request.method == "POST":
        workflow = request.form.get("workflow", "analysis").strip()

        if workflow == "preprocess":
            initial_tab = "preprocessing"
            previous_gene_name = str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)).strip().upper()
            requested_gene_name = (
                request.form.get("gene_name", "").strip() or DEFAULT_GENE_NAME
            ).upper()
            preprocess_state.update(
                {
                    "gene_name": requested_gene_name,
                    "region": request.form.get("preprocess_region", "").strip()
                    or str(preprocess_state.get("region", DEFAULT_REGION)),
                    "manifest_source": request.form.get("manifest_source", "").strip()
                    or str(preprocess_state.get("manifest_source", "")),
                }
            )
            if requested_gene_name != previous_gene_name:
                preprocess_state["region_ready"] = False
                preprocess_state["manifest_ready"] = False
                preprocess_state["analysis_ready"] = False
                preprocess_state["filtered_manifest"] = ""
                preprocess_state["probe_count"] = 0
                preprocess_state["selected_sources"] = []
                preprocess_state["region_candidates"] = []
                preprocess_state["logs"] = []
                preprocess_state["region_recently_updated"] = False
            preprocess_action = request.form.get("preprocess_action", "").strip()

            try:
                if preprocess_action == "find_region":
                    _append_preprocess_log(
                        preprocess_state,
                        f"Starting gene-region lookup for {preprocess_state['gene_name']}.",
                    )
                    lookup, captured_stdout, captured_stderr = _capture_preprocess_call(
                        find_gene_region,
                        preprocess_state["gene_name"],
                    )
                    _append_preprocess_log(preprocess_state, captured_stdout, stream="stdout")
                    _append_preprocess_log(preprocess_state, captured_stderr, stream="stderr")
                    preprocess_state["gene_name"] = str(lookup["gene_name"])
                    preprocess_state["region"] = str(lookup["selected_region"])
                    preprocess_state["selected_sources"] = list(lookup["selected_sources"])
                    preprocess_state["region_candidates"] = list(lookup["candidate_regions"])
                    preprocess_state["region_ready"] = True
                    preprocess_state["manifest_ready"] = False
                    preprocess_state["analysis_ready"] = False
                    preprocess_state["filtered_manifest"] = ""
                    preprocess_state["probe_count"] = 0
                    preprocess_state["region_recently_updated"] = True
                    _append_preprocess_log(
                        preprocess_state,
                        f"Resolved {preprocess_state['gene_name']} to {preprocess_state['region']} using {', '.join(preprocess_state['selected_sources']) or 'the available sources'}.",
                    )
                    preprocess_notice = (
                        f"Resolved {preprocess_state['gene_name']} to {preprocess_state['region']}."
                    )
                elif preprocess_action == "select_methylation":
                    if not preprocess_state.get("region"):
                        raise AnalysisError(
                            "Resolve the gene coordinates first, or enter a region before selecting methylation data."
                        )
                    if not preprocess_state.get("manifest_source"):
                        raise AnalysisError(
                            "Enter an EPIC manifest path before selecting methylation data."
                        )

                    manifest_source = str(_resolve_user_path(str(preprocess_state["manifest_source"])))
                    _append_preprocess_log(
                        preprocess_state,
                        f"Filtering EPIC manifest {manifest_source} for {preprocess_state['region']}.",
                    )
                    selection, captured_stdout, captured_stderr = _capture_preprocess_call(
                        save_filtered_manifest,
                        gene_name=str(preprocess_state["gene_name"]),
                        manifest_path=manifest_source,
                        region=str(preprocess_state["region"]),
                        genome_build=str(preprocess_state.get("build", "hg19")),
                        output_dir=PREPROCESSED_MANIFEST_DIR,
                    )
                    _append_preprocess_log(preprocess_state, captured_stdout, stream="stdout")
                    _append_preprocess_log(preprocess_state, captured_stderr, stream="stderr")
                    preprocess_state["filtered_manifest"] = _as_relative_display(selection["output_path"])
                    preprocess_state["manifest_ready"] = True
                    preprocess_state["analysis_ready"] = True
                    preprocess_state["probe_count"] = int(selection["probe_count"])
                    preprocess_state["region_recently_updated"] = False
                    _append_preprocess_log(
                        preprocess_state,
                        f"Saved {selection['probe_count']} probes to {preprocess_state['filtered_manifest']}.",
                    )
                    preprocess_notice = (
                        f"Saved {selection['probe_count']} probe(s) to {preprocess_state['filtered_manifest']}."
                    )
                else:
                    raise AnalysisError("Choose a preprocessing action before submitting the form.")
            except (AnalysisError, ValueError) as exc:
                preprocess_error = str(exc)
                _append_preprocess_log(preprocess_state, str(exc), stream="stderr")
            except Exception as exc:
                preprocess_error = str(exc)
                _append_preprocess_log(preprocess_state, traceback.format_exc(), stream="stderr")

            _store_preprocess_state(preprocess_state)
            analysis_unlocked = bool(preprocess_state.get("analysis_ready", False))
            form = _empty_form_state()
            _apply_preprocessing_defaults(form, preprocess_state)

        else:
            if not analysis_unlocked:
                preprocess_error = "Complete preprocessing before running the analysis workflow."
                initial_tab = "preprocessing"
            else:
                initial_tab = "analysis"
                form = _empty_form_state()
                _apply_preprocessing_defaults(form, preprocess_state)
                form.update(
                    {
                        "vcf": request.form.get("vcf", "").strip(),
                        "idat": request.form.get("idat", "").strip(),
                        "out": request.form.get("out", "").strip() or form["out"],
                        "region": request.form.get("region", "").strip() or form["region"],
                        "popstats": request.form.get("popstats", "").strip(),
                        "manifest_file": request.form.get("manifest_file", "").strip() or form["manifest_file"],
                    }
                )

                try:
                    analysis_result = run_analysis(
                        vcf_path=str(_resolve_user_path(form["vcf"])),
                        idat_base=str(_resolve_user_path(form["idat"])),
                        output_path=str(_resolve_user_path(form["out"])),
                        gene_name=str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)),
                        region=form["region"],
                        popstats_source=str(_resolve_user_path(form["popstats"])) if form["popstats"] else None,
                        manifest_filepath=(
                            str(_resolve_user_path(form["manifest_file"])) if form["manifest_file"] else None
                        ),
                    )

                    methylation_probe_preview = analysis_result.methylation_insights.get("probe_preview")
                    result = {
                        "report_path": _as_relative_display(analysis_result.report_path),
                        "methylation_output_path": _as_relative_display(analysis_result.methylation_output_path),
                        "variant_count": len(analysis_result.variants),
                        "methylation_count": len(analysis_result.methylation),
                        "variant_preview": _render_table(analysis_result.variants),
                        "methylation_preview": _render_table(analysis_result.methylation),
                        "popstats_present": analysis_result.popstats is not None,
                        "variant_interpretations": analysis_result.variant_interpretations,
                        "population_insights": analysis_result.population_insights,
                        "methylation_insights": {
                            **analysis_result.methylation_insights,
                            "probe_preview": (
                                _render_table(methylation_probe_preview, rows=10)
                                if isinstance(methylation_probe_preview, pd.DataFrame)
                                and not methylation_probe_preview.empty
                                else None
                            ),
                        },
                        "knowledge_base_name": analysis_result.knowledge_base.get(
                            "database_name", "Local interpretation database"
                        ),
                        "knowledge_base_version": analysis_result.knowledge_base.get("version", "curated"),
                        "population_database_name": analysis_result.population_database.get(
                            "database_name", "Local population database"
                        ),
                        "population_database_version": analysis_result.population_database.get(
                            "version", "curated"
                        ),
                    }
                except AnalysisError as exc:
                    analysis_error = str(exc)

    preprocess_result = _build_preprocess_result(preprocess_state)
    available_tabs = ["overview", "preprocessing", "proteins", "structure"]
    if analysis_unlocked:
        available_tabs.insert(2, "analysis")
    if initial_tab not in available_tabs:
        initial_tab = "preprocessing" if "preprocessing" in available_tabs else "overview"

    return render_template(
        "index.html",
        form=form,
        error=analysis_error,
        preprocess_error=preprocess_error,
        preprocess_notice=preprocess_notice,
        preprocess_state=preprocess_state,
        preprocess_result=preprocess_result,
        analysis_unlocked=analysis_unlocked,
        result=result,
        initial_tab=initial_tab,
        field_info=_build_field_info(
            form,
            preprocess_state=preprocess_state,
            vcf_files=vcf_files,
            idat_prefixes=idat_prefixes,
            popstats_files=popstats_files,
        ),
        preprocess_field_info=_build_preprocess_field_info(
            preprocess_state,
            manifest_files=manifest_files,
        ),
        data_dir=_as_relative_display(DATA_DIR),
        results_dir=_as_relative_display(RESULTS_DIR),
        vcf_files=vcf_files,
        idat_prefixes=idat_prefixes,
        popstats_files=popstats_files,
        manifest_files=manifest_files,
        featured_protein_queries=FEATURED_HUMAN_PROTEIN_QUERIES,
    )


def run_server(host: str = "0.0.0.0", port: int = 8000, debug: bool = False) -> None:
    """Start the web server."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PREPROCESSED_MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    app.run(host=host, port=port, debug=debug)
