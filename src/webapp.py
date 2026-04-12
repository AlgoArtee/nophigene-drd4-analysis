"""Flask-based web UI for the DRD4 analysis workflow."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from flask import Flask, render_template, request

try:
    from .analysis import AnalysisError, DEFAULT_REGION, DEFAULT_REPORT_NAME, run_analysis
except ImportError:
    from analysis import AnalysisError, DEFAULT_REGION, DEFAULT_REPORT_NAME, run_analysis

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

app = Flask(__name__, template_folder=str(Path(__file__).resolve().parent / "templates"))


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
    """Return the initial form state used on first page load."""
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


def _build_field_info(
    form: dict[str, str],
    *,
    vcf_files: list[str],
    idat_prefixes: list[str],
    popstats_files: list[str],
) -> dict[str, dict[str, str]]:
    """Build UI metadata for field examples and foldable biological explanations."""
    return {
        "vcf": {
            "example": form["vcf"] or (vcf_files[0] if vcf_files else "data/drd4.vcf.gz"),
            "details": (
                "A VCF, or Variant Call Format file, lists sequence differences observed "
                "between your sample and the reference genome. In this app it provides the "
                "DRD4-region SNPs and other variant calls that will be summarized in the report."
            ),
        },
        "idat": {
            "example": form["idat"] or (idat_prefixes[0] if idat_prefixes else "data/202277800037_R01C01"),
            "details": (
                "An IDAT base path points to the paired red and green Illumina methylation "
                "array intensity files. Those raw signal files are what methylprep uses to "
                "calculate probe-level methylation beta values for the DRD4-associated probes."
            ),
        },
        "out": {
            "example": form["out"] or f"results/{DEFAULT_REPORT_NAME}",
            "details": (
                "This is not a biological input, but it controls where the generated report "
                "artifact is written so you can review or share the analysis output after the run."
            ),
        },
        "region": {
            "example": form["region"] or DEFAULT_REGION,
            "details": (
                "The genomic region limits the analysis to a specific coordinate interval. "
                "Biologically, this defines the stretch of the genome around DRD4 whose variants "
                "you want to inspect in the current run."
            ),
        },
        "popstats": {
            "example": form["popstats"] or form["suggested_popstats"] or (popstats_files[0] if popstats_files else "data/gnomad.json"),
            "details": (
                "Population statistics files add cohort-level context such as allele frequency "
                "or prevalence summaries. They help interpret whether a DRD4 variant looks rare, "
                "common, or enriched in reference populations."
            ),
        },
        "manifest_file": {
            "example": form["manifest_file"] or "data/custom_manifest.csv.gz",
            "details": (
                "A manifest maps array probe identifiers to genomic coordinates and annotations. "
                "In methylation analysis, it is what turns raw probe signals into biologically "
                "located CpG measurements tied to genes, islands, enhancers, or other features."
            ),
        },
    }


@app.route("/", methods=["GET", "POST"])
def index() -> str:
    """Render the landing page and handle analysis submissions."""
    form = _empty_form_state()
    result = None
    error = None
    initial_tab = "overview"
    vcf_files = discover_vcf_files()
    idat_prefixes = discover_idat_prefixes()
    popstats_files = discover_population_stats_files()

    if request.method == "POST":
        initial_tab = "analysis"
        form.update(
            {
                "vcf": request.form.get("vcf", "").strip(),
                "idat": request.form.get("idat", "").strip(),
                "out": request.form.get("out", "").strip() or f"results/{DEFAULT_REPORT_NAME}",
                "region": request.form.get("region", "").strip() or DEFAULT_REGION,
                "popstats": request.form.get("popstats", "").strip(),
                "manifest_file": request.form.get("manifest_file", "").strip(),
            }
        )

        try:
            analysis_result = run_analysis(
                vcf_path=str(_resolve_user_path(form["vcf"])),
                idat_base=str(_resolve_user_path(form["idat"])),
                output_path=str(_resolve_user_path(form["out"])),
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
                    "database_name", "Local DRD4 interpretation database"
                ),
                "knowledge_base_version": analysis_result.knowledge_base.get("version", "curated"),
            }
        except AnalysisError as exc:
            error = str(exc)

    return render_template(
        "index.html",
        form=form,
        error=error,
        result=result,
        initial_tab=initial_tab,
        field_info=_build_field_info(
            form,
            vcf_files=vcf_files,
            idat_prefixes=idat_prefixes,
            popstats_files=popstats_files,
        ),
        data_dir=_as_relative_display(DATA_DIR),
        results_dir=_as_relative_display(RESULTS_DIR),
        vcf_files=vcf_files,
        idat_prefixes=idat_prefixes,
        popstats_files=popstats_files,
    )


def run_server(host: str = "0.0.0.0", port: int = 8000, debug: bool = False) -> None:
    """Start the web server."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    app.run(host=host, port=port, debug=debug)
