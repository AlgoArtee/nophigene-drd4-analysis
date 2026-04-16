"""Regression tests for non-DRD4 gene-specific analysis behavior."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.analysis import (
    build_population_insights,
    build_methylation_insights,
    build_variant_interpretations,
    load_gene_interpretation_database,
    load_gene_population_database,
    load_methylation,
)
from src.webapp import (
    _apply_preprocessing_defaults,
    _build_field_info,
    _build_population_context_status,
    discover_population_stats_files,
)


def test_igf1r_gene_databases_load_from_gene_data() -> None:
    """IGF1R should load dedicated interpretation and population databases."""
    knowledge_base = load_gene_interpretation_database("IGF1R")
    population_database = load_gene_population_database("igf1r")

    assert knowledge_base is not None
    assert population_database is not None
    assert knowledge_base["gene_context"]["gene_name"] == "IGF1R"
    assert knowledge_base["gene_context"]["gene_region"]["start"] == 99191768
    assert len(knowledge_base["variant_records"]) >= 2
    assert knowledge_base["gene_context"]["relevant_methylation_probe_ids"]
    assert population_database["database_name"].startswith("NophiGene IGF1R Population")
    assert len(population_database["variant_population_records"]) >= 2
    assert population_database["gene_population_patterns"]
    assert population_database["gene_population_patterns_intro"].startswith(
        "Broader population patterns curated from cohort-focused IGF1R"
    )


def test_igf1r_curated_copy_replaces_drd4_text() -> None:
    """Gene-specific IGF1R interpretation and methylation summaries should stay gene-specific."""
    knowledge_base = load_gene_interpretation_database("IGF1R")
    assert knowledge_base is not None

    variants = pd.DataFrame(
        [
            {
                "chrom": "15",
                "id": "rs2229765",
                "pos": 99300000,
                "ref": "G",
                "alt": "A",
                "qual": 61.2,
                "filter_pass": True,
            }
        ]
    )
    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg19620752",
                "beta": 0.34,
                "UCSC_RefGene_Group": "TSS1500",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "chr15:99190446-99194559",
            },
            {
                "probe_id": "cg07779120",
                "beta": 0.57,
                "UCSC_RefGene_Group": "Body",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "chr15:99190446-99194559",
            },
        ]
    )

    interpretation = build_variant_interpretations(
        variants,
        knowledge_base,
        region="15:99191768-99507759",
    )
    methylation_insights = build_methylation_insights(methylation, knowledge_base)

    assert interpretation["gene_name"] == "IGF1R"
    assert "IGF1R is located at" in interpretation["summary"]
    assert "DRD4" not in interpretation["summary"]
    assert interpretation["matched_records"][0]["variant"] == "rs2229765"
    assert interpretation["sample_highlights"]["summary"].startswith("This sample matched 1 curated IGF1R")
    assert interpretation["curated_named_markers_summary"].startswith(
        "The local IGF1R bundle seeds 2 curated named marker(s)."
    )
    assert len(interpretation["curated_named_markers"]) == 2

    matched_marker = next(
        item for item in interpretation["curated_named_markers"] if item["variant"] == "rs2229765"
    )
    assert matched_marker["observed_in_run"] is True
    assert matched_marker["observed_variants"] == ["rs2229765"]
    assert any(
        item["variant"] == "rs2016347" and item["observed_in_run"] is False
        for item in interpretation["curated_named_markers"]
    )

    assert methylation_insights["gene_name"] == "IGF1R"
    assert methylation_insights["observed_probe_count"] == 2
    assert "IGF1R" in methylation_insights["summary"]
    assert "DRD4" not in methylation_insights["summary"]


def test_igf1r_population_insights_use_curated_population_database() -> None:
    """IGF1R population summaries should expose the bundled curated frequency records."""
    knowledge_base = load_gene_interpretation_database("IGF1R")
    population_database = load_gene_population_database("IGF1R")

    assert knowledge_base is not None
    assert population_database is not None

    variants = pd.DataFrame(
        [
            {
                "chrom": "15",
                "id": "rs2229765",
                "pos": 99478225,
                "ref": "G",
                "alt": "A",
                "qual": 61.2,
                "filter_pass": True,
            }
        ]
    )

    insights = build_population_insights(variants, knowledge_base, population_database)

    assert "IGF1R" in insights["summary"]
    assert "DRD4" not in insights["gene_population_patterns_intro"]
    assert insights["gene_population_patterns_intro"].startswith(
        "Broader population patterns curated from cohort-focused IGF1R"
    )

    matched_record = next(
        item for item in insights["variant_population_records"] if item["variant"] == "rs2229765"
    )
    assert matched_record["observed_in_run"] is True
    assert matched_record["population_extremes"] is not None
    assert matched_record["top_level_location_frequencies"]


def test_load_methylation_retries_without_custom_manifest(monkeypatch, tmp_path: Path) -> None:
    """A bad custom manifest should retry with methylprep's default manifest."""
    sample_prefix = tmp_path / "sample"
    for suffix in ("_Grn.idat", "_Red.idat"):
        sample_prefix.with_name(sample_prefix.name + suffix).write_bytes(b"")

    region_manifest = tmp_path / "igf1r_subset.csv"
    pd.DataFrame(
        [
            {
                "IlmnID": "cg19620752",
                "CHR": "15",
                "MAPINFO": 99192520,
                "GencodeBasicV12_NAME": "IGF1R",
                "UCSC_RefGene_Name": "IGF1R",
                "UCSC_RefGene_Accession": "NM_000875",
                "UCSC_RefGene_Group": "TSS1500",
                "UCSC_CpG_Islands_Name": "chr15:99190446-99194559",
                "Relation_to_UCSC_CpG_Island": "Island",
            }
        ]
    ).to_csv(region_manifest, index=False)

    monkeypatch.setattr("src.analysis._prepare_gene_manifest_subset", lambda **_: region_manifest)

    pipeline_calls: list[str | None] = []

    def fake_run_pipeline(data_dir: str, *, manifest_filepath: str | None):
        pipeline_calls.append(manifest_filepath)
        if manifest_filepath:
            raise ValueError("Number of passed names did not match number of header fields in the file")
        return pd.DataFrame(
            {sample_prefix.name: [0.42]},
            index=pd.Index(["cg19620752"], name="probe_id"),
        )

    monkeypatch.setattr("src.analysis._run_methylprep_pipeline", fake_run_pipeline)

    methylation = load_methylation(
        str(sample_prefix),
        manifest_filepath="data/custom_manifest.csv",
        gene_name="IGF1R",
        region="15:99191768-99507759",
    )

    assert pipeline_calls == ["data/custom_manifest.csv", None]
    assert methylation["probe_id"].tolist() == ["cg19620752"]
    assert methylation["beta"].tolist() == [0.42]


def test_preprocessing_defaults_do_not_force_custom_manifest() -> None:
    """The analysis form should keep the custom manifest override empty by default."""
    form = {
        "vcf": "",
        "idat": "",
        "out": "results/drd4_report.html",
        "region": "",
        "popstats": "",
        "manifest_file": "",
        "suggested_popstats": "",
    }
    preprocess_state = {
        "gene_name": "IGF1R",
        "region": "15:99191768-99507759",
        "manifest_source": "data/infinium-methylationepic-v-1-0-b5-manifest-file.csv",
        "filtered_manifest": "src/gene_data/IGF1R_epigenetics_hg19.csv",
    }

    _apply_preprocessing_defaults(form, preprocess_state)

    assert form["region"] == "15:99191768-99507759"
    assert form["manifest_file"] == ""
    assert form["out"] == "results/igf1r_report.html"


def test_population_stats_suggestions_ignore_processed_sample_outputs(
    monkeypatch, tmp_path: Path
) -> None:
    """Population-stat suggestions should prefer reference files over sample outputs."""
    project_root = tmp_path
    data_dir = project_root / "data"
    processed_dir = data_dir / "202277800037"
    processed_dir.mkdir(parents=True)

    (processed_dir / "202277800037_R01C01_processed.csv").write_text("probe_id,beta\ncg1,0.42\n", encoding="utf-8")
    (data_dir / "reference_population.json").write_text("[]\n", encoding="utf-8")
    (data_dir / "gnomad_ancestry_frequencies.csv").write_text("variant,frequency\nrs1,0.12\n", encoding="utf-8")
    (data_dir / "sample_sheet.csv").write_text("sample\nx\n", encoding="utf-8")
    (data_dir / "epic_manifest.csv").write_text("IlmnID\ncg1\n", encoding="utf-8")

    monkeypatch.setattr("src.webapp.PROJECT_ROOT", project_root)
    monkeypatch.setattr("src.webapp.DATA_DIR", data_dir)

    popstats_files = discover_population_stats_files()

    assert popstats_files == [
        "data/gnomad_ancestry_frequencies.csv",
        "data/reference_population.json",
    ]

    field_info = _build_field_info(
        {
            "vcf": "",
            "idat": "",
            "out": "results/igf1r_report.html",
            "region": "",
            "popstats": "",
            "manifest_file": "",
            "suggested_popstats": "data/reference_population.json",
        },
        preprocess_state={
            "gene_name": "IGF1R",
            "region": "15:99191768-99507759",
            "manifest_source": "",
        },
        vcf_files=[],
        idat_prefixes=[],
        popstats_files=popstats_files,
    )

    assert field_info["popstats"]["example"] == "data/reference_population.json"


def test_population_context_status_reports_curated_database_without_sidecar() -> None:
    """The result viewer should treat the built-in population DB as loaded context."""
    status = _build_population_context_status(
        popstats=None,
        population_database={
            "database_name": "NophiGene IGF1R Population Database",
            "version": "2026-04-15",
        },
        population_insights={
            "variant_population_records": [{"variant": "rs2229765"}],
            "gene_population_patterns": [],
        },
    )

    assert status == "Database"
