"""Regression tests for the local DRD4 interpretation database."""

from __future__ import annotations

import pandas as pd

from src.analysis import (
    GENERAL_ANALYSIS_DATABASE_COLUMNS,
    build_methylation_insights,
    build_population_insights,
    build_predictive_theses,
    build_variant_interpretations,
    generate_report,
    load_gene_interpretation_database,
    load_interpretation_database,
    load_population_database,
    load_synthesis_database,
    update_general_analysis_database,
)


def test_interpretation_database_loads_expected_gene_context() -> None:
    """The bundled knowledge base should stay loadable and DRD4-specific."""
    knowledge_base = load_interpretation_database()

    assert knowledge_base["gene_context"]["gene_name"] == "DRD4"
    assert knowledge_base["gene_context"]["gene_region"]["start"] == 637269
    assert len(knowledge_base["variant_records"]) >= 5
    assert knowledge_base["gene_context"]["relevant_methylation_probe_ids"]
    assert knowledge_base["gene_context"]["variant_effect_overview"]
    assert knowledge_base["gene_context"]["methylation_effects"]
    rs1800955_record = next(
        item for item in knowledge_base["variant_records"] if item["variant"] == "rs1800955"
    )
    assert rs1800955_record["literature_findings"]
    assert rs1800955_record["literature_findings"][0]["genotypes"]


def test_coordinate_alias_can_match_curated_variant_record() -> None:
    """Curated coordinate aliases should still resolve to the intended variant."""
    knowledge_base = load_interpretation_database()
    variants = pd.DataFrame(
        [
            {
                "chrom": "11",
                "id": None,
                "pos": 636784,
                "ref": "C",
                "alt": "T",
                "qual": 48.5,
                "filter_pass": True,
            }
        ]
    )

    interpretation = build_variant_interpretations(
        variants, knowledge_base, region="11:636269-640706"
    )

    assert interpretation["matched_records"]
    assert interpretation["matched_records"][0]["variant"] == "rs1800955"
    assert interpretation["promoter_analysis"]["included"] is True
    assert interpretation["matched_records"][0]["functional_effects"]
    assert interpretation["matched_records"][0]["associated_conditions"]
    assert interpretation["matched_records"][0]["literature_findings"]
    assert interpretation["sample_highlights"]["highlight_items"]
    assert interpretation["sample_highlights"]["highlight_items"][0]["literature_findings"]
    assert interpretation["sample_highlights"]["result_table_rows"][0]["variant_label"] == "None"
    assert interpretation["sample_highlights"]["result_table_rows"][0]["change"] == "C -> T"
    assert interpretation["sample_highlights"]["result_table_rows"][0]["linked_to"]
    assert len(interpretation["region_recommendations"]) == 3


def test_gene_interval_without_promoter_is_reported_cleanly() -> None:
    """The locus audit should distinguish promoter coverage from gene-body coverage."""
    knowledge_base = load_interpretation_database()
    variants = pd.DataFrame(
        [
            {
                "chrom": "11",
                "id": None,
                "pos": 637933,
                "ref": "G",
                "alt": "A",
                "qual": 47.5,
                "filter_pass": True,
            }
        ]
    )

    interpretation = build_variant_interpretations(
        variants, knowledge_base, region="11:637269-640706"
    )

    assert interpretation["gene_analysis"]["included"] is True
    assert interpretation["promoter_analysis"]["included"] is False
    assert interpretation["gene_analysis"]["found_variant_count"] == 1
    assert interpretation["sample_highlights"]["summary"]


def test_methylation_insights_use_curated_probe_subset() -> None:
    """Gene-level methylation summaries should surface the curated DRD4 probes."""
    knowledge_base = load_interpretation_database()
    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg11335335",
                "beta": 0.72,
                "GencodeBasicV12_NAME": "DRD4",
                "UCSC_RefGene_Group": "Body",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "DRD4_CGI",
            },
            {
                "probe_id": "cg20411756",
                "beta": 0.41,
                "GencodeBasicV12_NAME": "DRD4;DRD4",
                "UCSC_RefGene_Group": "3'UTR",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "DRD4_CGI",
            },
        ]
    )

    insights = build_methylation_insights(methylation, knowledge_base)

    assert insights["gene_name"] == "DRD4"
    assert insights["observed_probe_count"] == 2
    assert insights["probe_preview"].shape[0] == 2
    assert insights["whitelist_mean_beta"] == 0.565
    assert insights["whitelist_mean_beta_probe_count"] == 2
    assert insights["gene_name_mean_beta"] == 0.565
    assert insights["gene_name_mean_beta_probe_count"] == 2
    assert insights["all_numeric_mean_beta"] == 0.565
    assert insights["all_numeric_mean_beta_probe_count"] == 2
    assert len(insights["summary_metric_rows"]) == 3
    assert insights["summary_metric_rows"][0]["metric"] == "Whitelist mean beta"
    assert insights["summary_metric_rows"][1]["metric"] == "DRD4-named row mean beta"
    assert insights["whitelist_probe_statuses"][0]["observed_in_run"] is True
    assert insights["gene_name_match_columns"] == ["GencodeBasicV12_NAME"]
    assert insights["methylation_effects"]
    assert insights["methylation_condition_research"]
    assert "relevant_methylation_probe_ids" in insights["whitelist_explanation"]


def test_methylation_insights_distinguish_curated_subset_mean_from_raw_table_mean() -> None:
    """Whitelist, gene-name, and raw-table means should stay explicit and numerically separate."""
    knowledge_base = load_interpretation_database()
    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg11335335",
                "beta": "0.91",
                "GencodeBasicV12_NAME": "DRD4",
                "UCSC_RefGene_Group": "Body",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "DRD4_CGI",
            },
            {
                "probe_id": "cg20411756",
                "beta": "0.93",
                "GencodeBasicV12_NAME": "DRD4;DRD4",
                "UCSC_RefGene_Group": "3'UTR",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "DRD4_CGI",
            },
            {
                "probe_id": "cg_not_curated",
                "beta": "0.12",
                "GencodeBasicV12_NAME": "",
                "UCSC_RefGene_Group": "Body",
                "Relation_to_UCSC_CpG_Island": "OpenSea",
                "UCSC_CpG_Islands_Name": "",
            },
        ]
    )

    insights = build_methylation_insights(methylation, knowledge_base)

    assert insights["whitelist_mean_beta"] == 0.92
    assert insights["whitelist_mean_beta_probe_count"] == 2
    assert insights["gene_name_mean_beta"] == 0.92
    assert insights["gene_name_mean_beta_probe_count"] == 2
    assert insights["raw_mean_beta"] == 0.653
    assert insights["raw_mean_beta_probe_count"] == 3
    assert "whitelist mean beta is 0.92 from 2 numeric whitelist probe(s)" in insights["summary"]


def test_methylation_insights_handle_repeated_gene_names_and_blank_annotations() -> None:
    """Gene-name means should match semicolon-delimited gene tokens while excluding NaN beta values."""
    knowledge_base = load_gene_interpretation_database("SIRT6")
    assert knowledge_base is not None

    methylation = pd.DataFrame(
        [
            {"probe_id": "cg01983454", "beta": 0.927, "GencodeBasicV12_NAME": "SIRT6;SIRT6"},
            {"probe_id": "cg09936839", "beta": 0.466, "GencodeBasicV12_NAME": "SIRT6"},
            {"probe_id": "cg10400916", "beta": 0.926, "GencodeBasicV12_NAME": "SIRT6"},
            {"probe_id": "cg13663667", "beta": None, "GencodeBasicV12_NAME": "SIRT6"},
            {"probe_id": "cg15635336", "beta": 0.02, "GencodeBasicV12_NAME": "SIRT6;SIRT6;SIRT6;SIRT6"},
            {"probe_id": "cg18065174", "beta": 0.901, "GencodeBasicV12_NAME": None},
            {"probe_id": "cg21192671", "beta": 0.891, "GencodeBasicV12_NAME": None},
            {"probe_id": "cg24152732", "beta": 0.316, "GencodeBasicV12_NAME": "SIRT6"},
        ]
    )

    insights = build_methylation_insights(methylation, knowledge_base)

    assert insights["whitelist_mean_beta"] == 0.243
    assert insights["whitelist_mean_beta_probe_count"] == 2
    assert insights["gene_name_mean_beta"] == 0.531
    assert insights["gene_name_mean_beta_probe_count"] == 5
    assert insights["gene_name_row_count"] == 6
    assert insights["all_numeric_mean_beta"] == 0.635
    assert insights["all_numeric_mean_beta_probe_count"] == 7
    assert insights["gene_name_match_columns"] == ["GencodeBasicV12_NAME"]
    assert "semicolon-delimited token" in insights["gene_name_match_rule"]


def test_methylation_insights_build_whitelist_probe_reference_rows() -> None:
    """Whitelist probe reference rows should expose linked variants, nearby loci, and papers."""
    knowledge_base = load_gene_interpretation_database("SIRT6")
    assert knowledge_base is not None

    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg15635336",
                "beta": 0.02,
                "chrom": "19",
                "pos": 4182521,
                "GencodeBasicV12_NAME": "SIRT6;SIRT6;SIRT6;SIRT6",
            },
            {
                "probe_id": "cg09936839",
                "beta": 0.466,
                "chrom": "19",
                "pos": 4181854,
                "GencodeBasicV12_NAME": "SIRT6",
            },
        ]
    )

    insights = build_methylation_insights(methylation, knowledge_base)
    assert insights["whitelist_probe_reference_rows"]
    assert insights["whitelist_probe_reference_summary"].startswith("Each whitelist probe is cross-referenced")

    probe_row = next(
        row for row in insights["whitelist_probe_reference_rows"] if row["probe_id"] == "cg15635336"
    )
    assert probe_row["observed_in_run"] is True
    assert probe_row["beta"] == 0.02
    assert probe_row["probe_locus"] == "chr19:4,182,521"
    assert any(item["label"] == "rs350846" for item in probe_row["linked_variants"])
    assert any("centSIRT6" in item["label"] for item in probe_row["linked_variants"])
    assert any(item["variant"] == "rs201182672" for item in probe_row["nearby_manifest_variants"])
    assert any("Li et al., 2016" in item["label"] for item in probe_row["papers"])
    assert any("Simon et al., 2022" in item["label"] for item in probe_row["papers"])


def test_methylation_probe_reference_rows_hide_when_no_curated_variant_was_observed() -> None:
    """Variant-linked methylation rows should disappear when no curated variant was matched in the run."""
    knowledge_base = load_gene_interpretation_database("SIRT6")
    assert knowledge_base is not None

    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg15635336",
                "beta": 0.02,
                "chrom": "19",
                "pos": 4182521,
                "GencodeBasicV12_NAME": "SIRT6;SIRT6;SIRT6;SIRT6",
            }
        ]
    )

    insights = build_methylation_insights(
        methylation,
        knowledge_base,
        matched_variant_ids=set(),
    )

    assert insights["whitelist_probe_reference_rows"] == []
    assert "table is hidden" in insights["whitelist_probe_reference_summary"]


def test_population_database_loads_expected_variant_records() -> None:
    """The bundled population database should expose DRD4 geography summaries."""
    population_database = load_population_database()

    assert population_database["database_name"].startswith("NophiGene DRD4 Population")
    assert len(population_database["variant_population_records"]) >= 3
    assert population_database["gene_population_patterns"]


def test_population_insights_flag_observed_curated_variant() -> None:
    """Observed curated SNPs should be marked in the population summary block."""
    knowledge_base = load_interpretation_database()
    population_database = load_population_database()
    variants = pd.DataFrame(
        [
            {
                "chrom": "11",
                "id": "rs1800955",
                "pos": 636784,
                "ref": "C",
                "alt": "T",
                "qual": 51.2,
                "filter_pass": True,
            }
        ]
    )

    insights = build_population_insights(variants, knowledge_base, population_database)

    assert insights["variant_population_records"]
    matched_record = next(
        item for item in insights["variant_population_records"] if item["variant"] == "rs1800955"
    )
    assert matched_record["observed_in_run"] is True
    assert matched_record["population_extremes"] is not None
    assert insights["location_groups"]


def test_generate_report_includes_variant_and_methylation_interpretation_sections(tmp_path) -> None:
    """HTML reports should include the richer interpretation/data sections now shown in the UI."""
    knowledge_base = load_interpretation_database()
    population_database = load_population_database()

    variants = pd.DataFrame(
        [
            {
                "chrom": "11",
                "id": None,
                "pos": 636784,
                "ref": "C",
                "alt": "T",
                "qual": 51.2,
                "filter_pass": True,
            }
        ]
    )
    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg11335335",
                "beta": 0.72,
                "chrom": "11",
                "pos": 637050,
                "GencodeBasicV12_NAME": "DRD4",
                "UCSC_RefGene_Group": "Body",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "DRD4_CGI",
            }
        ]
    )

    variant_interpretations = build_variant_interpretations(
        variants,
        knowledge_base,
        region="11:636269-640706",
    )
    methylation_insights = build_methylation_insights(methylation, knowledge_base)
    population_insights = build_population_insights(variants, knowledge_base, population_database)
    predictive_theses = build_predictive_theses(
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        knowledge_base=knowledge_base,
        synthesis_database=load_synthesis_database(),
    )

    report_path = generate_report(
        variants,
        methylation,
        None,
        str(tmp_path / "report.html"),
        gene_name="DRD4",
        region="11:636269-640706",
        methylation_output_path=tmp_path / "report_methylation.csv",
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        population_insights=population_insights,
        predictive_theses=predictive_theses,
    )

    report_html = report_path.read_text(encoding="utf-8")
    assert "Genetic Variant Results" in report_html
    assert "Sample Results" in report_html
    assert "Matched Variant Interpretations" in report_html
    assert "Predictive Theses" in report_html
    assert "Variant Prediction" in report_html
    assert "Methylation Prediction" in report_html
    assert "Synthesis" in report_html
    assert "Sample allele-change thesis" in report_html
    assert "report-table-shell" in report_html
    assert "table-layout: fixed" in report_html
    assert "width: min(98vw, 1800px)" in report_html
    assert "Methylation Summary Metrics" in report_html
    assert "Methylation Raw Results" in report_html


def test_general_analysis_database_adds_once_and_overwrites_by_gene(tmp_path) -> None:
    """The central database should keep one row per gene unless overwrite is requested."""
    database_path = tmp_path / "general_gene_analysis_database.csv"
    variant_interpretations = {
        "gene_region": {"display": "chr15:28,356,186-28,567,325"},
        "search_region": {"display": "chr15:28,356,000-28,567,325"},
    }
    methylation_insights = {
        "whitelist_mean_beta": 0.71,
        "gene_name_mean_beta": 0.62,
        "all_numeric_mean_beta": 0.53,
    }
    first_variants = pd.DataFrame(
        [
            {
                "chrom": "15",
                "id": "rs12913832",
                "pos": 28365618,
                "ref": "A",
                "alt": "G",
                "qual": 88.0,
                "filter_pass": True,
            }
        ]
    )
    second_variants = pd.DataFrame(
        [
            {
                "chrom": "15",
                "id": "rs12913832",
                "pos": 28365618,
                "ref": "G",
                "alt": "A",
                "qual": 91.25,
                "filter_pass": True,
            }
        ]
    )

    added = update_general_analysis_database(
        gene_name="HERC2",
        variants=first_variants,
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        database_path=database_path,
    )
    skipped = update_general_analysis_database(
        gene_name="herc2",
        variants=second_variants,
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        database_path=database_path,
    )

    database = pd.read_csv(database_path)
    assert added["action"] == "added"
    assert skipped["action"] == "skipped_existing"
    assert database.columns.tolist() == GENERAL_ANALYSIS_DATABASE_COLUMNS
    assert len(database) == 1
    assert database.loc[0, "gene"] == "HERC2"
    assert database.loc[0, "observed gene variant"] == "rs12913832"
    assert database.loc[0, "gene variant label"] == "rs12913832"
    assert database.loc[0, "change"] == "A -> G"
    assert database.loc[0, "gene location"] == "chr15:28,356,186-28,567,325"
    assert database.loc[0, "source"] == "VCF"
    assert database.loc[0, "(VCF) quality (qual)"] == 88.0
    assert database.loc[0, "mean beta whitelist"] == 0.71
    assert database.loc[0, "mean beta related to gene"] == 0.62
    assert database.loc[0, "mean beta on found probes in the area (numerical rows)"] == 0.53

    overwritten = update_general_analysis_database(
        gene_name="HERC2",
        variants=second_variants,
        variant_interpretations=variant_interpretations,
        methylation_insights=methylation_insights,
        overwrite=True,
        database_path=database_path,
    )

    database = pd.read_csv(database_path)
    assert overwritten["action"] == "overwritten"
    assert len(database) == 1
    assert database.loc[0, "change"] == "G -> A"
    assert database.loc[0, "(VCF) quality (qual)"] == 91.25
