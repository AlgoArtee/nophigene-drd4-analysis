"""Regression tests for the local DRD4 interpretation database."""

from __future__ import annotations

import pandas as pd

from src.analysis import (
    build_methylation_insights,
    build_population_insights,
    build_variant_interpretations,
    load_gene_interpretation_database,
    load_interpretation_database,
    load_population_database,
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
