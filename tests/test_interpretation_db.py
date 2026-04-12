"""Regression tests for the local DRD4 interpretation database."""

from __future__ import annotations

import pandas as pd

from src.analysis import (
    build_methylation_insights,
    build_population_insights,
    build_variant_interpretations,
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
                "UCSC_RefGene_Group": "Body",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "DRD4_CGI",
            },
            {
                "probe_id": "cg20411756",
                "beta": 0.41,
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
    assert insights["methylation_effects"]
    assert insights["methylation_condition_research"]


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
