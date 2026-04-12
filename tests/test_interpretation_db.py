"""Regression tests for the local DRD4 interpretation database."""

from __future__ import annotations

import pandas as pd

from src.analysis import (
    build_methylation_insights,
    build_variant_interpretations,
    load_interpretation_database,
)


def test_interpretation_database_loads_expected_gene_context() -> None:
    """The bundled knowledge base should stay loadable and DRD4-specific."""
    knowledge_base = load_interpretation_database()

    assert knowledge_base["gene_context"]["gene_name"] == "DRD4"
    assert len(knowledge_base["variant_records"]) >= 3
    assert knowledge_base["gene_context"]["relevant_methylation_probe_ids"]


def test_coordinate_alias_can_match_curated_variant_record() -> None:
    """Curated coordinate aliases should still resolve to the intended variant."""
    knowledge_base = load_interpretation_database()
    variants = pd.DataFrame(
        [
            {
                "chrom": "11",
                "id": None,
                "pos": 62636784,
                "ref": "C",
                "alt": "T",
                "qual": 48.5,
                "filter_pass": True,
            }
        ]
    )

    interpretation = build_variant_interpretations(variants, knowledge_base)

    assert interpretation["matched_records"]
    assert interpretation["matched_records"][0]["variant"] == "rs1800955"


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
