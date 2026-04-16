"""Regression coverage for the newly bundled curated gene knowledge bases."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src.analysis import (
    build_population_insights,
    build_variant_interpretations,
    load_gene_interpretation_database,
    load_gene_population_database,
)


CURATED_GENES = {
    "FOXO3": 108881028,
    "MTOR": 11166592,
    "RPS6": 19375713,
    "SIK3": 116714118,
    "FLCN": 17115526,
    "SIRT6": 4174106,
    "PRKAA1": 40759491,
    "NAMPT": 105888744,
    "CDKN2A": 21967751,
    "TERT": 1253282,
}


@pytest.mark.parametrize("gene_name,gene_start", CURATED_GENES.items())
def test_curated_gene_bundle_loads_with_manifest_subset(gene_name: str, gene_start: int) -> None:
    """Each requested curated gene should ship interpretation, population, and probe data."""
    knowledge_base = load_gene_interpretation_database(gene_name)
    population_database = load_gene_population_database(gene_name)

    assert knowledge_base is not None
    assert population_database is not None
    assert knowledge_base["gene_context"]["gene_name"] == gene_name
    assert knowledge_base["gene_context"]["gene_region"]["start"] == gene_start
    assert knowledge_base["gene_context"]["variant_effect_overview"]
    assert knowledge_base["gene_context"]["relevant_methylation_probe_ids"]
    assert population_database["database_name"].startswith(f"NophiGene {gene_name} Population")
    assert population_database["gene_population_patterns"]

    subset_path = Path("src/gene_data") / f"{gene_name}_epigenetics_hg19.csv"
    assert subset_path.exists()
    manifest_subset = pd.read_csv(subset_path)
    assert not manifest_subset.empty


def test_foxo3_curated_bundle_drives_interpretation_and_population_helpers() -> None:
    """FOXO3 should behave like a curated gene even with pattern-only population notes."""
    knowledge_base = load_gene_interpretation_database("FOXO3")
    population_database = load_gene_population_database("FOXO3")

    assert knowledge_base is not None
    assert population_database is not None

    variants = pd.DataFrame(
        [
            {
                "chrom": "6",
                "id": "rs2802292",
                "pos": 108900000,
                "ref": "T",
                "alt": "G",
                "qual": 61.2,
                "filter_pass": True,
            }
        ]
    )

    interpretation = build_variant_interpretations(
        variants,
        knowledge_base,
        region="6:108881028-109005977",
    )
    population_insights = build_population_insights(variants, knowledge_base, population_database)

    assert interpretation["gene_name"] == "FOXO3"
    assert "FOXO3 is located at" in interpretation["summary"]
    assert interpretation["matched_records"][0]["variant"] == "rs2802292"

    assert population_insights["variant_population_records"] == []
    assert population_insights["gene_population_patterns"]
    assert "No embedded allele-frequency panel is bundled for FOXO3 yet" in population_insights["summary"]
    assert population_insights["gene_population_patterns_intro"].startswith(
        "Broader population patterns curated from FOXO3"
    )
