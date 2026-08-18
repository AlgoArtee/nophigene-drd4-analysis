"""Genotype-aware VCF interpretation tests."""

from __future__ import annotations

import pandas as pd

from src.analysis import (
    _load_raw_vcf_sample_fields,
    build_canonical_genotype,
    build_variant_interpretations,
    load_gene_interpretation_database,
    load_variants,
)


def test_gt_decoding_distinguishes_reference_heterozygous_and_alt_states() -> None:
    """REF/ALT define site alleles; GT defines the sample genotype."""
    assert build_canonical_genotype({"ref": "A", "alt": "G", "gt_raw": "0/0"})["genotype"] == "A/A"
    assert build_canonical_genotype({"ref": "A", "alt": "G", "gt_raw": "0/1"})["genotype"] == "A/G"
    assert build_canonical_genotype({"ref": "A", "alt": "G", "gt_raw": "1/1"})["genotype"] == "G/G"


def test_gt_decoding_supports_multiallelic_and_missing_calls() -> None:
    """Multiallelic GT codes should map back to the correct ALT alleles."""
    multiallelic = build_canonical_genotype({"ref": "A", "alt": "G,T", "gt_raw": "1/2"})
    missing = build_canonical_genotype({"ref": "A", "alt": "G", "gt_raw": "./."})

    assert multiallelic["genotype"] == "G/T"
    assert multiallelic["zygosity"] == "compound_heterozygous"
    assert multiallelic["allele_dosage_per_alt"] == {"G": 1, "T": 1}
    assert missing["zygosity"] == "missing"
    assert "missing_gt" in missing["qc_flags"]


def test_raw_vcf_sample_overlay_preserves_phased_gt_and_multiple_samples(tmp_path) -> None:
    """Raw FORMAT parsing keeps exact GT text that numeric arrays can lose."""
    vcf_path = tmp_path / "phased.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
                "1\t100\trs1\tA\tG\t60\tPASS\t.\tGT:AD:DP:GQ\t0|1:7,6:13:45\t1/1:0,12:12:50",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    fields = _load_raw_vcf_sample_fields(str(vcf_path), "1:90-110")

    assert fields[("1", 100, "A", "G", "s1")]["gt_raw"] == "0|1"
    assert fields[("1", 100, "A", "G", "s2")]["gt_raw"] == "1/1"


def test_raw_vcf_sample_overlay_reads_low_coordinate_mt_variants(tmp_path) -> None:
    """MT-RNR1's standard scope should include low-coordinate mitochondrial calls."""
    vcf_path = tmp_path / "mt.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
                "MT\t64\t.\tC\tT\t.\tPASS\tDP=906\tGT:SQ:AD:AF:DP\t0/1:9.4:854,4:0.005:858",
                "MT\t73\t.\tA\tG\t.\tPASS\tDP=933\tGT:SQ:AD:AF:DP\t1/1:98.13:0,889:1:889",
                "MT\t143\t.\tG\tA\t.\tPASS\tDP=929\tGT:SQ:AD:AF:DP\t0/1:1.53:926,3:0.003:929",
                "MT\t146\t.\tT\tC\t.\tPASS\tDP=915\tGT:SQ:AD:AF:DP\t1/1:98.13:0,914:1:914",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    fields = _load_raw_vcf_sample_fields(str(vcf_path), "MT:1-1601")

    assert fields[("MT", 64, "C", "T", "s1")]["gt_raw"] == "0/1"
    assert fields[("MT", 73, "A", "G", "s1")]["gt_raw"] == "1/1"
    assert fields[("MT", 143, "G", "A", "s1")]["sample_af"] == "0.003"
    assert fields[("MT", 146, "T", "C", "s1")]["dp"] == "914"


def test_load_variants_returns_empty_table_for_valid_region_without_calls(monkeypatch, tmp_path) -> None:
    """An empty VCF interval should not abort gene analysis."""
    vcf_path = tmp_path / "empty_region.vcf.gz"
    vcf_path.write_bytes(b"placeholder")
    monkeypatch.setattr("src.analysis.allel.read_vcf", lambda *args, **kwargs: {})

    variants = load_variants(str(vcf_path), "15:21405401-21441499")

    assert variants.empty
    assert {"chrom", "pos", "id", "gt_raw", "zygosity", "confidence_score"} <= set(variants.columns)


def test_qc_keeps_imbalanced_heterozygous_call_as_heterozygous() -> None:
    """Allelic imbalance lowers confidence, but it must not rewrite GT."""
    genotype = build_canonical_genotype(
        {
            "ref": "A",
            "alt": "G",
            "gt_raw": "0/1",
            "ad": [13, 5],
            "dp": 18,
            "gq": 42,
            "qual": 88,
            "filter_status": "PASS",
        }
    )

    assert genotype["genotype"] == "A/G"
    assert genotype["zygosity"] == "heterozygous"
    assert genotype["allele_dosage_per_alt"] == {"G": 1}
    assert "heterozygous_allelic_imbalance_mild" in genotype["qc_flags"]


def test_qc_supports_homozygous_alt_and_penalizes_non_pass_filter() -> None:
    """Call support and FILTER status should affect confidence without changing GT."""
    pass_call = build_canonical_genotype(
        {
            "ref": "A",
            "alt": "G",
            "gt_raw": "1/1",
            "ad": [0, 13],
            "dp": 13,
            "gq": 50,
            "qual": 90,
            "filter_status": "PASS",
        }
    )
    non_pass_call = build_canonical_genotype(
        {
            "ref": "A",
            "alt": "G",
            "gt_raw": "1/1",
            "ad": [0, 13],
            "dp": 13,
            "gq": 50,
            "qual": 90,
            "filter_status": "LowQual",
        }
    )

    assert pass_call["genotype"] == "G/G"
    assert pass_call["zygosity"] == "homozygous_alternate"
    assert pass_call["confidence_score"] >= 0.7
    assert non_pass_call["confidence_score"] < pass_call["confidence_score"]
    assert "filter_non_pass" in non_pass_call["qc_flags"]
