"""Regression tests for non-DRD4 gene-specific analysis behavior."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.analysis import (
    annotate_known_variant_ids,
    build_predictive_theses,
    build_population_insights,
    build_methylation_insights,
    build_variant_interpretations,
    load_gene_interpretation_database,
    load_gene_population_database,
    load_gene_synthesis_database,
    load_methylation,
)
from src.webapp import (
    _apply_preprocessing_defaults,
    _build_analysis_scope_regions,
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


def test_herc2_gene_databases_load_from_gene_data() -> None:
    """HERC2 should load dedicated interpretation and population databases."""
    knowledge_base = load_gene_interpretation_database("HERC2")
    population_database = load_gene_population_database("herc2")
    synthesis_database = load_gene_synthesis_database("HERC2")

    assert knowledge_base is not None
    assert population_database is not None
    assert synthesis_database is not None
    assert knowledge_base["gene_context"]["gene_name"] == "HERC2"
    assert knowledge_base["gene_context"]["gene_region"]["start"] == 28356186
    assert len(knowledge_base["variant_records"]) >= 5
    assert knowledge_base["gene_context"]["relevant_methylation_probe_ids"]
    assert population_database["database_name"].startswith("NophiGene HERC2 Population")
    assert len(population_database["variant_population_records"]) >= 5
    assert population_database["gene_population_patterns"]
    assert population_database["gene_population_patterns_intro"].startswith(
        "Broader population patterns curated from HERC2/OCA2"
    )
    assert synthesis_database["database_name"].startswith("NophiGene HERC2 Predictive")
    assert synthesis_database["case_count"] == 10
    assert len(synthesis_database["cases"]) == 10
    assert "brown-eye tendency" in synthesis_database["concrete_variant_prediction"]
    assert "blue-eye tendency" in synthesis_database["concrete_variant_prediction"]
    rs12913832_rule = next(
        rule
        for rule in synthesis_database["variant_prediction_rules"]
        if rule["variant"] == "rs12913832"
    )
    assert "OCA2-driven iris melanin" in rs12913832_rule["prediction"]
    assert "{change}" in rs12913832_rule["sample_change_template"]
    assert any(
        rule.get("change") == "A -> G" and rule.get("alt_allele") == "G"
        for rule in rs12913832_rule["allele_change_rules"]
    )
    assert any(
        rule.get("change") == "G -> A" and rule.get("alt_allele") == "A"
        for rule in rs12913832_rule["allele_change_rules"]
    )


def test_reverse_strand_scope_regions_use_valid_promoter_gene_union() -> None:
    """Reverse-strand promoter+gene regions should cover both body and upstream promoter."""
    scope_regions = _build_analysis_scope_regions("SIRT6", "19:4174106-4182560")

    assert scope_regions["gene_only"] == "19:4174106-4182560"
    assert scope_regions["promoter_only"] == "19:4182561-4183560"
    assert scope_regions["promoter_plus_gene"] == "19:4174106-4183560"


def test_all_local_interpretation_databases_have_valid_combined_regions() -> None:
    """Every bundled promoter+gene recommendation should cover promoter and gene body."""
    gene_data_dir = Path(__file__).resolve().parents[1] / "src" / "gene_data"

    for database_path in gene_data_dir.glob("*_interpretation_db.json"):
        gene_name = database_path.name.removesuffix("_interpretation_db.json")
        knowledge_base = load_gene_interpretation_database(gene_name)
        assert knowledge_base is not None
        gene_context = knowledge_base["gene_context"]
        combined = gene_context["recommended_promoter_plus_gene_region"].replace("chr", "")
        combined_chrom, combined_span = combined.split(":")
        combined_start, combined_end = [int(value) for value in combined_span.split("-")]
        gene_region = gene_context["gene_region"]
        promoter_region = gene_context["promoter_review_region"]

        assert combined_start <= combined_end, database_path.name
        assert str(gene_context["chromosome"]).removeprefix("chr") == combined_chrom
        assert combined_start <= min(gene_region["start"], gene_region["end"])
        assert combined_end >= max(gene_region["start"], gene_region["end"])
        assert combined_start <= min(promoter_region["start"], promoter_region["end"])
        assert combined_end >= max(promoter_region["start"], promoter_region["end"])


def test_herc2_predictive_theses_use_concrete_eye_colour_variant_rules() -> None:
    """Matched HERC2 pigmentation markers should surface concrete eye-colour predictions."""
    knowledge_base = load_gene_interpretation_database("HERC2")
    synthesis_database = load_gene_synthesis_database("HERC2")

    assert knowledge_base is not None
    assert synthesis_database is not None

    variants = pd.DataFrame(
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
    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg14091419",
                "beta": 0.71,
                "GencodeBasicV12_NAME": "HERC2",
                "UCSC_RefGene_Group": "TSS1500",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "HERC2_CGI",
            }
        ]
    )

    interpretation = build_variant_interpretations(
        variants,
        knowledge_base,
        region="15:28356000-28567325",
    )
    methylation_insights = build_methylation_insights(
        methylation,
        knowledge_base,
        matched_variant_ids={"rs12913832"},
    )
    predictive_theses = build_predictive_theses(
        variant_interpretations=interpretation,
        methylation_insights=methylation_insights,
        knowledge_base=knowledge_base,
        synthesis_database=synthesis_database,
    )

    assert predictive_theses["variant_found"] is True
    assert "blue-eye tendency" in predictive_theses["variant_summary"]
    concrete_rows = [
        row
        for row in predictive_theses["variant_prediction_rows"]
        if row["source"] == "Sample allele-change thesis"
    ]
    assert concrete_rows
    assert "A -> G" in concrete_rows[0]["observed_signal"]
    assert "observed alternate allele is G" in concrete_rows[0]["prediction"]
    assert "blue or lighter-eye tendency" in concrete_rows[0]["prediction"]
    assert "reduced OCA2 expression" in concrete_rows[0]["research_focus"]


def test_herc2_predictive_theses_use_actual_reverse_eye_colour_change() -> None:
    """The reverse HERC2 REF -> ALT direction should produce the darker-eye thesis."""
    knowledge_base = load_gene_interpretation_database("HERC2")
    synthesis_database = load_gene_synthesis_database("HERC2")

    assert knowledge_base is not None
    assert synthesis_database is not None

    variants = pd.DataFrame(
        [
            {
                "chrom": "15",
                "id": "rs12913832",
                "pos": 28365618,
                "ref": "G",
                "alt": "A",
                "qual": 88.0,
                "filter_pass": True,
            }
        ]
    )
    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg14091419",
                "beta": 0.42,
                "GencodeBasicV12_NAME": "HERC2",
                "UCSC_RefGene_Group": "TSS1500",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "HERC2_CGI",
            }
        ]
    )

    interpretation = build_variant_interpretations(
        variants,
        knowledge_base,
        region="15:28356000-28567325",
    )
    methylation_insights = build_methylation_insights(
        methylation,
        knowledge_base,
        matched_variant_ids={"rs12913832"},
    )
    predictive_theses = build_predictive_theses(
        variant_interpretations=interpretation,
        methylation_insights=methylation_insights,
        knowledge_base=knowledge_base,
        synthesis_database=synthesis_database,
    )

    concrete_rows = [
        row
        for row in predictive_theses["variant_prediction_rows"]
        if row["source"] == "Sample allele-change thesis"
    ]
    assert concrete_rows
    assert "G -> A" in concrete_rows[0]["observed_signal"]
    assert "observed alternate allele is A" in concrete_rows[0]["prediction"]
    assert "brown or darker-eye tendency" in concrete_rows[0]["prediction"]


def test_predictive_theses_match_variant_and_all_three_methylation_views() -> None:
    """A curated variant plus three numeric methylation summaries should match four synthesis cases."""
    knowledge_base = load_gene_interpretation_database("IGF1R")
    synthesis_database = load_gene_synthesis_database("IGF1R")

    assert knowledge_base is not None
    assert synthesis_database is not None

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
    methylation = pd.DataFrame(
        [
            {
                "probe_id": "cg19620752",
                "beta": 0.82,
                "GencodeBasicV12_NAME": "IGF1R",
                "UCSC_RefGene_Group": "TSS1500",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "chr15:99190446-99194559",
            },
            {
                "probe_id": "cg07779120",
                "beta": 0.84,
                "GencodeBasicV12_NAME": "IGF1R",
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
    matched_variant_ids = {
        str(record.get("variant", "")).strip()
        for record in interpretation.get("matched_records", [])
        if str(record.get("variant", "")).strip()
    }
    methylation_insights = build_methylation_insights(
        methylation,
        knowledge_base,
        matched_variant_ids=matched_variant_ids,
    )
    predictive_theses = build_predictive_theses(
        variant_interpretations=interpretation,
        methylation_insights=methylation_insights,
        knowledge_base=knowledge_base,
        synthesis_database=synthesis_database,
    )

    assert predictive_theses["variant_found"] is True
    assert predictive_theses["matched_case_count"] == 4
    assert predictive_theses["case_catalog_size"] == 10
    assert predictive_theses["matched_cases"][0]["case_label"] == "Gene variant found"
    assert all(row["band_display"] == "High" for row in predictive_theses["methylation_prediction_rows"])
    assert all(row["matched"] is True for row in predictive_theses["methylation_prediction_rows"])
    assert "rs2229765" in predictive_theses["variant_summary"]


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


def test_curated_variant_match_survives_missing_source_id_when_coordinates_match() -> None:
    """Coordinate and allele matching should still identify curated variants when the VCF ID field is missing."""
    knowledge_base = load_gene_interpretation_database("DRD4")
    assert knowledge_base is not None

    variants = pd.DataFrame(
        [
            {
                "chrom": "11",
                "id": None,
                "pos": 636689,
                "ref": "G",
                "alt": "C",
                "qual": 75.35,
                "filter_pass": True,
            }
        ]
    )

    interpretation = build_variant_interpretations(
        variants,
        knowledge_base,
        region="11:637269-640706",
    )

    assert interpretation["matched_records"]
    assert interpretation["matched_records"][0]["variant"] == "rs747302"
    assert interpretation["matched_records"][0]["observed_variant"] == "11:636689 G>C"


def test_herc2_known_variants_are_labeled_when_source_vcf_ids_are_missing() -> None:
    """Curated HERC2 markers should fill visible IDs from coordinate matches when the source VCF leaves them blank."""
    knowledge_base = load_gene_interpretation_database("HERC2")
    population_database = load_gene_population_database("HERC2")

    assert knowledge_base is not None
    assert population_database is not None

    variants = pd.DataFrame(
        [
            {
                "chrom": "15",
                "id": None,
                "pos": 28356859,
                "ref": "C",
                "alt": "T",
                "qual": 50.0,
                "filter_pass": True,
            },
            {
                "chrom": "15",
                "id": None,
                "pos": 28365618,
                "ref": "A",
                "alt": "G",
                "qual": 42.0,
                "filter_pass": True,
            },
            {
                "chrom": "15",
                "id": None,
                "pos": 28427986,
                "ref": "T",
                "alt": "A",
                "qual": 86.13,
                "filter_pass": True,
            },
            {
                "chrom": "15",
                "id": None,
                "pos": 28513364,
                "ref": "T",
                "alt": "C",
                "qual": 92.15,
                "filter_pass": True,
            },
        ]
    )

    labeled_variants = annotate_known_variant_ids(variants, knowledge_base)
    interpretation = build_variant_interpretations(
        labeled_variants,
        knowledge_base,
        region="15:28356000-28567325",
    )
    population_insights = build_population_insights(
        labeled_variants,
        knowledge_base,
        population_database,
    )

    assert labeled_variants["id"].tolist() == [
        "rs1129038",
        "rs12913832",
        "rs7170852",
        "rs916977",
    ]
    assert labeled_variants["id_source"].tolist() == [
        "Knowledge base match",
        "Knowledge base match",
        "Knowledge base match",
        "Knowledge base match",
    ]
    assert interpretation["matched_records"]
    assert interpretation["matched_records"][0]["observed_variant"] == "rs1129038"
    assert any(record["variant"] == "rs12913832" for record in interpretation["matched_records"])
    assert any(record["variant"] == "rs916977" for record in interpretation["matched_records"])

    matched_population_record = next(
        item for item in population_insights["variant_population_records"] if item["variant"] == "rs12913832"
    )
    assert matched_population_record["observed_in_run"] is True
    assert matched_population_record["observed_variants"] == ["rs12913832"]
    assert matched_population_record["population_extremes"] is not None


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
        "analysis_scope": "",
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
    assert form["analysis_scope"] == "promoter_plus_gene"
    assert form["manifest_file"] == ""
    assert form["out"] == "results/igf1r_promoter_plus_gene_report.html"


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
