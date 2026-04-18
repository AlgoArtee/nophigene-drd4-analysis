"""Regression tests for preprocessing form behavior in the Flask UI."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from src.webapp import app


def test_preprocessing_template_preserves_clicked_submit_button() -> None:
    """The loading helper should preserve the clicked preprocessing action."""
    template_path = Path(__file__).resolve().parent.parent / "src" / "templates" / "index.html"
    template_text = template_path.read_text(encoding="utf-8")

    assert "const submitter = event.submitter instanceof HTMLElement ? event.submitter : null;" in template_text
    assert "formElement.requestSubmit(submitter);" in template_text
    assert "temporarySubmitterInput.name = submitter.name;" in template_text
    assert "Bundled named {{ result.variant_interpretations.gene_name }} markers" in template_text
    assert "{{ result.variant_interpretations.curated_named_markers_summary }}" in template_text
    assert "{{ result.methylation_insights.whitelist_mean_beta_label }}" in template_text
    assert "{{ result.methylation_insights.gene_name_mean_beta_label }}" in template_text
    assert "{{ result.methylation_insights.all_numeric_mean_beta_label }}" in template_text
    assert "{{ result.methylation_insights.whitelist_explanation }}" in template_text
    assert "{{ result.methylation_insights.gene_name_match_rule }}" in template_text
    assert "{{ result.methylation_insights.whitelist_probe_reference_summary }}" in template_text
    assert 'data-detail-target="qa"' in template_text
    assert "{{ app_structure_qa_items }}" not in template_text
    assert "{{ item.question }}" in template_text
    assert "App Structure Q&amp;A" in template_text
    assert "Curated probe-to-variant and paper links" in template_text
    assert "These are the exact observed whitelist probe rows used to build the whitelist mean beta shown above." in template_text
    assert "flex-wrap: wrap;" in template_text
    assert "window.setTimeout(() => formElement.submit(), 40);" not in template_text
    assert "const preprocessLoadingStepsByAction = {" in template_text
    assert 'data-preprocess-loading-steps' in template_text
    assert "Confirm the selected coordinates for the current gene." in template_text
    assert ".preprocess-detail-grid {" in template_text
    assert ".preprocess-panel {\n      padding: 18px;\n      display: grid;" in template_text
    assert "panel.hidden = !isActive;" in template_text
    assert "[hidden] {" in template_text
    assert 'style="display:{% if initial_tab == ' in template_text
    assert 'panel.style.display = isActive ? "block" : "none";' in template_text
    assert 'data-variant-raw-table' in template_text
    assert 'variant-raw-data' in template_text
    assert "function renderVariantRawPage()" in template_text


def test_preprocess_find_region_submission_updates_session(monkeypatch) -> None:
    """A valid preprocessing action should resolve the gene interval and persist it."""

    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr(
        "src.webapp.find_gene_region",
        lambda gene_symbol: {
            "gene_name": gene_symbol.upper(),
            "selected_region": "chr11:639677-643057",
            "selected_sources": ["NCBI RefSeq"],
            "candidate_regions": [{"source": "NCBI RefSeq", "region": "chr11:639677-643057"}],
        },
    )

    client = app.test_client()
    response = client.post(
        "/",
        data={
            "workflow": "preprocess",
            "gene_name": "drd4",
            "preprocess_action": "find_region",
        },
    )

    assert response.status_code == 200
    assert "Resolved DRD4 to chr11:639677-643057." in response.get_data(as_text=True)

    with client.session_transaction() as session_state:
        preprocess_state = session_state["preprocess_state"]

    assert preprocess_state["gene_name"] == "DRD4"
    assert preprocess_state["region"] == "chr11:639677-643057"
    assert preprocess_state["region_ready"] is True
    assert preprocess_state["manifest_ready"] is False
    assert preprocess_state["analysis_ready"] is False
    assert preprocess_state["selected_sources"] == ["NCBI RefSeq"]


def test_get_request_resets_preprocessing_workspace(monkeypatch) -> None:
    """A fresh page load should not keep the previous gene unlocked in the UI."""
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_report_history", lambda: [])

    client = app.test_client()
    with client.session_transaction() as session_state:
        session_state["preprocess_state"] = {
            "gene_name": "IGF1R",
            "region": "15:99191768-99507759",
            "manifest_source": "data/infinium-methylationepic-v-1-0-b5-manifest-file.csv",
            "filtered_manifest": "src/gene_data/IGF1R_epigenetics_hg19.csv",
            "region_candidates": [],
            "selected_sources": [],
            "region_ready": True,
            "manifest_ready": True,
            "analysis_ready": True,
            "probe_count": 208,
            "build": "hg19",
            "logs": ["[stdout] stale state"],
            "region_recently_updated": False,
            "overwrite_filtered_manifest": False,
        }

    response = client.get("/")
    page = response.get_data(as_text=True)

    assert response.status_code == 200
    assert 'data-tab-target="analysis"' not in page
    assert 'value="IGF1R"' not in page
    assert 'value="DRD4"' in page

    with client.session_transaction() as session_state:
        assert "preprocess_state" not in session_state


def test_app_structure_page_includes_general_probe_mapping_qa(monkeypatch) -> None:
    """The App Structure page should include the general literature-to-probe mapping explanation."""
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_report_history", lambda: [])

    client = app.test_client()
    response = client.get("/")
    page = response.get_data(as_text=True)

    assert response.status_code == 200
    assert "How are the local databases built from literature, and how does a probe get mapped to a locus or variant?" in page
    assert "the current whitelist probes are usually not stored as" in page
    assert "SIRT6 on the reverse strand" in page
    assert "That nearby-locus column is purely manifest-derived proximity annotation" in page


def test_history_tab_lists_saved_reports_and_serves_artifacts(monkeypatch, tmp_path: Path) -> None:
    """The History tab should surface prior reports with openable links."""
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    report_path = results_dir / "igf1r_report.html"
    report_path.write_text("<html><body>IGF1R report</body></html>", encoding="utf-8")
    methylation_path = results_dir / "igf1r_report_methylation.csv"
    methylation_path.write_text("probe_id,beta\ncg1,0.42\n", encoding="utf-8")

    monkeypatch.setattr("src.webapp.RESULTS_DIR", results_dir)
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])

    client = app.test_client()
    response = client.get("/")
    page = response.get_data(as_text=True)

    assert response.status_code == 200
    assert 'data-tab-target="history"' in page
    assert "igf1r_report.html" in page
    assert "/results/igf1r_report.html" in page
    assert "/results/igf1r_report_methylation.csv" in page

    artifact_response = client.get("/results/igf1r_report.html")
    assert artifact_response.status_code == 200
    assert "IGF1R report" in artifact_response.get_data(as_text=True)


def test_analysis_result_keeps_full_curated_methylation_probe_preview(monkeypatch, tmp_path: Path) -> None:
    """The curated methylation probe preview should include every row used for interpretation."""
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_report_history", lambda: [])

    captured_context: dict[str, object] = {}

    def fake_render_template(template_name: str, **context: object) -> str:
        captured_context["template_name"] = template_name
        captured_context["context"] = context
        return "ok"

    probe_preview = pd.DataFrame(
        [
            {
                "probe_id": f"cg{i:08d}",
                "beta": 0.8 + i / 1000,
                "UCSC_RefGene_Group": "Body",
                "Relation_to_UCSC_CpG_Island": "Island",
                "UCSC_CpG_Islands_Name": "TEST_CGI",
            }
            for i in range(15)
        ]
    )

    monkeypatch.setattr("src.webapp.render_template", fake_render_template)
    monkeypatch.setattr(
        "src.webapp.run_analysis",
        lambda **_: SimpleNamespace(
            report_path=tmp_path / "mock_report.html",
            methylation_output_path=tmp_path / "mock_report_methylation.csv",
            variants=pd.DataFrame(
                [{"chrom": "1", "id": "rs1", "pos": 1, "ref": "A", "alt": "G", "qual": 99.0, "filter_pass": True}]
            ),
            methylation=probe_preview.copy(),
            popstats=None,
            knowledge_base={"database_name": "Mock interpretation DB", "version": "test"},
            population_database={"database_name": "Mock population DB", "version": "test"},
            population_insights={"variant_population_records": [], "gene_population_patterns": []},
            variant_interpretations={"gene_name": "TEST", "matched_records": [], "sample_highlights": {"summary": ""}},
            methylation_insights={
                "gene_name": "TEST",
                "clinical_context": "",
                "summary": "",
                "mean_beta": 0.9,
                "mean_beta_label": "Whitelist mean beta",
                "mean_beta_probe_count": 15,
                "whitelist_mean_beta": 0.9,
                "whitelist_mean_beta_label": "Whitelist mean beta",
                "whitelist_mean_beta_probe_count": 15,
                "whitelist_probe_count": 15,
                "whitelist_observed_probe_count": 15,
                "whitelist_explanation": "Mock whitelist explanation",
                "whitelist_literature_context": "Mock literature context",
                "whitelist_probe_statuses": [
                    {"probe_id": probe_id, "observed_in_run": True}
                    for probe_id in probe_preview["probe_id"].tolist()
                ],
                "whitelist_probe_reference_rows": [
                    {
                        "probe_id": "cg00000000",
                        "observed_in_run": True,
                        "beta": 0.9,
                        "probe_locus": "chr1:1",
                        "linked_variants": [{"label": "rs1", "common_name": "mock variant", "locus": "chr1:1"}],
                        "nearby_manifest_variants": [{"variant": "rs1", "distance": "0"}],
                        "papers": [{"label": "Mock Paper", "url": "https://example.com/mock", "source_variant": "rs1"}],
                    }
                ],
                "whitelist_probe_reference_summary": "Mock whitelist reference summary",
                "gene_name_mean_beta": 0.9,
                "gene_name_mean_beta_label": "TEST-named row mean beta",
                "gene_name_mean_beta_probe_count": 15,
                "gene_name_row_count": 15,
                "gene_name_match_columns": ["GencodeBasicV12_NAME"],
                "gene_name_match_rule": "Mock gene-name rule",
                "raw_mean_beta": 0.9,
                "raw_mean_beta_label": "All numeric-row mean beta",
                "raw_probe_count": 15,
                "raw_mean_beta_probe_count": 15,
                "all_numeric_mean_beta": 0.9,
                "all_numeric_mean_beta_label": "All numeric-row mean beta",
                "all_numeric_mean_beta_probe_count": 15,
                "beta_band": "high",
                "beta_band_source_label": "Whitelist mean beta",
                "observed_probe_count": 15,
                "curated_probe_count": 15,
                "probe_ids": probe_preview["probe_id"].tolist(),
                "group_breakdown": {"Body": 15},
                "methylation_effects": [],
                "methylation_condition_research": [],
                "evidence": [],
                "probe_preview": probe_preview,
            },
        ),
    )

    client = app.test_client()
    with client.session_transaction() as session_state:
        session_state["preprocess_state"] = {
            "gene_name": "TEST",
            "region": "1:1-100",
            "manifest_source": "",
            "filtered_manifest": "",
            "region_candidates": [],
            "selected_sources": [],
            "region_ready": True,
            "manifest_ready": True,
            "analysis_ready": True,
            "probe_count": 15,
            "build": "hg19",
            "logs": [],
            "region_recently_updated": False,
            "overwrite_filtered_manifest": False,
        }

    response = client.post(
        "/",
        data={
            "workflow": "analysis",
            "vcf": "data/mock.vcf.gz",
            "idat": "data/mock_sample",
            "out": "results/mock_report.html",
            "region": "1:1-100",
            "popstats": "",
            "manifest_file": "",
        },
    )

    assert response.status_code == 200
    result = captured_context["context"]["result"]
    probe_preview_html = result["methylation_insights"]["probe_preview"]
    assert probe_preview_html.count("<tr") >= 16
    assert "cg00000000" in probe_preview_html
    assert "cg00000014" in probe_preview_html


def test_analysis_result_labels_missing_variant_ids_in_preview(monkeypatch, tmp_path: Path) -> None:
    """Variant previews should explain when the source VCF does not provide a named ID."""
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_report_history", lambda: [])

    captured_context: dict[str, object] = {}

    def fake_render_template(template_name: str, **context: object) -> str:
        captured_context["template_name"] = template_name
        captured_context["context"] = context
        return "ok"

    monkeypatch.setattr("src.webapp.render_template", fake_render_template)
    monkeypatch.setattr(
        "src.webapp.run_analysis",
        lambda **_: SimpleNamespace(
            report_path=tmp_path / "mock_report.html",
            methylation_output_path=tmp_path / "mock_report_methylation.csv",
            variants=pd.DataFrame(
                [{"chrom": "11", "id": None, "pos": 636689, "ref": "G", "alt": "C", "qual": 75.35, "filter_pass": True}]
            ),
            methylation=pd.DataFrame([{"probe_id": "cg1", "beta": 0.42}]),
            popstats=None,
            knowledge_base={"database_name": "Mock interpretation DB", "version": "test"},
            population_database={"database_name": "Mock population DB", "version": "test"},
            population_insights={"variant_population_records": [], "gene_population_patterns": []},
            variant_interpretations={"gene_name": "TEST", "matched_records": [], "sample_highlights": {"summary": ""}},
            methylation_insights={
                "gene_name": "TEST",
                "clinical_context": "",
                "summary": "",
                "mean_beta": 0.42,
                "mean_beta_label": "Whitelist mean beta",
                "mean_beta_probe_count": 1,
                "whitelist_mean_beta": 0.42,
                "whitelist_mean_beta_label": "Whitelist mean beta",
                "whitelist_mean_beta_probe_count": 1,
                "whitelist_probe_count": 1,
                "whitelist_observed_probe_count": 1,
                "whitelist_explanation": "",
                "whitelist_literature_context": "",
                "whitelist_probe_statuses": [{"probe_id": "cg1", "observed_in_run": True}],
                "whitelist_probe_reference_rows": [],
                "whitelist_probe_reference_summary": "",
                "gene_name_mean_beta": 0.42,
                "gene_name_mean_beta_label": "TEST-named row mean beta",
                "gene_name_mean_beta_probe_count": 1,
                "gene_name_row_count": 1,
                "gene_name_match_columns": [],
                "gene_name_match_rule": "",
                "raw_mean_beta": 0.42,
                "raw_mean_beta_label": "All numeric-row mean beta",
                "raw_probe_count": 1,
                "raw_mean_beta_probe_count": 1,
                "all_numeric_mean_beta": 0.42,
                "all_numeric_mean_beta_label": "All numeric-row mean beta",
                "all_numeric_mean_beta_probe_count": 1,
                "beta_band": "intermediate",
                "beta_band_source_label": "Whitelist mean beta",
                "observed_probe_count": 1,
                "curated_probe_count": 1,
                "probe_ids": ["cg1"],
                "group_breakdown": {},
                "methylation_effects": [],
                "methylation_condition_research": [],
                "evidence": [],
                "probe_preview": pd.DataFrame([{"probe_id": "cg1", "beta": 0.42}]),
            },
        ),
    )

    client = app.test_client()
    with client.session_transaction() as session_state:
        session_state["preprocess_state"] = {
            "gene_name": "TEST",
            "region": "11:636000-637000",
            "manifest_source": "",
            "filtered_manifest": "",
            "region_candidates": [],
            "selected_sources": [],
            "region_ready": True,
            "manifest_ready": True,
            "analysis_ready": True,
            "probe_count": 1,
            "build": "hg19",
            "logs": [],
            "region_recently_updated": False,
            "overwrite_filtered_manifest": False,
        }

    response = client.post(
        "/",
        data={
            "workflow": "analysis",
            "vcf": "data/mock.vcf.gz",
            "idat": "data/mock_sample",
            "out": "results/mock_report.html",
            "region": "11:636000-637000",
            "popstats": "",
            "manifest_file": "",
        },
    )

    assert response.status_code == 200
    result = captured_context["context"]["result"]
    assert "Unlabeled in source VCF" in result["variant_preview"]
    assert result["variant_rows"][0]["id"] == "Unlabeled in source VCF"
    assert result["variant_raw_page_size"] == 25


def test_analysis_result_keeps_full_variant_rows_for_pagination(monkeypatch, tmp_path: Path) -> None:
    """The raw variant panel should keep all rows for client-side pagination."""
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_report_history", lambda: [])

    captured_context: dict[str, object] = {}

    def fake_render_template(template_name: str, **context: object) -> str:
        captured_context["template_name"] = template_name
        captured_context["context"] = context
        return "ok"

    variant_rows = pd.DataFrame(
        [
            {
                "chrom": "15",
                "id": f"rs{i:05d}",
                "pos": i,
                "ref": "A",
                "alt": "G",
                "qual": 99.0,
                "filter_pass": True,
            }
            for i in range(1, 31)
        ]
    )

    monkeypatch.setattr("src.webapp.render_template", fake_render_template)
    monkeypatch.setattr(
        "src.webapp.run_analysis",
        lambda **_: SimpleNamespace(
            report_path=tmp_path / "mock_report.html",
            methylation_output_path=tmp_path / "mock_report_methylation.csv",
            variants=variant_rows.copy(),
            methylation=pd.DataFrame([{"probe_id": "cg1", "beta": 0.42}]),
            popstats=None,
            knowledge_base={"database_name": "Mock interpretation DB", "version": "test"},
            population_database={"database_name": "Mock population DB", "version": "test"},
            population_insights={"variant_population_records": [], "gene_population_patterns": []},
            variant_interpretations={"gene_name": "TEST", "matched_records": [], "sample_highlights": {"summary": ""}},
            methylation_insights={
                "gene_name": "TEST",
                "clinical_context": "",
                "summary": "",
                "mean_beta": 0.42,
                "mean_beta_label": "Whitelist mean beta",
                "mean_beta_probe_count": 1,
                "whitelist_mean_beta": 0.42,
                "whitelist_mean_beta_label": "Whitelist mean beta",
                "whitelist_mean_beta_probe_count": 1,
                "whitelist_probe_count": 1,
                "whitelist_observed_probe_count": 1,
                "whitelist_explanation": "",
                "whitelist_literature_context": "",
                "whitelist_probe_statuses": [{"probe_id": "cg1", "observed_in_run": True}],
                "whitelist_probe_reference_rows": [],
                "whitelist_probe_reference_summary": "",
                "gene_name_mean_beta": 0.42,
                "gene_name_mean_beta_label": "TEST-named row mean beta",
                "gene_name_mean_beta_probe_count": 1,
                "gene_name_row_count": 1,
                "gene_name_match_columns": [],
                "gene_name_match_rule": "",
                "raw_mean_beta": 0.42,
                "raw_mean_beta_label": "All numeric-row mean beta",
                "raw_probe_count": 1,
                "raw_mean_beta_probe_count": 1,
                "all_numeric_mean_beta": 0.42,
                "all_numeric_mean_beta_label": "All numeric-row mean beta",
                "all_numeric_mean_beta_probe_count": 1,
                "beta_band": "intermediate",
                "beta_band_source_label": "Whitelist mean beta",
                "observed_probe_count": 1,
                "curated_probe_count": 1,
                "probe_ids": ["cg1"],
                "group_breakdown": {},
                "methylation_effects": [],
                "methylation_condition_research": [],
                "evidence": [],
                "probe_preview": pd.DataFrame([{"probe_id": "cg1", "beta": 0.42}]),
            },
        ),
    )

    client = app.test_client()
    with client.session_transaction() as session_state:
        session_state["preprocess_state"] = {
            "gene_name": "TEST",
            "region": "15:1-1000",
            "manifest_source": "",
            "filtered_manifest": "",
            "region_candidates": [],
            "selected_sources": [],
            "region_ready": True,
            "manifest_ready": True,
            "analysis_ready": True,
            "probe_count": 1,
            "build": "hg19",
            "logs": [],
            "region_recently_updated": False,
            "overwrite_filtered_manifest": False,
        }

    response = client.post(
        "/",
        data={
            "workflow": "analysis",
            "vcf": "data/mock.vcf.gz",
            "idat": "data/mock_sample",
            "out": "results/mock_report.html",
            "region": "15:1-1000",
            "popstats": "",
            "manifest_file": "",
        },
    )

    assert response.status_code == 200
    result = captured_context["context"]["result"]
    assert len(result["variant_rows"]) == 30
    assert "rs00001" in result["variant_preview"]
    assert "rs00025" in result["variant_preview"]
    assert "rs00026" not in result["variant_preview"]
