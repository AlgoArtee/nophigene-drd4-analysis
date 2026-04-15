"""Regression tests for preprocessing form behavior in the Flask UI."""

from __future__ import annotations

from pathlib import Path

from src.webapp import app


def test_preprocessing_template_preserves_clicked_submit_button() -> None:
    """The loading helper should preserve the clicked preprocessing action."""
    template_path = Path(__file__).resolve().parent.parent / "src" / "templates" / "index.html"
    template_text = template_path.read_text(encoding="utf-8")

    assert "const submitter = event.submitter instanceof HTMLElement ? event.submitter : null;" in template_text
    assert "formElement.requestSubmit(submitter);" in template_text
    assert "temporarySubmitterInput.name = submitter.name;" in template_text
    assert "window.setTimeout(() => formElement.submit(), 40);" not in template_text


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
