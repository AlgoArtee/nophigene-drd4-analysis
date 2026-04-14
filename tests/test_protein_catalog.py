"""Regression tests for the UniProt-backed human protein catalog helpers."""

from __future__ import annotations

from src.human_protein_catalog import (
    build_human_protein_query,
    normalize_protein_record,
    parse_next_cursor,
)
from src.webapp import app


def test_build_human_protein_query_supports_symbol_search() -> None:
    """Gene-symbol searches should stay scoped to human UniProt entries."""
    query = build_human_protein_query("DRD4", reviewed_only=True)

    assert "organism_id:9606" in query
    assert "reviewed:true" in query
    assert "gene:DRD4" in query


def test_parse_next_cursor_extracts_uniprot_cursor_value() -> None:
    """UniProt pagination headers should yield the next cursor token."""
    link_header = '<https://rest.uniprot.org/uniprotkb/search?query=organism_id%3A9606&cursor=eyJvZmZzZXQiOjI0fQ%3D%3D&size=24>; rel="next"'

    assert parse_next_cursor(link_header) == "eyJvZmZzZXQiOjI0fQ=="


def test_normalize_protein_record_extracts_drd4_details() -> None:
    """The DRD4 UniProt example should normalize into the UI card shape."""
    record = {
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "primaryAccession": "P21917",
        "uniProtkbId": "DRD4_HUMAN",
        "annotationScore": 5.0,
        "proteinExistence": "1: Evidence at protein level",
        "organism": {"scientificName": "Homo sapiens"},
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": "D(4) dopamine receptor"}},
            "alternativeNames": [
                {"fullName": {"value": "Dopamine D4 receptor"}},
                {"fullName": {"value": "D(2C) dopamine receptor"}},
            ],
        },
        "genes": [{"geneName": {"value": "DRD4"}}],
        "comments": [
            {
                "commentType": "FUNCTION",
                "texts": [
                    {
                        "value": "Dopamine receptor responsible for neuronal signaling in the mesolimbic system of the brain."
                    }
                ],
            },
            {
                "commentType": "SUBCELLULAR LOCATION",
                "subcellularLocations": [
                    {
                        "location": {"value": "Cell membrane"},
                        "topology": {"value": "Multi-pass membrane protein"},
                    }
                ],
            },
        ],
        "uniProtKBCrossReferences": [
            {"database": "AlphaFoldDB", "id": "P21917"},
            {"database": "GeneID", "id": "1815"},
        ],
        "sequence": {"length": 419},
    }

    normalized = normalize_protein_record(record)

    assert normalized["gene_name"] == "DRD4"
    assert normalized["protein_name"] == "D(4) dopamine receptor"
    assert normalized["length"] == 419
    assert normalized["gene_id"] == "1815"
    assert normalized["alphafold_id"] == "P21917"
    assert normalized["subcellular_locations"] == ["Cell membrane; Multi-pass membrane protein"]
    assert normalized["uniprot_url"].endswith("/P21917/entry")


def test_human_proteins_api_returns_json_payload(monkeypatch) -> None:
    """The UI endpoint should expose the protein catalog as JSON."""
    sample_payload = {
        "query": "DRD4",
        "reviewed_only": True,
        "cursor": None,
        "next_cursor": None,
        "has_next_page": False,
        "page_size": 24,
        "records_returned": 1,
        "total_results": 1,
        "catalog_label": "UniProt human UniProtKB catalog",
        "catalog_scope": "Reviewed human proteins",
        "catalog_url": "https://rest.uniprot.org/uniprotkb/search",
        "featured_queries": ["DRD4"],
        "proteins": [
            {
                "gene_name": "DRD4",
                "protein_name": "D(4) dopamine receptor",
                "accession": "P21917",
                "entry_name": "DRD4_HUMAN",
                "entry_type": "Reviewed",
                "entry_type_label": "UniProtKB reviewed (Swiss-Prot)",
                "gene_synonyms": [],
                "alternative_names": ["Dopamine D4 receptor"],
                "organism": "Homo sapiens",
                "length": 419,
                "annotation_score": 5.0,
                "protein_existence": "1: Evidence at protein level",
                "function_summary": "Example function",
                "subcellular_locations": ["Cell membrane"],
                "gene_id": "1815",
                "alphafold_id": "P21917",
                "uniprot_url": "https://www.uniprot.org/uniprotkb/P21917/entry",
                "alphafold_url": "https://alphafold.ebi.ac.uk/entry/P21917",
                "ncbi_gene_url": "https://www.ncbi.nlm.nih.gov/gene/1815",
            }
        ],
        "error": None,
    }

    monkeypatch.setattr("src.webapp.get_human_protein_catalog", lambda **_: sample_payload)

    client = app.test_client()
    response = client.get("/api/human-proteins?q=DRD4")

    assert response.status_code == 200
    payload = response.get_json()
    assert payload["query"] == "DRD4"
    assert payload["proteins"][0]["gene_name"] == "DRD4"
