"""Helpers for browsing the human UniProt protein catalog from the web UI."""

from __future__ import annotations

import re
from typing import Any
from urllib.parse import unquote

import requests

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_PAGE_SIZE = 24
DEFAULT_TIMEOUT_SEC = 20
FEATURED_HUMAN_PROTEIN_QUERIES = ["DRD4", "TP53", "BRCA1", "EGFR", "APOE", "ACE2"]
HUMAN_PROTEOME_SCOPE = "organism_id:9606"
NEXT_CURSOR_PATTERN = re.compile(r"<[^>]*[?&]cursor=([^&>]+)[^>]*>;\s*rel=\"next\"")
REQUEST_HEADERS = {
    "User-Agent": "NophiGene/1.0 (+https://rest.uniprot.org)",
    "Accept": "application/json",
}
SEARCH_FIELDS = ",".join(
    [
        "accession",
        "id",
        "protein_name",
        "gene_names",
        "organism_name",
        "length",
        "annotation_score",
        "protein_existence",
        "cc_function",
        "cc_subcellular_location",
        "xref_geneid",
        "xref_alphafolddb",
    ]
)


class ProteinCatalogError(RuntimeError):
    """Raised when the live UniProt catalog cannot be reached or parsed."""


def build_human_protein_query(search_term: str = "", *, reviewed_only: bool = True) -> str:
    """Build the UniProt query used for the human protein catalog."""
    clauses = [HUMAN_PROTEOME_SCOPE]
    if reviewed_only:
        clauses.append("reviewed:true")

    normalized_term = " ".join(search_term.split())
    if normalized_term:
        protein_term = f"\"{normalized_term}\"" if " " in normalized_term else normalized_term
        clauses.append(
            f"(gene:{normalized_term} OR protein_name:{protein_term} OR id:{normalized_term} OR accession:{normalized_term})"
        )

    return " AND ".join(clauses)


def parse_next_cursor(link_header: str | None) -> str | None:
    """Extract the next-page cursor from a UniProt pagination header."""
    if not link_header:
        return None

    match = NEXT_CURSOR_PATTERN.search(link_header)
    if match is None:
        return None
    return unquote(match.group(1))


def _extract_primary_gene_name(record: dict[str, Any]) -> str:
    """Return the preferred gene symbol for a UniProt result record."""
    genes = record.get("genes") or []
    if not genes:
        return ""

    gene_name = genes[0].get("geneName") or {}
    return str(gene_name.get("value", "")).strip()


def _extract_gene_synonyms(record: dict[str, Any]) -> list[str]:
    """Return alternate gene symbols from a UniProt result record."""
    genes = record.get("genes") or []
    if not genes:
        return []

    synonyms = genes[0].get("synonyms") or []
    return [
        str(item.get("value", "")).strip()
        for item in synonyms
        if str(item.get("value", "")).strip()
    ]


def _extract_protein_name(description: dict[str, Any]) -> str:
    """Return the preferred protein name from the UniProt description block."""
    recommended_name = (description.get("recommendedName") or {}).get("fullName") or {}
    if recommended_name.get("value"):
        return str(recommended_name["value"]).strip()

    submission_names = description.get("submissionNames") or []
    if submission_names:
        submission_name = (submission_names[0].get("fullName") or {}).get("value")
        if submission_name:
            return str(submission_name).strip()

    return ""


def _extract_alternative_names(description: dict[str, Any]) -> list[str]:
    """Return alternative protein names from a UniProt description block."""
    alternative_names = description.get("alternativeNames") or []
    normalized: list[str] = []
    for item in alternative_names:
        value = str((item.get("fullName") or {}).get("value", "")).strip()
        if value:
            normalized.append(value)
    return normalized


def _extract_comment_text(record: dict[str, Any], *, comment_type: str) -> str:
    """Return the first free-text comment for a given UniProt comment type."""
    for comment in record.get("comments") or []:
        if comment.get("commentType") != comment_type:
            continue
        for text_block in comment.get("texts") or []:
            value = str(text_block.get("value", "")).strip()
            if value:
                return value
    return ""


def _extract_subcellular_locations(record: dict[str, Any]) -> list[str]:
    """Return flattened subcellular-location labels from a UniProt result record."""
    locations: list[str] = []
    for comment in record.get("comments") or []:
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for item in comment.get("subcellularLocations") or []:
            location = str((item.get("location") or {}).get("value", "")).strip()
            topology = str((item.get("topology") or {}).get("value", "")).strip()
            if location and topology:
                locations.append(f"{location}; {topology}")
            elif location:
                locations.append(location)
            elif topology:
                locations.append(topology)
    return locations


def _extract_cross_reference(record: dict[str, Any], *, database: str) -> str | None:
    """Return a cross-reference identifier for a requested linked database."""
    for item in record.get("uniProtKBCrossReferences") or []:
        if item.get("database") == database and item.get("id"):
            return str(item["id"]).strip()
    return None


def normalize_protein_record(record: dict[str, Any]) -> dict[str, Any]:
    """Convert a raw UniProt record into the UI-friendly protein card shape."""
    description = record.get("proteinDescription") or {}
    accession = str(record.get("primaryAccession", "")).strip()
    gene_name = _extract_primary_gene_name(record)
    gene_id = _extract_cross_reference(record, database="GeneID")
    alphafold_id = _extract_cross_reference(record, database="AlphaFoldDB")
    function_summary = _extract_comment_text(record, comment_type="FUNCTION")
    subcellular_locations = _extract_subcellular_locations(record)

    entry_type = str(record.get("entryType", "")).strip()
    entry_kind = "Reviewed" if "reviewed" in entry_type.lower() else "Unreviewed"

    return {
        "accession": accession,
        "entry_name": str(record.get("uniProtkbId", "")).strip(),
        "entry_type": entry_kind,
        "entry_type_label": entry_type or entry_kind,
        "gene_name": gene_name or accession,
        "gene_synonyms": _extract_gene_synonyms(record),
        "protein_name": _extract_protein_name(description) or accession,
        "alternative_names": _extract_alternative_names(description),
        "organism": str((record.get("organism") or {}).get("scientificName", "")).strip() or "Homo sapiens",
        "length": int((record.get("sequence") or {}).get("length") or 0),
        "annotation_score": record.get("annotationScore"),
        "protein_existence": str(record.get("proteinExistence", "")).strip() or "Not provided",
        "function_summary": function_summary or "UniProt did not return a function summary for this entry in the current response.",
        "subcellular_locations": subcellular_locations,
        "gene_id": gene_id,
        "alphafold_id": alphafold_id,
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry" if accession else "",
        "alphafold_url": f"https://alphafold.ebi.ac.uk/entry/{alphafold_id}" if alphafold_id else "",
        "ncbi_gene_url": f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}" if gene_id else "",
    }


def fetch_human_protein_catalog(
    *,
    query: str = "",
    reviewed_only: bool = True,
    cursor: str | None = None,
    page_size: int = DEFAULT_PAGE_SIZE,
) -> dict[str, Any]:
    """Fetch a page of the human UniProt protein catalog."""
    params = {
        "query": build_human_protein_query(query, reviewed_only=reviewed_only),
        "fields": SEARCH_FIELDS,
        "size": max(1, min(page_size, 50)),
        "format": "json",
    }
    if cursor:
        params["cursor"] = cursor

    try:
        response = requests.get(
            UNIPROT_SEARCH_URL,
            params=params,
            headers=REQUEST_HEADERS,
            timeout=DEFAULT_TIMEOUT_SEC,
        )
        response.raise_for_status()
    except requests.RequestException as exc:
        raise ProteinCatalogError(f"UniProt request failed: {exc}") from exc

    payload = response.json()
    proteins = [normalize_protein_record(record) for record in payload.get("results", [])]
    total_results = response.headers.get("X-Total-Results")
    next_cursor = parse_next_cursor(response.headers.get("Link"))

    return {
        "query": query,
        "reviewed_only": reviewed_only,
        "cursor": cursor,
        "next_cursor": next_cursor,
        "has_next_page": next_cursor is not None,
        "page_size": params["size"],
        "records_returned": len(proteins),
        "total_results": int(total_results) if total_results and str(total_results).isdigit() else None,
        "catalog_label": "UniProt human UniProtKB catalog",
        "catalog_scope": "Reviewed human proteins" if reviewed_only else "Reviewed and unreviewed human proteins",
        "catalog_url": "https://rest.uniprot.org/uniprotkb/search",
        "featured_queries": FEATURED_HUMAN_PROTEIN_QUERIES,
        "proteins": proteins,
        "error": None,
    }


def get_human_protein_catalog(
    *,
    query: str = "",
    reviewed_only: bool = True,
    cursor: str | None = None,
    page_size: int = DEFAULT_PAGE_SIZE,
) -> dict[str, Any]:
    """Return a UI-safe protein catalog payload, even when the live fetch fails."""
    try:
        return fetch_human_protein_catalog(
            query=query,
            reviewed_only=reviewed_only,
            cursor=cursor,
            page_size=page_size,
        )
    except ProteinCatalogError as exc:
        return {
            "query": query,
            "reviewed_only": reviewed_only,
            "cursor": cursor,
            "next_cursor": None,
            "has_next_page": False,
            "page_size": page_size,
            "records_returned": 0,
            "total_results": None,
            "catalog_label": "UniProt human UniProtKB catalog",
            "catalog_scope": "Reviewed human proteins" if reviewed_only else "Reviewed and unreviewed human proteins",
            "catalog_url": "https://rest.uniprot.org/uniprotkb/search",
            "featured_queries": FEATURED_HUMAN_PROTEIN_QUERIES,
            "proteins": [],
            "error": str(exc),
        }
