"""Helpers for browsing the human UniProt protein catalog from the web UI."""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.parse import unquote

import requests

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_PAGE_SIZE = 24
DEFAULT_TIMEOUT_SEC = 20
FEATURED_HUMAN_PROTEIN_QUERIES = ["DRD4", "TP53", "BRCA1", "EGFR", "APOE", "ACE2"]
HUMAN_PROTEOME_SCOPE = "organism_id:9606"
NEXT_CURSOR_PATTERN = re.compile(r"<[^>]*[?&]cursor=([^&>]+)[^>]*>;\s*rel=\"next\"")
LONGEVITY_DB_PATH = Path(__file__).resolve().parent / "gene_data" / "longevity_genes_hagr.json"
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


@lru_cache(maxsize=1)
def load_longevity_gene_database(database_path: str | Path = LONGEVITY_DB_PATH) -> dict[str, Any]:
    """Load the curated local longevity-gene symbol set derived from HAGR."""
    payload_path = Path(database_path)
    if not payload_path.exists():
        raise ProteinCatalogError(f"Longevity gene database not found: {payload_path}")

    try:
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ProteinCatalogError(f"Longevity gene database is not valid JSON: {payload_path}") from exc

    gene_symbols = payload.get("gene_symbols")
    if not isinstance(gene_symbols, list) or not gene_symbols:
        raise ProteinCatalogError(f"Longevity gene database has no gene_symbols list: {payload_path}")
    return payload


def _get_longevity_gene_symbols() -> list[str]:
    """Return the ordered longevity-gene symbols used by the longevity filter."""
    payload = load_longevity_gene_database()
    return [str(symbol).strip() for symbol in payload.get("gene_symbols", []) if str(symbol).strip()]


def _get_longevity_gene_symbol_set() -> set[str]:
    """Return the longevity-gene symbols as an uppercase set for fast matching."""
    return {symbol.upper() for symbol in _get_longevity_gene_symbols()}


def parse_next_cursor(link_header: str | None) -> str | None:
    """Extract the next-page cursor from a UniProt pagination header."""
    if not link_header:
        return None

    match = NEXT_CURSOR_PATTERN.search(link_header)
    if match is None:
        return None
    return unquote(match.group(1))


def _quote_gene_symbol_for_uniprot(symbol: str) -> str:
    """Format a gene symbol safely for a UniProt gene query clause."""
    cleaned_symbol = symbol.replace('"', "").strip()
    if re.fullmatch(r"[A-Za-z0-9]+", cleaned_symbol):
        return f"gene:{cleaned_symbol}"
    return f'gene:"{cleaned_symbol}"'


def _build_gene_symbol_subset_query(
    gene_symbols: list[str],
    *,
    reviewed_only: bool,
) -> str:
    """Build a human UniProt query restricted to a given symbol subset."""
    clauses = [HUMAN_PROTEOME_SCOPE]
    if reviewed_only:
        clauses.append("reviewed:true")

    symbol_terms = [_quote_gene_symbol_for_uniprot(symbol) for symbol in gene_symbols if symbol]
    if symbol_terms:
        clauses.append(f"({' OR '.join(symbol_terms)})")
    return " AND ".join(clauses)


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


def _protein_matches_gene_symbol_set(protein: dict[str, Any], gene_symbols: set[str]) -> bool:
    """Return whether a normalized protein card belongs to the requested gene set."""
    candidates = {str(protein.get("gene_name", "")).upper()}
    candidates.update(str(item).upper() for item in protein.get("gene_synonyms", []))
    return any(candidate and candidate in gene_symbols for candidate in candidates)


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


def _execute_uniprot_search(
    *,
    query_string: str,
    cursor: str | None = None,
    page_size: int = DEFAULT_PAGE_SIZE,
) -> tuple[list[dict[str, Any]], requests.Response]:
    """Execute a UniProt search query and normalize the returned records."""
    params = {
        "query": query_string,
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
    return proteins, response


def _fetch_longevity_protein_catalog(
    *,
    query: str = "",
    reviewed_only: bool = True,
    page_index: int = 1,
    page_size: int = DEFAULT_PAGE_SIZE,
) -> dict[str, Any]:
    """Fetch the protein catalog restricted to HAGR longevity genes."""
    longevity_db = load_longevity_gene_database()
    all_symbols = _get_longevity_gene_symbols()
    normalized_query = query.strip().lower()

    filtered_symbols = [
        symbol for symbol in all_symbols if not normalized_query or normalized_query in symbol.lower()
    ]

    safe_page_index = max(1, page_index)
    start_index = (safe_page_index - 1) * page_size
    end_index = start_index + page_size
    page_symbols = filtered_symbols[start_index:end_index]

    proteins: list[dict[str, Any]] = []
    error_message = None
    if page_symbols:
        proteins, _ = _execute_uniprot_search(
            query_string=_build_gene_symbol_subset_query(page_symbols, reviewed_only=reviewed_only),
            page_size=page_size,
        )
        page_symbol_index = {symbol.upper(): index for index, symbol in enumerate(page_symbols)}
        proteins = [
            protein
            for protein in proteins
            if _protein_matches_gene_symbol_set(protein, {symbol.upper() for symbol in page_symbols})
        ]
        proteins.sort(
            key=lambda protein: page_symbol_index.get(str(protein.get("gene_name", "")).upper(), len(page_symbol_index))
        )
    elif normalized_query:
        error_message = (
            "No HAGR longevity genes matched that symbol filter. Try a direct gene symbol such as APOE, FOXO3, TERT, or DRD4."
        )

    return {
        "query": query,
        "reviewed_only": reviewed_only,
        "cursor": None,
        "next_cursor": None,
        "has_next_page": end_index < len(filtered_symbols),
        "has_previous_page": safe_page_index > 1,
        "page_index": safe_page_index,
        "next_page_index": safe_page_index + 1 if end_index < len(filtered_symbols) else None,
        "previous_page_index": safe_page_index - 1 if safe_page_index > 1 else None,
        "pagination_mode": "page",
        "page_size": page_size,
        "records_returned": len(proteins),
        "total_results": len(filtered_symbols),
        "catalog_label": "UniProt human proteins filtered by HAGR LongevityMap genes",
        "catalog_scope": "Reviewed longevity-associated human proteins" if reviewed_only else "Longevity-associated human proteins",
        "catalog_url": "https://rest.uniprot.org/uniprotkb/search",
        "featured_queries": FEATURED_HUMAN_PROTEIN_QUERIES,
        "longevity_only": True,
        "longevity_source": {
            "database_name": longevity_db.get("database_name"),
            "source": longevity_db.get("source"),
            "build": longevity_db.get("build"),
            "release_date": longevity_db.get("release_date"),
            "source_url": longevity_db.get("source_url"),
            "download_url": longevity_db.get("download_url"),
            "selection_rule": longevity_db.get("selection_rule"),
            "total_significant_gene_symbols": longevity_db.get("total_significant_gene_symbols"),
        },
        "proteins": proteins,
        "error": error_message,
    }


def fetch_human_protein_catalog(
    *,
    query: str = "",
    reviewed_only: bool = True,
    cursor: str | None = None,
    page_size: int = DEFAULT_PAGE_SIZE,
    longevity_only: bool = False,
    longevity_page: int = 1,
) -> dict[str, Any]:
    """Fetch a page of the human UniProt protein catalog."""
    if longevity_only:
        return _fetch_longevity_protein_catalog(
            query=query,
            reviewed_only=reviewed_only,
            page_index=longevity_page,
            page_size=page_size,
        )

    proteins, response = _execute_uniprot_search(
        query_string=build_human_protein_query(query, reviewed_only=reviewed_only),
        cursor=cursor,
        page_size=page_size,
    )
    total_results = response.headers.get("X-Total-Results")
    next_cursor = parse_next_cursor(response.headers.get("Link"))

    return {
        "query": query,
        "reviewed_only": reviewed_only,
        "cursor": cursor,
        "next_cursor": next_cursor,
        "has_next_page": next_cursor is not None,
        "has_previous_page": cursor is not None,
        "page_index": None,
        "next_page_index": None,
        "previous_page_index": None,
        "pagination_mode": "cursor",
        "page_size": max(1, min(page_size, 50)),
        "records_returned": len(proteins),
        "total_results": int(total_results) if total_results and str(total_results).isdigit() else None,
        "catalog_label": "UniProt human UniProtKB catalog",
        "catalog_scope": "Reviewed human proteins" if reviewed_only else "Reviewed and unreviewed human proteins",
        "catalog_url": "https://rest.uniprot.org/uniprotkb/search",
        "featured_queries": FEATURED_HUMAN_PROTEIN_QUERIES,
        "longevity_only": False,
        "longevity_source": None,
        "proteins": proteins,
        "error": None,
    }


def get_human_protein_catalog(
    *,
    query: str = "",
    reviewed_only: bool = True,
    cursor: str | None = None,
    page_size: int = DEFAULT_PAGE_SIZE,
    longevity_only: bool = False,
    longevity_page: int = 1,
) -> dict[str, Any]:
    """Return a UI-safe protein catalog payload, even when the live fetch fails."""
    longevity_source = None
    if longevity_only:
        try:
            longevity_source = load_longevity_gene_database()
        except ProteinCatalogError:
            longevity_source = None

    try:
        return fetch_human_protein_catalog(
            query=query,
            reviewed_only=reviewed_only,
            cursor=cursor,
            page_size=page_size,
            longevity_only=longevity_only,
            longevity_page=longevity_page,
        )
    except ProteinCatalogError as exc:
        return {
            "query": query,
            "reviewed_only": reviewed_only,
            "cursor": cursor,
            "next_cursor": None,
            "has_next_page": False,
            "has_previous_page": False,
            "page_index": longevity_page if longevity_only else None,
            "next_page_index": None,
            "previous_page_index": None,
            "pagination_mode": "page" if longevity_only else "cursor",
            "page_size": page_size,
            "records_returned": 0,
            "total_results": None,
            "catalog_label": (
                "UniProt human proteins filtered by HAGR LongevityMap genes"
                if longevity_only
                else "UniProt human UniProtKB catalog"
            ),
            "catalog_scope": (
                "Reviewed longevity-associated human proteins"
                if longevity_only and reviewed_only
                else "Longevity-associated human proteins"
                if longevity_only
                else "Reviewed human proteins"
                if reviewed_only
                else "Reviewed and unreviewed human proteins"
            ),
            "catalog_url": "https://rest.uniprot.org/uniprotkb/search",
            "featured_queries": FEATURED_HUMAN_PROTEIN_QUERIES,
            "longevity_only": longevity_only,
            "longevity_source": longevity_source,
            "proteins": [],
            "error": str(exc),
        }
