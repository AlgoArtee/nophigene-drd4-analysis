"""Fetch candidate DRD4 genomic intervals from multiple public annotation APIs."""

from typing import Iterable

import requests

DEFAULT_TIMEOUT_SECONDS = 30


def fetch_refseq_region(gene_symbol: str = "DRD4") -> str | None:
    """Look up a gene interval with the NCBI Entrez Gene APIs.

    The function first resolves the gene symbol to an Entrez Gene identifier and
    then fetches the gene summary payload to extract the genomic coordinates.
    Returning ``None`` instead of raising on lookup problems lets the caller
    combine several public sources opportunistically and decide later which
    interval to trust.

    Parameters
    ----------
    gene_symbol : str, optional
        HGNC-style gene symbol to search for. Defaults to ``"DRD4"``.

    Returns
    -------
    str | None
        Region string in ``chrom:start-end`` format, or ``None`` when the API
        does not return a usable interval.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    try:
        search_response = requests.get(
            base_url + "esearch.fcgi",
            params={
                "db": "gene",
                "term": f"{gene_symbol}[sym] AND Homo sapiens[orgn]",
                "retmode": "json",
            },
            timeout=DEFAULT_TIMEOUT_SECONDS,
        )
        search_response.raise_for_status()
        gene_ids = search_response.json().get("esearchresult", {}).get("idlist", [])
        if not gene_ids:
            return None

        gene_id = gene_ids[0]
        summary_response = requests.get(
            base_url + "esummary.fcgi",
            params={"db": "gene", "id": gene_id, "retmode": "json"},
            timeout=DEFAULT_TIMEOUT_SECONDS,
        )
        summary_response.raise_for_status()
        genomic_info = summary_response.json().get("result", {}).get(gene_id, {}).get("genomicinfo", [])
        if not genomic_info:
            return None

        location = genomic_info[0]
        chr_num = location.get("chr")
        start = location.get("chrstart")
        end = location.get("chrstop")
        if chr_num is None or start is None or end is None:
            return None
    except (requests.RequestException, ValueError):
        return None

    return f"{chr_num}:{min(start, end)}-{max(start, end)}"


def fetch_ensembl_region(
    gene_symbol: str = "DRD4",
    server: str = "https://grch37.rest.ensembl.org",
) -> str | None:
    """Resolve a gene interval from the Ensembl REST API.

    Parameters
    ----------
    gene_symbol : str, optional
        Gene symbol to resolve. Defaults to ``"DRD4"``.
    server : str, optional
        Base URL for the Ensembl REST server. The default points to the GRCh37
        archive because the rest of the repository is largely hg19/GRCh37-based.

    Returns
    -------
    str | None
        Region string in ``chrom:start-end`` format, or ``None`` when Ensembl
        does not return a successful lookup payload.
    """
    try:
        response = requests.get(
            f"{server}/lookup/symbol/homo_sapiens/{gene_symbol}",
            params={"expand": 1},
            headers={"Content-Type": "application/json"},
            timeout=DEFAULT_TIMEOUT_SECONDS,
        )
        response.raise_for_status()
        data = response.json()
    except (requests.RequestException, ValueError):
        return None

    seq_region_name = data.get("seq_region_name")
    start = data.get("start")
    end = data.get("end")
    if seq_region_name is None or start is None or end is None:
        return None
    return f"{seq_region_name}:{start}-{end}"


def fetch_ucsc_region(gene_symbol: str = "DRD4") -> str | None:
    """Look up a gene interval from the UCSC knownGene track.

    Notes
    -----
    This helper currently queries chromosome 11 directly because it was created
    for DRD4-focused exploratory work. It still filters by ``gene_symbol`` once
    the UCSC payload is downloaded, but it should be generalized before being
    reused for genes outside chromosome 11.

    Parameters
    ----------
    gene_symbol : str, optional
        Gene symbol to match against the UCSC ``name2`` field.

    Returns
    -------
    str | None
        Widest transcript interval returned by the UCSC payload, or ``None`` if
        no matching gene entries are found.
    """
    base_url = "https://api.genome.ucsc.edu/getData/track?genome=hg19;track=knownGene;chrom=chr11"

    try:
        response = requests.get(base_url, timeout=DEFAULT_TIMEOUT_SECONDS)
        response.raise_for_status()
        entries = response.json().get("knownGene", [])
    except (requests.RequestException, ValueError):
        return None

    regions: list[tuple[int, int]] = []
    for entry in entries:
        if entry.get("name2") == gene_symbol:
            regions.append((entry["txStart"], entry["txEnd"]))

    if not regions:
        return None

    start = min(start for start, _ in regions)
    end = max(end for _, end in regions)
    return f"11:{start}-{end}"


def get_widest_region(regions: Iterable[str]) -> str | None:
    """Return the widest interval from an iterable of ``chrom:start-end`` strings.

    Parameters
    ----------
    regions : Iterable[str]
        Candidate genomic intervals collected from one or more sources.

    Returns
    -------
    str | None
        The interval with the greatest span, or ``None`` when the iterable is
        empty.
    """
    max_span = -1
    widest = None

    for region in regions:
        _chrom, genomic_range = region.split(":")
        start, end = map(int, genomic_range.split("-"))
        span = end - start
        if span > max_span:
            max_span = span
            widest = region

    return widest


def main() -> None:
    """Fetch DRD4 intervals from all configured sources and print a summary."""
    refseq_region = fetch_refseq_region()
    ensembl_region = fetch_ensembl_region()
    ucsc_region = fetch_ucsc_region()

    all_regions = [region for region in [refseq_region, ensembl_region, ucsc_region] if region]
    widest = get_widest_region(all_regions)

    print("RefSeq region:", refseq_region)
    print("Ensembl region:", ensembl_region)
    print("UCSC region:", ucsc_region)
    print("Widest region:", widest)


if __name__ == "__main__":
    main()
