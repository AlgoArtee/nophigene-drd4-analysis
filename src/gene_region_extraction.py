import requests
from xml.etree import ElementTree

def fetch_refseq_region(gene_symbol="DRD4"):
    """
    Fetch DRD4 region using NCBI Entrez eUtils (RefSeq).
    Output format: 'chr:start-end'
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # 1. Get Gene ID from gene symbol
    query = f"{base_url}esearch.fcgi?db=gene&term={gene_symbol}[sym]+AND+Homo+sapiens[orgn]&retmode=json"
    gene_id = requests.get(query).json()['esearchresult']['idlist'][0]

    # 2. Get gene summary to extract location
    summary_url = f"{base_url}esummary.fcgi?db=gene&id={gene_id}&retmode=json"
    summary = requests.get(summary_url).json()
    location = summary['result'][gene_id]['genomicinfo'][0]
    chr_num = location['chr']
    start = location['chrstart']
    end = location['chrstop']
    return f"{chr_num}:{min(start, end)}-{max(start, end)}"


def fetch_ensembl_region(gene_symbol="DRD4", server="https://grch37.rest.ensembl.org"):
    """
    Fetch DRD4 region using Ensembl REST API (GRCh37 build).
    Output format: 'chr:start-end'
    """
    headers = {"Content-Type": "application/json"}
    ext = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    r = requests.get(server + ext, headers=headers)
    if not r.ok:
        return None
    data = r.json()
    return f"{data['seq_region_name']}:{data['start']}-{data['end']}"


def fetch_ucsc_region(gene_symbol="DRD4"):
    """
    Fetch DRD4 region using UCSC Table Browser JSON endpoint.
    Output format: 'chr:start-end'
    Note: Uses UCSC knownGene track and xref to gene symbol.
    """
    base_url = "https://api.genome.ucsc.edu/getData/track?genome=hg19;track=knownGene;chrom=chr11"
    r = requests.get(base_url)
    if not r.ok:
        return None

    entries = r.json().get('knownGene', [])
    regions = []
    for entry in entries:
        if 'name2' in entry and entry['name2'] == gene_symbol:
            regions.append((entry['txStart'], entry['txEnd']))
    if not regions:
        return None

    start = min(r[0] for r in regions)
    end = max(r[1] for r in regions)
    return f"11:{start}-{end}"


def get_widest_region(regions):
    """
    Determine the widest genomic region among all.
    """
    max_span = 0
    widest = None
    for r in regions:
        chrom, rng = r.split(":")
        start, end = map(int, rng.split("-"))
        span = end - start
        if span > max_span:
            max_span = span
            widest = r
    return widest


# Run fetch functions
refseq_region = fetch_refseq_region()
ensembl_region = fetch_ensembl_region()
ucsc_region = fetch_ucsc_region()

all_regions = [r for r in [refseq_region, ensembl_region, ucsc_region] if r]
widest = get_widest_region(all_regions)

# Output results
print("RefSeq region:", refseq_region)
print("Ensembl region:", ensembl_region)
print("UCSC region:", ucsc_region)
print("Widest region:", widest)
