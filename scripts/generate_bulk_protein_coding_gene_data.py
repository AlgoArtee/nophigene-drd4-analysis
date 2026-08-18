#!/usr/bin/env python3
"""Generate sharded curated-lite knowledge bases for HGNC protein-coding genes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


GENE_DATA_DIR = PROJECT_ROOT / "src" / "gene_data"
REFERENCE_CACHE_DIR = PROJECT_ROOT / "data" / "reference_cache"
HGNC_COMPLETE_SET_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
GENCODE_V19_GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
HGNC_COMPLETE_SET_PATH = REFERENCE_CACHE_DIR / "hgnc_complete_set.txt"
GENCODE_V19_GTF_PATH = REFERENCE_CACHE_DIR / "gencode.v19.annotation.gtf.gz"
MANIFEST_PATH = PROJECT_ROOT / "data" / "infinium-methylationepic-v-1-0-b5-manifest-file.csv"
CURATED_BUNDLE_PATH = GENE_DATA_DIR / "gene_data_bundle.zip"
BULK_SHARD_DIR = GENE_DATA_DIR / "bulk_gene_data_shards"
INDEX_PATH = GENE_DATA_DIR / "gene_data_index.json"
FORMAT_VERSION = "bulk-protein-coding-v2"
DATA_VERSION = "2026-08-17"
SHARD_COUNT = 64
ZIP_TIMESTAMP = (2026, 1, 1, 0, 0, 0)
GENERATED_SUFFIXES = (
    "_interpretation_db.json",
    "_population_db.json",
    "_epigenetics_hg19.csv",
)
MANIFEST_COLUMNS = [
    "IlmnID",
    "CHR",
    "MAPINFO",
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Accession",
    "UCSC_RefGene_Group",
    "UCSC_CpG_Islands_Name",
    "Relation_to_UCSC_CpG_Island",
    "GencodeBasicV12_NAME",
    "GencodeBasicV12_Accession",
    "GencodeBasicV12_Group",
    "SNP_ID",
    "SNP_DISTANCE",
    "CHR_hg38",
    "Start_hg38",
    "End_hg38",
]
def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _ensure_reference_file(path: Path, url: str, *, allow_download: bool) -> None:
    if path.exists():
        return
    if not allow_download:
        raise FileNotFoundError(f"Missing reference file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as response, path.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def _split_tokens(value: str | None) -> list[str]:
    if not value:
        return []
    return [token.strip() for token in re.split(r"[|;]", value) if token.strip()]


def sanitize_gene_name_for_filename(gene_name: str) -> str:
    """Convert a gene symbol into the filename-safe stem used by the app."""
    cleaned = gene_name.strip()
    if not cleaned:
        raise ValueError("Gene symbol is required.")
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", cleaned)
    return sanitized.strip("_") or "gene"


def _parse_gtf_attributes(attribute_text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]*)"', attribute_text))


def _load_hgnc_protein_coding(path: Path) -> dict[str, dict[str, str]]:
    genes: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") != "Approved":
                continue
            if row.get("locus_group") != "protein-coding gene":
                continue
            symbol = row.get("symbol", "").strip()
            if symbol:
                genes[symbol.upper()] = row
    return genes


def _prefer_coordinate_record(existing: dict[str, Any] | None, candidate: dict[str, Any]) -> dict[str, Any]:
    if existing is None:
        return candidate
    if existing.get("gene_type") != "protein_coding" and candidate.get("gene_type") == "protein_coding":
        return candidate
    return existing


def _load_gencode_gene_coordinates(path: Path) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    by_ensembl: dict[str, dict[str, Any]] = {}
    by_symbol: dict[str, dict[str, Any]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id", "").split(".")[0]
            gene_name = attrs.get("gene_name", "")
            record = {
                "chromosome": parts[0].removeprefix("chr"),
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "ensembl_gene_id": gene_id,
                "gencode_gene_name": gene_name,
                "gencode_gene_type": attrs.get("gene_type", ""),
                "gencode_gene_status": attrs.get("gene_status", ""),
            }
            if gene_id:
                by_ensembl[gene_id] = _prefer_coordinate_record(by_ensembl.get(gene_id), record)
            if gene_name:
                key = gene_name.upper()
                by_symbol[key] = _prefer_coordinate_record(by_symbol.get(key), record)
    return by_ensembl, by_symbol


def _load_curated_gene_symbols(bundle_path: Path) -> set[str]:
    if not bundle_path.exists():
        return set()
    curated: set[str] = set()
    with zipfile.ZipFile(bundle_path) as bundle:
        for name in bundle.namelist():
            if not name.endswith("_interpretation_db.json"):
                continue
            payload = json.loads(bundle.read(name).decode("utf-8"))
            symbol = str(payload.get("gene_context", {}).get("gene_name", "")).strip()
            if symbol:
                curated.add(symbol.upper())
    return curated


def _coordinate_for_gene(
    symbol: str,
    hgnc_record: dict[str, str],
    by_ensembl: dict[str, dict[str, Any]],
    by_symbol: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    for ensembl_id in _split_tokens(hgnc_record.get("ensembl_gene_id")):
        coordinate = by_ensembl.get(ensembl_id.split(".")[0])
        if coordinate:
            return coordinate
    return by_symbol.get(symbol.upper())


def _csv_line(values: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(values)
    return buffer.getvalue()


def _manifest_row_priority(row: dict[str, str], gene_symbol: str) -> tuple[int, int, int, str]:
    groups = f"{row.get('UCSC_RefGene_Group', '')};{row.get('GencodeBasicV12_Group', '')}".lower()
    relation = row.get("Relation_to_UCSC_CpG_Island", "").lower()
    if "tss200" in groups:
        group_priority = 0
    elif "tss1500" in groups:
        group_priority = 1
    elif "1stexon" in groups:
        group_priority = 2
    elif "5'utr" in groups:
        group_priority = 3
    elif "island" in relation:
        group_priority = 4
    elif "body" in groups:
        group_priority = 5
    else:
        group_priority = 6
    try:
        position = int(float(row.get("MAPINFO", "") or "0"))
    except ValueError:
        position = 0
    return group_priority, position, 0 if gene_symbol in row.get("UCSC_RefGene_Name", "").split(";") else 1, row.get("IlmnID", "")


def _load_manifest_rows_by_gene(
    manifest_path: Path,
    target_symbols: set[str],
) -> tuple[dict[str, list[tuple[tuple[int, int, int, str], str, str]]], list[str]]:
    by_gene: dict[str, list[tuple[tuple[int, int, int, str], str, str]]] = defaultdict(list)
    with manifest_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        header: list[str] | None = None
        for row in reader:
            if row and row[0] == "IlmnID":
                header = row
                break
        if header is None:
            raise ValueError(f"Could not locate IlmnID header in {manifest_path}")

        available_columns = [column for column in MANIFEST_COLUMNS if column in header]
        index = {column: header.index(column) for column in available_columns}
        gene_columns = [column for column in ("UCSC_RefGene_Name", "GencodeBasicV12_NAME") if column in index]
        for row in reader:
            if not row:
                continue
            row_values = {column: row[index[column]] if index[column] < len(row) else "" for column in available_columns}
            symbols_for_row: set[str] = set()
            for column in gene_columns:
                for token in _split_tokens(row_values.get(column)):
                    token_key = token.upper()
                    if token_key in target_symbols:
                        symbols_for_row.add(token_key)
            if not symbols_for_row:
                continue
            line = _csv_line([row_values.get(column, "") for column in available_columns])
            probe_id = row_values.get("IlmnID", "")
            for symbol in symbols_for_row:
                by_gene[symbol].append((_manifest_row_priority(row_values, symbol), probe_id, line))
    return by_gene, available_columns


def _promoter_region(start: int, end: int, strand: str) -> tuple[int, int]:
    if strand == "-":
        return end + 1, end + 1000
    return max(1, start - 1000), max(1, start - 1)


def _region_text(chromosome: str, start: int, end: int) -> str:
    return f"{chromosome}:{min(start, end)}-{max(start, end)}"


def _source_links(hgnc: dict[str, str], coordinate: dict[str, Any]) -> list[dict[str, str]]:
    links = [
        {
            "label": f"HGNC {hgnc.get('hgnc_id', '')}: {hgnc.get('symbol', '')}",
            "url": f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{hgnc.get('hgnc_id', '')}",
        }
    ]
    ensembl_id = coordinate.get("ensembl_gene_id") or hgnc.get("ensembl_gene_id", "")
    if ensembl_id:
        links.append(
            {
                "label": f"Ensembl GRCh37 {ensembl_id}",
                "url": f"https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl_id}",
            }
        )
    if hgnc.get("entrez_id"):
        links.append(
            {
                "label": f"NCBI Gene {hgnc['entrez_id']}",
                "url": f"https://www.ncbi.nlm.nih.gov/gene/{hgnc['entrez_id']}",
            }
        )
    return links


def _build_interpretation_database(
    symbol: str,
    hgnc: dict[str, str],
    coordinate: dict[str, Any],
    probe_ids: list[str],
    probe_count: int,
) -> dict[str, Any]:
    start = int(coordinate["start"])
    end = int(coordinate["end"])
    promoter_start, promoter_end = _promoter_region(start, end, str(coordinate.get("strand", "+")))
    combined_start = min(start, end, promoter_start, promoter_end)
    combined_end = max(start, end, promoter_start, promoter_end)
    name = hgnc.get("name", symbol)
    return {
        "database_name": f"NophiGene {symbol} Curated-Lite Interpretation Database",
        "version": DATA_VERSION,
        "curation_level": "curated-lite",
        "gene_context": {
            "gene_name": symbol,
            "hgnc_id": hgnc.get("hgnc_id", ""),
            "approved_name": name,
            "assembly": "GRCh37 / hg19",
            "cytoband": hgnc.get("location", ""),
            "chromosome": coordinate["chromosome"],
            "strand": coordinate.get("strand", ""),
            "gene_region": {
                "label": f"{symbol} transcribed interval",
                "start": start,
                "end": end,
                "definition": "GENCODE v19 GRCh37 gene interval joined to the HGNC-approved protein-coding symbol.",
            },
            "promoter_review_region": {
                "label": "Operational promoter review window",
                "start": promoter_start,
                "end": promoter_end,
                "definition": "A strand-aware 1 kb upstream window used by the app for promoter-adjacent review.",
            },
            "recommended_promoter_plus_gene_region": _region_text(coordinate["chromosome"], combined_start, combined_end),
            "gene_summary": (
                f"{symbol} ({name}) is an HGNC-approved protein-coding gene. "
                "This curated-lite record provides coordinate, nomenclature, and methylation-manifest context for exploratory NophiGene analysis."
            ),
            "clinical_context": (
                "This automatically generated record is metadata-backed and intended for research triage. "
                "It does not assert pathogenicity, diagnosis, treatment response, or curated variant-level clinical interpretation."
            ),
            "variant_effect_overview": [
                f"Variants in the {symbol} interval should be interpreted through decoded genotype dosage, call quality, and external clinical or literature evidence.",
                "No hand-curated high-impact variant assertions are bundled in this curated-lite layer.",
                "Existing richer curated NophiGene bundles take precedence over this bulk fallback whenever both exist.",
            ],
            "condition_research_overview": [
                "HGNC, Ensembl/GENCODE, and NCBI identifiers provide the starting point for gene-specific literature and database review.",
                "Use matched variants as locus context unless a separate curated source supports a stronger interpretation.",
            ],
            "relevant_methylation_probe_ids": probe_ids,
            "manifest_probe_count": probe_count,
            "methylation_interpretation": (
                f"{symbol} methylation is represented by compact EPIC manifest rows annotated to the gene when available. "
                "Treat beta values as regulatory context rather than as a standalone biomarker."
            ),
            "methylation_effects": [
                "Promoter-proximal probes can suggest local regulatory state, but tissue, cell mixture, and assay context matter.",
                "Gene-body probes are descriptive methylation context and should not be over-interpreted without study-specific evidence.",
            ],
            "methylation_condition_research": [
                "Use the retained manifest probe IDs to connect sample beta values with gene-region regulatory context.",
            ],
            "source_links": _source_links(hgnc, coordinate),
            "bulk_generated": True,
            "bulk_format_version": FORMAT_VERSION,
        },
        "variant_records": [],
    }


def _build_population_database(symbol: str, hgnc: dict[str, str], coordinate: dict[str, Any]) -> dict[str, Any]:
    return {
        "database_name": f"NophiGene {symbol} Curated-Lite Population Database",
        "version": DATA_VERSION,
        "curation_level": "curated-lite",
        "gene_name": symbol,
        "variant_population_records": [],
        "gene_population_patterns_intro": (
            f"Curated-lite population context for {symbol} is limited to stable gene identifiers and coordinate metadata."
        ),
        "gene_population_patterns": [
            {
                "category": "Reference metadata",
                "summary": (
                    f"{symbol} is represented as HGNC {hgnc.get('hgnc_id', 'unknown')} at "
                    f"GRCh37 {coordinate['chromosome']}:{coordinate['start']}-{coordinate['end']}."
                ),
                "evidence": _source_links(hgnc, coordinate),
            }
        ],
    }


def _json_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode("utf-8")


def _csv_bytes(rows: list[tuple[tuple[int, int, int, str], str, str]], columns: list[str]) -> bytes:
    lines = [_csv_line(columns)]
    lines.extend(line for _, _, line in sorted(rows, key=lambda item: item[0]))
    return "".join(lines).encode("utf-8")


def _shard_id_for_gene(symbol: str, shard_count: int = SHARD_COUNT) -> int:
    digest = hashlib.sha256(symbol.upper().encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % shard_count


def _zip_write_member(bundle: zipfile.ZipFile, filename: str, payload: bytes) -> None:
    info = zipfile.ZipInfo(filename, ZIP_TIMESTAMP)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    bundle.writestr(info, payload)


def build_bulk_gene_data(
    *,
    hgnc_path: Path = HGNC_COMPLETE_SET_PATH,
    gencode_path: Path = GENCODE_V19_GTF_PATH,
    manifest_path: Path = MANIFEST_PATH,
    curated_bundle_path: Path = CURATED_BUNDLE_PATH,
    output_dir: Path = GENE_DATA_DIR,
    shard_count: int = SHARD_COUNT,
    limit: int | None = None,
    allow_download: bool = False,
) -> dict[str, Any]:
    _ensure_reference_file(hgnc_path, HGNC_COMPLETE_SET_URL, allow_download=allow_download)
    _ensure_reference_file(gencode_path, GENCODE_V19_GTF_URL, allow_download=allow_download)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing EPIC manifest file: {manifest_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    shard_dir = output_dir / "bulk_gene_data_shards"
    tmp_shard_dir = output_dir / "bulk_gene_data_shards.tmp"
    if tmp_shard_dir.exists():
        shutil.rmtree(tmp_shard_dir)
    tmp_shard_dir.mkdir(parents=True, exist_ok=True)

    hgnc_genes = _load_hgnc_protein_coding(hgnc_path)
    by_ensembl, by_symbol = _load_gencode_gene_coordinates(gencode_path)
    curated_symbols = _load_curated_gene_symbols(curated_bundle_path)

    selected: dict[str, tuple[dict[str, str], dict[str, Any]]] = {}
    skipped_missing_coordinates: list[str] = []
    for symbol in sorted(hgnc_genes):
        if symbol in curated_symbols:
            continue
        coordinate = _coordinate_for_gene(symbol, hgnc_genes[symbol], by_ensembl, by_symbol)
        if coordinate is None:
            skipped_missing_coordinates.append(symbol)
            continue
        selected[symbol] = (hgnc_genes[symbol], coordinate)
        if limit is not None and len(selected) >= limit:
            break

    manifest_rows_by_gene, manifest_columns = _load_manifest_rows_by_gene(manifest_path, set(selected))
    shard_payloads: dict[int, dict[str, bytes]] = defaultdict(dict)
    index_genes: dict[str, Any] = {}
    index_files: dict[str, str] = {}

    for symbol, (hgnc, coordinate) in selected.items():
        rows = manifest_rows_by_gene.get(symbol, [])
        probe_ids = [probe_id for _, probe_id, _ in sorted(rows, key=lambda item: item[0])[:10] if probe_id]
        file_stem = sanitize_gene_name_for_filename(symbol)
        filenames = {
            "interpretation": f"{file_stem.lower()}_interpretation_db.json",
            "population": f"{file_stem.lower()}_population_db.json",
            "epigenetics": f"{file_stem}_epigenetics_hg19.csv",
        }
        payloads = {
            filenames["interpretation"]: _json_bytes(
                _build_interpretation_database(symbol, hgnc, coordinate, probe_ids, len(rows))
            ),
            filenames["population"]: _json_bytes(_build_population_database(symbol, hgnc, coordinate)),
            filenames["epigenetics"]: _csv_bytes(rows, manifest_columns),
        }
        shard_id = _shard_id_for_gene(symbol, shard_count)
        shard_name = f"gene_data_bulk_shard_{shard_id:02d}.zip"
        for filename, payload in payloads.items():
            shard_payloads[shard_id][filename] = payload
            index_files[filename] = shard_name
        index_genes[symbol] = {
            "symbol": symbol,
            "hgnc_id": hgnc.get("hgnc_id", ""),
            "name": hgnc.get("name", ""),
            "ensembl_gene_id": coordinate.get("ensembl_gene_id", hgnc.get("ensembl_gene_id", "")),
            "entrez_id": hgnc.get("entrez_id", ""),
            "chromosome": coordinate["chromosome"],
            "start": coordinate["start"],
            "end": coordinate["end"],
            "strand": coordinate.get("strand", ""),
            "manifest_probe_count": len(rows),
            "shard": shard_name,
            "files": filenames,
        }

    shard_index: dict[str, Any] = {}
    for shard_id in range(shard_count):
        shard_name = f"gene_data_bulk_shard_{shard_id:02d}.zip"
        shard_path = tmp_shard_dir / shard_name
        payloads = shard_payloads.get(shard_id, {})
        with zipfile.ZipFile(shard_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
            for filename in sorted(payloads):
                _zip_write_member(bundle, filename, payloads[filename])
        shard_index[shard_name] = {
            "path": f"bulk_gene_data_shards/{shard_name}",
            "gene_count": sum(1 for gene in index_genes.values() if gene["shard"] == shard_name),
            "file_count": len(payloads),
            "size_bytes": shard_path.stat().st_size,
            "sha256": _sha256_file(shard_path),
        }

    if shard_dir.exists():
        shutil.rmtree(shard_dir)
    tmp_shard_dir.replace(shard_dir)

    index = {
        "format_version": FORMAT_VERSION,
        "data_version": DATA_VERSION,
        "source_snapshot": {
            "hgnc_complete_set": {
                "url": HGNC_COMPLETE_SET_URL,
                "path": str(hgnc_path.relative_to(PROJECT_ROOT)) if hgnc_path.is_relative_to(PROJECT_ROOT) else str(hgnc_path),
                "sha256": _sha256_file(hgnc_path),
                "size_bytes": hgnc_path.stat().st_size,
            },
            "gencode_v19_gtf": {
                "url": GENCODE_V19_GTF_URL,
                "path": str(gencode_path.relative_to(PROJECT_ROOT)) if gencode_path.is_relative_to(PROJECT_ROOT) else str(gencode_path),
                "sha256": _sha256_file(gencode_path),
                "size_bytes": gencode_path.stat().st_size,
            },
            "manifest": {
                "path": str(manifest_path.relative_to(PROJECT_ROOT)) if manifest_path.is_relative_to(PROJECT_ROOT) else str(manifest_path),
                "sha256": _sha256_file(manifest_path),
                "size_bytes": manifest_path.stat().st_size,
            },
        },
        "shard_count": shard_count,
        "hgnc_protein_coding_count": len(hgnc_genes),
        "curated_excluded_count": len(curated_symbols & set(hgnc_genes)),
        "bulk_gene_count": len(index_genes),
        "bulk_file_count": len(index_files),
        "skipped_missing_grch37_coordinate_count": len(skipped_missing_coordinates),
        "skipped_missing_grch37_coordinate_symbols": skipped_missing_coordinates,
        "manifest_columns": manifest_columns,
        "genes": dict(sorted(index_genes.items())),
        "files": dict(sorted(index_files.items())),
        "shards": dict(sorted(shard_index.items())),
    }
    index_path = output_dir / "gene_data_index.json"
    index_path.write_text(json.dumps(index, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n", encoding="utf-8")
    return index


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate sharded curated-lite protein-coding gene data.")
    parser.add_argument("--download", action="store_true", help="Download missing HGNC/GENCODE reference files.")
    parser.add_argument("--limit", type=int, default=None, help="Generate only the first N genes for smoke testing.")
    args = parser.parse_args()

    index = build_bulk_gene_data(allow_download=args.download, limit=args.limit)
    compressed_size = sum(shard["size_bytes"] for shard in index["shards"].values())
    print(
        "Generated {genes:,} bulk protein-coding gene bundle(s), {files:,} files, {size:,} compressed bytes across {shards} shards.".format(
            genes=index["bulk_gene_count"],
            files=index["bulk_file_count"],
            size=compressed_size,
            shards=index["shard_count"],
        )
    )
    if index["skipped_missing_grch37_coordinate_count"]:
        print(
            "Skipped {count:,} HGNC protein-coding symbols without a GENCODE v19 GRCh37 gene coordinate.".format(
                count=index["skipped_missing_grch37_coordinate_count"]
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
