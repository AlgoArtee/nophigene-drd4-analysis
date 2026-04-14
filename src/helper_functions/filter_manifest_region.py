"""Filter Illumina methylation manifest rows down to a target genomic interval."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "src" / "gene_data"


def load_manifest(filepath: str) -> pd.DataFrame:
    """Load an Illumina manifest while skipping the vendor metadata preamble.

    Illumina manifest files often begin with a variable-length block of metadata
    before the real CSV header. This helper scans for the row starting with
    ``IlmnID`` and then hands the file off to ``pandas.read_csv`` from that
    point forward. The scan works for both plain-text and ``.gz`` inputs so the
    rest of the workflow can treat manifests uniformly.

    Parameters
    ----------
    filepath : str
        Path to a manifest file in CSV or gzipped CSV form.

    Returns
    -------
    pd.DataFrame
        Parsed manifest with whitespace-trimmed headers and numeric coordinate
        columns coerced to numbers where available.

    Raises
    ------
    ValueError
        Raised when the file does not contain a detectable ``IlmnID`` header
        row.
    """
    compression = "gzip" if filepath.endswith(".gz") else None
    open_fn = gzip.open if compression == "gzip" else open
    skiprows = None

    # Find the first real header row instead of assuming a fixed number of
    # vendor metadata lines at the top of the file.
    with open_fn(filepath, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("IlmnID"):
                skiprows = line_number
                break

    if skiprows is None:
        raise ValueError(f"Could not locate an 'IlmnID' header row in {filepath}")

    df = pd.read_csv(filepath, skiprows=skiprows, compression=compression, low_memory=False)

    # Normalize column names first so follow-up selection logic does not have to
    # worry about stray whitespace in vendor-provided headers.
    df.columns = [column.strip() for column in df.columns]

    # Coordinate columns are stored as strings in some manifest versions. Coerce
    # them once so downstream region comparisons can stay numeric.
    for column in ("MAPINFO", "Start_hg38", "End_hg38"):
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    return df


def filter_probes_by_region(
    df: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    genome_build: str = "hg19",
) -> pd.DataFrame:
    """Filter manifest probes to a chromosome interval for a chosen genome build.

    Parameters
    ----------
    df : pd.DataFrame
        Full Illumina manifest DataFrame produced by :func:`load_manifest`.
    chrom : str
        Chromosome identifier, such as ``"11"`` or ``"chr11"``.
    start : int
        Inclusive interval start coordinate.
    end : int
        Inclusive interval end coordinate.
    genome_build : str, optional
        Either ``"hg19"`` or ``"hg38"``. The helper maps this value to the
        correct chromosome and coordinate columns before filtering.

    Returns
    -------
    pd.DataFrame
        Copy of the manifest rows that fall inside the requested interval.

    Raises
    ------
    ValueError
        Raised when an unsupported genome build is requested.
    """
    if genome_build == "hg19":
        chrom_col = "CHR"
        pos_col = "MAPINFO"
        chrom = chrom.lstrip("chr")
    elif genome_build == "hg38":
        chrom_col = "CHR_hg38"
        pos_col = "Start_hg38"
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
    else:
        raise ValueError("Unsupported genome build. Use 'hg19' or 'hg38'.")

    # Cast chromosome labels to strings because manifests sometimes mix numeric
    # and string chromosome representations within the same column.
    filtered = df[
        (df[chrom_col].astype(str) == chrom)
        & (df[pos_col] >= start)
        & (df[pos_col] <= end)
    ].copy()

    return filtered


def parse_region_string(region: str) -> tuple[str, int, int]:
    """Parse a ``chr:start-end`` style region string into its coordinate parts."""
    cleaned_region = region.strip().replace(",", "")
    match = re.fullmatch(r"(?:chr)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)", cleaned_region)
    if match is None:
        raise ValueError(
            f"Unsupported region format '{region}'. Use chr:start-end, for example 11:637269-640706."
        )

    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        raise ValueError(f"Invalid region '{region}': start must be <= end.")

    return match.group("chrom"), start, end


def sanitize_gene_name_for_filename(gene_name: str) -> str:
    """Convert a gene symbol into a filesystem-safe filename stem."""
    cleaned = gene_name.strip()
    if not cleaned:
        raise ValueError("Enter a gene symbol before saving a filtered manifest.")
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", cleaned)
    return sanitized.strip("_") or "gene"


def save_filtered_manifest(
    *,
    gene_name: str,
    manifest_path: str,
    region: str,
    genome_build: str = "hg19",
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, object]:
    """Filter a manifest to one gene interval and save the resulting CSV subset.

    Parameters
    ----------
    gene_name : str
        Gene symbol used to name the output subset file.
    manifest_path : str
        Filesystem path to the source Illumina manifest.
    region : str
        Gene interval in ``chr:start-end`` form.
    genome_build : str, optional
        Genome build used for manifest coordinate matching.
    output_dir : str | Path, optional
        Destination directory for the filtered CSV subset.

    Returns
    -------
    dict[str, object]
        Summary containing the saved path, the filtered probe count, and a
        preview DataFrame for UI display.
    """
    manifest_source = Path(manifest_path)
    if not manifest_source.exists():
        raise ValueError(f"Manifest file not found: {manifest_source}")

    chrom, start, end = parse_region_string(region)
    manifest_df = load_manifest(str(manifest_source))
    filtered = filter_probes_by_region(manifest_df, chrom, start, end, genome_build)
    if filtered.empty:
        raise ValueError(
            f"No manifest probes were found in {chrom}:{start}-{end} for genome build {genome_build}."
        )

    destination_dir = Path(output_dir)
    destination_dir.mkdir(parents=True, exist_ok=True)
    output_path = destination_dir / (
        f"{sanitize_gene_name_for_filename(gene_name)}_epigenetics_{genome_build}.csv"
    )
    filtered.to_csv(output_path, index=False)

    preview_columns = [
        "IlmnID",
        "CHR",
        "MAPINFO",
        "CHR_hg38",
        "Start_hg38",
        "End_hg38",
        "UCSC_RefGene_Name",
    ]
    available_preview_columns = [column for column in preview_columns if column in filtered.columns]
    preview = filtered[available_preview_columns].head(12).copy() if available_preview_columns else filtered.head(12).copy()

    return {
        "gene_name": gene_name,
        "chrom": chrom,
        "start": start,
        "end": end,
        "build": genome_build,
        "probe_count": int(len(filtered)),
        "output_path": output_path,
        "preview": preview,
    }


def main() -> None:
    """Parse CLI arguments, filter the manifest, and save the region subset."""
    parser = argparse.ArgumentParser(
        description="Filter Illumina manifest probes in a specific genomic region."
    )
    parser.add_argument("--gene_name", required=True, help="Gene name used to name the output file, for example DRD4.")
    parser.add_argument("--manifest", required=True, help="Path to an Illumina manifest file (CSV or CSV.GZ).")
    parser.add_argument("--chrom", required=True, help="Chromosome number or label, for example 11 or chr11.")
    parser.add_argument("--start", type=int, required=True, help="Start coordinate.")
    parser.add_argument("--end", type=int, required=True, help="End coordinate.")
    parser.add_argument(
        "--build",
        choices=["hg19", "hg38"],
        default="hg19",
        help="Genome build used for coordinate matching (default: hg19).",
    )

    args = parser.parse_args()

    result = save_filtered_manifest(
        gene_name=args.gene_name,
        manifest_path=args.manifest,
        region=f"{args.chrom}:{args.start}-{args.end}",
        genome_build=args.build,
        output_dir=PROJECT_ROOT / "src" / "gene_data",
    )

    print("Saving to:", result["output_path"])
    print(f"Found {result['probe_count']} probes in {args.chrom}:{args.start}-{args.end} ({args.build})")
    print(result["preview"])


if __name__ == "__main__":
    main()

# Example use:
# py src/helper_functions/filter_manifest_region.py \
#   --gene_name DRD4 \
#   --manifest data/infinium-methylationepic-v-1-0-b5-manifest-file.csv \
#   --chrom 11 \
#   --start 637269 \
#   --end 640706 \
#   --build hg19
