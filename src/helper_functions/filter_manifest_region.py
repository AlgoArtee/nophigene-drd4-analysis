"""Filter Illumina methylation manifest rows down to a target genomic interval."""

import argparse
import gzip
from pathlib import Path

import pandas as pd


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

    df = load_manifest(args.manifest)
    filtered = filter_probes_by_region(df, args.chrom, args.start, args.end, args.build)

    out_filename = Path("src/gene_data") / f"{args.gene_name}_epigenetics_{args.build}.csv"
    out_filename = out_filename.resolve()
    print("Saving to:", out_filename)
    filtered.to_csv(out_filename, index=False)

    print(f"Found {len(filtered)} probes in {args.chrom}:{args.start}-{args.end} ({args.build})")

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
    if available_preview_columns:
        print(filtered[available_preview_columns].head(10))


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
