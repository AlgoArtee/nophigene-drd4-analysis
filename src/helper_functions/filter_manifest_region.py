import os
import pandas as pd
import argparse

def load_manifest(filepath: str) -> pd.DataFrame:
    """
    Loads an Illumina manifest file, skipping the header metadata.
    """
    compression = 'gzip' if filepath.endswith('.gz') else None

    # Find the line where the actual CSV header starts
    with open(filepath, 'rt') as f:
        for i, line in enumerate(f):
            if line.startswith('IlmnID'):
                skiprows = i
                break

    df = pd.read_csv(filepath, skiprows=skiprows, compression=compression, low_memory=False)

    # Normalize and cast
    df.columns = [col.strip() for col in df.columns]
    df['MAPINFO'] = pd.to_numeric(df.get('MAPINFO'), errors='coerce')
    df['Start_hg38'] = pd.to_numeric(df.get('Start_hg38'), errors='coerce')
    df['End_hg38'] = pd.to_numeric(df.get('End_hg38'), errors='coerce')

    return df



def filter_probes_by_region(df: pd.DataFrame, chrom: str, start: int, end: int, genome_build: str = 'hg19') -> pd.DataFrame:
    """
    Filters the manifest DataFrame for probes located within a specified chromosomal region.

    Parameters:
    - df: manifest DataFrame
    - chrom: chromosome (e.g., '11')
    - start: start coordinate (inclusive)
    - end: end coordinate (inclusive)
    - genome_build: 'hg19' or 'hg38'

    Returns:
    - DataFrame of matching probes
    """
    if genome_build == 'hg19':
        chrom_col = 'CHR'
        pos_col = 'MAPINFO'
        chrom = chrom.lstrip("chr")
    elif genome_build == 'hg38':
        chrom_col = 'CHR_hg38'
        pos_col = 'Start_hg38'
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
    else:
        raise ValueError("Unsupported genome build. Use 'hg19' or 'hg38'.")

    filtered = df[
        (df[chrom_col].astype(str) == chrom) &
        (df[pos_col] >= start) &
        (df[pos_col] <= end)
    ].copy()

    return filtered


def main():
    parser = argparse.ArgumentParser(description="Filter Illumina manifest probes in a specific genomic region.")
    parser.add_argument("--gene_name", required=True, help="Gene Name (e.g.: DRD4")
    parser.add_argument("--manifest", required=True, help="Path to Illumina manifest file (CSV or CSV.GZ)")
    parser.add_argument("--chrom", required=True, help="Chromosome number (e.g.: 11)")
    parser.add_argument("--start", type=int, required=True, help="Start coordinate")
    parser.add_argument("--end", type=int, required=True, help="End coordinate")
    parser.add_argument("--build", choices=['hg19', 'hg38'], default='hg19', help="Genome build (default: hg19)")

    args = parser.parse_args()

    df = load_manifest(args.manifest)
    filtered = filter_probes_by_region(df, args.chrom, args.start, args.end, args.build)

    # Output file name
    out_filename = os.path.abspath(f"src/gene_data/{args.gene_name}_epigenetics_{args.build}.csv")
    print("Saving to:", out_filename)
    filtered.to_csv(out_filename, index=False)

    print(f"Found {len(filtered)} probes in {args.chrom}:{args.start}-{args.end} ({args.build})")
    print(filtered[['IlmnID', 'CHR', 'MAPINFO', 'CHR_hg38', 'Start_hg38', 'End_hg38', 'UCSC_RefGene_Name']].head(10))


if __name__ == "__main__":
    main()

# Example Use:
# python src/helper_functions/filter_manifest_region.py \
#  --gene_name DRD4
#  --manifest data/infinium-methylationepic-v-1-0-b5-manifest-file.csv \
#  --chrom 11 \
#  --start 637269 \
#  --end 640706 \
#  --build hg19


# python src/helper_functions/filter_manifest_region.py   --manifest data/infinium-methylationepic-v-1-0-b5-manifest-file.csv   --chrom 11   --start 637269   --end 640706   --build hg19