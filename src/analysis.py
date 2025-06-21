#!/usr/bin/env python3
"""
DRD4 Gene Analysis Pipeline

Usage example:
    python src/analysis.py \
        --vcf data/GFXC926398.filtered.snp.vcf.gz \
        --idat data/202277800037_R01C01 \
        --out results/drd4_report.html \
        --region 11:63671737-63677367
"""

import argparse
import os
import sys
import allel
import pandas as pd
from methylprep import run_pipeline


def parse_args():
    parser = argparse.ArgumentParser(
        description="DRD4 gene analysis: variants, methylation, population stats, and report generation."
    )
    parser.add_argument("--vcf", required=True,
                        help="Tabix-indexed VCF file with filtered SNPs (e.g. *.vcf.gz).")
    parser.add_argument("--idat", required=True,
                        help="Path prefix to your Illumina IDATs (e.g. data/R01C01, without _Grn/_Red.suffix).")
    parser.add_argument("--out", required=True,
                        help="Path for the output report (HTML, CSV, or JSON).")
    parser.add_argument("--region", default="11:63671737-63677367",
                        help="Genomic region (chr:start-end), default DRD4 GRCh37.")
    parser.add_argument("--popstats", default=None,
                        help="Optional JSON or CSV with population-frequency data.")

    return parser.parse_args()


def load_variants(vcf_path: str, region: str) -> pd.DataFrame:
    """Load SNP variants from a tabix-indexed VCF into a pandas DataFrame."""
    if not os.path.exists(vcf_path):
        sys.exit(f"ERROR: VCF not found: {vcf_path}")

    callset = allel.read_vcf(
        vcf_path,
        region=region,
        fields=[
            'variants/CHROM',
            'variants/POS',
            'variants/REF',
            'variants/ALT',
            'variants/QUAL',
            'variants/FILTER_PASS'
        ]
    )
    df = pd.DataFrame({
        'chrom': callset['variants/CHROM'],
        'pos': callset['variants/POS'],
        'ref': callset['variants/REF'],
        'alt': [alts[0] if len(alts) > 0 else None
                for alts in callset['variants/ALT']],
        'qual': callset['variants/QUAL'],
        'filter_pass': callset['variants/FILTER_PASS']
    })

    if df.empty:
        sys.exit(f"ERROR: No PASS variants found in {region} for {vcf_path}")

    print(f"Loaded {len(df)} variants from {vcf_path} in region {region}")
    print(df.head(), "\n")

    return df


def load_methylation(idat_base: str, region: str, manifest_filepath: str = None) -> pd.DataFrame:
    """
    Parse Illumina IDATs (no sample sheet needed), compute beta values
    in-memory, and return a DataFrame of probes that map to the specified region.

    Parameters
    ----------
    idat_base : str
        Path prefix (without _Grn/_Red) to your IDAT files, e.g. 'data/R01C01'.
        Expects to find 'R01C01_Grn.idat' and 'R01C01_Red.idat' in the same folder.
    region : str
        Genomic interval in 'chr:start-end' or '11:start-end' format, e.g. '11:63671737-63677367'.

    Returns
    -------
    pd.DataFrame
        Long-form DataFrame with columns:
          - probe_id: Illumina probe identifier
          - beta: methylation beta value (0–1)
          - chrom: chromosome of the probe
          - pos: genomic coordinate of the probe
          - pval: detection p-value for the probe
    """
    # 1. Determine data directory and sample name
    data_dir    = os.path.dirname(idat_base)
    sample_name = os.path.basename(idat_base)

    # 2. Confirm the IDAT files exist
    grn = os.path.join(data_dir, f"{sample_name}_Grn.idat")
    red = os.path.join(data_dir, f"{sample_name}_Red.idat")
    if not os.path.exists(grn) or not os.path.exists(red):
        sys.exit(f"ERROR: Missing IDAT files: {grn}, {red}")

    # Single-sample pipeline run with automatic manifest download
    betas_wide = run_pipeline(
        data_dir,
        export=True,            
        betas=True,
        make_sample_sheet=True,
        
    )

    processed_dir_path = sample_name.split("_")[0]
    print("Processed files folder: ", processed_dir_path)
    csv_path = os.path.join(data_dir, processed_dir_path,  f"{sample_name}_processed.csv")

    # 3) Read and rename
    df = pd.read_csv(csv_path, index_col=0)
    df = df.rename(
        columns={
            'beta_value':    'beta',
            'Detection_Pval':'pval',
            'CHR':           'chrom',
            'MAPINFO':       'pos'
        }
    )
    df = df.reset_index().rename(columns={'index':'probe_id'})

    print(df.head)

    
    # 7. Filter to the specified region
    chrom, coords = region.replace('chr','').split(':')
    start, end    = map(int, coords.split('-'))
    df_region = df[
        (df['chrom'] == chrom) &
        (df['pos'] >= start) &
        (df['pos'] <= end)
    ].copy()
    if df_region.empty:
        sys.exit(f"ERROR: No methylation probes found in region {region}.")

    # 8. Return the filtered DataFrame
    print(f"Loaded {len(df_region)} probes in region {region}")
    return df_region[['probe_id','beta','chrom','pos','pval']]


def fetch_population_stats(popstats_source, variants):
    """Annotate variants with population frequencies."""
    # TODO: implement in Step 4
    raise NotImplementedError("fetch_population_stats() not yet implemented")


def generate_report(variants, methylation, popstats, output_path: str):
    """Compile variants, methylation, and popstats into a report."""
    # TODO: implement Step 5 (e.g., Jinja2 template → HTML)
    print(f"Report would be written to {output_path}")


def main():
    args = parse_args()

    variants = load_variants(args.vcf, args.region)
    meth = load_methylation(args.idat, args.region)
    popstats = None
    if args.popstats:
        popstats = fetch_population_stats(args.popstats, variants)

    generate_report(variants, meth, popstats, args.out)


if __name__ == "__main__":
    print(__doc__)
    main()


