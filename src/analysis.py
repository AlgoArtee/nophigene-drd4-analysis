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


def load_methylation(idat_base: str, region: str) -> pd.DataFrame:
    """
    Parse Illumina IDATs, compute beta values, and return
    a DataFrame of probes that map to the DRD4 region.
    
    Parameters
    ----------
    idat_base : str
        Path prefix (without _Grn/_Red) to your IDATs, e.g. 'data/R01C01'
    region : str
        Genomic interval 'chr:start-end' or '11:start-end'
    """

    # Build paths
    data_dir   = os.path.dirname(idat_base)
    print(data_dir)
    sample_name = os.path.basename(idat_base)
    print(sample_name)
    grn = os.path.join(data_dir, sample_name + "_Grn.idat")
    red = os.path.join(data_dir, sample_name + "_Red.idat")

    # Sanity checks
    if not os.path.exists(grn) or not os.path.exists(red):
        sys.exit(f"ERROR: Missing IDAT files: {grn}, {red}")


    # run_pipeline will create a cache and export a <sample>_betas.csv
    containers = run_pipeline(
        data_dir,
        export=True,
        export_format=['csv']
    )

    # run_pipeline returns a list of DataContainer; ours is containers[0]
    # and it writes out a CSV named '<sample_name>_betas.csv' in its export folder.
    out_dir     = containers[0].export_folder
    csv_path    = os.path.join(out_dir, f"{sample_name}_betas.csv")
    betas       = pd.read_csv(csv_path, index_col=0)

    # Select & rename columns
    df = betas[['CHR','MAPINFO','Beta_value','Detection_Pval']].rename(
        columns={'CHR':'chrom','MAPINFO':'pos',
                 'Beta_value':'beta','Detection_Pval':'pval'}
    ).reset_index().rename(columns={'index':'probe_id'})

    # Filter to the requested region
    chrom, coords = region.split(':')
    start, end    = map(int, coords.split('-'))
    df_region = df[
        (df['chrom'] == chrom) &
        (df['pos'] >= start) &
        (df['pos'] <= end)
    ].copy()

    if df_region.empty:
        sys.exit(f"ERROR: No methylation probes found in {region}")

    print(f"Loaded {len(df_region)} probes in region {region}")
    print(df_region.head(), "\n")

    return df_region


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


