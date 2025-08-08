#!/usr/bin/env python3
"""
DRD4 Gene Analysis Pipeline

Usage example:
    python src/analysis.py \
        --vcf data/GFXC926398.filtered.snp.vcf.gz \
        --idat data/202277800037_R01C01 \
        --out results/drd4_report.html \
        --region 11:637269-640706 // source: https://genome-euro.ucsc.edu/cgi-bin/hgGene?db=hg19&hgg_gene=DRD4  old:11:63671737-63677367
"""

import argparse
import os
import sys
import logging
import allel
import pandas as pd
from methylprep import Manifest, run_pipeline
from pathlib import Path

# -------------------------------------------------
# Configure root logger for stdout at DEBUG level
# -------------------------------------------------
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


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

    # === 1. Setup paths and validate input ===
    logger.info("Starting methylation loading for sample")

    data_dir    = os.path.dirname(idat_base) # directory containing the IDATs
    sample_name = os.path.basename(idat_base) # ex: 202277800037_R01C01
    #sample_id   = sample_name.split("_")[0]

    logger.debug(f"Data directory:  {data_dir}")
    logger.debug(f"Sample name:     {sample_name}")
    #logger.debug(f"Sample id:       {sample_id}")


    # 2a) === 2. Check presence of both Red and Green Channel IDAT files ===
    for suffix in ("_Grn.idat", "_Red.idat"):
        path = os.path.join(data_dir, sample_name + suffix)
        logger.debug(f"Checking IDAT: {path}")
        if not os.path.isfile(path):
            logger.error(f"Missing IDAT file: {path}")
            sys.exit(1)


    # === 2. Run pipeline to get wide beta matrix ===
    try:
        logger.info("Running methylprep pipeline with betas=True …")
        beta_values = run_pipeline(
            data_dir,
            export=True,
            betas=True,
            manifest_filepath=manifest_filepath
        )
        logger.debug(f"β-values DataFrame loaded, shape = {beta_values.shape}")
        # Print the first 5 lines for inspection
        logger.debug("First 5 rows of wide β‑values DataFrame:\n%s", beta_values.head(5))
    except Exception:
        logger.exception("run_pipeline failed")
        sys.exit(1)

    if sample_name not in beta_values.columns:
        logger.error(f"No β column named '{sample_name}' in output")
        sys.exit(1)


    # Result is DataFrame: rows = probes, columns = samples


    # 3. Extract the beta series directly (no melt needed)
    beta_series = beta_values[sample_name]
    beta_df = beta_series.reset_index().rename(columns={'IlmnID': 'probe_id', sample_name: 'beta'})
    logger.debug("First 5 rows of β‑values DataFrame with Renamed Columns:\n%s", beta_values.head(5))



    # === 4. Load region manifest and validate ===
    region_manifest_file = "src/gene_data/drd4_epigenetics_hg19.csv"
    try:
        manifest_region = pd.read_csv(region_manifest_file)
    except Exception:
        logger.exception("Failed to read region_manifest_file")
        sys.exit(1)

    
    # === Rename columns from Illumina format to expected schema ===
    rename_cols = {
        'IlmnID': 'probe_id',
        'CHR': 'chrom',
        'MAPINFO': 'pos',
        'UCSC_RefGene_Name': 'gene'
    }
    manifest_region.rename(columns=rename_cols, inplace=True)
    logger.debug("First 5 rows of renamed filtered manifest:\n%s", manifest_region.head(5))

    # === Validate required columns ===
    required = {'probe_id', 'chrom', 'pos', 'gene'}
    missing = required - set(manifest_region.columns)
    if missing:
        logger.error(f"Region manifest is missing required columns: {missing}")
        sys.exit(1)

    # 5. Merge and filter
    # Inner join keeps only probes present in both the sample and the manifest
    merged = pd.merge(beta_df, manifest_region, on='probe_id', how='inner')
    logger.info("After merging: %d probes remain", len(merged))

    # === Keep only desired columns ===
    keep_columns = [
        'probe_id',
        'beta',
        'chrom',
        'pos',
        'GencodeBasicV12_NAME',
        'UCSC_RefGene_Name',
        'UCSC_RefGene_Accession',
        'UCSC_RefGene_Group',
        'UCSC_CpG_Islands_Name',
        'Relation_to_UCSC_CpG_Island',
        'Phantom4_Enhancers',
        'Phantom5_Enhancers,DMR,450k_Enhancer,HMM_Island',
        'DNase_Hypersensitivity_NAME',
        'DNase_Hypersensitivity_Evidence_Count'
    ]

    available = set(merged.columns)
    missing = [col for col in keep_columns if col not in available]
    if missing:
        logger.warning(f"Some expected columns are missing: {missing}")

    final_df = merged[[col for col in keep_columns if col in available]]
    logger.debug("First 5 rows of renamed filtered epigenetic data table:\n%s", final_df.head(5))
    return final_df


# # === 5. Find the processed output CSV file ===
    # sample_folder = os.path.join(data_dir, sample_id)
    # candidates = [
    #     f"{sample_name}_processed.csv",
    #     f"{sample_name}_betas.csv",
    #     "beta_values.csv"
    # ]
    # for fname in candidates:
    #     csv_path = os.path.join(sample_folder, fname)
    #     if os.path.isfile(csv_path):
    #         logger.info(f"Found processed CSV file: {csv_path}")
    #         break
    # else:
    #     logger.error("No valid processed CSV found.")
    #     sys.exit(1)


    # # === 6. Load methylation data into a DataFrame ===
    # try:
    #     df = pd.read_csv(csv_path, index_col=0)
    #     logger.debug(f"Loaded methylation DataFrame with shape {df.shape}")
    # except Exception as e:
    #     logger.exception("Failed to read processed CSV.")
    #     sys.exit(1)



    # # === 7. Rename columns to standardized schema ===
    # df.rename(columns={
    #     'beta_value': 'beta',
    #     'Detection_Pval': 'pval',
    #     'CHR': 'chrom',
    #     'MAPINFO': 'pos'
    # }, inplace=True)

    # # Move probe IDs from index to a column
    # df.reset_index(inplace=True)
    # df.rename(columns={'index': 'probe_id'}, inplace=True)

    # logger.debug(f"Data columns after rename: {list(df.columns)}")

    # # === 8. Platform inference (based on probe count) ===
    # n_probes = df.shape[0]
    # if n_probes > 900_000:
    #     array_type = "epic"
    # elif n_probes > 400_000:
    #     array_type = "450k"
    # elif n_probes < 30_000:
    #     array_type = "27k"
    # else:
    #     array_type = "custom"
    # logger.info(f"Inferred array type: {array_type}")


    # # === 9. Load manifest and annotate ===
    # try:
    #     manifest = Manifest(array_type)
    #     manifest_df = manifest.data_frame
    # except Exception as e:
    #     logger.exception("Failed to load manifest using methylprep.")
    #     sys.exit(1)

    # # Try to retrieve gene symbols from full manifest file if not in data_frame
    # if 'UCSC_RefGene_Name' not in manifest_df.columns:
    #     cache_dir = Path.home() / ".methylprep_manifest_files"
    #     manifest_files = list(cache_dir.glob(f"*{array_type}*manifest*csv*"))
    #     if not manifest_files:
    #         logger.error("No cached manifest file found.")
    #         sys.exit(1)
    #     full_manifest = pd.read_csv(manifest_files[0], comment='[', header=0, low_memory=False)
    #     full_manifest.set_index('IlmnID', inplace=True)
    #     manifest_df = manifest_df.join(full_manifest[['UCSC_RefGene_Name']], how='left')

    
    # # === 10. Standardize annotation columns ===
    # manifest_df = manifest_df.rename(columns={
    #     'CHR': 'chrom',
    #     'MAPINFO': 'pos',
    #     'UCSC_RefGene_Name': 'gene'
    # })
    # annot = manifest_df[['chrom', 'pos', 'gene']]
    
    # # === 11. Merge annotation with beta data ===
    # annotated_df = pd.merge(df, annot, left_on='probe_id', right_index=True, how='left')

    # # === 12. Parse region string ===
    # try:
    #     chrom_str, pos_str = region.replace("chr", "").split(":")
    #     start, end = map(int, pos_str.split("-"))
    #     chrom_str = str(chrom_str)
    # except Exception:
    #     logger.error("Region format should be like 'chr11:123456-124000'")
    #     sys.exit(1)
    
    
    # # === 13. Filter to region ===
    # annotated_df['chrom'] = annotated_df['chrom'].astype(str)
    # region_df = annotated_df[
    #     (annotated_df['chrom'] == chrom_str) &
    #     (annotated_df['pos'] >= start) &
    #     (annotated_df['pos'] <= end)
    # ].copy()
    
    
    # logger.info(f"Filtered to region chr{chrom_str}:{start}-{end} — {len(region_df)} probes found")
    # if region_df.empty:
    #     logger.warning("No probes found in region.")
    #     sys.exit(1)

    # # Final output
    # return region_df[['probe_id', 'beta', 'chrom', 'pos', 'gene', 'pval']]
    


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
    
    
    # Save Methylation Data to CSV
    output_path = "results/methylation_output.csv"
    meth.to_csv(output_path, index=False)
    print(f"Saved methylation data to {output_path}")
    
    
    
    
    popstats = None
    if args.popstats:
        popstats = fetch_population_stats(args.popstats, variants)

    generate_report(variants, meth, popstats, args.out)


if __name__ == "__main__":
    print(__doc__)
    main()


