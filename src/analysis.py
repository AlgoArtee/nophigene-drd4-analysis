#!/usr/bin/env python3
"""
DRD4 Gene Analysis Pipeline

Usage example:
    python src/analysis.py \
        --vcf data/GFXC926398.filtered.snp.vcf.gz \
        --idat data/202277800037_R01C01 \
        --out results/drd4_report.html \
        --region 11:637269-640706
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any

import allel
import pandas as pd
from methylprep import run_pipeline

# Configure the root logger once so CLI runs stream progress to stdout.
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the DRD4 analysis workflow.

    The CLI is intentionally small: it expects a regional VCF, a paired set of
    Illumina IDAT files expressed as a shared path prefix, and an output path
    for the final report. Optional population statistics can be provided for
    later enrichment once that stage is implemented.

    Returns
    -------
    argparse.Namespace
        Parsed arguments ready to be consumed by ``main()``.
    """
    parser = argparse.ArgumentParser(
        description="DRD4 gene analysis: variants, methylation, population stats, and report generation."
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Tabix-indexed VCF file with filtered SNPs (for example, *.vcf.gz).",
    )
    parser.add_argument(
        "--idat",
        required=True,
        help="Path prefix to Illumina IDATs, without the _Grn/_Red suffixes.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Path for the output report (HTML, CSV, or JSON).",
    )
    parser.add_argument(
        "--region",
        default="11:63671737-63677367",
        help="Genomic region in chr:start-end format. Defaults to the DRD4 GRCh37 interval.",
    )
    parser.add_argument(
        "--popstats",
        default=None,
        help="Optional JSON or CSV with population-frequency data.",
    )

    return parser.parse_args()


def load_variants(vcf_path: str, region: str) -> pd.DataFrame:
    """Load PASS variants for the requested region from a tabix-indexed VCF.

    Parameters
    ----------
    vcf_path : str
        Filesystem path to a bgzip-compressed and tabix-indexed VCF.
    region : str
        Region string understood by ``scikit-allel`` (for example,
        ``11:63671737-63677367``).

    Returns
    -------
    pd.DataFrame
        Variant table with one row per PASS site and the columns ``chrom``,
        ``pos``, ``ref``, ``alt``, ``qual``, and ``filter_pass``.

    Notes
    -----
    ``variants/ALT`` is modeled as an array in the VCF reader. The pipeline
    currently keeps only the first alternate allele so downstream tables stay
    scalar and easy to serialize.
    """
    if not os.path.exists(vcf_path):
        sys.exit(f"ERROR: VCF not found: {vcf_path}")

    callset = allel.read_vcf(
        vcf_path,
        region=region,
        fields=[
            "variants/CHROM",
            "variants/POS",
            "variants/REF",
            "variants/ALT",
            "variants/QUAL",
            "variants/FILTER_PASS",
        ],
    )
    if not callset:
        sys.exit(f"ERROR: No variants could be read from {vcf_path} in region {region}")

    df = pd.DataFrame(
        {
            "chrom": callset["variants/CHROM"],
            "pos": callset["variants/POS"],
            "ref": callset["variants/REF"],
            "alt": [alts[0] if len(alts) > 0 else None for alts in callset["variants/ALT"]],
            "qual": callset["variants/QUAL"],
            "filter_pass": callset["variants/FILTER_PASS"],
        }
    )

    # Keep only variants that passed upstream filtering so the status column and
    # the error message below are consistent with the actual output table.
    df = df[df["filter_pass"].fillna(False)].reset_index(drop=True)

    if df.empty:
        sys.exit(f"ERROR: No PASS variants found in {region} for {vcf_path}")

    print(f"Loaded {len(df)} PASS variants from {vcf_path} in region {region}")
    print(df.head(), "\n")
    return df


def load_methylation(idat_base: str, manifest_filepath: str | None = None) -> pd.DataFrame:
    """Load methylation beta values and annotate them with DRD4 region metadata.

    Parameters
    ----------
    idat_base : str
        Path prefix shared by the two Illumina IDAT files. For a sample stored as
        ``data/R01C01_Grn.idat`` and ``data/R01C01_Red.idat``, pass
        ``data/R01C01``.
    manifest_filepath : str | None, optional
        Optional path to a custom manifest file consumed directly by
        ``methylprep.run_pipeline``. When omitted, methylprep chooses its
        default manifest behavior.

    Returns
    -------
    pd.DataFrame
        Probe-level methylation table limited to probes present in the curated
        DRD4 region manifest. The returned DataFrame keeps the beta values plus
        the most useful annotation columns that already exist in the local CSV.

    Notes
    -----
    This function intentionally joins against
    ``src/gene_data/drd4_epigenetics_hg19.csv`` instead of filtering a whole
    manifest on the fly. That keeps the runtime path deterministic and makes the
    output consistent with the repository's curated DRD4 probe selection.
    """
    logger.info("Starting methylation loading for sample")

    # ``idat_base`` is a prefix, not a literal file. Split it once so we can
    # validate the paired files and later match the sample column produced by
    # methylprep.
    data_dir = os.path.dirname(idat_base) or "."
    sample_name = os.path.basename(idat_base)

    logger.debug("Data directory: %s", data_dir)
    logger.debug("Sample name: %s", sample_name)

    # The assay is only valid when both channels are present. Failing early here
    # produces a clearer message than letting methylprep error deeper inside.
    for suffix in ("_Grn.idat", "_Red.idat"):
        path = os.path.join(data_dir, sample_name + suffix)
        logger.debug("Checking IDAT: %s", path)
        if not os.path.isfile(path):
            logger.error("Missing IDAT file: %s", path)
            sys.exit(1)

    try:
        logger.info("Running methylprep pipeline with betas=True")
        beta_values = run_pipeline(
            data_dir,
            export=True,
            betas=True,
            manifest_filepath=manifest_filepath,
        )
        logger.debug("Beta-values DataFrame loaded, shape = %s", beta_values.shape)
        logger.debug("First 5 rows of wide beta-values DataFrame:\n%s", beta_values.head(5))
    except Exception:
        logger.exception("run_pipeline failed")
        sys.exit(1)

    if sample_name not in beta_values.columns:
        logger.error("No beta column named '%s' in output", sample_name)
        sys.exit(1)

    # The pipeline returns one column per sample. Reset the index so probe IDs
    # become an explicit merge key instead of staying hidden in the index.
    beta_df = beta_values[sample_name].rename("beta").reset_index()
    probe_column = beta_df.columns[0]
    beta_df = beta_df.rename(columns={probe_column: "probe_id"})
    logger.debug("First 5 rows of beta-values DataFrame after renaming:\n%s", beta_df.head(5))

    region_manifest_file = Path(__file__).resolve().parent / "gene_data" / "drd4_epigenetics_hg19.csv"
    try:
        manifest_region = pd.read_csv(region_manifest_file)
    except Exception:
        logger.exception("Failed to read region manifest file: %s", region_manifest_file)
        sys.exit(1)

    # Normalize the key columns so the join below can stay simple and the output
    # schema is easier to reason about in notebooks or reports.
    rename_cols = {
        "IlmnID": "probe_id",
        "CHR": "chrom",
        "MAPINFO": "pos",
        "UCSC_RefGene_Name": "gene",
    }
    manifest_region = manifest_region.rename(columns=rename_cols)
    logger.debug("First 5 rows of renamed filtered manifest:\n%s", manifest_region.head(5))

    required = {"probe_id", "chrom", "pos", "gene"}
    missing = required - set(manifest_region.columns)
    if missing:
        logger.error("Region manifest is missing required columns: %s", missing)
        sys.exit(1)

    # An inner join keeps only probes observed in both the sample output and the
    # curated DRD4 manifest, which is the narrow result set expected downstream.
    merged = pd.merge(beta_df, manifest_region, on="probe_id", how="inner")
    logger.info("After merging, %d probes remain", len(merged))

    keep_columns = [
        "probe_id",
        "beta",
        "chrom",
        "pos",
        "GencodeBasicV12_NAME",
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Accession",
        "UCSC_RefGene_Group",
        "UCSC_CpG_Islands_Name",
        "Relation_to_UCSC_CpG_Island",
        "Phantom4_Enhancers",
        "Phantom5_Enhancers,DMR,450k_Enhancer,HMM_Island",
        "DNase_Hypersensitivity_NAME",
        "DNase_Hypersensitivity_Evidence_Count",
    ]

    available = set(merged.columns)
    missing = [column for column in keep_columns if column not in available]
    if missing:
        logger.warning("Some expected columns are missing: %s", missing)

    final_df = merged[[column for column in keep_columns if column in available]]
    logger.debug("First 5 rows of renamed filtered epigenetic data table:\n%s", final_df.head(5))
    return final_df


def fetch_population_stats(popstats_source: str, variants: pd.DataFrame) -> Any:
    """Load or derive population statistics for the supplied variants.

    Parameters
    ----------
    popstats_source : str
        Path or identifier for the external population statistics source.
    variants : pd.DataFrame
        Variant table returned by :func:`load_variants`.

    Returns
    -------
    Any
        Placeholder return type until this stage is implemented.

    Notes
    -----
    The function is kept as an explicit stub so the pipeline shape is visible
    while the population enrichment step is still under development.
    """
    raise NotImplementedError("fetch_population_stats() not yet implemented")


def generate_report(
    variants: pd.DataFrame,
    methylation: pd.DataFrame,
    popstats: Any,
    output_path: str,
) -> None:
    """Generate the final analysis report from the assembled data tables.

    Parameters
    ----------
    variants : pd.DataFrame
        PASS variant table for the requested DRD4 interval.
    methylation : pd.DataFrame
        Annotated probe-level beta-value table returned by
        :func:`load_methylation`.
    popstats : Any
        Optional population statistics payload. The current pipeline leaves this
        as ``None`` unless the enrichment step is implemented.
    output_path : str
        Destination path for the final report artifact.

    Notes
    -----
    The implementation is intentionally still a placeholder. Keeping the
    function boundary documented now makes it easier to replace the print with a
    real templating layer later.
    """
    _ = (variants, methylation, popstats)
    print(f"Report would be written to {output_path}")


def main() -> None:
    """Run the end-to-end DRD4 analysis workflow from the command line.

    The current pipeline performs three visible steps:

    1. load regional variants from the VCF,
    2. process methylation IDATs and save the annotated probe table to CSV,
    3. call the report generator placeholder.

    Population statistics are wired in as an optional future step and therefore
    only execute when the corresponding flag is provided.
    """
    args = parse_args()

    variants = load_variants(args.vcf, args.region)
    meth = load_methylation(args.idat)

    # Persist the intermediate methylation table so notebooks and manual review
    # can inspect the processed probe-level output independently of report
    # generation.
    output_path = Path("results/methylation_output.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    meth.to_csv(output_path, index=False)
    print(f"Saved methylation data to {output_path}")

    popstats = None
    if args.popstats:
        popstats = fetch_population_stats(args.popstats, variants)

    generate_report(variants, meth, popstats, args.out)


if __name__ == "__main__":
    print(__doc__)
    main()
