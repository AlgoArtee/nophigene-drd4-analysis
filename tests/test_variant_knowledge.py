import json
from pathlib import Path

import pandas as pd
import pytest
import requests

from src.variant_knowledge import client as request_client_module
from src.variant_knowledge.client import KnowledgeRequestError, RequestClient
from src.variant_knowledge.connectors import connector_for
from src.variant_knowledge.credentials import ResolvedCredential, credential_status_for_specs
from src.variant_knowledge.imports import parse_source_import
from src.variant_knowledge.local_articles import (
    LOCAL_ARTICLE_SOURCE_KEY,
    LOCAL_ARTICLE_WORKFLOW_KEY,
    extract_local_article_evidence,
)
from src.variant_knowledge.merger import merge_dynamic_knowledge_base
from src.variant_knowledge.models import KnowledgeQuery, QueryVariant
from src.variant_knowledge.orchestrator import build_dynamic_knowledge_base
from src.variant_knowledge.registry import RESOURCES_PATH, get_source_spec, list_source_cards, list_source_specs
from src.variant_knowledge.workflows import CORE_SAFETY_WORKFLOW_KEYS, list_workflow_specs


class FakeClient:
    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        if "esearch.fcgi" in url:
            return {"esearchresult": {"idlist": ["123"]}}
        if "esummary.fcgi" in url:
            return {
                "result": {
                    "123": {
                        "uid": "123",
                        "title": "ClinVar rs123 pathogenic record",
                        "description": "Pathogenic clinical assertion for rs123",
                    }
                }
            }
        if "variation/v0/refsnp" in url:
            return {}
        raise AssertionError(f"Unexpected fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class CountingFakeClient(FakeClient):
    def __init__(self) -> None:
        self.get_counts: dict[str, int] = {}

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        if "esearch.fcgi" in url:
            term = str((params or {}).get("term", ""))
            self.get_counts[term] = self.get_counts.get(term, 0) + 1
        return super().get_json(
            url,
            params=params,
            headers=headers,
            rate_limit_per_second=rate_limit_per_second,
            timeout=timeout,
        )


class RecordingGnomadClient:
    def __init__(self, *, graphql_errors: bool = False) -> None:
        self.graphql_errors = graphql_errors
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        raise AssertionError(f"Unexpected fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        self.calls.append(
            {
                "url": url,
                "json_payload": json_payload or {},
                "headers": headers or {},
                "rate_limit_per_second": rate_limit_per_second,
            }
        )
        if self.graphql_errors:
            return {"errors": [{"message": "Cannot query field \"badField\" on type \"VariantDetails\"."}]}
        return {
            "data": {
                "meta": {"clinvar_release_date": "2026-06-06"},
                "variant": {
                    "variantId": "11-637293-C-T",
                    "variant_id": "11-637293-C-T",
                    "reference_genome": "GRCh38",
                    "chrom": "11",
                    "pos": 637293,
                    "ref": "C",
                    "alt": "T",
                    "caid": "CA216942642",
                    "rsid": "rs927984495",
                    "rsids": ["rs927984495"],
                    "flags": [],
                    "exome": {
                        "ac": 14,
                        "an": 1044702,
                        "af": 1.3400950701731212e-05,
                        "homozygote_count": 0,
                        "hemizygote_count": 0,
                        "filters": [],
                        "flags": [],
                        "faf95": {"popmax": 9.05e-06, "popmax_population": "nfe"},
                        "faf99": {"popmax": None, "popmax_population": None},
                        "populations": [
                            {"id": "nfe", "ac": 14, "an": 897410, "homozygote_count": 0, "hemizygote_count": 0},
                            {"id": "afr", "ac": 0, "an": 20982, "homozygote_count": 0, "hemizygote_count": 0},
                        ],
                    },
                    "genome": {
                        "ac": 6,
                        "an": 151146,
                        "af": 3.969671708149736e-05,
                        "homozygote_count": 0,
                        "hemizygote_count": 0,
                        "filters": [],
                        "flags": [],
                        "faf95": {"popmax": 3.776e-05, "popmax_population": "nfe"},
                        "faf99": {"popmax": None, "popmax_population": None},
                        "populations": [
                            {"id": "nfe", "ac": 6, "an": 67718, "homozygote_count": 0, "hemizygote_count": 0},
                            {"id": "XX", "ac": 3, "an": 77334, "homozygote_count": 0, "hemizygote_count": 0},
                        ],
                    },
                    "joint": {
                        "ac": 20,
                        "an": 1195848,
                        "homozygote_count": 0,
                        "hemizygote_count": 0,
                        "filters": ["discrepant_frequencies"],
                        "faf95": {"popmax": 1.353e-05, "popmax_population": "nfe"},
                        "faf99": None,
                        "populations": [
                            {"id": "nfe", "ac": 20, "an": 965128, "homozygote_count": 0, "hemizygote_count": 0},
                            {"id": "nfe_XX", "ac": 8, "an": 512460, "homozygote_count": 0, "hemizygote_count": 0},
                        ],
                    },
                    "transcript_consequences": [
                        {
                            "gene_id": "ENSG00000069696",
                            "gene_symbol": "DRD4",
                            "transcript_id": "ENST00000176183",
                            "transcript_version": "5",
                            "is_canonical": True,
                            "is_mane_select": True,
                            "major_consequence": "5_prime_UTR_variant",
                            "consequence_terms": ["5_prime_UTR_variant"],
                            "hgvsc": "ENST00000176183.5:c.-11C>T",
                            "hgvsp": None,
                            "hgvs": "ENST00000176183.5:c.-11C>T",
                            "lof": None,
                            "lof_filter": None,
                            "lof_flags": None,
                            "polyphen_prediction": None,
                            "sift_prediction": None,
                            "domains": [],
                        }
                    ],
                    "in_silico_predictors": [{"id": "CADD", "value": "3.2", "flags": []}],
                    "lof_curations": [],
                    "non_coding_constraint": None,
                },
                "gene": {
                    "gene_id": "ENSG00000069696",
                    "gene_version": "17",
                    "symbol": "DRD4",
                    "gencode_symbol": "DRD4",
                    "hgnc_id": "HGNC:3025",
                    "ncbi_id": "1815",
                    "omim_id": "126452",
                    "name": "dopamine receptor D4",
                    "reference_genome": "GRCh38",
                    "chrom": "11",
                    "start": 637269,
                    "stop": 640706,
                    "strand": "+",
                    "canonical_transcript_id": "ENST00000176183",
                    "mane_select_transcript": {
                        "ensembl_id": "ENST00000176183",
                        "ensembl_version": "5",
                        "refseq_id": "NM_000797",
                        "refseq_version": "4",
                    },
                    "flags": [],
                    "gnomad_constraint": {
                        "exp_lof": 3.405,
                        "exp_mis": 79.211,
                        "exp_syn": 31.2,
                        "obs_lof": 4,
                        "obs_mis": 114,
                        "obs_syn": 31,
                        "oe_lof": 1.1746694306580123,
                        "oe_lof_lower": 0.512,
                        "oe_lof_upper": 2.211,
                        "oe_lof_percentile": 81,
                        "oe_mis": 1.4437142640429352,
                        "oe_mis_lower": 1.18,
                        "oe_mis_upper": 1.76,
                        "oe_syn": 0.99,
                        "oe_syn_lower": 0.7,
                        "oe_syn_upper": 1.2,
                        "lof_z": -0.6971243887353076,
                        "mis_z": -4.75889174564733,
                        "syn_z": 0.1,
                        "pli": 1.0280908139937256e-10,
                        "pLI": 1.0280908139937256e-10,
                        "flags": ["outlier_mis", "outlier_syn"],
                    },
                },
            }
        }


class RecordingGwasCatalogClient:
    METADATA = {
        "title": "GWAS Catalog Rest API 2.0",
        "version": "2.0",
        "data_release_date": "2026-06-22",
        "api_release_date": "2025-08-01",
        "dbsnp_build": "156",
        "gene_build": "GRCh38.p14",
        "efo_version": "v3.91.0",
    }
    ASSOCIATION = {
        "association_id": 216548355,
        "risk_frequency": "0.2688",
        "pvalue_description": "",
        "pvalue_mantissa": 5,
        "pvalue_exponent": -18,
        "multi_snp_haplotype": False,
        "snp_interaction": False,
        "range": "[0.027-0.043]",
        "beta": "0.0352 unit decrease",
        "p_value": 5e-18,
        "efo_traits": [{"efo_id": "MONDO_0005420", "efo_trait": "hypothyroidism"}],
        "reported_trait": ["Hypothyroidism"],
        "accession_id": "GCST90627750",
        "locations": ["11:640349"],
        "mapped_genes": ["DRD4"],
        "bg_efo_traits": [],
        "pubmed_id": "41644669",
        "first_author": "White SL",
        "ci_lower": 0.027,
        "ci_upper": 0.043,
        "_links": {
            "self": {"href": "https://www.ebi.ac.uk/gwas/rest/api/v2/associations/216548355"},
            "loci": {"href": "https://www.ebi.ac.uk/gwas/rest/api/v2/associations/216548355/loci"},
            "snp": {
                "href": "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms/rs1870723"
            },
        },
        "snp_effect_allele": ["rs1870723-A"],
        "snp_allele": [{"rs_id": "rs1870723", "effect_allele": "A"}],
    }
    STUDY = {
        "initial_sample_size": (
            "257,365 African ancestry, admixed American, East Asian ancestry, European ancestry, "
            "South Asian ancestry, Greater Middle Eastern cases, 2,186,763 controls"
        ),
        "replication_sample_size": "NA",
        "snp_count": 28645100,
        "imputed": True,
        "pooled": False,
        "accession_id": "GCST90627750",
        "full_summary_stats_available": True,
        "pubmed_id": 41644669,
        "platforms": "NR [28645100] (imputed)",
        "disease_trait": "Hypothyroidism",
        "genotyping_technologies": [
            "Genome-wide sequencing",
            "Exome genotyping array",
            "Genome-wide genotyping array",
        ],
        "discovery_ancestry": [
            "2444128 African unspecified,Hispanic or Latin American,East Asian,European,South Asian,Greater Middle Eastern"
        ],
        "replication_ancestry": [],
        "full_summary_stats": (
            "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90627001-GCST90628000/GCST90627750"
        ),
        "terms_of_license": "https://creativecommons.org/publicdomain/zero/1.0/",
        "cohort": ["AllofUs", "BBJ", "BioVU", "FinnGen", "UKB"],
    }
    LOCI = {
        "_embedded": {
            "loci": [
                {
                    "description": "Single variant",
                    "strongest_risk_alleles": [
                        {"risk_allele_name": "rs1870723-A", "genome_wide": False, "limited_list": False}
                    ],
                    "author_reported_genes": ["DRD4"],
                    "_links": {
                        "self": {
                            "href": "https://www.ebi.ac.uk/gwas/rest/api/v2/associations/216548355/loci/216548354"
                        }
                    },
                }
            ]
        }
    }

    def __init__(self, *, rsid_has_records: bool = False, fail_optional: bool = False) -> None:
        self.rsid_has_records = rsid_has_records
        self.fail_optional = fail_optional
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if str(url).endswith("/metadata"):
            return dict(self.METADATA)
        if str(url).endswith("/associations"):
            if params.get("rs_id") and not self.rsid_has_records:
                return {"_embedded": {"associations": []}, "page": {"totalElements": 0, "size": 5, "number": 0}}
            return {
                "_embedded": {"associations": [json.loads(json.dumps(self.ASSOCIATION))]},
                "page": {"totalElements": 1, "size": 5, "number": 0},
            }
        if str(url).endswith("/studies/GCST90627750"):
            if self.fail_optional:
                raise KnowledgeRequestError("GET GWAS study failed: timed out")
            return json.loads(json.dumps(self.STUDY))
        if str(url).endswith("/associations/216548355/loci"):
            if self.fail_optional:
                raise KnowledgeRequestError("GET GWAS loci failed: timed out")
            return json.loads(json.dumps(self.LOCI))
        raise AssertionError(f"Unexpected GWAS Catalog fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingPgsCatalogClient:
    VARIANT = {"id": "rs123", "associated_pgs_ids": ["PGS000001"]}
    SCORE = {
        "id": "PGS000001",
        "name": "PRS77_BC",
        "ftp_scoring_file": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/PGS000001.txt.gz",
        "ftp_harmonized_scoring_files": {
            "GRCh37": {
                "positions": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/Harmonized/PGS000001_hmPOS_GRCh37.txt.gz"
            },
            "GRCh38": {
                "positions": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/Harmonized/PGS000001_hmPOS_GRCh38.txt.gz"
            },
        },
        "publication": {
            "id": "PGP000001",
            "title": "Prediction of breast cancer risk based on profiling with common genetic variants.",
            "doi": "10.1093/jnci/djv036",
            "PMID": 25855707,
            "journal": "J Natl Cancer Inst",
            "firstauthor": "Mavaddat N",
            "date_publication": "2015-04-08",
        },
        "matches_publication": True,
        "samples_variants": [
            {
                "sample_number": 22627,
                "ancestry_broad": "European",
                "ancestry_country": "Finland, Sweden, U.S., Australia, Netherlands, Germany, U.K.",
                "source_GWAS_catalog": "GCST001937",
                "source_PMID": 23535729,
                "cohorts": [],
            }
        ],
        "samples_training": [],
        "trait_reported": "Breast cancer",
        "trait_additional": None,
        "trait_efo": [
            {
                "id": "MONDO_0004989",
                "label": "breast carcinoma",
                "description": "A carcinoma that arises from epithelial cells of the breast",
                "url": "http://purl.obolibrary.org/obo/MONDO_0004989",
            }
        ],
        "method_name": "SNPs passing genome-wide significance",
        "method_params": "P<5x10-8",
        "variants_number": 77,
        "variants_interactions": 0,
        "variants_genomebuild": "NR",
        "weight_type": "beta",
        "ancestry_distribution": {
            "eval": {"dist": {"EUR": 80, "NR": 20}, "count": 10},
            "gwas": {"dist": {"EUR": 100}, "count": 22627},
        },
        "date_release": "2019-10-14",
        "license": "PGS obtained from the Catalog should be cited appropriately.",
    }
    PERFORMANCE = {
        "size": 1,
        "count": 1,
        "results": [
            {
                "id": "PPM000001",
                "associated_pgs_id": "PGS000001",
                "phenotyping_reported": "All breast cancer",
                "publication": SCORE["publication"],
                "sampleset": {
                    "id": "PSS000001",
                    "samples": [
                        {
                            "sample_number": 67054,
                            "sample_cases": 33673,
                            "sample_controls": 33381,
                            "phenotyping_free": "All breast cancer",
                            "ancestry_broad": "European",
                            "ancestry_country": "Australia, Canada, Denmark, Finland, France, Germany, UK, USA",
                            "cohorts": [{"name_short": "ABCFS", "name_full": "Australian Breast Cancer Family Study"}],
                        }
                    ],
                },
                "performance_metrics": {
                    "effect_sizes": [
                        {
                            "name_long": "Odds Ratio",
                            "name_short": "OR",
                            "estimate": 1.55,
                            "ci_lower": 1.52,
                            "ci_upper": 1.58,
                        }
                    ],
                    "class_acc": [
                        {
                            "name_long": "Concordance Statistic",
                            "name_short": "C-index",
                            "estimate": 0.622,
                            "ci_lower": 0.619,
                            "ci_upper": 0.627,
                        }
                    ],
                    "othermetrics": [],
                },
                "covariates": None,
                "performance_comments": None,
            }
        ],
    }

    def __init__(self, *, variant_not_found: bool = False, fail_performance: bool = False) -> None:
        self.variant_not_found = variant_not_found
        self.fail_performance = fail_performance
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if "/rest/variant/" in str(url):
            if self.variant_not_found:
                raise KnowledgeRequestError("GET PGS variant failed: 404 Client Error: Not Found")
            return json.loads(json.dumps(self.VARIANT))
        if str(url).endswith("/rest/score/PGS000001"):
            return json.loads(json.dumps(self.SCORE))
        if str(url).endswith("/rest/performance/search"):
            if self.fail_performance:
                raise KnowledgeRequestError("GET PGS performance failed: timed out")
            return json.loads(json.dumps(self.PERFORMANCE))
        raise AssertionError(f"Unexpected PGS Catalog fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingIgsrClient:
    HIGH_LISTING = """
<html><body><table>
<tr><td><img></td><td><a href="CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz">CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz</a></td><td align="right">2020-10-29 16:36  </td><td align="right">1.6G</td><td></td></tr>
<tr><td><img></td><td><a href="CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz.tbi">CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz.tbi</a></td><td align="right">2020-10-29 16:36  </td><td align="right">124K</td><td></td></tr>
</table></body></html>
"""
    PHASE3_LISTING = """
<html><body><table>
<tr><td><img></td><td><a href="ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz">ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz</a></td><td align="right">2021-03-16 15:54  </td><td align="right">701M</td><td></td></tr>
<tr><td><img></td><td><a href="ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi">ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi</a></td><td align="right">2021-03-16 16:00  </td><td align="right">130K</td><td></td></tr>
<tr><td><img></td><td><a href="ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz">ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz</a></td><td align="right">2021-03-16 15:53  </td><td align="right">1.4G</td><td></td></tr>
<tr><td><img></td><td><a href="ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz.tbi">ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz.tbi</a></td><td align="right">2021-03-16 15:56  </td><td align="right">2.3M</td><td></td></tr>
<tr><td><img></td><td><a href="integrated_call_samples_v3.20130502.ALL.panel">integrated_call_samples_v3.20130502.ALL.panel</a></td><td align="right">2014-09-09 12:00  </td><td align="right">54K</td><td></td></tr>
</table></body></html>
"""

    def __init__(self, *, fail_high: bool = False, fail_phase3: bool = False) -> None:
        self.fail_high = fail_high
        self.fail_phase3 = fail_phase3
        self.calls: list[dict[str, object]] = []

    def get_text(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        self.calls.append(
            {
                "url": str(url),
                "params": dict(params or {}),
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if str(url).endswith("/working/20201028_3202_phased/"):
            if self.fail_high:
                raise KnowledgeRequestError("GET IGSR high-coverage listing failed: timed out")
            return self.HIGH_LISTING
        if str(url).endswith("/release/20130502/"):
            if self.fail_phase3:
                raise KnowledgeRequestError("GET IGSR Phase 3 listing failed: timed out")
            return self.PHASE3_LISTING
        raise AssertionError(f"Unexpected IGSR fake GET URL: {url}")

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        raise AssertionError(f"Unexpected fake JSON GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingCivicClient:
    def __init__(self, *, graphql_errors: bool = False, empty_gene: bool = False) -> None:
        self.graphql_errors = graphql_errors
        self.empty_gene = empty_gene
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        raise AssertionError(f"Unexpected fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        self.calls.append(
            {
                "url": url,
                "json_payload": json_payload or {},
                "headers": headers or {},
                "rate_limit_per_second": rate_limit_per_second,
            }
        )
        if self.graphql_errors:
            return {"errors": [{"message": "Cannot query field \"gene\" on type \"Query\"."}]}
        if self.empty_gene:
            return {"data": {"gene": None}}
        return {
            "data": {
                "gene": {
                    "id": "1815",
                    "name": "DRD4",
                    "link": "/features/1815",
                    "description": "Dopamine receptor D4",
                    "featureAliases": ["D4DR"],
                    "variants": {
                        "totalCount": 1,
                        "nodes": [
                            {
                                "id": "101",
                                "name": "DRD4 V194G",
                                "link": "/variants/101",
                                "variantAliases": ["rs1800955"],
                                "variantTypes": [{"id": "1", "name": "missense_variant"}],
                                "hgvsDescriptions": ["NM_000797.4:c.581T>G"],
                                "clinvarIds": ["123"],
                                "alleleRegistryId": "CA123",
                                "maneSelectTranscript": "NM_000797.4",
                                "singleVariantMolecularProfile": {
                                    "id": "201",
                                    "name": "DRD4 V194G",
                                    "link": "/molecular-profiles/201",
                                    "description": "Single variant molecular profile.",
                                    "evidenceItems": {
                                        "totalCount": 2,
                                        "nodes": [
                                            {
                                                "id": "301",
                                                "name": "EID301",
                                                "link": "/evidence/301",
                                                "description": "Example predictive evidence description.",
                                                "descriptionWithNames": (
                                                    "DRD4 V194G predicts response in Example cancer."
                                                ),
                                                "status": "ACCEPTED",
                                                "evidenceType": "PREDICTIVE",
                                                "evidenceLevel": "B",
                                                "evidenceRating": 4,
                                                "evidenceDirection": "SUPPORTS",
                                                "significance": "SENSITIVITYRESPONSE",
                                                "variantOrigin": "SOMATIC",
                                                "variantHgvs": "NM_000797.4:c.581T>G",
                                                "therapyInteractionType": "COMBINATION",
                                                "disease": {
                                                    "id": "401",
                                                    "name": "Example cancer",
                                                    "link": "/diseases/401",
                                                },
                                                "therapies": [
                                                    {
                                                        "id": "501",
                                                        "name": "Example therapy",
                                                        "link": "/therapies/501",
                                                    }
                                                ],
                                                "source": {
                                                    "id": "601",
                                                    "citation": "Example et al. 2026",
                                                    "citationId": "123456",
                                                    "sourceType": "PUBMED",
                                                    "link": "/sources/601",
                                                },
                                                "phenotypes": [
                                                    {
                                                        "id": "701",
                                                        "name": "response",
                                                        "link": "/phenotypes/701",
                                                    }
                                                ],
                                            }
                                        ],
                                    },
                                },
                            }
                        ],
                    },
                }
            }
        }


class RecordingPanelAppClient:
    EXACT_ENTRY = {
        "gene_data": {
            "alias": ["RNF53", "FANCS"],
            "biotype": "protein_coding",
            "hgnc_id": "HGNC:1100",
            "gene_name": "BRCA1, DNA repair associated",
            "omim_gene": ["113705"],
            "gene_symbol": "BRCA1",
            "hgnc_symbol": "BRCA1",
            "ensembl_genes": {
                "GRch37": {"82": {"location": "17:41196312-41277500", "ensembl_id": "ENSG00000012048"}},
                "GRch38": {"90": {"location": "17:43044295-43170245", "ensembl_id": "ENSG00000012048"}},
            },
        },
        "entity_type": "gene",
        "entity_name": "BRCA1",
        "confidence_level": "3",
        "penetrance": "Complete",
        "mode_of_pathogenicity": "",
        "publications": ["29661970"],
        "evidence": ["NHS GMS", "Expert Review Green", "Expert list"],
        "phenotypes": ["{Breast-ovarian cancer, familial, 1}, OMIM:604370"],
        "mode_of_inheritance": "BOTH monoallelic and biallelic, autosomal or pseudoautosomal",
        "tags": ["watchlist_moi"],
        "panel": {
            "id": 143,
            "name": "Inherited ovarian cancer (without breast cancer)",
            "disease_group": "Inherited cancer",
            "disease_sub_group": "",
            "status": "public",
            "version": "5.1",
            "version_created": "2026-05-06T16:02:21.284533Z",
            "relevant_disorders": ["Familial ovarian cancer", "R207"],
            "stats": {"number_of_genes": 27, "number_of_strs": 0, "number_of_regions": 0},
            "types": [
                {"name": "Rare Disease 100K", "slug": "rare-disease-100k"},
                {"name": "GMS signed-off", "slug": "gms-signed-off"},
            ],
        },
        "transcript": ["ENST00000357654.8", "NM_007294.3"],
    }

    def __init__(self, *, empty: bool = False, non_exact: bool = False) -> None:
        self.empty = empty
        self.non_exact = non_exact
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if "panelapp.genomicsengland.co.uk/api/v1/genes" not in str(url):
            raise AssertionError(f"Unexpected PanelApp fake GET URL: {url}")
        if self.empty:
            return {"count": 0, "next": None, "previous": None, "results": []}
        if self.non_exact:
            row = dict(self.EXACT_ENTRY)
            row["entity_name"] = "BRIP1"
            row["gene_data"] = {**self.EXACT_ENTRY["gene_data"], "gene_symbol": "BRIP1", "hgnc_symbol": "BRIP1"}
            return {"count": 1, "next": None, "previous": None, "results": [row]}
        return {"count": 1, "next": None, "previous": None, "results": [self.EXACT_ENTRY]}

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingMaveDbClient:
    BRCA1_RESPONSE = {
        "symbol": "BRCA1",
        "name": "BRCA1 DNA repair associated",
        "hgncId": "HGNC:1100",
        "locusGroup": "protein-coding gene",
        "location": "17q21.31",
        "omimId": "113705",
        "limit": 5,
        "offset": 0,
        "total": 1,
        "totalScoredVariants": 2271,
        "scoreSets": [
            {
                "urn": "urn:mavedb:00001222-b-2",
                "title": "Scores from multiplexed functional assay of BRCA1 variants",
                "shortDescription": (
                    "Multiplexed assay of BRCA1 variants measuring homology directed repair activity"
                ),
                "publishedDate": "2025-10-22",
                "numVariants": 2271,
                "experiment": {
                    "urn": "urn:mavedb:00001222-b",
                    "title": "Multiplexed functional assay of BRCA1 variants",
                    "shortDescription": (
                        "Multiplexed assay of BRCA1 variants measuring homology directed repair activity"
                    ),
                },
                "primaryPublicationIdentifiers": [
                    {"identifier": "39999999", "dbName": "PubMed", "title": "BRCA1 MAVE paper"}
                ],
                "secondaryPublicationIdentifiers": [],
                "license": {
                    "shortName": "CC BY",
                    "longName": "Creative Commons Attribution",
                    "version": "4.0",
                },
                "targetGenes": [
                    {
                        "name": "BRCA1",
                        "category": "protein_coding",
                        "mappedHgncName": "BRCA1",
                        "uniprotIdFromMappedMetadata": "P38398",
                        "externalIdentifiers": [
                            {"identifier": {"dbName": "HGNC", "identifier": "HGNC:1100"}, "offset": 0}
                        ],
                    }
                ],
                "private": False,
                "recordType": "ScoreSet",
            }
        ],
    }
    DRD4_RESPONSE = {
        "symbol": "DRD4",
        "name": "dopamine receptor D4",
        "hgncId": "HGNC:3025",
        "locusGroup": "protein-coding gene",
        "location": "11p15.5",
        "omimId": "126452",
        "scoreSets": [],
        "limit": 5,
        "offset": 0,
        "total": 0,
        "totalScoredVariants": 0,
    }

    def __init__(self, *, empty: bool = False) -> None:
        self.empty = empty
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        self.calls.append(
            {
                "url": str(url),
                "params": dict(params or {}),
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if "api.mavedb.org/api/v1/genes/" not in str(url):
            raise AssertionError(f"Unexpected MaveDB fake GET URL: {url}")
        payload = self.DRD4_RESPONSE if self.empty else self.BRCA1_RESPONSE
        return json.loads(json.dumps(payload))

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingClinVarClient:
    def __init__(
        self,
        *,
        ids_by_term: dict[str, list[str]] | None = None,
        summaries: dict[str, dict[str, object]] | None = None,
        clinical_tables_payload: list[object] | None = None,
    ) -> None:
        self.ids_by_term = ids_by_term or {}
        self.summaries = summaries or {}
        self.clinical_tables_payload = clinical_tables_payload or [0, [], {}, []]
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None):
        params = dict(params or {})
        self.calls.append({"url": url, "params": params})
        if "esearch.fcgi" in url:
            return {"esearchresult": {"idlist": self.ids_by_term.get(str(params.get("term", "")), [])}}
        if "esummary.fcgi" in url:
            return {
                "result": {
                    item_id: self.summaries[item_id]
                    for item_id in str(params.get("id", "")).split(",")
                    if item_id in self.summaries
                }
            }
        if "clinicaltables.nlm.nih.gov" in url:
            return self.clinical_tables_payload
        raise AssertionError(f"Unexpected ClinVar fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingClinGenClient:
    VALIDITY_CSV = (
        '"CLINGEN GENE DISEASE VALIDITY CURATIONS","",""\n'
        '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
        '"+++++++++++","++++","++++","++++","++++","++++","++++","++++","++++","++++"\n'
        '"GENE1","HGNC:1","Example syndrome","MONDO:0000001","AD","SOP10","Definitive","https://search.clinicalgenome.org/kb/gene-validity/GENE1","2026-01-01T00:00:00.000Z","Example GCEP"\n'
    )
    DOSAGE_CSV = (
        '"CLINGEN DOSAGE SENSITIVITY CURATIONS","",""\n'
        '"GENE SYMBOL","HGNC ID","HAPLOINSUFFICIENCY","TRIPLOSENSITIVITY","ONLINE REPORT","DATE"\n'
        '"+++++++++++","++++","++++","++++","++++","++++"\n'
        '"GENE1","HGNC:1","Sufficient Evidence for Haploinsufficiency","No Evidence for Triplosensitivity","https://search.clinicalgenome.org/kb/gene-dosage/HGNC:1","2026-01-02T00:00:00+00:00"\n'
    )
    SUMMARY_CSV = (
        '"README","",""\n'
        '"gene_symbol","hgnc_id","gene_url","disease_label","mondo_id","disease_url","mode_of_inheritance","dosage_haploinsufficiency_assertion","dosage_triplosensitivity_assertion","dosage_report","dosage_group","gene_disease_validity_assertion_classifications","gene_disease_validity_assertion_reports","gene_disease_validity_gceps","actionability_assertion_classifications","actionability_assertion_reports","actionability_groups"\n'
        '"GENE1","HGNC:1","https://search.clinicalgenome.org/kb/genes/HGNC:1","Example syndrome","MONDO:0000001","https://search.clinicalgenome.org/kb/conditions/MONDO:0000001","Autosomal dominant inheritance","Sufficient Evidence for Haploinsufficiency","No Evidence for Triplosensitivity","https://search.clinicalgenome.org/kb/gene-dosage/HGNC:1","Dosage Working Group","Definitive","https://search.clinicalgenome.org/kb/gene-validity/GENE1","Example GCEP","Actionable","https://actionability.clinicalgenome.org/ac/Adult/ui/summ","Example Actionability Group"\n'
    )
    ACTIONABILITY_JSON = {
        "columns": [
            "docId",
            "contextIri",
            "context",
            "releaseDate",
            "geneOrVariant",
            "disease",
            "status-overall",
            "outcome",
            "intervention",
            "severity",
            "likelihood",
            "natureOfIntervention",
            "effectiveness",
            "overall",
        ],
        "rows": [
            [
                "AC1",
                "https://actionability.clinicalgenome.org/ac/Adult/api/sepio/doc/AC1",
                "Adult",
                "Thu, 01 Jan 2026 00:00:00 -0000",
                "GENE1",
                "Example syndrome",
                "Released",
                "Preventable morbidity",
                "Surveillance",
                "2",
                "3C",
                "3",
                "2N",
                "10CN",
            ]
        ],
    }

    def __init__(
        self,
        *,
        empty: bool = False,
        malformed_validity: bool = False,
        fail_all: bool = False,
        summary_timeout: bool = False,
    ) -> None:
        self.empty = empty
        self.malformed_validity = malformed_validity
        self.fail_all = fail_all
        self.summary_timeout = summary_timeout
        self.calls: list[str] = []
        self.timeouts: dict[str, int | float | None] = {}

    def get_text(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        self.calls.append(str(url))
        self.timeouts[str(url)] = timeout
        if self.fail_all:
            raise KnowledgeRequestError(
                "GET ClinGen failed: TLS certificate verification failed.",
                code="tls_certificate_verification_failed",
                remediation="Configure NOPHIGENE_CA_BUNDLE.",
            )
        if "gene-validity/download" in url:
            return "not,a,valid,clingen,file\n" if self.malformed_validity else self._maybe_empty(self.VALIDITY_CSV)
        if "gene-dosage/download" in url:
            return self._maybe_empty(self.DOSAGE_CSV)
        if "curation-activity-summary-report" in url:
            if self.summary_timeout:
                raise KnowledgeRequestError(
                    "GET https://search.clinicalgenome.org/kb/reports/curation-activity-summary-report "
                    "failed: HTTPSConnectionPool(host='search.clinicalgenome.org', port=443): "
                    "Read timed out. (read timeout=20)"
                )
            return self._maybe_empty(self.SUMMARY_CSV)
        raise AssertionError(f"Unexpected ClinGen fake text URL: {url}")

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        self.calls.append(str(url))
        if self.fail_all:
            raise KnowledgeRequestError(
                "GET ClinGen failed: TLS certificate verification failed.",
                code="tls_certificate_verification_failed",
                remediation="Configure NOPHIGENE_CA_BUNDLE.",
            )
        return {"columns": self.ACTIONABILITY_JSON["columns"], "rows": []} if self.empty else self.ACTIONABILITY_JSON

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None, timeout=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")

    def _maybe_empty(self, csv_text: str) -> str:
        if not self.empty:
            return csv_text
        line_count = 2 if csv_text.startswith('"README"') else 3
        return "\n".join(csv_text.splitlines()[:line_count]) + "\n"


class RecordingMedGenClient:
    def __init__(
        self,
        *,
        ids_by_term: dict[str, list[str]] | None = None,
        summaries: dict[str, dict[str, object]] | None = None,
        malformed_concept_meta: bool = False,
        fail_all: bool = False,
    ) -> None:
        self.ids_by_term = ids_by_term or {
            "GENE1[gene]": ["100", "101"],
            '"medgen gtr tests clinical"[Filter] AND GENE1[gene]': ["101"],
        }
        concept_meta = (
            "<GeneSymbol>GENE1</GeneSymbol> OMIM:123456 HP:0000001 Orphanet:12345 "
            "ClinVar Genetic Testing Registry"
        )
        if malformed_concept_meta:
            concept_meta = "<ConceptMeta><Broken>"
        self.summaries = summaries or {
            "100": {
                "uid": "100",
                "title": "Example syndrome",
                "conceptid": "C0000001",
                "definition": "Example syndrome associated with GENE1.",
                "semanticid": "T047",
                "semantictype": "Disease or Syndrome",
                "modificationdate": "2026/01/01",
                "conceptmeta": concept_meta,
            },
            "101": {
                "uid": "101",
                "title": "GENE1 phenotype",
                "conceptid": "CN000002",
                "definition": "",
                "semanticid": "T033",
                "semantictype": "Finding",
                "modificationdate": "2026/01/02",
                "conceptmeta": "GTR",
            },
        }
        self.fail_all = fail_all
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None):
        params = dict(params or {})
        self.calls.append({"url": url, "params": params})
        if self.fail_all:
            raise KnowledgeRequestError(
                "GET MedGen failed: TLS certificate verification failed.",
                code="tls_certificate_verification_failed",
                remediation="Configure NOPHIGENE_CA_BUNDLE.",
            )
        if "esearch.fcgi" in url:
            return {"esearchresult": {"idlist": self.ids_by_term.get(str(params.get("term", "")), [])}}
        if "esummary.fcgi" in url:
            return {
                "result": {
                    item_id: self.summaries[item_id]
                    for item_id in str(params.get("id", "")).split(",")
                    if item_id in self.summaries
                }
            }
        raise AssertionError(f"Unexpected MedGen fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingDbSnpClient:
    ESUMMARY_ITEM = {
        "uid": "2533154733",
        "snp_id": "2533154733",
        "chr": "11",
        "chrpos": "11:637373",
        "spdi": "NC_000011.10:637372:T:",
        "snp_class": "del",
        "fxn_class": "frameshift_variant,coding_sequence_variant",
        "validated": "by-frequency",
        "handle": "GNOMAD",
        "createdate": "2024/01/01",
        "updatedate": "2026/01/01",
        "global_mafs": [{"study": "GnomAD_exomes", "freq": "-=0.000001/1"}],
        "genes": [{"name": "DRD4", "gene_id": "1815"}],
        "docsum": "HGVS=NC_000011.10:g.637373del|SEQ=[T/]|LEN=1|GENE=DRD4:1815",
    }
    REFSNP_PAYLOAD = {
        "refsnp_id": "2533154733",
        "create_date": "2024-01-01T00:00Z",
        "last_update_date": "2026-01-01T00:00Z",
        "citations": [{"pmid": 1}],
        "primary_snapshot_data": {
            "variant_type": "del",
            "mane_select_ids": ["NM_000797.4"],
            "placements_with_allele": [
                {
                    "placement_annot": {
                        "is_ptlp": True,
                        "seq_id_traits_by_assembly": [{"assembly_name": "GRCh38.p14"}],
                    },
                    "alleles": [
                        {
                            "hgvs": "NC_000011.10:g.637373del",
                            "allele": {
                                "spdi": {
                                    "seq_id": "NC_000011.10",
                                    "position": 637372,
                                    "deleted_sequence": "T",
                                    "inserted_sequence": "",
                                }
                            },
                        }
                    ],
                }
            ],
            "allele_annotations": [
                {
                    "frequency": [
                        {
                            "study_name": "GnomAD_exomes",
                            "allele_count": 1,
                            "total_count": 998848,
                        }
                    ],
                    "clinical": [{"clinical_significance": "not provided"}],
                    "assembly_annotation": [
                        {
                            "genes": [
                                {
                                    "name": "DRD4",
                                    "id": 1815,
                                    "rnas": [
                                        {
                                            "accession_version": "NM_000797.4",
                                            "hgvs": "NM_000797.4:c.69del",
                                            "sequence_ontology": [
                                                {"name": "coding_sequence_variant"},
                                                {"name": "frameshift_variant"},
                                            ],
                                            "protein": {
                                                "accession_version": "NP_000788.2",
                                                "hgvs": "NP_000788.2:p.Ala24fs",
                                                "sequence_ontology": [{"name": "frameshift_variant"}],
                                            },
                                        }
                                    ],
                                }
                            ]
                        }
                    ],
                }
            ],
            "support": [{"submitter_handle": "GNOMAD"}],
        },
    }

    def __init__(self, *, fail_refsnp: bool = False) -> None:
        self.fail_refsnp = fail_refsnp
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "timeout": timeout,
            }
        )
        if "esearch.fcgi" in url:
            return {"esearchresult": {"idlist": ["2533154733"]}}
        if "esummary.fcgi" in url:
            return {"result": {"2533154733": dict(self.ESUMMARY_ITEM)}}
        if "variation/v0/refsnp/2533154733" in url:
            if self.fail_refsnp:
                raise KnowledgeRequestError("GET dbSNP RefSNP failed: timed out")
            return self.REFSNP_PAYLOAD
        raise AssertionError(f"Unexpected dbSNP fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingEnsemblClient:
    def __init__(self, *, fail_variation: bool = False, fail_vep: bool = False) -> None:
        self.fail_variation = fail_variation
        self.fail_vep = fail_vep
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None):
        params = dict(params or {})
        self.calls.append({"url": str(url), "params": params, "headers": dict(headers or {})})
        if "lookup/symbol/homo_sapiens/DRD4" in url:
            return {
                "id": "ENSG00000069696",
                "display_name": "DRD4",
                "seq_region_name": "11",
                "start": 637293,
                "end": 640706,
                "assembly_name": "GRCh37",
                "canonical_transcript": "ENST00000176183.5",
                "biotype": "protein_coding",
            }
        if "overlap/region/homo_sapiens/11:637293-637293" in url:
            return [
                {
                    "id": "rs927984495",
                    "start": 637293,
                    "end": 637293,
                    "seq_region_name": "11",
                    "assembly_name": "GRCh37",
                    "alleles": ["C", "T"],
                    "source": "dbSNP",
                    "consequence_type": "intergenic_variant",
                    "feature_type": "variation",
                },
                {
                    "id": "rs927984495",
                    "start": 637293,
                    "end": 637293,
                    "seq_region_name": "11",
                    "assembly_name": "GRCh37",
                    "alleles": ["C", "T"],
                    "source": "dbSNP",
                    "consequence_type": "intergenic_variant",
                    "feature_type": "variation",
                },
            ]
        if "variation/homo_sapiens/rs927984495" in url:
            if self.fail_variation:
                raise KnowledgeRequestError("GET Ensembl variation failed: timed out")
            return {
                "name": "rs927984495",
                "var_class": "SNP",
                "most_severe_consequence": "5_prime_UTR_variant",
                "mappings": [
                    {
                        "assembly_name": "GRCh37",
                        "location": "11:637293-637293",
                        "allele_string": "C/T",
                        "seq_region_name": "11",
                        "start": 637293,
                        "end": 637293,
                    }
                ],
                "clinical_significance": [],
                "evidence": ["Frequency", "TOPMed", "gnomAD"],
                "synonyms": ["example-synonym"],
                "minor_allele": None,
                "MAF": None,
                "phenotypes": [],
                "source": "Variants (including SNPs and indels) imported from dbSNP",
            }
        if "vep/homo_sapiens/region/11:637293:637293/T" in url:
            if self.fail_vep:
                raise KnowledgeRequestError("GET Ensembl VEP failed: timed out")
            return [
                {
                    "id": "11_637293_C/T",
                    "assembly_name": "GRCh37",
                    "most_severe_consequence": "5_prime_UTR_variant",
                    "variant_class": "SNV",
                    "colocated_variants": [
                        {
                            "id": "rs927984495",
                            "allele_string": "C/T",
                            "seq_region_name": "11",
                            "start": 637293,
                            "end": 637293,
                        }
                    ],
                    "transcript_consequences": [
                        {
                            "gene_symbol": "DRD4",
                            "gene_id": "ENSG00000069696",
                            "transcript_id": "ENST00000176183",
                            "canonical": 1,
                            "exon": "1/4",
                            "cdna_start": 1,
                            "cdna_end": 1,
                            "consequence_terms": ["5_prime_UTR_variant"],
                            "impact": "MODIFIER",
                            "biotype": "protein_coding",
                        }
                    ],
                }
            ]
        raise AssertionError(f"Unexpected Ensembl fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingUcscClient:
    def __init__(self, *, fail_track: str = "") -> None:
        self.fail_track = fail_track
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if str(url).endswith("/getData/sequence"):
            assert params["start"] == 637292
            assert params["end"] == 637293
            return {"downloadTime": "2026:06:23T10:00:00Z", "dna": "C"}
        if str(url).endswith("/search"):
            return {
                "positionMatches": [
                    {
                        "matches": [
                            {
                                "posName": "DRD4",
                                "position": "chr11:637293-640706",
                            }
                        ]
                    }
                ]
            }
        if not str(url).endswith("/getData/track"):
            raise AssertionError(f"Unexpected UCSC fake GET URL: {url}")

        track = params.get("track")
        if self.fail_track and track == self.fail_track:
            raise KnowledgeRequestError(f"GET UCSC {track} failed: timed out")
        if track in {"ncbiRefSeq", "cpgIslandExt", "encodeCcreCombined", "encRegTfbsClustered", "rmsk"}:
            assert params["start"] == 637292
            assert params["end"] == 640706
        if track == "ncbiRefSeq":
            return {
                "ncbiRefSeq": [
                    {
                        "name": "NM_000797.4",
                        "name2": "DRD4",
                        "chrom": "chr11",
                        "strand": "+",
                        "txStart": 637268,
                        "txEnd": 640706,
                        "cdsStart": 637304,
                        "cdsEnd": 640603,
                        "exonCount": 4,
                        "exonStarts": "637268,639432,639647,640400,",
                        "exonEnds": "637589,639545,640306,640706,",
                        "cdsStartStat": "cmpl",
                        "cdsEndStat": "cmpl",
                    }
                ]
            }
        if track == "cpgIslandExt":
            return {
                "cpgIslandExt": [
                    {
                        "name": "CpG: 313",
                        "chrom": "chr11",
                        "chromStart": 636435,
                        "chromEnd": 640628,
                        "length": 4193,
                        "cpgNum": 313,
                        "gcNum": 2277,
                        "perCpg": 14.9,
                        "perGc": 54.3,
                        "obsExp": 1.01,
                    }
                ]
            }
        if track == "encodeCcreCombined":
            return {
                "encodeCcreCombined": [
                    {
                        "name": "EH38E1513769",
                        "chrom": "chr11",
                        "chromStart": 637042,
                        "chromEnd": 637385,
                        "ccre": "PLS,CTCF-bound",
                        "encodeLabel": "PLS",
                        "ucscLabel": "prom",
                        "accessionLabel": "E1513769",
                        "description": "EH38E1513769 promoter-like signature",
                        "zScore": 2.69406,
                    },
                    {
                        "name": "EH38E1513770",
                        "chrom": "chr11",
                        "chromStart": 637540,
                        "chromEnd": 637693,
                        "ccre": "pELS,CTCF-bound",
                        "encodeLabel": "pELS",
                        "ucscLabel": "enhP",
                        "accessionLabel": "E1513770",
                        "description": "EH38E1513770 proximal enhancer-like signature",
                        "zScore": 2.08296,
                    },
                ]
            }
        if track == "encRegTfbsClustered":
            return {
                "maxItemsLimit": True,
                "encRegTfbsClustered": [
                    {"name": "POLR2G", "chrom": "chr11", "chromStart": 636191, "chromEnd": 637340, "score": 184, "sourceCount": 2},
                    {"name": "EZH2", "chrom": "chr11", "chromStart": 636719, "chromEnd": 637361, "score": 564, "sourceCount": 2},
                    {"name": "KDM4A", "chrom": "chr11", "chromStart": 636925, "chromEnd": 637642, "score": 857, "sourceCount": 1},
                    {"name": "EZH2", "chrom": "chr11", "chromStart": 637387, "chromEnd": 639433, "score": 1000, "sourceCount": 4},
                ],
            }
        if track == "rmsk":
            return {
                "maxItemsLimit": True,
                "rmsk": [
                    {"repName": "G-rich", "chrom": "chr11", "chromStart": 637321, "chromEnd": 637431, "repClass": "Low_complexity", "repFamily": "Low_complexity"},
                    {"repName": "(AC)n", "chrom": "chr11", "chromStart": 637986, "chromEnd": 638032, "repClass": "Simple_repeat", "repFamily": "Simple_repeat"},
                    {"repName": "MIR", "chrom": "chr11", "chromStart": 638516, "chromEnd": 638693, "repClass": "SINE", "repFamily": "MIR"},
                ],
            }
        if track == "snp151Common":
            if params["start"] == 637292 and params["end"] == 637293:
                return {"snp151Common": []}
            assert params["start"] == 637292
            assert params["end"] == 640706
            return {
                "snp151Common": [
                    {
                        "name": "rs146680769",
                        "chrom": "chr11",
                        "chromStart": 637293,
                        "chromEnd": 637294,
                        "observed": "C/T",
                        "class": "single",
                        "valid": "by-frequency,by-1000genomes",
                    }
                ]
            }
        raise AssertionError(f"Unexpected UCSC fake track: {track}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingEncodeClient:
    def __init__(self, *, empty: bool = False, fail_region: bool = False) -> None:
        self.empty = empty
        self.fail_region = fail_region
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if str(url).endswith("/region-search/"):
            if self.fail_region:
                raise KnowledgeRequestError("GET ENCODE region search failed: timed out")
            if self.empty:
                return {
                    "@graph": [],
                    "total": 0,
                    "notification": "Error during search",
                    "assembly": "GRCh38",
                }
            return {
                "@graph": [
                    {
                        "@id": "/experiments/ENCSRREGION1/",
                        "accession": "ENCSRREGION1",
                        "assay_title": "DNase-seq",
                        "biosample_summary": "Homo sapiens K562 cell line",
                        "status": "released",
                    }
                ],
                "total": 1,
                "notification": "Success",
                "assembly": "GRCh38",
            }
        if not str(url).endswith("/search/"):
            raise AssertionError(f"Unexpected ENCODE fake GET URL: {url}")
        if params.get("target.genes.symbol") == "DRD4":
            if self.empty:
                return {"@graph": [], "total": 0, "notification": "No results found"}
            return {
                "@graph": [
                    {
                        "@id": "/experiments/ENCSRDRD4TF/",
                        "accession": "ENCSRDRD4TF",
                        "assay_title": "TF ChIP-seq",
                        "assay_term_name": "ChIP-seq",
                        "target": {
                            "label": "DRD4",
                            "genes": [{"symbol": "DRD4"}],
                        },
                        "biosample_summary": "Homo sapiens K562 cell line",
                        "biosample_ontology": {"term_name": "K562", "classification": "cell line"},
                        "replicates": [
                            {
                                "biological_replicate_number": 1,
                                "technical_replicate_number": 1,
                                "library": {
                                    "biosample": {
                                        "organism": {"scientific_name": "Homo sapiens"},
                                    }
                                },
                            },
                            {
                                "biological_replicate_number": 2,
                                "technical_replicate_number": 1,
                                "library": {
                                    "biosample": {
                                        "organism": {"scientific_name": "Homo sapiens"},
                                    }
                                },
                            },
                        ],
                        "files": [
                            {
                                "status": "released",
                                "assembly": "GRCh38",
                                "output_type": "optimal IDR thresholded peaks",
                                "file_format": "bed narrowPeak",
                            },
                            {
                                "status": "released",
                                "assembly": "GRCh38",
                                "output_type": "signal p-value",
                                "file_format": "bigWig",
                            },
                        ],
                        "lab": {"title": "ENCODE Processing Pipeline"},
                        "award": {"project": "ENCODE"},
                        "status": "released",
                        "date_released": "2026-01-15",
                        "replication_type": "isogenic",
                    }
                ],
                "total": 1,
                "notification": "Success",
            }
        if params.get("searchTerm") == "DRD4":
            if self.empty:
                return {"@graph": [], "total": 0, "notification": "No results found"}
            return {
                "@graph": [
                    {
                        "@id": "/experiments/ENCSRDRD4TF/",
                        "accession": "ENCSRDRD4TF",
                        "assay_title": "TF ChIP-seq",
                    },
                    {
                        "@id": "/experiments/ENCSRDRD4RNA/",
                        "accession": "ENCSRDRD4RNA",
                        "assay_title": "total RNA-seq",
                        "biosample_summary": "Homo sapiens neural progenitor cell",
                        "biosample_ontology": {"term_name": "neural progenitor cell", "classification": "cell line"},
                        "replicates": [
                            {
                                "biological_replicate_number": 1,
                                "library": {
                                    "biosample": {
                                        "organism": {"scientific_name": "Homo sapiens"},
                                    }
                                },
                            }
                        ],
                        "files": [{"status": "released", "assembly": "GRCh38", "output_type": "gene quantifications"}],
                        "lab": {"title": "Thomas Gingeras, CSHL"},
                        "award": {"project": "ENCODE"},
                        "status": "released",
                    },
                ],
                "total": 2,
                "notification": "Success",
            }
        raise AssertionError(f"Unexpected ENCODE fake params: {params}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingScreenClient:
    def __init__(self, *, fail_ccre: bool = False, graphql_errors: bool = False) -> None:
        self.fail_ccre = fail_ccre
        self.graphql_errors = graphql_errors
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        raise AssertionError(f"Unexpected fake GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        payload = dict(json_payload or {})
        self.calls.append(
            {
                "url": str(url),
                "json_payload": payload,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
            }
        )
        if not str(url).endswith("/api/screen-graphql"):
            raise AssertionError(f"Unexpected SCREEN fake POST URL: {url}")
        if self.fail_ccre:
            raise KnowledgeRequestError("POST SCREEN cCRE query failed: 500 Internal Server Error")
        if self.graphql_errors:
            return {
                "errors": [{"message": "temporary SCREEN GraphQL validation issue"}],
                "data": {"cCREQuery": []},
            }
        return {
            "data": {
                "cCREQuery": [
                    {
                        "accession": "EH38E2937824",
                        "assembly": "grch38",
                        "rDHS": "EH38D4573084",
                        "group": "PLS",
                        "ctcf_bound": True,
                        "dnaseMax": 3.523491859436035,
                        "h3k4me3Max": 3.9978396892547607,
                        "h3k27acMax": 2.0321149826049805,
                        "ctcfMax": 1.881090521812439,
                        "coordinates": {"chromosome": "chr11", "start": 637025, "end": 637375},
                        "nearby_genes": {
                            "intersecting_genes": [
                                {
                                    "id": "ENSG00000069696.7",
                                    "name": "DRD4",
                                    "gene_type": "protein_coding",
                                    "strand": "+",
                                },
                                {
                                    "id": "ENSG00000069696.7",
                                    "name": "DRD4",
                                    "gene_type": "protein_coding",
                                    "strand": "+",
                                },
                            ]
                        },
                    },
                    {
                        "accession": "EH38E2937825",
                        "assembly": "grch38",
                        "rDHS": "EH38D4573085",
                        "group": "pELS",
                        "ctcf_bound": True,
                        "dnaseMax": 2.406550407409668,
                        "h3k4me3Max": 4.400539875030518,
                        "h3k27acMax": 3.127086877822876,
                        "ctcfMax": 1.8611128330230713,
                        "coordinates": {"chromosome": "chr11", "start": 637523, "end": 637724},
                        "nearby_genes": {
                            "intersecting_genes": [
                                {
                                    "id": "ENSG00000069696.7",
                                    "name": "DRD4",
                                    "gene_type": "protein_coding",
                                    "strand": "+",
                                }
                            ]
                        },
                    },
                ]
            }
        }


class RecordingEwasCatalogClient:
    def __init__(self, *, fail_tsv: bool = False) -> None:
        self.fail_tsv = fail_tsv
        self.calls: list[dict[str, object]] = []

    def get_text(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if str(url) == "https://www.ewascatalog.org/":
            assert params == {"gene": "DRD4"}
            return (
                '<html><body><a href="/static//tmp/DRD4_test.tsv" download="DRD4.tsv">'
                "Download</a></body></html>"
            )
        if str(url).endswith("/static//tmp/DRD4_test.tsv"):
            if self.fail_tsv:
                raise KnowledgeRequestError("GET EWAS Catalog TSV failed: timed out")
            return (
                "Author\tPMID\tOutcome\tExposure\tTissue\tAnalysis\tN\tCpG\tLocation\tGene\tBeta\tP\tStudyID\n"
                "Mulder RH\t33450751\tDNA methylation\tage\tWhole blood\tModel 1 with age as a fixed effect\t"
                "2338\tcg02762115\tchr11:640446\tDRD4\t0.0021\t0E+00\tEWAS0001\n"
                "Example AB\t12345678\tDNA methylation\tcognitive function\tBrain\tAdjusted linear model\t"
                "415\tcg15861585\tchr11:637038\tDRD4\t-0.032\t9.8E-05\tEWAS0002\n"
            )
        raise AssertionError(f"Unexpected EWAS Catalog fake GET URL: {url}")

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        raise AssertionError(f"Unexpected fake JSON GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


class RecordingEwasAtlasClient:
    def __init__(self, *, empty: bool = False) -> None:
        self.empty = empty
        self.calls: list[dict[str, object]] = []

    def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        params = dict(params or {})
        self.calls.append(
            {
                "url": str(url),
                "params": params,
                "headers": dict(headers or {}),
                "rate_limit_per_second": rate_limit_per_second,
                "timeout": timeout,
            }
        )
        if not str(url).endswith("/ewas/rest/pos"):
            raise AssertionError(f"Unexpected EWAS Atlas fake GET URL: {url}")
        assert params == {"chr": "11", "start": 637293, "end": 640706}
        if self.empty:
            return {"code": 0, "msg": "success", "data": []}
        return {
            "code": 0,
            "msg": "success",
            "data": [
                {
                    "probeId": "cg02762115",
                    "chrHg19": 11,
                    "posHg19": 640446,
                    "cpgIsland": "Shelf",
                    "relatedTranscription": [
                        {
                            "geneName": "DRD4",
                            "posToTss": 3153,
                            "ensemblTranscriptId": "ENST00000176183.5",
                        }
                    ],
                    "associationList": [
                        {
                            "studyId": "ES00218",
                            "trait": "maternal alcohol consumption",
                            "correlation": "neg",
                            "rank": 196,
                            "pmid": 29172695,
                        },
                        {
                            "studyId": "ES00743",
                            "trait": "cognitive function",
                            "correlation": "neg",
                            "rank": 4,
                            "pmid": 29311653,
                        },
                    ],
                },
                {
                    "probeId": "cg03909863",
                    "chrHg19": 11,
                    "posHg19": 638404,
                    "cpgIsland": "Shore",
                    "relatedTranscription": [
                        {
                            "geneName": "DRD4",
                            "posToTss": 1111,
                            "ensemblTranscriptId": "ENST00000176183.5",
                        }
                    ],
                    "associationList": [
                        {
                            "studyId": "ES00161",
                            "trait": "ascending aortic dissection (AD)",
                            "correlation": "pos",
                            "rank": 98,
                            "pmid": 28444195,
                        }
                    ],
                },
            ],
        }

    def get_text(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
        raise AssertionError(f"Unexpected fake text GET URL: {url}")

    def post_json(self, url, *, json_payload=None, headers=None, rate_limit_per_second=None):
        raise AssertionError(f"Unexpected fake POST URL: {url}")


def _write_simple_pdf(path: Path, text: str) -> None:
    escaped = text.replace("\\", "\\\\").replace("(", r"\(").replace(")", r"\)")
    stream = f"BT /F1 12 Tf 72 720 Td ({escaped}) Tj ET"
    payload = (
        "%PDF-1.4\n"
        "1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\n"
        "2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\n"
        "3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Contents 4 0 R >> endobj\n"
        f"4 0 obj << /Length {len(stream)} >> stream\n"
        f"{stream}\n"
        "endstream endobj\n"
        "trailer << /Root 1 0 R >>\n"
        "%%EOF\n"
    )
    path.write_bytes(payload.encode("latin-1", errors="ignore"))


def test_registry_has_one_card_per_resource_entry():
    active_resource_lines = [
        line
        for line in RESOURCES_PATH.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    specs = list_source_specs()

    assert len(specs) == len(active_resource_lines)
    assert len({spec.key for spec in specs}) == len(active_resource_lines)
    assert all(spec.name for spec in specs)
    assert all(spec.connector_kind for spec in specs)


def test_non_open_sources_advertise_compliant_ingestion_modes_only():
    target_keys = {
        "google_scholar",
        "hgmd",
        "genecards",
        "varsome",
        "franklin",
        "mastermind",
        "embase",
        "scopus",
        "web_of_science",
    }
    specs = {spec.key: spec for spec in list_source_specs()}
    cards = {card["key"]: card for card in list_source_cards()}

    for key in target_keys:
        spec = specs[key]
        assert "scrap" not in spec.connector_kind.lower()
        assert "user_export" in spec.ingestion_modes
        assert "linkout_only" in spec.ingestion_modes
        assert set(spec.ingestion_modes) <= {"official_api", "user_export", "linkout_only"}
        assert cards[key]["import_schema"]
        assert cards[key]["accepted_import_formats"] == ["csv", "json"]


def test_workflow_registry_references_valid_sources_and_core_defaults():
    source_keys = {spec.key for spec in list_source_specs()}
    synthetic_local_sources = {LOCAL_ARTICLE_SOURCE_KEY}
    workflows = list_workflow_specs()
    medgen_spec = get_source_spec("medgen")
    string_spec = get_source_spec("string")

    assert [workflow.key for workflow in workflows if workflow.default_selected] == list(CORE_SAFETY_WORKFLOW_KEYS)
    assert medgen_spec is not None
    assert medgen_spec.access_type == "open_api"
    assert medgen_spec.connector_kind == "medgen"
    assert medgen_spec.ingestion_modes == ("official_api", "linkout_only")
    assert string_spec is not None
    assert string_spec.access_type == "open_api"
    assert string_spec.connector_kind == "string"
    assert string_spec.lane == "interactions"
    assert workflows
    for workflow in workflows:
        assert workflow.label
        assert workflow.purpose
        assert workflow.report_section
        assert set(workflow.ordered_source_keys) <= source_keys | synthetic_local_sources
        if workflow.key == "clinical_variant_triage":
            assert workflow.ordered_source_keys[:4] == ("clinvar", "clingen", "medgen", "ensembl")
        if workflow.key == "gene_interaction_network":
            assert workflow.ordered_source_keys == ("string",)
            assert workflow.default_selected
            assert not workflow.requires_vcf
            assert not workflow.requires_manifest
        if workflow.key == LOCAL_ARTICLE_WORKFLOW_KEY:
            assert workflow.ordered_source_keys == (LOCAL_ARTICLE_SOURCE_KEY,)
            assert not workflow.default_selected
        if workflow.key == "licensed_aggregator_review":
            joined_notes = " ".join(workflow.licensed_notes).lower()
            assert "scraping" in joined_notes
            assert "captcha" in joined_notes


def test_string_connector_returns_only_direct_query_gene_functional_associations():
    class RecordingStringClient:
        def __init__(self) -> None:
            self.calls = []

        def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None, timeout=None):
            self.calls.append(
                {
                    "url": url,
                    "params": params,
                    "headers": headers,
                    "rate_limit_per_second": rate_limit_per_second,
                }
            )
            return [
                {
                    "stringId_A": "9606.ENSP_DRD4",
                    "stringId_B": "9606.ENSP_SLC6A4",
                    "preferredName_A": "DRD4",
                    "preferredName_B": "SLC6A4",
                    "score": 0.978,
                    "ascore": 0.054,
                    "escore": 0.0,
                    "dscore": 0.0,
                    "tscore": 0.978,
                },
                {
                    "stringId_A": "9606.ENSP_DRD3",
                    "stringId_B": "9606.ENSP_DRD4",
                    "preferredName_A": "DRD3",
                    "preferredName_B": "DRD4",
                    "score": 0.943,
                    "dscore": 0.9,
                    "tscore": 0.454,
                },
                {
                    "stringId_A": "9606.ENSP_DRD3",
                    "stringId_B": "9606.ENSP_SLC6A4",
                    "preferredName_A": "DRD3",
                    "preferredName_B": "SLC6A4",
                    "score": 0.7,
                },
            ]

    spec = get_source_spec("string")
    assert spec is not None
    client = RecordingStringClient()
    result = connector_for(spec, client, ResolvedCredential("string")).query(
        KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg38")
    )

    assert result.status == "ok"
    assert [record["partner_gene"] for record in result.records] == ["SLC6A4", "DRD3"]
    assert all(record["edge_type"] == "functional_association" for record in result.records)
    assert all(record["source_release"] == "12.0" for record in result.records)
    assert result.records[0]["evidence_channels"]["text_mining"] == 0.978
    assert "not necessarily direct physical binding" in result.records[0]["association_scope"]
    assert client.calls == [
        {
            "url": "https://version-12-0.string-db.org/api/json/network",
            "params": {
                "identifiers": "DRD4",
                "species": 9606,
                "required_score": 400,
                "add_nodes": 10,
                "caller_identity": "NophiGene",
            },
            "headers": {"Accept": "application/json"},
            "rate_limit_per_second": 1.0,
        }
    ]


def test_request_client_prefers_explicit_ca_bundle(monkeypatch, tmp_path: Path):
    ca_bundle = tmp_path / "custom-ca.pem"
    ca_bundle.write_text("", encoding="ascii")
    monkeypatch.setenv("NOPHIGENE_CA_BUNDLE", str(ca_bundle))
    monkeypatch.setenv("REQUESTS_CA_BUNDLE", str(tmp_path / "other-ca.pem"))

    client = RequestClient()

    assert client.verify == str(ca_bundle)


def test_request_client_uses_windows_merged_ca_bundle(monkeypatch, tmp_path: Path):
    ca_bundle = tmp_path / "windows-merged-ca.pem"
    ca_bundle.write_text("", encoding="ascii")
    for env_var in request_client_module.EXPLICIT_CA_ENV_VARS:
        monkeypatch.delenv(env_var, raising=False)
    monkeypatch.setattr(request_client_module.os, "name", "nt")
    monkeypatch.setattr(request_client_module, "_windows_merged_ca_bundle", lambda: str(ca_bundle))

    client = RequestClient()

    assert client.verify == str(ca_bundle)


def test_request_client_ssl_errors_are_redacted_and_actionable():
    client = RequestClient()

    def fail_get(url, *, params=None, headers=None, timeout=None, verify=None):
        raise requests.exceptions.SSLError(
            "certificate verify failed for https://example.test/?api_key=SECRET&token=TOKEN"
        )

    client.session.get = fail_get

    with pytest.raises(KnowledgeRequestError) as exc_info:
        client.get_json("https://example.test/endpoint", params={"api_key": "SECRET"})

    message = str(exc_info.value)
    assert exc_info.value.code == "tls_certificate_verification_failed"
    assert "NOPHIGENE_CA_BUNDLE" in message
    assert "SECRET" not in message
    assert "TOKEN" not in message
    assert "api_key=[redacted]" in message
    assert "token=[redacted]" in message


def test_clinvar_connector_uses_fielded_gene_query_and_ncbi_params(monkeypatch):
    monkeypatch.setenv("NOPHIGENE_NCBI_EMAIL", "researcher@example.org")
    monkeypatch.setenv("NOPHIGENE_NCBI_API_KEY", "secret-ncbi-key")
    client = RecordingClinVarClient(
        ids_by_term={"DRD4[gene] AND single_gene[prop]": ["123"]},
        summaries={
            "123": {
                "uid": "123",
                "title": "DRD4 variant",
                "clinical_significance": {"description": "Pathogenic"},
                "trait_set": [{"trait_name": "Example condition"}],
                "variation_set": [
                    {
                        "variation_name": "DRD4 c.1A>G",
                        "variation_xrefs": [{"db_source": "dbSNP", "db_id": "1800955"}],
                    }
                ],
            }
        },
    )
    connector = connector_for(get_source_spec("clinvar"), client, ResolvedCredential("clinvar"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:1-10", genome_build="hg19"))

    search_calls = [call for call in client.calls if "esearch.fcgi" in str(call["url"])]
    assert search_calls[0]["params"]["term"] == "DRD4[gene] AND single_gene[prop]"
    assert search_calls[0]["params"]["tool"] == "NophiGeneDynamicKB"
    assert search_calls[0]["params"]["email"] == "researcher@example.org"
    assert search_calls[0]["params"]["api_key"] == "secret-ncbi-key"
    assert result.status == "ok"
    assert result.records[0]["clinical_significance"] == "Pathogenic"
    assert result.records[0]["rsid"] == "rs1800955"
    assert result.records[0]["url"] == "https://www.ncbi.nlm.nih.gov/clinvar/variation/123/"
    serialized = json.dumps(result.to_status(get_source_spec("clinvar")), sort_keys=True)
    assert "secret-ncbi-key" not in serialized


def test_clinvar_connector_queries_rsid_first_and_deduplicates_gene_results():
    client = RecordingClinVarClient(
        ids_by_term={
            "rs123": ["123"],
            "DRD4[gene] AND single_gene[prop]": ["123"],
        },
        summaries={
            "123": {
                "uid": "123",
                "title": "ClinVar rs123 record",
                "description": "Clinical assertion",
            }
        },
    )
    connector = connector_for(get_source_spec("clinvar"), client, ResolvedCredential("clinvar"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:1-10",
        genome_build="hg19",
        variants=(QueryVariant(chrom="11", pos=1, rsid="rs123"),),
    )

    result = connector.query(query)

    search_terms = [call["params"]["term"] for call in client.calls if "esearch.fcgi" in str(call["url"])]
    assert search_terms[:2] == ["rs123", "DRD4[gene] AND single_gene[prop]"]
    assert len(result.records) == 1
    assert result.records[0]["variant"] == "rs123"


def test_clinvar_connector_falls_back_from_single_gene_to_gene_query():
    client = RecordingClinVarClient(
        ids_by_term={
            "DRD4[gene] AND single_gene[prop]": [],
            "DRD4[gene]": ["456"],
        },
        summaries={"456": {"uid": "456", "title": "ClinVar broad gene record"}},
    )
    connector = connector_for(get_source_spec("clinvar"), client, ResolvedCredential("clinvar"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:1-10", genome_build="hg19"))

    search_terms = [call["params"]["term"] for call in client.calls if "esearch.fcgi" in str(call["url"])]
    assert search_terms[:2] == ["DRD4[gene] AND single_gene[prop]", "DRD4[gene]"]
    assert result.status == "ok"
    assert result.records[0]["source_id"] == "456"


def test_clinvar_connector_returns_ok_for_zero_results():
    client = RecordingClinVarClient()
    connector = connector_for(get_source_spec("clinvar"), client, ResolvedCredential("clinvar"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:1-10", genome_build="hg19"))

    assert result.status == "ok"
    assert result.records == []
    assert "0 clinical variant record" in result.message


def test_clinvar_connector_reports_tls_error_without_traceback():
    class FailingClient:
        def get_json(self, url, *, params=None, headers=None, rate_limit_per_second=None):
            raise KnowledgeRequestError(
                "GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi failed: "
                "TLS certificate verification failed. Configure NOPHIGENE_CA_BUNDLE.",
                code="tls_certificate_verification_failed",
                remediation="Configure NOPHIGENE_CA_BUNDLE.",
            )

    connector = connector_for(get_source_spec("clinvar"), FailingClient(), ResolvedCredential("clinvar"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:1-10", genome_build="hg19"))
    status = result.to_status(get_source_spec("clinvar"))

    assert result.status == "failed"
    assert status["error_code"] == "tls_certificate_verification_failed"
    assert status["remediation"] == "Configure NOPHIGENE_CA_BUNDLE."
    assert "Traceback" not in json.dumps(status)


def test_medgen_connector_queries_gene_and_gtr_filter_with_ncbi_params(monkeypatch):
    monkeypatch.setenv("NOPHIGENE_NCBI_EMAIL", "researcher@example.org")
    monkeypatch.setenv("NOPHIGENE_NCBI_API_KEY", "secret-ncbi-key")
    client = RecordingMedGenClient()
    connector = connector_for(get_source_spec("medgen"), client, ResolvedCredential("medgen"))

    result = connector.query(KnowledgeQuery(gene="GENE1", region="1:1-10", genome_build="hg19"))

    search_calls = [call for call in client.calls if "esearch.fcgi" in str(call["url"])]
    summary_calls = [call for call in client.calls if "esummary.fcgi" in str(call["url"])]
    assert [call["params"]["term"] for call in search_calls] == [
        "GENE1[gene]",
        '"medgen gtr tests clinical"[Filter] AND GENE1[gene]',
    ]
    assert search_calls[0]["params"]["tool"] == "NophiGeneDynamicKB"
    assert search_calls[0]["params"]["email"] == "researcher@example.org"
    assert search_calls[0]["params"]["api_key"] == "secret-ncbi-key"
    assert len(summary_calls) == 1
    assert summary_calls[0]["params"]["id"] == "100,101"
    assert result.status == "ok"
    assert len(result.records) == 2
    first = result.records[0]
    assert first["category"] == "clinical_condition"
    assert first["concept_id"] == "C0000001"
    assert first["semantic_type"] == "Disease or Syndrome"
    assert first["url"] == "https://www.ncbi.nlm.nih.gov/medgen/C0000001"
    assert first["omim_ids"] == ["123456"]
    assert first["hpo_ids"] == ["HP:0000001"]
    assert first["orphanet_ids"] == ["orphanet_12345"]
    assert first["related_genes"] == ["GENE1"]
    assert first["has_clinvar"] is True
    assert first["has_gtr"] is True
    assert any(link["url"].endswith("#Additional_description") for link in first["research_links"])
    assert result.records[1]["query_contexts"] == ["gene", "gtr_clinical_tests"]
    serialized = json.dumps(result.to_status(get_source_spec("medgen")), sort_keys=True)
    assert "secret-ncbi-key" not in serialized


def test_medgen_connector_returns_ok_for_zero_results():
    client = RecordingMedGenClient(ids_by_term={})
    connector = connector_for(get_source_spec("medgen"), client, ResolvedCredential("medgen"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:1-10", genome_build="hg19"))

    assert result.status == "ok"
    assert result.records == []
    assert "no MedGen records found for DRD4" in result.message


def test_medgen_connector_keeps_records_when_conceptmeta_is_malformed():
    client = RecordingMedGenClient(malformed_concept_meta=True)
    connector = connector_for(get_source_spec("medgen"), client, ResolvedCredential("medgen"))

    result = connector.query(KnowledgeQuery(gene="GENE1", region="1:1-10", genome_build="hg19"))

    assert result.status == "ok"
    assert result.records[0]["title"] == "Example syndrome"
    assert any("MedGen ConceptMeta for UID 100 could not be fully parsed" in warning for warning in result.warnings)


def test_medgen_connector_reports_request_failure_without_traceback():
    client = RecordingMedGenClient(fail_all=True)
    connector = connector_for(get_source_spec("medgen"), client, ResolvedCredential("medgen"))

    result = connector.query(KnowledgeQuery(gene="GENE1", region="1:1-10", genome_build="hg19"))
    status = result.to_status(get_source_spec("medgen"))

    assert result.status == "failed"
    assert status["error_code"] == "tls_certificate_verification_failed"
    assert status["remediation"] == "Configure NOPHIGENE_CA_BUNDLE."
    assert "Traceback" not in json.dumps(status)


def test_gnomad_connector_queries_variant_frequency_and_gene_constraint_details():
    client = RecordingGnomadClient()
    connector = connector_for(get_source_spec("gnomad"), client, ResolvedCredential("gnomad"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert len(result.records) == 2
    call = client.calls[0]
    assert call["url"] == "https://gnomad.broadinstitute.org/api/"
    assert call["headers"]["Content-Type"] == "application/json"
    assert call["json_payload"]["variables"] == {
        "dataset": "gnomad_r4",
        "geneSymbol": "DRD4",
        "referenceGenome": "GRCh38",
        "variantId": "11-637293-C-T",
        "rsid": None,
    }
    assert "variant(variantId: $variantId, rsid: $rsid, dataset: $dataset)" in call["json_payload"]["query"]
    assert "gnomad_constraint" in call["json_payload"]["query"]

    variant_record = result.records[0]
    assert variant_record["category"] == "population_frequency"
    assert variant_record["label"] == "rs927984495"
    assert variant_record["source_id"] == "11-637293-C-T"
    assert variant_record["url"] == "https://gnomad.broadinstitute.org/variant/11-637293-C-T?dataset=gnomad_r4"
    assert variant_record["caid"] == "CA216942642"
    assert variant_record["rsids"] == ["rs927984495"]
    assert variant_record["reference_genome"] == "GRCh38"
    assert variant_record["frequencies"]["joint"]["ac"] == 20
    assert variant_record["frequencies"]["joint"]["an"] == 1195848
    assert variant_record["frequencies"]["joint"]["af"] == pytest.approx(20 / 1195848)
    assert variant_record["frequencies"]["genome"]["af"] == pytest.approx(3.969671708149736e-05)
    assert variant_record["frequencies"]["exome"]["faf95"]["popmax_population"] == "nfe"
    assert variant_record["top_populations"][0]["dataset"] == "genome"
    assert variant_record["top_populations"][0]["id"] == "nfe"
    assert variant_record["filters"] == ["discrepant_frequencies"]
    assert variant_record["transcript_consequence"]["transcript_id"] == "ENST00000176183"
    assert variant_record["transcript_consequence"]["consequence_terms"] == ["5_prime_UTR_variant"]
    assert variant_record["in_silico_predictors"] == [{"id": "CADD", "value": "3.2"}]
    assert variant_record["gene_constraint"]["pLI"] == pytest.approx(1.0280908139937256e-10)
    assert variant_record["clinvar_release_date"] == "2026-06-06"
    assert "gnomAD v4 rs927984495 (11-637293-C-T) at GRCh38 11:637293 C>T" in variant_record[
        "summary"
    ]
    assert "joint AF 1.67e-05 (AC 20/AN 1,195,848; hom 0; FAF95 nfe 1.35e-05)" in variant_record[
        "summary"
    ]
    assert "genomes AF 3.97e-05" in variant_record["summary"]
    assert "exomes AF 1.34e-05" in variant_record["summary"]
    assert "highest observed population genome nfe AF 8.86e-05" in variant_record["summary"]
    assert "consequence 5_prime_UTR_variant in DRD4 transcript ENST00000176183" in variant_record["summary"]
    assert "HGVS ENST00000176183.5:c.-11C>T" in variant_record["summary"]
    assert "predictors: CADD 3.2" in variant_record["summary"]
    assert "filters: discrepant_frequencies" in variant_record["summary"]

    gene_record = result.records[1]
    assert gene_record["category"] == "gene_constraint"
    assert gene_record["label"] == "DRD4 gnomAD constraint"
    assert gene_record["source_id"] == "ENSG00000069696"
    assert gene_record["url"] == "https://gnomad.broadinstitute.org/gene/DRD4?dataset=gnomad_r4"
    assert gene_record["location"] == "GRCh38 11:637269-640706"
    assert gene_record["constraint"]["oe_lof"] == pytest.approx(1.1746694306580123)
    assert gene_record["constraint"]["flags"] == ["outlier_mis", "outlier_syn"]
    assert gene_record["mane_select_transcript"]["refseq_id"] == "NM_000797"
    assert "gnomAD v4 DRD4 gene constraint (ENSG00000069696)" in gene_record["summary"]
    assert "pLI 1.03e-10" in gene_record["summary"]
    assert "LoF O/E 1.175 (0.512-2.211), Z -0.6971" in gene_record["summary"]
    assert "missense O/E 1.444 (1.18-1.76), Z -4.759" in gene_record["summary"]
    assert "constraint flags: outlier_mis, outlier_syn" in gene_record["summary"]


def test_gnomad_connector_reports_graphql_error_without_traceback():
    client = RecordingGnomadClient(graphql_errors=True)
    connector = connector_for(get_source_spec("gnomad"), client, ResolvedCredential("gnomad"))

    result = connector.query(
        KnowledgeQuery(
            gene="DRD4",
            region="11:637293-640706",
            genome_build="hg38",
            variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
        )
    )
    status = result.to_status(get_source_spec("gnomad"))

    assert result.status == "failed"
    assert "gnomAD GraphQL query failed: Cannot query field" in result.message
    assert "Traceback" not in json.dumps(status)


def test_gwas_catalog_connector_uses_v2_rsid_then_gene_fallback_and_enriches_associations():
    client = RecordingGwasCatalogClient()
    connector = connector_for(get_source_spec("gwas_catalog"), client, ResolvedCredential("gwas_catalog"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert len(result.records) == 1
    assert [call["url"] for call in client.calls[:3]] == [
        "https://www.ebi.ac.uk/gwas/rest/api/v2/metadata",
        "https://www.ebi.ac.uk/gwas/rest/api/v2/associations",
        "https://www.ebi.ac.uk/gwas/rest/api/v2/associations",
    ]
    assert client.calls[1]["params"] == {"rs_id": "rs927984495", "size": 5}
    assert client.calls[2]["params"] == {"mapped_gene": "DRD4", "size": 5}
    assert "no rsID-specific associations were found for rs927984495" in result.message

    record = result.records[0]
    assert record["category"] == "population_association"
    assert record["label"] == "rs1870723-A - Hypothyroidism"
    assert record["source_id"] == "216548355"
    assert record["url"] == "https://www.ebi.ac.uk/gwas/rest/api/v2/associations/216548355"
    assert record["browser_url"] == "https://www.ebi.ac.uk/gwas/associations/216548355"
    assert record["variant"] == "rs927984495"
    assert record["query_context"] == "gene_fallback"
    assert record["association_id"] == "216548355"
    assert record["accession_id"] == "GCST90627750"
    assert record["study_url"] == "https://www.ebi.ac.uk/gwas/studies/GCST90627750"
    assert record["snp_url"] == "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms/rs1870723"
    assert record["rsids"] == ["rs1870723"]
    assert record["snp_effect_alleles"] == ["rs1870723-A"]
    assert record["snp_alleles"] == [{"rsid": "rs1870723", "effect_allele": "A"}]
    assert record["risk_frequency"] == "0.2688"
    assert record["p_value"] == 5e-18
    assert record["beta"] == "0.0352 unit decrease"
    assert record["ci_lower"] == 0.027
    assert record["ci_upper"] == 0.043
    assert record["reported_traits"] == ["Hypothyroidism"]
    assert record["efo_traits"] == [{"id": "MONDO_0005420", "trait": "hypothyroidism"}]
    assert record["mapped_genes"] == ["DRD4"]
    assert record["locations"] == ["11:640349"]
    assert record["pubmed_id"] == "41644669"
    assert record["first_author"] == "White SL"
    assert record["study"]["full_summary_stats_available"] is True
    assert record["study"]["snp_count"] == 28645100
    assert record["study"]["cohort"] == ["AllofUs", "BBJ", "BioVU", "FinnGen", "UKB"]
    assert record["loci"][0]["strongest_risk_alleles"][0]["risk_allele_name"] == "rs1870723-A"
    assert record["data_release_date"] == "2026-06-22"
    assert record["dbsnp_build"] == "156"
    assert "GWAS Catalog association rs1870723-A with Hypothyroidism / hypothyroidism" in record["summary"]
    assert "mapped to DRD4 at 11:640349" in record["summary"]
    assert "p=5e-18" in record["summary"]
    assert "beta 0.0352 unit decrease; CI 0.027-0.043" in record["summary"]
    assert "risk frequency 0.2688" in record["summary"]
    assert "study GCST90627750, PMID 41644669, first author White SL" in record["summary"]
    assert "full summary statistics available" in record["summary"]
    assert "strongest risk allele rs1870723-A" in record["summary"]
    assert "data release 2026-06-22" in record["summary"]


def test_gwas_catalog_connector_keeps_association_when_optional_enrichment_fails():
    client = RecordingGwasCatalogClient(rsid_has_records=True, fail_optional=True)
    connector = connector_for(get_source_spec("gwas_catalog"), client, ResolvedCredential("gwas_catalog"))

    result = connector.query(
        KnowledgeQuery(
            gene="DRD4",
            region="11:637293-640706",
            genome_build="hg38",
            variants=(QueryVariant(chrom="11", pos=640349, ref="G", alt="A", rsid="rs1870723"),),
        )
    )

    assert result.status == "ok"
    assert len(result.records) == 1
    assert result.records[0]["query_context"] == "rsid"
    assert "study" not in result.records[0]
    assert "loci" not in result.records[0]
    assert result.warnings == [
        "Optional GWAS Catalog study detail lookup failed for GCST90627750; association details were still used.",
        "Optional GWAS Catalog loci lookup failed for association 216548355; association details were still used.",
    ]
    assert "Traceback" not in json.dumps(result.to_status(get_source_spec("gwas_catalog")))


def test_pgs_catalog_connector_enriches_linked_score_and_performance_details():
    client = RecordingPgsCatalogClient()
    connector = connector_for(get_source_spec("pgs_catalog"), client, ResolvedCredential("pgs_catalog"))
    query = KnowledgeQuery(
        gene="GENE1",
        region="1:1-10",
        genome_build="hg38",
        variants=(QueryVariant(chrom="1", pos=1, ref="A", alt="G", rsid="rs123"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert len(result.records) == 1
    assert [call["url"] for call in client.calls] == [
        "https://www.pgscatalog.org/rest/variant/rs123",
        "https://www.pgscatalog.org/rest/score/PGS000001",
        "https://www.pgscatalog.org/rest/performance/search",
    ]
    assert client.calls[2]["params"] == {"pgs_id": "PGS000001"}
    record = result.records[0]
    assert record["category"] == "polygenic_score"
    assert record["label"] == "PGS000001 - PRS77_BC - Breast cancer"
    assert record["source_id"] == "PGS000001"
    assert record["url"] == "https://www.pgscatalog.org/score/PGS000001/"
    assert record["variant_url"] == "https://www.pgscatalog.org/variant/rs123/"
    assert record["variant"] == "rs123"
    assert record["pgs_id"] == "PGS000001"
    assert record["trait_reported"] == "Breast cancer"
    assert record["trait_efo"] == [
        {
            "id": "MONDO_0004989",
            "label": "breast carcinoma",
            "description": "A carcinoma that arises from epithelial cells of the breast",
            "url": "http://purl.obolibrary.org/obo/MONDO_0004989",
        }
    ]
    assert record["method_name"] == "SNPs passing genome-wide significance"
    assert record["method_params"] == "P<5x10-8"
    assert record["variants_number"] == 77
    assert record["weight_type"] == "beta"
    assert record["publication"]["pmid"] == "25855707"
    assert record["samples_variants"][0]["sample_number"] == 22627
    assert record["samples_variants"][0]["source_gwas_catalog"] == "GCST001937"
    assert record["performance"][0]["id"] == "PPM000001"
    assert record["performance"][0]["samples"][0]["sample_cases"] == 33673
    assert record["performance"][0]["metrics"][0]["name_short"] == "OR"
    assert record["performance"][0]["metrics"][1]["name_short"] == "C-index"
    assert set(record["ftp_harmonized_scoring_files"]) == {"GRCh37", "GRCh38"}
    assert record["associated_pgs_ids"] == ["PGS000001"]
    assert "PGS Catalog score PGS000001 (PRS77_BC) includes variant rs123 and predicts Breast cancer" in record[
        "summary"
    ]
    assert "mapped traits: breast carcinoma" in record["summary"]
    assert "77 variants" in record["summary"]
    assert "method SNPs passing genome-wide significance (P<5x10-8)" in record["summary"]
    assert "weight type beta" in record["summary"]
    assert "variant-source sample 22,627 European" in record["summary"]
    assert "GWAS ancestry n=22,627 (EUR 100%)" in record["summary"]
    assert "evaluation ancestry n=10 (EUR 80%, NR 20%)" in record["summary"]
    assert "performance All breast cancer: OR 1.55 (1.52-1.58), C-index 0.622 (0.619-0.627)" in record[
        "summary"
    ]
    assert "evaluated sample 67,054 European" in record["summary"]
    assert "publication Mavaddat N, 2015-04-08, PMID 25855707, DOI 10.1093/jnci/djv036" in record[
        "summary"
    ]
    assert "harmonized scoring files: GRCh37, GRCh38" in record["summary"]
    assert "released 2019-10-14" in record["summary"]


def test_pgs_catalog_connector_returns_ok_when_variant_not_in_catalog():
    client = RecordingPgsCatalogClient(variant_not_found=True)
    connector = connector_for(get_source_spec("pgs_catalog"), client, ResolvedCredential("pgs_catalog"))

    result = connector.query(
        KnowledgeQuery(
            gene="DRD4",
            region="11:637293-640706",
            genome_build="hg38",
            variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
        )
    )

    assert result.status == "ok"
    assert result.records == []
    assert "no deposited score variant record found for rs927984495" in result.message
    assert "metadata" not in json.dumps(result.to_status(get_source_spec("pgs_catalog"))).lower()


def test_igsr_connector_returns_high_coverage_and_phase3_file_context():
    client = RecordingIgsrClient()
    connector = connector_for(get_source_spec("igsr"), client, ResolvedCredential("igsr"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert [call["url"] for call in client.calls] == [
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/",
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    ]
    assert len(result.records) == 1

    record = result.records[0]
    assert record["category"] == "population_reference_panel"
    assert record["source_id"] == "1000G_2504_high_coverage_GRCh38_chr11"
    assert record["variant"] == "rs927984495"
    assert record["rsid"] == "rs927984495"
    assert record["assembly"] == "GRCh38"
    assert record["samples"] == 2504
    assert record["populations"] == 26
    assert record["file_name"] == "CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz"
    assert record["file_size"] == "1.6G"
    assert record["index_url"].endswith(".vcf.gz.tbi")
    assert record["population_groups"] == ["AFR", "AMR", "EAS", "EUR", "SAS"]
    assert record["related_files"][1]["name"] == "ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    assert record["related_files"][2]["name"] == "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"
    assert record["related_files"][3]["name"] == "integrated_call_samples_v3.20130502.ALL.panel"
    assert "high-coverage GRCh38 data-access context for rs927984495 at 11:637293 C>T" in record["summary"]
    assert "30x NYGC 2504-sample Phase 3 panel" in record["summary"]
    assert "extract exact genotypes, AC/AN, or allele counts" in record["summary"]
    assert "legacy Phase 3 GRCh37 integrated data" in record["summary"]
    assert "global and superpopulation AF tags" in record["summary"]
    assert "rsIDs were removed from the Phase 3 v5b VCF" in record["summary"]


def test_igsr_connector_falls_back_to_phase3_when_high_coverage_listing_fails():
    client = RecordingIgsrClient(fail_high=True)
    connector = connector_for(get_source_spec("igsr"), client, ResolvedCredential("igsr"))

    result = connector.query(
        KnowledgeQuery(
            gene="DRD4",
            region="11:637293-640706",
            genome_build="hg38",
            variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
        )
    )

    assert result.status == "ok"
    assert len(result.records) == 1
    assert result.records[0]["source_id"] == "1000G_phase3_20130502_GRCh37_chr11"
    assert result.records[0]["assembly"] == "GRCh37"
    assert "input build is GRCh38" in result.records[0]["summary"]
    assert result.warnings == [
        "Optional IGSR high-coverage GRCh38 phased VCF listing lookup failed; available release details were still used."
    ]
    assert "Traceback" not in json.dumps(result.to_status(get_source_spec("igsr")))


def test_civic_connector_queries_current_graphql_schema_and_summarizes_evidence():
    client = RecordingCivicClient()
    connector = connector_for(get_source_spec("civic"), client, ResolvedCredential("civic"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg19",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs1800955"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert len(result.records) == 1
    call = client.calls[0]
    graphql_query = call["json_payload"]["query"]
    assert call["url"] == "https://civicdb.org/api/graphql"
    assert "gene(entrezSymbol: $entrezSymbol)" in graphql_query
    assert "singleVariantMolecularProfile" in graphql_query
    assert "evidenceItems(first: $evidenceFirst, includeRejected: false)" in graphql_query
    assert call["json_payload"]["variables"] == {
        "entrezSymbol": "DRD4",
        "variantFirst": 5,
        "evidenceFirst": 3,
    }
    assert call["headers"]["Content-Type"] == "application/json"

    record = result.records[0]
    assert record["label"] == "DRD4 V194G"
    assert record["source_id"] == "101"
    assert record["url"] == "https://civicdb.org/variants/101"
    assert record["variant"] == "DRD4 V194G"
    assert record["variant_aliases"] == ["rs1800955"]
    assert record["variant_types"] == ["missense_variant"]
    assert record["hgvs_descriptions"] == ["NM_000797.4:c.581T>G"]
    assert record["clinvar_ids"] == ["123"]
    assert record["molecular_profile"] == "DRD4 V194G"
    assert record["evidence_count"] == 2
    assert record["evidence_type"] == "PREDICTIVE"
    assert record["evidence_level"] == "B"
    assert record["evidence_rating"] == "4"
    assert record["evidence_direction"] == "SUPPORTS"
    assert record["significance"] == "SENSITIVITYRESPONSE"
    assert record["disease"] == "Example cancer"
    assert record["therapies"] == ["Example therapy"]
    assert record["citation_id"] == "123456"
    assert record["source_type"] == "PUBMED"
    assert record["evidence_items"][0]["url"] == "https://civicdb.org/evidence/301"
    assert record["evidence_items"][0]["source_url"] == "https://civicdb.org/sources/601"
    assert "CIViC DRD4 variant DRD4 V194G" in record["summary"]
    assert "variant type missense_variant" in record["summary"]
    assert "aliases: rs1800955" in record["summary"]
    assert "PREDICTIVE evidence, level B, rating 4, SUPPORTS" in record["summary"]
    assert "significance SENSITIVITYRESPONSE" in record["summary"]
    assert "disease Example cancer" in record["summary"]
    assert "therapies Example therapy" in record["summary"]
    assert "source PMID 123456" in record["summary"]
    assert "status ACCEPTED" in record["summary"]
    assert "2 total accepted evidence item(s)" in record["summary"]
    assert "DRD4 V194G predicts response in Example cancer." in record["summary"]


def test_civic_connector_reports_graphql_error_without_traceback():
    client = RecordingCivicClient(graphql_errors=True)
    connector = connector_for(get_source_spec("civic"), client, ResolvedCredential("civic"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:1-10", genome_build="hg19"))
    status = result.to_status(get_source_spec("civic"))

    assert result.status == "failed"
    assert result.records == []
    assert "CIViC GraphQL query failed: Cannot query field" in result.message
    assert "Traceback" not in json.dumps(status)


def test_panelapp_connector_uses_exact_gene_lookup_and_summarizes_panel_entry():
    client = RecordingPanelAppClient()
    connector = connector_for(get_source_spec("panelapp"), client, ResolvedCredential("panelapp"))

    result = connector.query(KnowledgeQuery(gene="BRCA1", region="17:41196312-41277500", genome_build="hg19"))

    assert result.status == "ok"
    assert len(result.records) == 1
    call = client.calls[0]
    assert call["url"] == "https://panelapp.genomicsengland.co.uk/api/v1/genes/"
    assert call["params"] == {"entity_name": "BRCA1"}
    assert call["headers"]["Accept"] == "application/json"
    record = result.records[0]
    assert record["category"] == "gene_panel"
    assert record["label"] == "BRCA1 in Inherited ovarian cancer (without breast cancer)"
    assert record["source_id"] == "143:BRCA1"
    assert record["url"] == "https://panelapp.genomicsengland.co.uk/panels/143/gene/BRCA1/"
    assert record["confidence_level"] == "3"
    assert record["confidence_label"] == "green/high evidence"
    assert record["panel_status"] == "public"
    assert record["panel_version"] == "5.1"
    assert record["panel_types"] == ["Rare Disease 100K", "GMS signed-off"]
    assert record["relevant_disorders"] == ["Familial ovarian cancer", "R207"]
    assert record["phenotypes"] == ["{Breast-ovarian cancer, familial, 1}, OMIM:604370"]
    assert record["evidence"] == ["NHS GMS", "Expert Review Green", "Expert list"]
    assert record["publications"] == ["29661970"]
    assert record["genomic_locations"][0]["assembly"] == "GRch37"
    assert "BRCA1 is listed on PanelApp panel Inherited ovarian cancer (without breast cancer)" in record["summary"]
    assert "confidence 3 (green/high evidence)" in record["summary"]
    assert "phenotypes: {Breast-ovarian cancer, familial, 1}, OMIM:604370" in record["summary"]
    assert "inheritance: BOTH monoallelic and biallelic, autosomal or pseudoautosomal" in record["summary"]
    assert "evidence: NHS GMS, Expert Review Green, Expert list" in record["summary"]
    assert "publications: 29661970" in record["summary"]
    assert "panel type: Rare Disease 100K, GMS signed-off" in record["summary"]
    assert "gene location: GRch37 17:41196312-41277500 (ENSG00000012048)" in record["summary"]


def test_panelapp_connector_returns_no_records_for_absent_exact_gene():
    client = RecordingPanelAppClient(empty=True)
    connector = connector_for(get_source_spec("panelapp"), client, ResolvedCredential("panelapp"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg19"))

    assert result.status == "ok"
    assert result.records == []
    assert "no exact gene panel entries found for DRD4" in result.message
    assert "metadata" not in json.dumps(result.to_status(get_source_spec("panelapp"))).lower()


def test_panelapp_connector_filters_non_exact_api_rows():
    client = RecordingPanelAppClient(non_exact=True)
    connector = connector_for(get_source_spec("panelapp"), client, ResolvedCredential("panelapp"))

    result = connector.query(KnowledgeQuery(gene="BRCA1", region="17:41196312-41277500", genome_build="hg19"))

    assert result.status == "ok"
    assert result.records == []
    assert "no exact gene panel entries found for BRCA1" in result.message


def test_mavedb_connector_fetches_gene_score_sets_and_summarizes_records():
    client = RecordingMaveDbClient()
    connector = connector_for(get_source_spec("mavedb"), client, ResolvedCredential("mavedb"))

    result = connector.query(KnowledgeQuery(gene="BRCA1", region="17:41196312-41277500", genome_build="hg19"))

    assert result.status == "ok"
    assert len(result.records) == 1
    call = client.calls[0]
    assert call["url"] == "https://api.mavedb.org/api/v1/genes/BRCA1"
    assert call["params"] == {"limit": 5, "offset": 0}
    assert call["headers"]["Accept"] == "application/json"
    record = result.records[0]
    assert record["category"] == "functional_assay_score_set"
    assert record["label"] == "Scores from multiplexed functional assay of BRCA1 variants"
    assert record["source_id"] == "urn:mavedb:00001222-b-2"
    assert record["url"] == "https://www.mavedb.org/score-sets/urn:mavedb:00001222-b-2"
    assert record["gene"] == "BRCA1"
    assert record["gene_name"] == "BRCA1 DNA repair associated"
    assert record["hgnc_id"] == "HGNC:1100"
    assert record["omim_id"] == "113705"
    assert record["gene_location"] == "17q21.31"
    assert record["num_variants"] == 2271
    assert record["total_gene_scored_variants"] == 2271
    assert record["experiment_urn"] == "urn:mavedb:00001222-b"
    assert record["experiment_title"] == "Multiplexed functional assay of BRCA1 variants"
    assert record["target_genes"] == [
        {
            "name": "BRCA1",
            "category": "protein_coding",
            "mapped_hgnc_name": "BRCA1",
            "uniprot_id": "P38398",
            "external_identifiers": ["HGNC:1100"],
        }
    ]
    assert record["publications"] == ["PMID 39999999"]
    assert record["license"] == "CC BY 4.0"
    assert "MaveDB score set urn:mavedb:00001222-b-2 for BRCA1" in record["summary"]
    assert "Scores from multiplexed functional assay of BRCA1 variants" in record["summary"]
    assert "2,271 scored variants" in record["summary"]
    assert "published 2025-10-22" in record["summary"]
    assert "experiment: Multiplexed functional assay of BRCA1 variants" in record["summary"]
    assert "assay summary: Multiplexed assay of BRCA1 variants measuring homology directed repair activity" in record[
        "summary"
    ]
    assert "target genes: BRCA1" in record["summary"]
    assert "publications: PMID 39999999" in record["summary"]
    assert "license CC BY 4.0" in record["summary"]


def test_mavedb_connector_returns_no_records_when_gene_has_no_score_sets():
    client = RecordingMaveDbClient(empty=True)
    connector = connector_for(get_source_spec("mavedb"), client, ResolvedCredential("mavedb"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg19"))

    assert result.status == "ok"
    assert result.records == []
    assert "DRD4 has no published MAVE score sets" in result.message
    assert "metadata" not in json.dumps(result.to_status(get_source_spec("mavedb"))).lower()


def test_dbsnp_connector_enriches_esummary_with_refsnp_detail():
    client = RecordingDbSnpClient()
    connector = connector_for(get_source_spec("dbsnp"), client, ResolvedCredential("dbsnp"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637373, ref="T", alt="", rsid="rs2533154733"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert len(result.records) == 1
    record = result.records[0]
    assert record["label"] == "rs2533154733"
    assert record["rsid"] == "rs2533154733"
    assert record["source_id"] == "2533154733"
    assert record["snp_class"] == "del"
    assert record["variant_type"] == "del"
    assert record["chromosome"] == "11"
    assert record["position"] == "637373"
    assert record["spdi"] == "NC_000011.10:637372:T:"
    assert record["genes"] == [{"name": "DRD4", "gene_id": "1815"}]
    assert record["sequence_ontology"] == ["coding_sequence_variant", "frameshift_variant"]
    assert record["transcripts"][0]["accession"] == "NM_000797.4"
    assert record["proteins"][0]["accession"] == "NP_000788.2"
    assert record["frequencies"][0]["allele_count"] == "1"
    assert record["frequencies"][0]["total_count"] == "998848"
    assert "rs2533154733: deletion at 11:637373" in record["summary"]
    assert "DRD4 coding_sequence_variant / frameshift_variant" in record["summary"]
    assert "HGVS NC_000011.10:g.637373del, NM_000797.4:c.69del, NP_000788.2:p.Ala24fs" in record[
        "summary"
    ]
    assert "frequency GnomAD_exomes 1/998848" in record["summary"]
    assert "validated by-frequency" in record["summary"]
    assert "submitted by GNOMAD" in record["summary"]
    refsnp_calls = [call for call in client.calls if "variation/v0/refsnp/2533154733" in str(call["url"])]
    assert len(refsnp_calls) == 1
    assert refsnp_calls[0]["headers"]["Accept"] == "application/json"
    assert refsnp_calls[0]["timeout"] == 8


def test_dbsnp_connector_keeps_esummary_details_when_optional_refsnp_fails():
    client = RecordingDbSnpClient(fail_refsnp=True)
    connector = connector_for(get_source_spec("dbsnp"), client, ResolvedCredential("dbsnp"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637373, ref="T", alt="", rsid="rs2533154733"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == [
        "Optional dbSNP RefSNP detail lookup failed for rs2533154733; ESummary details were still used."
    ]
    assert "Traceback" not in json.dumps(result.warnings)
    assert "timed out" not in json.dumps(result.warnings)
    record = result.records[0]
    assert record["label"] == "rs2533154733"
    assert record["summary"].startswith(
        "rs2533154733: deletion at 11:637373 (NC_000011.10:637372:T:)"
    )
    assert "DRD4 frameshift_variant / coding_sequence_variant" in record["summary"]
    assert "frequency GnomAD_exomes -=0.000001/1" in record["summary"]
    assert "validated by-frequency" in record["summary"]
    assert "submitted by GNOMAD" in record["summary"]


def test_ensembl_connector_enriches_overlap_with_variation_and_vep_details():
    client = RecordingEnsemblClient()
    connector = connector_for(get_source_spec("ensembl"), client, ResolvedCredential("ensembl"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg19",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert "ENSG00000069696" in result.records[0]["summary"]
    variant_records = [record for record in result.records if record["category"] == "variant_annotation"]
    assert len(variant_records) == 1
    record = variant_records[0]
    assert record["source_id"] == "rs927984495"
    assert record["rsid"] == "rs927984495"
    assert record["location"] == "GRCh37 11:637293-637293"
    assert record["alleles"] == "C/T"
    assert record["consequence"] == "5_prime_UTR_variant"
    assert record["variant_class"] == "SNP"
    assert record["evidence"] == ["Frequency", "TOPMed", "gnomAD"]
    assert record["synonyms"] == ["example-synonym"]
    assert record["transcript_consequence"]["transcript_id"] == "ENST00000176183"
    assert record["colocated_variant"]["id"] == "rs927984495"
    assert "rs927984495 (C/T) at GRCh37 11:637293-637293" in record["summary"]
    assert "5_prime_UTR_variant in DRD4 canonical transcript ENST00000176183" in record["summary"]
    assert "exon 1/4" in record["summary"]
    assert "cDNA position 1" in record["summary"]
    assert "evidence: Frequency, TOPMed, gnomAD" in record["summary"]
    assert any("variation/homo_sapiens/rs927984495" in call["url"] for call in client.calls)
    assert any("vep/homo_sapiens/region/11:637293:637293/T" in call["url"] for call in client.calls)


def test_ensembl_connector_keeps_overlap_record_when_optional_enrichment_fails():
    client = RecordingEnsemblClient(fail_variation=True, fail_vep=True)
    connector = connector_for(get_source_spec("ensembl"), client, ResolvedCredential("ensembl"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg19",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    variant_records = [record for record in result.records if record["category"] == "variant_annotation"]
    assert len(variant_records) == 1
    record = variant_records[0]
    assert record["source_id"] == "rs927984495"
    assert record["summary"].startswith("rs927984495 (C/T) at GRCh37 11:637293-637293: intergenic_variant")
    assert "source dbSNP" in record["summary"]
    assert "input variant 11:637293:C>T" in record["summary"]
    assert result.warnings == [
        "Optional Ensembl VEP annotation failed for 11:637293:C>T; overlap details were still used.",
        "Optional Ensembl variation detail failed for rs927984495; overlap details were still used.",
    ]
    assert "Traceback" not in json.dumps(result.warnings)
    assert "timed out" not in json.dumps(result.warnings)


def test_ucsc_connector_returns_gene_sequence_regulatory_and_track_context():
    client = RecordingUcscClient()
    connector = connector_for(get_source_spec("ucsc"), client, ResolvedCredential("ucsc"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert "7 compact annotation record" in result.message
    categories = {record["category"] for record in result.records}
    assert categories == {
        "gene_model",
        "reference_sequence",
        "cpg_island",
        "regulatory_element",
        "transcription_factor_binding",
        "repeat_annotation",
        "common_variant",
    }
    gene_record = next(record for record in result.records if record["category"] == "gene_model")
    assert gene_record["source_id"] == "NM_000797.4"
    assert gene_record["transcript_interval"] == "chr11:637269-640706"
    assert gene_record["cds_interval"] == "chr11:637305-640603"
    assert gene_record["exon_count"] == 4
    assert "UCSC ncbiRefSeq NM_000797.4 for DRD4 on hg38 chr11:637269-640706 (+): 4 exons" in gene_record[
        "summary"
    ]
    assert "query window hg38 chr11:637293-640706" in gene_record["summary"]

    sequence_record = next(record for record in result.records if record["category"] == "reference_sequence")
    assert sequence_record["reference_base"] == "C"
    assert sequence_record["position"] == 637293
    assert "reference base at rs927984495 / chr11:637293 is C" in sequence_record["summary"]
    assert "sample allele C>T" in sequence_record["summary"]

    cpg_record = next(record for record in result.records if record["category"] == "cpg_island")
    assert "CpG count 313" in cpg_record["summary"]
    assert "GC 54.3%" in cpg_record["summary"]
    assert "observed/expected CpG 1.01" in cpg_record["summary"]

    ccre_record = next(record for record in result.records if record["category"] == "regulatory_element")
    assert ccre_record["source_id"] == "EH38E1513769"
    assert "EH38E1513769 promoter-like signature" in ccre_record["summary"]
    assert "labels include PLS, pELS" in ccre_record["summary"]

    tfbs_record = next(record for record in result.records if record["category"] == "transcription_factor_binding")
    assert "capped by maxItemsOutput" in tfbs_record["summary"]
    assert "EZH2 score 1000 (4 sources)" in tfbs_record["summary"]

    repeat_record = next(record for record in result.records if record["category"] == "repeat_annotation")
    assert "G-rich Low_complexity/Low_complexity at hg38 chr11:637322-637431" in repeat_record["summary"]
    assert "(AC)n Simple_repeat/Simple_repeat" in repeat_record["summary"]

    snp_record = next(record for record in result.records if record["category"] == "common_variant")
    assert snp_record["scope"] == "query window"
    assert "rs146680769 at hg38 chr11:637294-637294 C/T single validated by-frequency/by-1000genomes" in snp_record[
        "summary"
    ]

    sequence_calls = [call for call in client.calls if str(call["url"]).endswith("/getData/sequence")]
    assert sequence_calls[0]["params"]["start"] == 637292
    assert sequence_calls[0]["params"]["end"] == 637293
    assert sequence_calls[0]["headers"]["Accept"] == "application/json"
    assert sequence_calls[0]["rate_limit_per_second"] == 1.0
    track_calls = [call for call in client.calls if str(call["url"]).endswith("/getData/track")]
    assert all(call["params"].get("maxItemsOutput") == 5 for call in track_calls)


def test_ucsc_connector_keeps_partial_annotations_when_optional_track_fails():
    client = RecordingUcscClient(fail_track="cpgIslandExt")
    connector = connector_for(get_source_spec("ucsc"), client, ResolvedCredential("ucsc"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == [
        "Optional UCSC cpgIslandExt track lookup failed; other UCSC annotations were still used."
    ]
    assert "Traceback" not in json.dumps(result.warnings)
    assert "timed out" not in json.dumps(result.warnings)
    categories = {record["category"] for record in result.records}
    assert "gene_model" in categories
    assert "reference_sequence" in categories
    assert "regulatory_element" in categories
    assert "cpg_island" not in categories


def test_encode_connector_returns_experiment_and_region_context():
    client = RecordingEncodeClient()
    connector = connector_for(get_source_spec("encode"), client, ResolvedCredential("encode"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert "2 experiment record(s), 1 region-search record(s)" in result.message
    categories = [record["category"] for record in result.records]
    assert categories == ["regulatory_experiment", "regulatory_experiment", "regulatory_region_hit"]

    target_record = result.records[0]
    assert target_record["label"] == "ENCSRDRD4TF - TF ChIP-seq"
    assert target_record["match_type"] == "target_gene"
    assert target_record["target"] == "DRD4"
    assert target_record["target_genes"] == ["DRD4"]
    assert target_record["biological_replicates"] == 2
    assert target_record["files"]["total"] == 2
    assert target_record["files"]["released"] == 2
    assert target_record["files"]["assemblies"] == ["GRCh38"]
    assert target_record["files"]["output_types"] == ["optimal IDR thresholded peaks", "signal p-value"]
    assert "ENCODE DCC ENCSRDRD4TF TF ChIP-seq matched DRD4 through target gene metadata" in target_record[
        "summary"
    ]
    assert "target DRD4 (DRD4)" in target_record["summary"]
    assert "biosample Homo sapiens K562 cell line" in target_record["summary"]
    assert "2 biological replicate(s)" in target_record["summary"]
    assert "2/2 released file(s); outputs optimal IDR thresholded peaks, signal p-value; assemblies GRCh38" in target_record[
        "summary"
    ]
    assert "lab ENCODE Processing Pipeline" in target_record["summary"]
    assert "released 2026-01-15" in target_record["summary"]

    text_record = result.records[1]
    assert text_record["source_id"] == "ENCSRDRD4RNA"
    assert text_record["match_type"] == "gene_text"
    assert "total RNA-seq matched DRD4 through portal text search" in text_record["summary"]
    assert "gene quantifications" in text_record["summary"]

    region_record = result.records[2]
    assert region_record["category"] == "regulatory_region_hit"
    assert region_record["source_id"] == "ENCSRREGION1"
    assert region_record["region"] == "chr11:637293-640706"
    assert "ENCODE DCC region search returned ENCSRREGION1 overlapping chr11:637293-640706 for DRD4" in region_record[
        "summary"
    ]
    assert "assay DNase-seq" in region_record["summary"]

    target_calls = [call for call in client.calls if call["params"].get("target.genes.symbol") == "DRD4"]
    assert len(target_calls) == 1
    assert target_calls[0]["params"]["frame"] == "object"
    assert target_calls[0]["params"]["format"] == "json"
    assert target_calls[0]["params"]["limit"] == 5
    assert target_calls[0]["headers"]["Accept"] == "application/json"
    assert target_calls[0]["rate_limit_per_second"] == 5.0
    assert any(call["params"].get("searchTerm") == "DRD4" for call in client.calls)
    assert any(str(call["url"]).endswith("/region-search/") and call["params"].get("region") == "chr11:637293-640706" for call in client.calls)


def test_encode_connector_records_clear_context_when_no_direct_records_exist():
    client = RecordingEncodeClient(empty=True)
    connector = connector_for(get_source_spec("encode"), client, ResolvedCredential("encode"))
    query = KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg38")

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == [
        "Optional ENCODE region search did not return compact DCC rows; experiment searches were still used."
    ]
    assert "query context recorded" in result.message
    assert len(result.records) == 1
    record = result.records[0]
    assert record["category"] == "encode_query_context"
    assert record["target_search_total"] == 0
    assert record["text_search_total"] == 0
    assert "no direct experiment records were returned" in record["summary"]
    assert "target-gene search total 0; gene-text search total 0" in record["summary"]
    assert "region search note: Error during search" in record["summary"]


def test_encode_connector_keeps_experiment_records_when_region_search_fails():
    client = RecordingEncodeClient(fail_region=True)
    connector = connector_for(get_source_spec("encode"), client, ResolvedCredential("encode"))
    query = KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg38")

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == ["Optional ENCODE region search failed; other ENCODE Portal results were still used."]
    assert "Traceback" not in json.dumps(result.warnings)
    categories = {record["category"] for record in result.records}
    assert categories == {"regulatory_experiment"}
    assert len(result.records) == 2


def test_screen_connector_returns_candidate_regulatory_element_details():
    client = RecordingScreenClient()
    connector = connector_for(get_source_spec("screen"), client, ResolvedCredential("screen"))
    query = KnowledgeQuery(
        gene="DRD4",
        region="11:637293-640706",
        genome_build="hg38",
        variants=(QueryVariant(chrom="11", pos=637293, ref="C", alt="T", rsid="rs927984495"),),
    )

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert "2 overlapping cCRE record(s)" in result.message
    assert all("www.encodeproject.org/search" not in str(call["url"]) for call in client.calls)
    assert len(client.calls) == 1
    call = client.calls[0]
    assert call["url"] == "https://screen.encodeproject.org/api/screen-graphql"
    assert "cCREQuery" in call["json_payload"]["query"]
    assert call["json_payload"]["variables"]["assembly"] == "GRCh38"
    assert call["json_payload"]["variables"]["coordinates"] == [
        {"chromosome": "chr11", "start": 637292, "end": 640706}
    ]
    assert call["headers"]["Accept"] == "application/json"
    assert call["headers"]["Content-Type"] == "application/json"
    assert call["rate_limit_per_second"] == 5.0

    record = result.records[0]
    assert record["category"] == "candidate_regulatory_element"
    assert record["label"] == "EH38E2937824 - PLS"
    assert record["source_id"] == "EH38E2937824"
    assert record["location"] == "GRCh38 chr11:637026-637375"
    assert record["nearby_genes"] == [
        {
            "name": "DRD4",
            "id": "ENSG00000069696.7",
            "gene_type": "protein_coding",
            "strand": "+",
        }
    ]
    assert record["max_z_scores"]["DNase"] == 3.523491859436035
    assert record["url"] == "https://screen.encodeproject.org/search/?q=EH38E2937824&assembly=GRCh38"
    assert "SCREEN cCRE EH38E2937824 overlaps DRD4 query window at GRCh38 chr11:637026-637375" in record[
        "summary"
    ]
    assert "promoter-like signature (PLS)" in record["summary"]
    assert "CTCF-bound" in record["summary"]
    assert "intersecting genes: DRD4 (protein_coding)" in record["summary"]
    assert "max assay Z-scores DNase 3.52, H3K4me3 4.00, H3K27ac 2.03, CTCF 1.88" in record["summary"]
    assert "rDHS EH38D4573084" in record["summary"]
    assert "query window chr11:637293-640706" in record["summary"]


def test_screen_connector_records_registry_context_when_ccre_query_fails():
    client = RecordingScreenClient(fail_ccre=True)
    connector = connector_for(get_source_spec("screen"), client, ResolvedCredential("screen"))
    query = KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg38")

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == ["Optional SCREEN cCRE GraphQL query failed; registry context was still recorded."]
    assert "Traceback" not in json.dumps(result.warnings)
    assert "500 Internal Server Error" not in json.dumps(result.warnings)
    assert "SCREEN cCRE Registry context recorded" in result.message
    assert len(result.records) == 1
    record = result.records[0]
    assert record["category"] == "screen_registry_context"
    assert "candidate cis-regulatory element registry for human GRCh38" in record["summary"]
    assert "PLS, pELS, dELS, DNase-H3K4me3, and CTCF-only" in record["summary"]
    assert record["url"] == "https://screen.encodeproject.org/search/?q=chr11%3A637293-640706&assembly=GRCh38"


def test_ewas_catalog_connector_parses_query_specific_tsv_records():
    client = RecordingEwasCatalogClient()
    connector = connector_for(get_source_spec("ewas_catalog"), client, ResolvedCredential("ewas_catalog"))
    query = KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg19")

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert "2 compact association record(s)" in result.message
    assert len(client.calls) == 2
    assert client.calls[0]["url"] == "https://www.ewascatalog.org/"
    assert client.calls[0]["params"] == {"gene": "DRD4"}
    assert client.calls[0]["rate_limit_per_second"] == 1.0
    assert client.calls[1]["url"] == "https://www.ewascatalog.org/static//tmp/DRD4_test.tsv"

    record = result.records[0]
    assert record["category"] == "ewas_association"
    assert record["label"] == "cg02762115 - age"
    assert record["source_id"] == "EWAS0001"
    assert record["cpg"] == "cg02762115"
    assert record["location"] == "chr11:640446"
    assert record["beta"] == "0.0021"
    assert record["p_value"] == "0E+00"
    assert record["n"] == "2338"
    assert record["url"] == "https://www.ewascatalog.org/?cpg=cg02762115"
    assert "EWAS Catalog DRD4 CpG cg02762115 at chr11:640446" in record["summary"]
    assert "exposure/trait age" in record["summary"]
    assert "tissue Whole blood" in record["summary"]
    assert "N=2,338" in record["summary"]
    assert "PMID 33450751, Mulder RH" in record["summary"]


def test_ewas_catalog_connector_records_context_when_tsv_download_fails():
    client = RecordingEwasCatalogClient(fail_tsv=True)
    connector = connector_for(get_source_spec("ewas_catalog"), client, ResolvedCredential("ewas_catalog"))
    query = KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg38")

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == [
        "EWAS Catalog reports CpG locations in hg19 coordinates; input build hg38 was used as provided.",
        "Optional EWAS Catalog query-specific TSV download failed; query context was still recorded.",
    ]
    assert "Traceback" not in json.dumps(result.warnings)
    assert "timed out" not in json.dumps(result.warnings)
    assert len(result.records) == 1
    record = result.records[0]
    assert record["category"] == "ewas_catalog_context"
    assert "Catalog rows include CpG, hg19 location, gene, beta" in record["summary"]


def test_ewas_atlas_connector_returns_position_association_details():
    client = RecordingEwasAtlasClient()
    connector = connector_for(get_source_spec("ewas_atlas"), client, ResolvedCredential("ewas_atlas"))
    query = KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg19")

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == []
    assert "position query hg19 chr11:637293-640706" in result.message
    assert "3 compact association record(s)" in result.message
    assert len(client.calls) == 1
    assert client.calls[0]["url"] == "https://ngdc.cncb.ac.cn/ewas/rest/pos"
    assert client.calls[0]["params"] == {"chr": "11", "start": 637293, "end": 640706}
    assert client.calls[0]["headers"]["Accept"] == "application/json"
    assert client.calls[0]["rate_limit_per_second"] == 1.0

    record = result.records[0]
    assert record["category"] == "ewas_association"
    assert record["label"] == "cg02762115 - cognitive function"
    assert record["source_id"] == "ES00743:cg02762115"
    assert record["location"] == "hg19 chr11:640446"
    assert record["rank"] == "4"
    assert record["trait"] == "cognitive function"
    assert record["correlation"] == "neg"
    assert record["pmid"] == "29311653"
    assert record["url"] == "https://ngdc.cncb.ac.cn/ewas/search?item=cg02762115&term=Probe+Id"
    assert "EWAS Atlas DRD4 probe cg02762115 at hg19 chr11:640446" in record["summary"]
    assert "DRD4 transcript ENST00000176183.5, 3153 bp from TSS" in record["summary"]
    assert "negative methylation-trait correlation" in record["summary"]
    assert "rank 4" in record["summary"]
    assert "PMID 29311653" in record["summary"]


def test_ewas_atlas_connector_records_context_when_no_associations_return():
    client = RecordingEwasAtlasClient(empty=True)
    connector = connector_for(get_source_spec("ewas_atlas"), client, ResolvedCredential("ewas_atlas"))
    query = KnowledgeQuery(gene="DRD4", region="11:637293-640706", genome_build="hg38")

    result = connector.query(query)

    assert result.status == "ok"
    assert result.warnings == [
        "EWAS Atlas REST reports CpG locations as chrHg19/posHg19; input build hg38 was used as provided."
    ]
    assert len(result.records) == 1
    record = result.records[0]
    assert record["category"] == "ewas_atlas_context"
    assert "probe, gene, position, study, and publication objects" in record["summary"]


def test_clingen_connector_returns_gene_centered_curations():
    client = RecordingClinGenClient()
    connector = connector_for(get_source_spec("clingen"), client, ResolvedCredential("clingen"))

    result = connector.query(KnowledgeQuery(gene="GENE1", region="1:1-10", genome_build="hg19"))

    assert result.status == "ok"
    assert "5 gene-centered curation record" in result.message
    categories = {record["category"] for record in result.records}
    assert categories == {
        "gene_disease_validity",
        "dosage_sensitivity",
        "clinical_gene_curation_summary",
        "clinical_actionability",
    }
    validity = next(record for record in result.records if record["category"] == "gene_disease_validity")
    assert validity["classification"] == "Definitive"
    assert validity["disease"] == "Example syndrome"
    assert validity["mondo_id"] == "MONDO:0000001"
    dosage = next(record for record in result.records if record["category"] == "dosage_sensitivity")
    assert dosage["haploinsufficiency"] == "Sufficient Evidence for Haploinsufficiency"
    actionability_records = [record for record in result.records if record["category"] == "clinical_actionability"]
    assert len(actionability_records) == 2
    actionability = actionability_records[0]
    assert actionability["actionability_score"] == "10CN"
    assert actionability["intervention"] == "Surveillance"
    assert all(record["source"] == get_source_spec("clingen").name for record in result.records)


def test_clingen_connector_reports_empty_live_query_without_metadata_record():
    client = RecordingClinGenClient(empty=True)
    connector = connector_for(get_source_spec("clingen"), client, ResolvedCredential("clingen"))

    result = connector.query(KnowledgeQuery(gene="DRD4", region="11:1-10", genome_build="hg19"))

    assert result.status == "ok"
    assert result.records == []
    assert "no ClinGen curation records found for DRD4" in result.message


def test_clingen_connector_keeps_partial_results_when_one_feed_is_malformed():
    client = RecordingClinGenClient(malformed_validity=True)
    connector = connector_for(get_source_spec("clingen"), client, ResolvedCredential("clingen"))

    result = connector.query(KnowledgeQuery(gene="GENE1", region="1:1-10", genome_build="hg19"))

    assert result.status == "ok"
    assert any("Gene-Disease Validity response could not be parsed" in warning for warning in result.warnings)
    assert {record["category"] for record in result.records} == {
        "dosage_sensitivity",
        "clinical_gene_curation_summary",
        "clinical_actionability",
    }


def test_clingen_connector_treats_summary_timeout_as_optional():
    client = RecordingClinGenClient(summary_timeout=True)
    connector = connector_for(get_source_spec("clingen"), client, ResolvedCredential("clingen"))

    result = connector.query(KnowledgeQuery(gene="GENE1", region="1:1-10", genome_build="hg19"))

    assert result.status == "ok"
    assert "4 gene-centered curation record" in result.message
    assert {record["category"] for record in result.records} == {
        "gene_disease_validity",
        "dosage_sensitivity",
        "clinical_actionability",
    }
    assert result.warnings == [
        "Optional ClinGen curation activity summary timed out; primary ClinGen feeds were still used."
    ]
    assert "Traceback" not in json.dumps(result.warnings)
    assert "HTTPSConnectionPool" not in json.dumps(result.warnings)
    summary_url = next(url for url in client.timeouts if "curation-activity-summary-report" in url)
    assert client.timeouts[summary_url] >= 15


def test_clingen_connector_reports_request_failure_without_traceback():
    client = RecordingClinGenClient(fail_all=True)
    connector = connector_for(get_source_spec("clingen"), client, ResolvedCredential("clingen"))

    result = connector.query(KnowledgeQuery(gene="GENE1", region="1:1-10", genome_build="hg19"))
    status = result.to_status(get_source_spec("clingen"))

    assert result.status == "failed"
    assert status["error_code"] == "tls_certificate_verification_failed"
    assert status["remediation"] == "Configure NOPHIGENE_CA_BUNDLE."
    assert "Traceback" not in json.dumps(status)


def test_local_article_extractor_finds_gene_snippets_without_full_text(tmp_path: Path):
    article_dir = tmp_path / "articles"
    article_dir.mkdir()
    _write_simple_pdf(
        article_dir / "gene1_results.pdf",
        (
            "Abstract. GENE1 rs123 showed a significant association with Example syndrome in the "
            "reported cohort. DO_NOT_SERIALIZE_PRIVATE_APPENDIX"
        ),
    )

    extraction = extract_local_article_evidence(
        gene="GENE1",
        pdf_folder=article_dir,
        generated_at="2026-06-17T00:00:00Z",
    )

    assert extraction["status"] == "ok"
    assert extraction["provenance"]["pdf_count"] == 1
    assert extraction["provenance"]["record_count"] == 1
    record = extraction["records"][0]
    assert record["source_key"] == LOCAL_ARTICLE_SOURCE_KEY
    assert record["gene"] == "GENE1"
    assert record["rsid"] == "rs123"
    assert record["claim_type"] in {"clinical_variant", "population_association"}
    serialized = json.dumps(extraction, sort_keys=True)
    assert "DO_NOT_SERIALIZE_PRIVATE_APPENDIX" not in serialized
    assert str(article_dir) not in serialized
    assert record["source_file_sha256"]
    assert record["pdf_path_hash"]


def test_dynamic_builder_adds_local_article_evidence_workflow(tmp_path: Path):
    article_dir = tmp_path / "articles"
    article_dir.mkdir()
    _write_simple_pdf(
        article_dir / "gene1_function.pdf",
        "Results. GENE1 expression was reduced after knockdown and the assay showed altered pathway activity.",
    )
    variants = pd.DataFrame([{"chrom": "1", "id": "rs123", "pos": 101, "ref": "A", "alt": "G"}])

    payload = build_dynamic_knowledge_base(
        gene="GENE1",
        region="1:90-110",
        genome_build="hg19",
        variants=variants,
        selected_sources=["clinvar"],
        use_local_article_evidence=True,
        article_pdf_folder=article_dir,
        output_dir=tmp_path,
        request_client=FakeClient(),
        generated_at="2026-06-17T00:00:00Z",
    )

    statuses = {status["source_key"]: status["status"] for status in payload["provider_statuses"]}
    assert statuses[LOCAL_ARTICLE_SOURCE_KEY] == "ok"
    assert payload["local_article_evidence"]["provenance"]["record_count"] == 1
    assert payload["local_article_evidence_artifacts"]["article_evidence_json"].endswith("article_evidence.json")
    assert any(record["source_key"] == LOCAL_ARTICLE_SOURCE_KEY for record in payload["literature_records"])
    assert payload["workflow_runs"][-1]["workflow_key"] == LOCAL_ARTICLE_WORKFLOW_KEY
    assert payload["workflow_runs"][-1]["status"] == "ok"
    assert payload["workflow_source_matrix"][LOCAL_ARTICLE_SOURCE_KEY] == [LOCAL_ARTICLE_WORKFLOW_KEY]
    assert (tmp_path / "local_article_evidence" / "article_evidence_summary.csv").is_file()


def test_dynamic_builder_reports_missing_local_article_folder_as_input_needed():
    payload = build_dynamic_knowledge_base(
        gene="GENE1",
        region="1:90-110",
        genome_build="hg19",
        selected_sources=[],
        selected_workflows=[LOCAL_ARTICLE_WORKFLOW_KEY],
        use_local_article_evidence=True,
        generated_at="2026-06-17T00:00:00Z",
    )

    statuses = {status["source_key"]: status["status"] for status in payload["provider_statuses"]}
    assert statuses[LOCAL_ARTICLE_SOURCE_KEY] == "needs_folder"
    assert payload["workflow_runs"][0]["workflow_key"] == LOCAL_ARTICLE_WORKFLOW_KEY
    assert payload["workflow_runs"][0]["status"] == "needs_input"
    assert payload["local_article_evidence"]["provenance"]["record_count"] == 0


def test_credential_status_redacts_session_secret():
    spec = get_source_spec("omim")
    statuses = credential_status_for_specs([spec], {"omim": "super-secret-token"})

    assert statuses["omim"] == "session:[redacted]"
    assert "super-secret-token" not in json.dumps(statuses)


def test_dynamic_builder_writes_partial_results_without_serializing_secret(tmp_path: Path):
    variants = pd.DataFrame(
        [
            {
                "sample": "sample-one",
                "chrom": "1",
                "id": "rs123",
                "pos": 101,
                "ref": "A",
                "alt": "G",
                "gt_raw": "0/1",
                "zygosity": "heterozygous",
            }
        ]
    )
    manifest = pd.DataFrame(
        [
            {
                "IlmnID": "cg00000001",
                "CHR": "1",
                "MAPINFO": 99,
                "UCSC_RefGene_Name": "GENE1",
                "Relation_to_UCSC_CpG_Island": "Island",
            }
        ]
    )

    payload = build_dynamic_knowledge_base(
        gene="GENE1",
        region="1:90-110",
        genome_build="hg19",
        variants=variants,
        manifest_subset=manifest,
        selected_sources=["clinvar", "omim", "hgmd"],
        credentials={"omim": "super-secret-token"},
        output_dir=tmp_path,
        request_client=FakeClient(),
        generated_at="2026-06-17T00:00:00Z",
    )

    artifact = tmp_path / "variant_kb.json"
    assert artifact.is_file()
    assert payload["variant_records"][0]["variant"] == "rs123"
    assert payload["epigenetic_locus_records"][0]["probe_id"] == "cg00000001"
    assert {status["source_key"]: status["status"] for status in payload["provider_statuses"]} == {
        "clinvar": "ok",
        "omim": "metadata_only",
        "hgmd": "needs_export",
    }
    serialized = json.dumps(payload, sort_keys=True)
    assert "super-secret-token" not in serialized
    assert "session:[redacted]" in serialized
    assert payload["workflow_runs"][0]["workflow_key"] == "clinical_variant_triage"
    assert payload["workflow_runs"][0]["status"] == "partial"
    assert payload["workflow_source_matrix"]["clinvar"] == ["clinical_variant_triage"]


def test_dynamic_builder_runs_workflows_sequentially_and_deduplicates_sources():
    variants = pd.DataFrame([{"chrom": "1", "id": "rs123", "pos": 101, "ref": "A", "alt": "G"}])
    client = CountingFakeClient()

    payload = build_dynamic_knowledge_base(
        gene="GENE1",
        region="1:90-110",
        genome_build="hg19",
        variants=variants,
        selected_workflows=["clinical_variant_triage", "population_frequency_association"],
        selected_sources=["dbsnp"],
        request_client=client,
        generated_at="2026-06-17T00:00:00Z",
    )

    assert len(payload["provider_statuses"]) == 1
    assert payload["provider_statuses"][0]["source_key"] == "dbsnp"
    assert client.get_counts["rs123"] == 1
    assert [run["workflow_key"] for run in payload["workflow_runs"]] == [
        "clinical_variant_triage",
        "population_frequency_association",
    ]
    assert all(run["selected_source_keys"] == ["dbsnp"] for run in payload["workflow_runs"])
    assert payload["workflow_source_matrix"]["dbsnp"] == [
        "clinical_variant_triage",
        "population_frequency_association",
    ]
    assert payload["source_records"][0]["evidence_id"] == "dbsnp:001"


def test_dynamic_builder_merges_clingen_curations_into_workflow_records():
    payload = build_dynamic_knowledge_base(
        gene="GENE1",
        region="1:1-10",
        genome_build="hg19",
        selected_sources=["clingen"],
        request_client=RecordingClinGenClient(),
        generated_at="2026-06-17T00:00:00Z",
    )

    assert payload["provider_statuses"][0]["source_key"] == "clingen"
    assert payload["provider_statuses"][0]["status"] == "ok"
    assert payload["workflow_source_matrix"]["clingen"] == ["clinical_variant_triage"]
    assert {record["category"] for record in payload["source_records"]} == {
        "gene_disease_validity",
        "dosage_sensitivity",
        "clinical_gene_curation_summary",
        "clinical_actionability",
    }
    assert any(record["classification"] == "Definitive" for record in payload["source_records"])
    assert payload["workflow_runs"][0]["record_counts"]["source_records"] == 5


def test_dynamic_builder_merges_medgen_conditions_into_workflow_records():
    payload = build_dynamic_knowledge_base(
        gene="GENE1",
        region="1:1-10",
        genome_build="hg19",
        selected_sources=["medgen"],
        request_client=RecordingMedGenClient(),
        generated_at="2026-06-17T00:00:00Z",
    )

    assert payload["provider_statuses"][0]["source_key"] == "medgen"
    assert payload["provider_statuses"][0]["status"] == "ok"
    assert payload["workflow_source_matrix"]["medgen"] == ["clinical_variant_triage"]
    assert {record["category"] for record in payload["source_records"]} == {"clinical_condition"}
    assert any(record["concept_id"] == "C0000001" for record in payload["source_records"])
    assert payload["workflow_runs"][0]["record_counts"]["source_records"] == 2


def test_source_import_parser_normalizes_json_csv_and_discards_raw_columns(tmp_path: Path):
    csv_path = tmp_path / "hgmd.csv"
    csv_path.write_text(
        "gene,rsid,mutation,classification,disease,pmid,private_payload\n"
        "GENE1,rs123,c.1A>G,Pathogenic,Example syndrome,PMID:1,do not serialize\n"
        "GENE1,rs123,c.1A>G,Pathogenic,Example syndrome,PMID:1,duplicate raw\n",
        encoding="utf-8",
    )
    bundle = parse_source_import(
        "hgmd",
        csv_path,
        spec=get_source_spec("hgmd"),
        imported_at="2026-06-17T00:00:00Z",
    )

    assert bundle.row_count == 2
    assert bundle.normalized_record_count == 1
    record = bundle.records[0]
    assert record["source_key"] == "hgmd"
    assert record["clinical_significance"] == "Pathogenic"
    assert record["category"] == "clinical_variant"
    serialized = json.dumps(bundle.to_provenance()) + json.dumps(bundle.records)
    assert "do not serialize" not in serialized
    assert "duplicate raw" not in serialized

    json_path = tmp_path / "mastermind.json"
    json_path.write_text(
        json.dumps(
            {
                "records": [
                    {
                        "gene": "GENE1",
                        "variant": "rs123",
                        "article_title": "GENE1 rs123 evidence",
                        "snippet": "Literature summary",
                        "pmid": "PMID:2",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    literature_bundle = parse_source_import(
        "mastermind",
        json_path,
        spec=get_source_spec("mastermind"),
        imported_at="2026-06-17T00:00:00Z",
    )
    assert literature_bundle.records[0]["category"] == "literature"

    bad_json = tmp_path / "bad.json"
    bad_json.write_text(json.dumps({"not_records": {"x": 1}}), encoding="utf-8")
    with pytest.raises(ValueError):
        parse_source_import("hgmd", bad_json, spec=get_source_spec("hgmd"))
    with pytest.raises(ValueError):
        parse_source_import("", csv_path, spec=get_source_spec("hgmd"))


def test_dynamic_builder_uses_user_export_for_licensed_source(tmp_path: Path):
    import_path = tmp_path / "hgmd.csv"
    import_path.write_text(
        "gene,rsid,mutation,classification,disease,pmid,raw_payload\n"
        "GENE1,rs123,c.1A>G,Pathogenic,Example syndrome,PMID:1,raw licensed text\n",
        encoding="utf-8",
    )
    variants = pd.DataFrame([{"chrom": "1", "id": "rs123", "pos": 101, "ref": "A", "alt": "G"}])

    payload = build_dynamic_knowledge_base(
        gene="GENE1",
        region="1:90-110",
        genome_build="hg19",
        variants=variants,
        selected_sources=["hgmd"],
        source_imports={"hgmd": import_path},
        request_client=FakeClient(),
        generated_at="2026-06-17T00:00:00Z",
    )

    assert payload["provider_statuses"][0]["status"] == "imported"
    assert payload["provenance"]["source_imports"][0]["row_count"] == 1
    assert payload["source_records"][0]["clinical_significance"] == "Pathogenic"
    assert payload["variant_records"][0]["evidence"][0]["source_key"] == "hgmd"
    serialized = json.dumps(payload, sort_keys=True)
    assert "raw licensed text" not in serialized


def test_dynamic_builder_is_deterministic_with_fixed_timestamp():
    variants = pd.DataFrame([{"chrom": "2", "id": "rs999", "pos": 200, "ref": "C", "alt": "T"}])

    first = build_dynamic_knowledge_base(
        gene="GENE2",
        region="2:190-210",
        genome_build="hg19",
        variants=variants,
        selected_sources=["clinvar"],
        request_client=FakeClient(),
        generated_at="2026-06-17T00:00:00Z",
    )
    second = build_dynamic_knowledge_base(
        gene="GENE2",
        region="2:190-210",
        genome_build="hg19",
        variants=variants,
        selected_sources=["clinvar"],
        request_client=FakeClient(),
        generated_at="2026-06-17T00:00:00Z",
    )

    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)


def test_dynamic_merger_creates_minimal_bundle_for_uncurated_gene():
    dynamic_payload = {
        "database_name": "Dynamic GENE3 KB",
        "generated_at": "2026-06-17T00:00:00Z",
        "genome_build": "hg19",
        "provider_statuses": [{"source_key": "clinvar", "status": "ok"}],
        "variant_records": [
            {
                "variant": "rs77",
                "display_name": "rs77",
                "lookup_keys": ["rs77", "3:77"],
            }
        ],
        "epigenetic_locus_records": [{"probe_id": "cg77"}],
    }

    merged = merge_dynamic_knowledge_base(
        None,
        dynamic_payload,
        gene_name="GENE3",
        region="3:1-100",
    )

    assert merged["gene_context"]["gene_name"] == "GENE3"
    assert merged["variant_records"][0]["variant"] == "rs77"
    assert merged["gene_context"]["dynamic_knowledge_base"]["provider_count"] == 1
