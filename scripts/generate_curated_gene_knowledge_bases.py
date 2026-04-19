#!/usr/bin/env python3
"""Generate bundled gene knowledge bases and methylation manifest subsets."""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.helper_functions.filter_manifest_region import save_filtered_manifest

GENE_DATA_DIR = PROJECT_ROOT / "src" / "gene_data"
MANIFEST_PATH = PROJECT_ROOT / "data" / "infinium-methylationepic-v-1-0-b5-manifest-file.csv"
ASSEMBLY = "GRCh37 / hg19"
VERSION = "2026-04-15"

COMMON_POPULATION_CATEGORIES = [
    {"key": "Global pattern", "label": "Global pattern"},
    {"key": "Longevity cohorts", "label": "Longevity and healthy-aging cohorts"},
    {"key": "Disease cohorts", "label": "Disease-focused cohorts"},
    {"key": "Cancer cohorts", "label": "Cancer-focused cohorts"},
    {"key": "Metabolic cohorts", "label": "Metabolic and endocrine cohorts"},
    {"key": "Rare disease families", "label": "Rare-disease and familial series"},
]


def _evidence(label: str, url: str) -> dict[str, str]:
    return {"label": label, "url": url}


GENE_DEFINITIONS: list[dict[str, Any]] = [
    {
        "gene_name": "FOXO3",
        "cytoband": "6q21",
        "chromosome": "6",
        "start": 108881028,
        "end": 109005977,
        "strand": "+",
        "gene_summary": (
            "FOXO3 encodes a forkhead-box transcription factor that coordinates stress resistance, "
            "autophagy, apoptosis, stem-cell maintenance, and metabolic adaptation downstream of "
            "insulin-IGF and nutrient-signaling pathways."
        ),
        "clinical_context": (
            "The bundled FOXO3 database is intended for pathway and longevity research context rather "
            "than for monogenic disease reporting. FOXO3 is one of the most replicated human longevity "
            "genes, but common variants are interpreted as modest effect-size healthy-aging modifiers."
        ),
        "variant_effect_overview": [
            "The strongest common FOXO3 signals involve intronic enhancer or haplotype markers rather than protein-truncating alleles.",
            "Published longevity associations are ancestry- and cohort-dependent and should be treated as probabilistic healthy-aging signals, not deterministic predictions.",
            "Observed PASS variants in the FOXO3 interval are therefore best read as locus context unless they directly match a curated longevity-associated marker or have outside clinical evidence.",
        ],
        "condition_research_overview": [
            "FOXO3 is repeatedly linked to exceptional longevity and survival to advanced age across multiple cohort designs.",
            "The gene is also studied in cardiometabolic resilience, oxidative-stress biology, immune regulation, and stem-cell maintenance.",
            "FOXO3 sits in the broader insulin-IGF-AMPK-mTOR longevity network, so interpretation benefits from pathway-level context.",
        ],
        "methylation_interpretation": (
            "FOXO3 methylation is best interpreted as chromatin and regulatory context around a stress-response transcription factor. "
            "Promoter-proximal methylation or broader locus methylation can inform whether the sampled tissue appears permissive or restrained, "
            "but it should not be treated as a direct readout of any one longevity SNP."
        ),
        "methylation_effects": [
            "Promoter-focused FOXO3 methylation is usually interpreted as regulatory context for transcriptional accessibility rather than as a stand-alone biomarker.",
            "Because FOXO3 activity is highly context-dependent, methylation readouts should be integrated with pathway state, tissue type, and age-related biology.",
            "When FOXO3 variants and methylation are discussed together, the most defensible interpretation is pathway tuning rather than one-to-one deterministic coupling.",
        ],
        "methylation_condition_research": [
            "Healthy-aging and lifespan studies that frame FOXO3 as a stress-response resilience factor.",
            "Cardiometabolic and inflammatory research that interprets FOXO3 as a node linking nutrient signaling to cellular maintenance.",
            "Epigenetic aging and chromatin-regulation studies that examine FOXO-family transcriptional accessibility.",
        ],
        "evidence": [
            _evidence("NCBI Gene 2309: FOXO3 gene summary", "https://www.ncbi.nlm.nih.gov/gene/2309"),
            _evidence("PMCID PMC5403515: FOXO3 as a major human longevity gene", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5403515/"),
            _evidence("PubMed 29497356: rs2802292 creates a stress-responsive FOXO3 enhancer", "https://pubmed.ncbi.nlm.nih.gov/29497356/"),
        ],
        "variants": [
            {
                "variant": "rs2802292",
                "common_name": "longevity-associated intronic enhancer SNP",
                "region_class": "gene_body",
                "interpretation_scope": "Research association / healthy-aging modifier",
                "clinical_interpretation": (
                    "rs2802292 is the best-known common FOXO3 longevity marker. The G allele has been associated "
                    "with survival to advanced age in multiple cohorts and is interpreted here as a research-grade "
                    "healthy-aging modifier rather than as a pathogenic or diagnostic allele."
                ),
                "clinical_significance": "Research-level longevity association; not a monogenic disease variant.",
                "functional_effects": [
                    "The longevity-associated allele has been linked to stress-responsive enhancer behavior and higher FOXO3 expression under cellular stress.",
                    "Any effect is best viewed as pathway tuning across stress resistance and metabolic adaptation rather than as deterministic causality.",
                ],
                "associated_conditions": [
                    "Exceptional longevity and survival to advanced age",
                    "Cardiometabolic resilience and lower age-related disease burden in some cohorts",
                    "Stress-response and cellular-maintenance biology",
                ],
                "research_context": [
                    "This is the canonical common FOXO3 longevity SNP used in centenarian and healthy-aging studies.",
                    "Effect size is modest and cohort dependent, so replication and ancestry matching matter.",
                ],
                "usual_variant_note": "Most cited common FOXO3 longevity marker.",
                "methylation_interpretation": (
                    "Treat rs2802292 alongside FOXO3 methylation as shared regulatory context rather than as a direct methylation biomarker."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PMCID PMC5403515: FOXO3 longevity review", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5403515/"),
                    _evidence("PubMed 29497356: rs2802292 enhancer mechanism", "https://pubmed.ncbi.nlm.nih.gov/29497356/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Morris et al., 2015 (PMCID PMC5403515)",
                        "genotypes": "GG, GT, and TT",
                        "phenotype": "Human longevity across replicated cohort studies",
                        "finding": "The G allele is repeatedly described as enriched in long-lived individuals, supporting a modest but reproducible healthy-aging association.",
                        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5403515/",
                    },
                    {
                        "paper": "Flachsbart et al., 2017/2018 (PMID 29497356)",
                        "genotypes": "G-allele carriers versus non-carriers",
                        "phenotype": "Stress-responsive FOXO3 expression control",
                        "finding": "Mechanistic follow-up work proposed that the rs2802292 G allele helps create or stabilize a stress-responsive enhancer element that promotes FOXO3 upregulation.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/29497356/",
                    },
                ],
            }
        ],
        "population_intro": "Broader population patterns curated from FOXO3 healthy-aging, longevity, and stress-response literature.",
        "population_coverage_note": (
            "This FOXO3 population database emphasizes cohort-level longevity patterns and ancestry-aware interpretation. "
            "It does not yet ship a full embedded allele-frequency panel for every reported FOXO3 marker."
        ),
        "population_sources": [
            _evidence("NCBI Gene 2309: FOXO3 gene summary", "https://www.ncbi.nlm.nih.gov/gene/2309"),
            _evidence("PMCID PMC5403515: FOXO3 longevity review", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5403515/"),
            _evidence("PubMed 29497356: rs2802292 enhancer mechanism", "https://pubmed.ncbi.nlm.nih.gov/29497356/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "rs2802292",
                "location_group": "Longevity cohorts",
                "summary": "FOXO3 longevity signals are most often discussed in centenarian or long-lived cohort designs, where rs2802292 is treated as a modest survival-enrichment marker rather than a deterministic lifespan allele.",
                "evidence": [
                    _evidence("PMCID PMC5403515", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5403515/"),
                ],
            },
            {
                "variant": "FOXO3 longevity haplotypes",
                "location_group": "Global pattern",
                "summary": "Different populations tag the FOXO3 longevity signal with different SNPs or haplotypes, so cross-cohort interpretation should focus on the locus-level signal rather than one ancestry-specific marker alone.",
                "evidence": [
                    _evidence("PMCID PMC5403515", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5403515/"),
                ],
            },
        ],
    },
    {
        "gene_name": "MTOR",
        "cytoband": "1p36.22",
        "chromosome": "1",
        "start": 11166592,
        "end": 11322608,
        "strand": "-",
        "gene_summary": (
            "MTOR encodes the catalytic kinase shared by mTORC1 and mTORC2, integrating nutrient availability, growth-factor signaling, "
            "stress cues, and cellular energy state to regulate protein synthesis, autophagy, survival, and cytoskeletal programs."
        ),
        "clinical_context": (
            "The local MTOR knowledge base is designed for signaling and research interpretation. Pathogenic MTOR variants do exist in developmental "
            "and focal cortical dysplasia syndromes, but the bundled records here focus on commonly discussed regulatory polymorphisms and pathway context."
        ),
        "variant_effect_overview": [
            "Common MTOR polymorphisms are mainly interpreted as expression or signaling modifiers with context-dependent effect sizes.",
            "The strongest disease literature around MTOR centers on pathway dysregulation, drug response, and rare activating or inactivating variants rather than on one universal common SNP.",
            "Observed variants in this region should therefore be read as pathway-context markers unless they directly overlap a curated research SNP or have outside clinical curation.",
        ],
        "condition_research_overview": [
            "MTOR sits at the center of growth, anabolic signaling, and autophagy control, making it a core aging and cancer pathway gene.",
            "Clinical research links MTOR to immunosuppression, targeted oncology, developmental syndromes, and neurodevelopmental disease.",
            "Common regulatory SNP work is mostly cancer-risk or outcome oriented rather than diagnostic on its own.",
        ],
        "methylation_interpretation": (
            "MTOR methylation is best treated as pathway-regulatory context around a very large signaling locus. "
            "Promoter-proximal methylation may reflect transcriptional accessibility, but most interpretation remains pathway level rather than variant specific."
        ),
        "methylation_effects": [
            "Promoter-focused methylation may reflect altered MTOR expression potential, but it is not a stand-alone disease classifier.",
            "Given the size of the locus and multiple transcripts, methylation should be interpreted as a broad regulatory view rather than a one-site assay.",
            "MTOR methylation is most meaningful when read alongside AMPK, nutrient-sensing, and growth-signaling biology.",
        ],
        "methylation_condition_research": [
            "Cancer and targeted-therapy studies that frame mTOR activity as a proliferation and drug-response axis.",
            "Aging and nutrient-signaling research focused on mTORC1 inhibition, autophagy, and proteostasis.",
            "Neurodevelopmental and cortical dysplasia research in which mTOR pathway dysregulation is a major mechanistic theme.",
        ],
        "evidence": [
            _evidence("NCBI Gene 2475: MTOR gene summary", "https://www.ncbi.nlm.nih.gov/gene/2475"),
            _evidence("PubMed 24816861: meta-analysis of MTOR polymorphisms and cancer risk", "https://pubmed.ncbi.nlm.nih.gov/24816861/"),
            _evidence("PubMed 29978580: rs2536 and gastric-cancer survival", "https://pubmed.ncbi.nlm.nih.gov/29978580/"),
        ],
        "variants": [
            {
                "variant": "rs2295080",
                "common_name": "promoter regulatory MTOR SNP",
                "region_class": "promoter",
                "interpretation_scope": "Research association / regulatory marker",
                "clinical_interpretation": (
                    "rs2295080 is a promoter-proximal MTOR polymorphism studied as a transcriptional and cancer-risk modifier. "
                    "The local database treats it as a regulatory research marker rather than as a pathogenic clinical allele."
                ),
                "clinical_significance": "Research-level regulatory association; not seeded as pathogenic.",
                "functional_effects": [
                    "Promoter variation at rs2295080 has been discussed as a potential modifier of MTOR transcriptional output.",
                    "Reported effects are strongest in cancer-association studies and remain cohort dependent.",
                ],
                "associated_conditions": [
                    "Cancer susceptibility and clinical outcome studies",
                    "Pathway-level growth and survival signaling research",
                ],
                "research_context": [
                    "This SNP is commonly bundled with other MTOR polymorphisms in association meta-analyses.",
                    "Interpret it as a modest signaling modifier rather than a stand-alone disease determinant.",
                ],
                "usual_variant_note": "Common promoter-facing MTOR regulatory SNP.",
                "methylation_interpretation": (
                    "Promoter methylation and rs2295080 should be read together only as shared transcriptional context."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PubMed 24816861: MTOR polymorphism meta-analysis", "https://pubmed.ncbi.nlm.nih.gov/24816861/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Xu et al., 2014 (PMID 24816861)",
                        "genotypes": "TT versus GT/GG",
                        "phenotype": "Cancer risk across pooled association studies",
                        "finding": "The meta-analysis reported that the rs2295080 TT genotype was associated with higher cancer risk in the Chinese studies included, supporting a modest context-specific regulatory effect.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/24816861/",
                    }
                ],
            },
            {
                "variant": "rs2536",
                "common_name": "3'UTR functional MTOR SNP",
                "region_class": "gene_body",
                "interpretation_scope": "Research association / outcome modifier",
                "clinical_interpretation": (
                    "rs2536 is a commonly discussed MTOR 3'UTR polymorphism curated here as a research association marker for cancer susceptibility and outcome, "
                    "not as a diagnostic or pathogenic allele."
                ),
                "clinical_significance": "Research-level association; best interpreted in cohort context.",
                "functional_effects": [
                    "The variant is discussed as a post-transcriptional regulatory marker in the MTOR pathway.",
                    "Outcome studies suggest it may stratify prognosis in selected cancer cohorts.",
                ],
                "associated_conditions": [
                    "Gastric cancer risk and survival studies",
                    "Broader oncology association meta-analyses",
                ],
                "research_context": [
                    "This is one of the most frequently cited common MTOR SNPs in outcome-oriented cancer literature.",
                ],
                "usual_variant_note": "Common MTOR 3'UTR association SNP.",
                "methylation_interpretation": (
                    "No one-to-one methylation consequence is bundled for rs2536; treat methylation as broader MTOR locus context."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PubMed 29978580: rs2536 and gastric-cancer survival", "https://pubmed.ncbi.nlm.nih.gov/29978580/"),
                    _evidence("PubMed 24816861: MTOR polymorphism meta-analysis", "https://pubmed.ncbi.nlm.nih.gov/24816861/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Cheng et al., 2019 (PMID 29978580)",
                        "genotypes": "Common rs2536 genotypes",
                        "phenotype": "Survival of Chinese gastric-cancer patients",
                        "finding": "The study reported that rs2536 was functionally relevant in survival modeling for gastric-cancer patients, reinforcing its use as a regulatory and prognostic research marker.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/29978580/",
                    }
                ],
            },
        ],
        "population_intro": "Broader population patterns curated from MTOR cancer-association, therapeutic-response, and signaling literature.",
        "population_coverage_note": (
            "The bundled MTOR population database is literature oriented. It currently summarizes cohort patterns and pathway interpretation rather than shipping a full allele-frequency reference panel."
        ),
        "population_sources": [
            _evidence("NCBI Gene 2475: MTOR gene summary", "https://www.ncbi.nlm.nih.gov/gene/2475"),
            _evidence("PubMed 24816861: MTOR polymorphism meta-analysis", "https://pubmed.ncbi.nlm.nih.gov/24816861/"),
            _evidence("PubMed 29978580: rs2536 survival study", "https://pubmed.ncbi.nlm.nih.gov/29978580/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "rs2295080 and rs2536",
                "location_group": "Cancer cohorts",
                "summary": "Most common-SNP MTOR literature comes from cancer cohorts in East Asian populations, where promoter and 3'UTR markers are interpreted as modest expression or prognosis modifiers.",
            },
            {
                "variant": "Common MTOR polymorphisms",
                "location_group": "Global pattern",
                "summary": "Because mTOR biology is heavily pathway driven, cross-population interpretation should emphasize signaling context and cohort design rather than assuming a universal effect from one common SNP.",
            },
        ],
    },
    {
        "gene_name": "RPS6",
        "cytoband": "9p22.1",
        "chromosome": "9",
        "start": 19375713,
        "end": 19380234,
        "strand": "-",
        "gene_summary": (
            "RPS6 encodes the ribosomal protein S6, a core 40S subunit protein and major substrate of S6 kinases downstream of mTOR signaling. "
            "Its phosphorylation state is widely used as a readout of growth-factor and nutrient-responsive translation control."
        ),
        "clinical_context": (
            "The bundled RPS6 knowledge base is pathway oriented. RPS6 is typically interpreted as a translational-control and mTOR-output node rather than as a common-SNP clinical gene in routine VCF review."
        ),
        "variant_effect_overview": [
            "RPS6 is most commonly discussed as a signaling substrate and biomarker of pathway activity rather than as a locus with a broad common-variant literature.",
            "Observed variants in the compact RPS6 interval should therefore be treated as locus context unless outside evidence supports a specific interpretation.",
        ],
        "condition_research_overview": [
            "RPS6 phosphorylation is a canonical downstream readout of mTORC1-S6K activity in cancer and growth biology.",
            "The gene is relevant to translational control, cell growth, proliferation, and treatment-response biomarker work.",
        ],
        "methylation_interpretation": (
            "RPS6 methylation should be read as local chromatin context around a compact ribosomal-protein locus. "
            "Because the biology of interest is often protein phosphorylation rather than gene methylation, these values are mainly supportive context."
        ),
        "methylation_effects": [
            "Promoter methylation can provide a rough view of transcriptional accessibility, but most RPS6 functional work centers on protein phosphorylation and translational output.",
            "Interpret methylation alongside mTOR-S6K pathway activity rather than as a stand-alone biomarker.",
        ],
        "methylation_condition_research": [
            "Cancer biomarker studies using phospho-S6 as an mTOR-pathway activity readout.",
            "Translational-control and growth-signaling research.",
        ],
        "evidence": [
            _evidence("NCBI Gene 6194: RPS6 gene summary", "https://www.ncbi.nlm.nih.gov/gene/6194"),
            _evidence("PubMed 18157089: phospho-S6 as a marker for targeted mTOR therapy", "https://pubmed.ncbi.nlm.nih.gov/18157089/"),
        ],
        "variants": [],
        "population_intro": "Broader population patterns curated from RPS6 pathway-biomarker and translational-control literature.",
        "population_coverage_note": (
            "Clinically relevant RPS6 interpretation is usually pathway and biomarker driven rather than based on common allele-frequency panels."
        ),
        "population_sources": [
            _evidence("NCBI Gene 6194: RPS6 gene summary", "https://www.ncbi.nlm.nih.gov/gene/6194"),
            _evidence("PubMed 18157089: phospho-S6 biomarker paper", "https://pubmed.ncbi.nlm.nih.gov/18157089/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "RPS6 pathway context",
                "location_group": "Global pattern",
                "summary": "RPS6 is generally interpreted through pathway activation and phospho-S6 biomarker studies rather than through a dense catalog of common inherited markers.",
            }
        ],
    },
    {
        "gene_name": "SIK3",
        "cytoband": "11q23.3",
        "chromosome": "11",
        "start": 116714118,
        "end": 116969144,
        "strand": "-",
        "gene_summary": (
            "SIK3 encodes a salt-inducible serine-threonine kinase in the AMPK-related kinase family that participates in metabolic control, "
            "sleep regulation, transcriptional programs, and TOR-linked signaling."
        ),
        "clinical_context": (
            "The local SIK3 knowledge base focuses on emerging functional and sleep-trait biology. The locus has a growing GWAS and mechanistic literature, "
            "but interpretation is still mainly research grade for most variants."
        ),
        "variant_effect_overview": [
            "SIK3 is better characterized through pathway studies, sleep phenotypes, and selected rare functional variants than through a mature clinical common-SNP catalog.",
            "Observed variants in the SIK3 interval should be interpreted as research context unless they match a seeded functional record or have outside clinical curation.",
        ],
        "condition_research_overview": [
            "Human and model-organism studies link SIK3 to sleep need, circadian or sleep-duration biology, and neuronal regulation.",
            "Additional research connects SIK3 to metabolism, adiposity, and oncogenic or transcriptional signaling programs.",
        ],
        "methylation_interpretation": (
            "SIK3 methylation is best treated as regulatory context around an emerging signaling gene. "
            "Interpret any promoter-proximal patterns as supportive evidence rather than as a direct proxy for sleep or metabolic phenotype."
        ),
        "methylation_effects": [
            "Methylation may help contextualize transcriptional accessibility in tissues where SIK3 is active, but functional interpretation remains pathway level.",
            "Sleep-trait and signaling phenotypes should not be inferred from methylation alone.",
        ],
        "methylation_condition_research": [
            "Sleep and circadian-biology studies.",
            "Metabolic-signaling and adiposity research.",
            "Transcriptional-control and leukemia-dependency studies involving SIK kinases.",
        ],
        "evidence": [
            _evidence("NCBI Gene 23387: SIK3 gene summary", "https://www.ncbi.nlm.nih.gov/gene/23387"),
            _evidence("DOI 10.1073/pnas.2500356122: SIK3-N783Y and natural short sleep", "https://doi.org/10.1073/pnas.2500356122"),
            _evidence("PubMed 32126566: SIK3 haplotypes and noise-induced hearing loss", "https://pubmed.ncbi.nlm.nih.gov/32126566/"),
        ],
        "variants": [
            {
                "variant": "SIK3 p.Asn783Tyr",
                "display_name": "SIK3 p.Asn783Tyr (N783Y)",
                "common_name": "natural short sleep-associated missense variant",
                "region_class": "gene_body",
                "interpretation_scope": "Rare functional sleep-trait variant",
                "clinical_interpretation": (
                    "SIK3 p.Asn783Tyr is a rare functional missense variant reported in association with the human natural short sleep trait. "
                    "It is curated here as a mechanistic sleep-biology variant, not as a routine pathogenic diagnostic allele."
                ),
                "clinical_significance": "Rare functional research variant with emerging human sleep-trait evidence.",
                "functional_effects": [
                    "The variant is discussed as a gain-of-function or altered-function SIK3 allele in sleep regulation.",
                    "Its interpretation is primarily mechanistic and research facing at present.",
                ],
                "associated_conditions": [
                    "Human natural short sleep trait",
                    "Sleep-need and neuronal signaling research",
                ],
                "research_context": [
                    "This record captures a rare high-interest functional variant rather than a common population marker.",
                ],
                "usual_variant_note": "Emerging rare SIK3 human sleep variant.",
                "methylation_interpretation": (
                    "SIK3 methylation provides broad locus context but does not directly read out the functional consequence of p.Asn783Tyr."
                ),
                "is_assayable_in_snp_vcf": False,
                "evidence": [
                    _evidence("DOI 10.1073/pnas.2500356122: SIK3-N783Y and natural short sleep", "https://doi.org/10.1073/pnas.2500356122"),
                ],
                "literature_findings": [
                    {
                        "paper": "PNAS, 2025 (DOI 10.1073/pnas.2500356122)",
                        "genotypes": "Rare missense carrier state",
                        "phenotype": "Human natural short sleep",
                        "finding": "The report linked SIK3 N783Y to the human natural short sleep trait and framed the variant as a direct functional probe of sleep-need biology.",
                        "url": "https://doi.org/10.1073/pnas.2500356122",
                    }
                ],
            }
        ],
        "population_intro": "Broader population patterns curated from SIK3 sleep, hearing, and metabolic-signaling literature.",
        "population_coverage_note": (
            "SIK3 population interpretation currently depends more on cohort-specific association studies and rare-variant reports than on a mature embedded allele-frequency panel."
        ),
        "population_sources": [
            _evidence("NCBI Gene 23387: SIK3 gene summary", "https://www.ncbi.nlm.nih.gov/gene/23387"),
            _evidence("PubMed 32126566: SIK3 haplotypes and hearing-loss risk", "https://pubmed.ncbi.nlm.nih.gov/32126566/"),
            _evidence("DOI 10.1073/pnas.2500356122: SIK3-N783Y and natural short sleep", "https://doi.org/10.1073/pnas.2500356122"),
        ],
        "gene_population_patterns": [
            {
                "variant": "SIK3 haplotypes",
                "location_group": "Disease cohorts",
                "summary": "Published common-variant SIK3 literature is still cohort specific, including hearing-loss and adiposity studies, so the locus is best interpreted as emerging rather than fully settled population genetics.",
            },
            {
                "variant": "SIK3 p.Asn783Tyr",
                "location_group": "Rare disease families",
                "summary": "The most striking human SIK3 result to date is a rare functional sleep-trait variant, underscoring that this gene is currently driven by mechanistic and rare-variant evidence more than by broad common-SNP panels.",
            },
        ],
    },
    {
        "gene_name": "FLCN",
        "cytoband": "17p11.2",
        "chromosome": "17",
        "start": 17115526,
        "end": 17140482,
        "strand": "-",
        "gene_summary": (
            "FLCN encodes folliculin, a tumor suppressor involved in lysosomal signaling, AMPK-mTOR regulation, TFEB-TFE3 control, and renal and pulmonary homeostasis."
        ),
        "clinical_context": (
            "The strongest clinical interpretation for FLCN is rare-disease oriented. Germline loss-of-function variants cause Birt-Hogg-Dube syndrome, "
            "so the local database emphasizes rare pathogenic and tumor-suppressor biology rather than common-SNP association signals."
        ),
        "variant_effect_overview": [
            "The clinically important FLCN variants are usually rare truncating, splice, or deletion events rather than common polymorphisms.",
            "Observed PASS variants in the FLCN interval should therefore be treated as locus context unless they overlap a known rare pathogenic event or have independent clinical curation.",
        ],
        "condition_research_overview": [
            "FLCN is the canonical gene for Birt-Hogg-Dube syndrome with renal-tumor, pulmonary-cyst, and pneumothorax risk.",
            "The gene is also studied in lysosome signaling, TFEB/TFE3 regulation, and AMPK-mTOR pathway cross-talk.",
        ],
        "methylation_interpretation": (
            "FLCN methylation is best read as tumor-suppressor regulatory context. It can help frame whether the locus appears transcriptionally restrained or accessible, "
            "but it does not replace rare-variant interpretation in Birt-Hogg-Dube syndrome."
        ),
        "methylation_effects": [
            "For FLCN, rare coding or structural variants usually carry more clinical weight than methylation patterns.",
            "Methylation can still provide useful context in renal-tumor and tumor-suppressor regulation studies.",
        ],
        "methylation_condition_research": [
            "Birt-Hogg-Dube syndrome tumor-suppressor biology.",
            "Renal neoplasia and pneumothorax predisposition research.",
        ],
        "evidence": [
            _evidence("NCBI Gene 201163: FLCN gene summary", "https://www.ncbi.nlm.nih.gov/gene/201163"),
            _evidence("PubMed 21538689: BHD-associated FLCN mutations disrupt protein stability", "https://pubmed.ncbi.nlm.nih.gov/21538689/"),
            _evidence("PubMed 22146830: renal cancer and pneumothorax risk in FLCN mutation carriers", "https://pubmed.ncbi.nlm.nih.gov/22146830/"),
        ],
        "variants": [
            {
                "variant": "FLCN c.1285dupC",
                "display_name": "FLCN c.1285dupC",
                "common_name": "Birt-Hogg-Dube frameshift hotspot",
                "region_class": "gene_body",
                "interpretation_scope": "Rare pathogenic tumor-suppressor variant",
                "clinical_interpretation": (
                    "c.1285dupC is one of the best-known recurrent FLCN frameshift variants in Birt-Hogg-Dube syndrome and is interpreted as a pathogenic rare-disease allele. "
                    "It is bundled here as a named clinical hotspot for context, although the SNP-oriented VCF preview is not optimized to call it directly."
                ),
                "clinical_significance": "Pathogenic loss-of-function hotspot in Birt-Hogg-Dube syndrome.",
                "functional_effects": [
                    "Frameshift loss of folliculin function with tumor-suppressor consequences.",
                    "Expected to disrupt lysosomal and AMPK-mTOR regulatory roles of FLCN.",
                ],
                "associated_conditions": [
                    "Birt-Hogg-Dube syndrome",
                    "Renal tumors, lung cysts, and spontaneous pneumothorax",
                ],
                "research_context": [
                    "This is a rare pathogenic hotspot rather than a common association marker.",
                ],
                "usual_variant_note": "Classic recurrent BHD loss-of-function variant.",
                "methylation_interpretation": (
                    "FLCN methylation can add tumor-suppressor context, but it does not substitute for direct rare-variant testing when Birt-Hogg-Dube syndrome is suspected."
                ),
                "is_assayable_in_snp_vcf": False,
                "evidence": [
                    _evidence("PubMed 21538689: BHD-associated FLCN mutations", "https://pubmed.ncbi.nlm.nih.gov/21538689/"),
                    _evidence("PubMed 22146830: phenotype risk in FLCN carriers", "https://pubmed.ncbi.nlm.nih.gov/22146830/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Nahorski et al., 2011 (PMID 21538689)",
                        "genotypes": "Rare germline loss-of-function carrier state",
                        "phenotype": "Birt-Hogg-Dube syndrome protein-stability consequences",
                        "finding": "The study framed recurrent FLCN truncating mutations as protein-destabilizing tumor-suppressor events central to BHD pathogenesis.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/21538689/",
                    }
                ],
            }
        ],
        "population_intro": "Broader population patterns curated from FLCN rare-disease, renal-tumor, and pneumothorax literature.",
        "population_coverage_note": (
            "FLCN population interpretation is primarily family and rare-variant driven. The current database therefore emphasizes rare-disease patterns rather than common allele-frequency panels."
        ),
        "population_sources": [
            _evidence("NCBI Gene 201163: FLCN gene summary", "https://www.ncbi.nlm.nih.gov/gene/201163"),
            _evidence("PubMed 22146830: FLCN carrier phenotype analysis", "https://pubmed.ncbi.nlm.nih.gov/22146830/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "Rare FLCN loss-of-function variants",
                "location_group": "Rare disease families",
                "summary": "FLCN clinical interpretation is dominated by rare germline truncating or deletion events that segregate in Birt-Hogg-Dube families rather than by common population markers.",
            },
            {
                "variant": "Birt-Hogg-Dube syndrome",
                "location_group": "Global pattern",
                "summary": "Across populations, the major interpretation question for FLCN is whether a rare loss-of-function event is present and whether renal and pulmonary surveillance is warranted, not whether a common SNP modestly shifts risk.",
            },
        ],
    },
    {
        "gene_name": "SIRT6",
        "cytoband": "19p13.3",
        "chromosome": "19",
        "start": 4174106,
        "end": 4182560,
        "strand": "-",
        "gene_summary": (
            "SIRT6 encodes a nuclear NAD-dependent deacylase and mono-ADP-ribosyltransferase that supports DNA repair, telomere maintenance, transposon suppression, inflammation control, and energy homeostasis."
        ),
        "clinical_context": (
            "The bundled SIRT6 knowledge base is centered on aging and genome-maintenance biology. Most currently relevant human variants are interpreted as research-grade longevity or regulatory signals rather than as routine pathogenic clinical alleles."
        ),
        "variant_effect_overview": [
            "SIRT6 is a core longevity and genome-stability gene, so common and rare human variants are mainly interpreted through aging and cellular-maintenance phenotypes.",
            "Many effects are subtle or cohort specific; the strongest current evidence pairs genetic variation with mechanistic functional follow-up.",
        ],
        "condition_research_overview": [
            "SIRT6 is heavily studied in aging, genome stability, inflammation, lipid and glucose metabolism, and cardiovascular disease.",
            "Rare centenarian-enriched alleles and selected common SNPs are being used to connect human genetics with SIRT6 functional biochemistry.",
        ],
        "methylation_interpretation": (
            "SIRT6 methylation should be interpreted as chromatin and transcriptional context around an NAD-dependent longevity enzyme. "
            "Because SIRT6 function is closely tied to chromatin regulation, methylation can be informative, but not as a one-to-one biomarker."
        ),
        "methylation_effects": [
            "Promoter methylation may reflect altered SIRT6 transcriptional accessibility and should be integrated with broader NAD and stress-response biology.",
            "Functional inference still depends more on enzyme activity and downstream chromatin effects than on methylation alone.",
        ],
        "methylation_condition_research": [
            "Aging and longevity research.",
            "Genome-stability and DNA-repair studies.",
            "Cardiometabolic and inflammatory disease biology.",
        ],
        "evidence": [
            _evidence("NCBI Gene 51548: SIRT6 gene summary", "https://www.ncbi.nlm.nih.gov/gene/51548"),
            _evidence("PubMed 28032059: rs350846 and human longevity", "https://pubmed.ncbi.nlm.nih.gov/28032059/"),
            _evidence("PubMed 36215696: centenarian SIRT6 variant", "https://pubmed.ncbi.nlm.nih.gov/36215696/"),
            _evidence("PubMed 36465178: SIRT6 review in aging and metabolism", "https://pubmed.ncbi.nlm.nih.gov/36465178/"),
        ],
        "variants": [
            {
                "variant": "rs350846",
                "common_name": "longevity-associated SIRT6 SNP",
                "region_class": "gene_body",
                "interpretation_scope": "Research association / longevity marker",
                "clinical_interpretation": (
                    "rs350846 is a common SIRT6 polymorphism curated here as a research-grade human longevity marker. "
                    "The local database does not treat it as a diagnostic allele, but as a cohort-dependent aging association."
                ),
                "clinical_significance": "Research-level longevity association.",
                "functional_effects": [
                    "Reported as a human longevity-associated marker at the SIRT6 locus.",
                    "Likely reflects pathway-level modulation of SIRT6 expression or linked regulatory state rather than a deterministic coding lesion.",
                ],
                "associated_conditions": [
                    "Human longevity and healthy aging",
                    "Genome-maintenance and metabolic-aging biology",
                ],
                "research_context": [
                    "This SNP is best interpreted through cohort and population context rather than as a stand-alone predictive marker.",
                ],
                "usual_variant_note": "Most cited common SIRT6 longevity SNP.",
                "methylation_interpretation": (
                    "SIRT6 methylation provides broader chromatin context and should not be treated as a direct biomarker of rs350846 status."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PubMed 28032059: rs350846 and human longevity", "https://pubmed.ncbi.nlm.nih.gov/28032059/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Li et al., 2016 (PMID 28032059)",
                        "genotypes": "CC, CG, and GG",
                        "phenotype": "Human longevity in long-lived versus control cohorts",
                        "finding": "The study identified rs350846 as a SIRT6 locus polymorphism associated with human longevity, supporting its use as a research-grade aging marker.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/28032059/",
                    }
                ],
            },
            {
                "variant": "centSIRT6 (N308K/A313S)",
                "display_name": "centSIRT6 (N308K/A313S)",
                "common_name": "centenarian-enriched SIRT6 linked missense allele",
                "region_class": "gene_body",
                "interpretation_scope": "Rare functional longevity variant",
                "clinical_interpretation": (
                    "The linked SIRT6 missense substitutions N308K and A313S define a rare centenarian-enriched allele with enhanced genome-maintenance functions in experimental follow-up. "
                    "It is curated as a rare functional longevity variant, not as a routine clinical diagnostic allele."
                ),
                "clinical_significance": "Rare functional longevity-associated research allele.",
                "functional_effects": [
                    "Enhanced mono-ADP-ribosylase activity and stronger interaction with Lamin A/C were reported in functional studies.",
                    "The allele improved LINE1 suppression and DNA double-strand break repair in experimental systems.",
                ],
                "associated_conditions": [
                    "Exceptional longevity",
                    "Genome stability and stress-resistance biology",
                ],
                "research_context": [
                    "This is a rare mechanistic allele of high interest because it links centenarian enrichment to direct SIRT6 functional changes.",
                ],
                "usual_variant_note": "Rare centenarian SIRT6 functional allele.",
                "methylation_interpretation": (
                    "Methylation provides broad chromatin context for SIRT6 but does not directly substitute for sequencing or functional assessment of the centSIRT6 allele."
                ),
                "is_assayable_in_snp_vcf": False,
                "evidence": [
                    _evidence("PubMed 36215696: centenarian SIRT6 variant", "https://pubmed.ncbi.nlm.nih.gov/36215696/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Simon et al., 2022 (PMID 36215696)",
                        "genotypes": "Rare linked missense carrier state",
                        "phenotype": "Centenarian enrichment and genome-stability phenotypes",
                        "finding": "The centSIRT6 allele was enriched in Ashkenazi Jewish centenarians and displayed enhanced DNA repair, LINE1 suppression, and cancer-cell killing in functional assays.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/36215696/",
                    }
                ],
            },
        ],
        "population_intro": "Broader population patterns curated from SIRT6 longevity, genome-stability, and metabolic-aging literature.",
        "population_coverage_note": (
            "The bundled SIRT6 population database emphasizes longevity cohorts and rare centenarian alleles instead of a full general-population allele-frequency reference panel."
        ),
        "population_sources": [
            _evidence("NCBI Gene 51548: SIRT6 gene summary", "https://www.ncbi.nlm.nih.gov/gene/51548"),
            _evidence("PubMed 28032059: rs350846 and longevity", "https://pubmed.ncbi.nlm.nih.gov/28032059/"),
            _evidence("PubMed 36215696: centenarian SIRT6 variant", "https://pubmed.ncbi.nlm.nih.gov/36215696/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "rs350846",
                "location_group": "Longevity cohorts",
                "summary": "Common SIRT6 longevity signals such as rs350846 are cohort dependent and are best interpreted as probabilistic aging markers rather than universally penetrant alleles.",
            },
            {
                "variant": "centSIRT6 (N308K/A313S)",
                "location_group": "Rare disease families",
                "summary": "The most compelling current SIRT6 human-genetics result is a rare centenarian-enriched functional allele, underscoring that the locus is shaped by both common regulatory and rare mechanistic variation.",
            },
        ],
    },
    {
        "gene_name": "PRKAA1",
        "cytoband": "5p13.1",
        "chromosome": "5",
        "start": 40759491,
        "end": 40798476,
        "strand": "-",
        "gene_summary": (
            "PRKAA1 encodes the alpha-1 catalytic subunit of AMPK, the major cellular energy sensor that restrains anabolic signaling, supports metabolic adaptation, and suppresses mTORC1 under low-energy conditions."
        ),
        "clinical_context": (
            "The local PRKAA1 knowledge base is built around metabolic and cancer-association literature. "
            "Common PRKAA1 variants are best interpreted as pathway modifiers rather than as highly penetrant clinical alleles."
        ),
        "variant_effect_overview": [
            "Common PRKAA1 variants are typically interpreted as low-penetrance metabolic or cancer-risk modifiers.",
            "Because PRKAA1 acts inside the AMPK-mTOR energy-stress axis, published associations often reflect pathway state and environmental context.",
        ],
        "condition_research_overview": [
            "PRKAA1 is central to AMPK signaling, nutrient sensing, mitochondrial stress responses, and mTOR restraint.",
            "The strongest common-variant literature centers on gastric-cancer susceptibility and metabolic-endocrine phenotypes.",
        ],
        "methylation_interpretation": (
            "PRKAA1 methylation should be read as energy-sensing pathway context. Promoter-proximal methylation may suggest altered catalytic-subunit expression, "
            "but interpretation remains broader than any one common SNP."
        ),
        "methylation_effects": [
            "Promoter methylation may help contextualize PRKAA1 transcriptional tone in metabolism-focused tissues.",
            "Integrate methylation with AMPK-mTOR signaling, energy stress, and nutritional context.",
        ],
        "methylation_condition_research": [
            "Energy-stress and nutrient-sensing biology.",
            "Cancer, especially gastric-cancer association work.",
            "Metabolic and endocrine signaling research.",
        ],
        "evidence": [
            _evidence("NCBI Gene 5562: PRKAA1 gene summary", "https://www.ncbi.nlm.nih.gov/gene/5562"),
            _evidence("PubMed 26485766: PRKAA1 rs13361707 and gastric-cancer risk", "https://pubmed.ncbi.nlm.nih.gov/26485766/"),
        ],
        "variants": [
            {
                "variant": "rs13361707",
                "common_name": "gastric-cancer-associated PRKAA1 SNP",
                "region_class": "gene_body",
                "interpretation_scope": "Research association / low-penetrance risk marker",
                "clinical_interpretation": (
                    "rs13361707 is the best-established common PRKAA1 marker in gastric-cancer association literature. "
                    "The local database treats it as a low-penetrance pathway-associated risk marker rather than as a pathogenic clinical allele."
                ),
                "clinical_significance": "Research-level low-penetrance association marker.",
                "functional_effects": [
                    "Likely acts through regulatory or linked-locus effects on AMPK-pathway signaling rather than as a direct protein-disrupting lesion.",
                    "Effect sizes reported in cancer cohorts are modest and ancestry specific.",
                ],
                "associated_conditions": [
                    "Gastric-cancer susceptibility",
                    "Metabolic and stress-signaling pathway context",
                ],
                "research_context": [
                    "This is the canonical common PRKAA1 association SNP in East Asian gastric-cancer studies.",
                ],
                "usual_variant_note": "Most cited common PRKAA1 association SNP.",
                "methylation_interpretation": (
                    "PRKAA1 methylation should be read as broader energy-sensing locus context rather than as a direct biomarker of rs13361707."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PubMed 26485766: PRKAA1 rs13361707 and gastric-cancer risk", "https://pubmed.ncbi.nlm.nih.gov/26485766/"),
                    _evidence("PubMed 37670284: PRKAA1 polymorphisms and gastric cancer", "https://pubmed.ncbi.nlm.nih.gov/37670284/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Chen et al., 2015 (PMID 26485766)",
                        "genotypes": "TT, CT, and CC",
                        "phenotype": "Gastric-cancer susceptibility in an eastern Chinese population",
                        "finding": "The study reported that the rs13361707 C allele increased gastric-cancer risk, consistent with its use as a low-penetrance PRKAA1 research marker.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/26485766/",
                    }
                ],
            }
        ],
        "population_intro": "Broader population patterns curated from PRKAA1 gastric-cancer, metabolic, and AMPK-pathway literature.",
        "population_coverage_note": (
            "This PRKAA1 population database emphasizes cohort-level association patterns instead of a full embedded allele-frequency panel."
        ),
        "population_sources": [
            _evidence("NCBI Gene 5562: PRKAA1 gene summary", "https://www.ncbi.nlm.nih.gov/gene/5562"),
            _evidence("PubMed 26485766: PRKAA1 rs13361707 and gastric-cancer risk", "https://pubmed.ncbi.nlm.nih.gov/26485766/"),
            _evidence("PubMed 37670284: PRKAA1 polymorphism follow-up", "https://pubmed.ncbi.nlm.nih.gov/37670284/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "rs13361707",
                "location_group": "Cancer cohorts",
                "summary": "PRKAA1 common-variant literature is strongest in East Asian gastric-cancer cohorts, where rs13361707 is interpreted as a low-penetrance risk marker at the AMPK pathway locus.",
            },
            {
                "variant": "PRKAA1 pathway context",
                "location_group": "Metabolic cohorts",
                "summary": "Across populations, PRKAA1 interpretation benefits from energy-stress and nutrient-signaling context rather than from assuming a uniform common-SNP effect across phenotypes.",
            },
        ],
    },
    {
        "gene_name": "NAMPT",
        "cytoband": "7q22.3",
        "chromosome": "7",
        "start": 105888744,
        "end": 105925429,
        "strand": "-",
        "gene_summary": (
            "NAMPT encodes nicotinamide phosphoribosyltransferase, the rate-limiting enzyme in the major NAD salvage pathway and a key regulator of metabolic stress, inflammation, and cellular resilience."
        ),
        "clinical_context": (
            "The bundled NAMPT database is focused on metabolism, inflammation, and disease-association research. "
            "Most current human variant interpretation at this locus is research grade and cohort dependent."
        ),
        "variant_effect_overview": [
            "Common NAMPT variants are usually interpreted as regulatory or biomarker-linked modifiers rather than as high-penetrance alleles.",
            "Because NAMPT biology spans NAD metabolism, inflammation, and secreted eNAMPT signaling, effect sizes and directions can vary substantially by phenotype.",
        ],
        "condition_research_overview": [
            "NAMPT is studied in metabolism, obesity, inflammation, cardiovascular disease, cancer, and aging.",
            "The gene also links NAD salvage to sirtuin biology, stress responses, and systemic inflammatory signaling.",
        ],
        "methylation_interpretation": (
            "NAMPT methylation is best interpreted as regulatory context around NAD-salvage capacity and inflammatory-metabolic signaling. "
            "It should be integrated with broader pathway state rather than read as a direct biomarker for any one common SNP."
        ),
        "methylation_effects": [
            "Promoter methylation may help contextualize NAMPT transcriptional potential and systemic metabolic stress responses.",
            "Most phenotype interpretation still depends more on pathway activity, circulating NAMPT, and cohort design than on methylation alone.",
        ],
        "methylation_condition_research": [
            "NAD metabolism and healthy-aging biology.",
            "Inflammation and cardiovascular-risk studies.",
            "Cancer and systemic stress-response research.",
        ],
        "evidence": [
            _evidence("NCBI Gene 10135: NAMPT gene summary", "https://www.ncbi.nlm.nih.gov/gene/10135"),
            _evidence("PubMed 25896907: NAMPT polymorphisms in esophageal carcinoma", "https://pubmed.ncbi.nlm.nih.gov/25896907/"),
            _evidence("PubMed 33446046: NAMPT rs61330082 and cardiovascular risk after HCV clearance", "https://pubmed.ncbi.nlm.nih.gov/33446046/"),
        ],
        "variants": [
            {
                "variant": "rs61330082",
                "common_name": "regulatory NAMPT SNP",
                "region_class": "promoter",
                "interpretation_scope": "Research association / metabolic-inflammatory marker",
                "clinical_interpretation": (
                    "rs61330082 is a promoter-facing NAMPT variant studied in cancer, lipid, and cardiovascular-risk cohorts. "
                    "The local database treats it as a research association marker for inflammatory-metabolic signaling rather than as a pathogenic clinical allele."
                ),
                "clinical_significance": "Research-level metabolic and inflammatory association marker.",
                "functional_effects": [
                    "Likely regulatory or linked-locus effects on NAMPT expression and pathway tone.",
                    "Phenotypic interpretation varies by inflammatory and metabolic context.",
                ],
                "associated_conditions": [
                    "Esophageal squamous-cell carcinoma susceptibility",
                    "Cardiovascular-event risk after chronic viral disease treatment",
                    "Lipid and metabolic biomarker studies",
                ],
                "research_context": [
                    "This SNP is one of the better characterized common NAMPT regulatory markers in disease-association literature.",
                ],
                "usual_variant_note": "Most cited common regulatory NAMPT SNP.",
                "methylation_interpretation": (
                    "Treat NAMPT methylation as broader regulatory context and not as a direct surrogate for rs61330082 genotype."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PubMed 25896907: NAMPT polymorphisms in esophageal carcinoma", "https://pubmed.ncbi.nlm.nih.gov/25896907/"),
                    _evidence("PubMed 33446046: NAMPT rs61330082 and cardiovascular risk", "https://pubmed.ncbi.nlm.nih.gov/33446046/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Chen et al., 2015 (PMID 25896907)",
                        "genotypes": "CC, CT, and TT",
                        "phenotype": "Esophageal squamous-cell carcinoma susceptibility",
                        "finding": "The study examined rs61330082 among the main common NAMPT polymorphisms contributing to esophageal-squamous-cell-carcinoma risk modeling.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/25896907/",
                    },
                    {
                        "paper": "Huang et al., 2021 (PMID 33446046)",
                        "genotypes": "TT versus non-TT",
                        "phenotype": "Long-term cardiovascular events after HCV clearance",
                        "finding": "The rs61330082 TT genotype was associated with a higher cumulative incidence of cardiovascular events in the treated cohort, reinforcing its value as a context-dependent metabolic-inflammatory marker.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/33446046/",
                    },
                ],
            }
        ],
        "population_intro": "Broader population patterns curated from NAMPT inflammatory, cardiovascular, cancer, and NAD-metabolism literature.",
        "population_coverage_note": (
            "This NAMPT population database emphasizes cohort-specific metabolic and inflammatory patterns instead of a full embedded allele-frequency panel."
        ),
        "population_sources": [
            _evidence("NCBI Gene 10135: NAMPT gene summary", "https://www.ncbi.nlm.nih.gov/gene/10135"),
            _evidence("PubMed 25896907: NAMPT polymorphism study", "https://pubmed.ncbi.nlm.nih.gov/25896907/"),
            _evidence("PubMed 33446046: NAMPT rs61330082 cardiovascular study", "https://pubmed.ncbi.nlm.nih.gov/33446046/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "rs61330082",
                "location_group": "Metabolic cohorts",
                "summary": "NAMPT common-variant interpretation is highly phenotype specific, spanning lipid traits, inflammatory tone, and long-term cardiovascular outcomes rather than a single stable population pattern.",
            },
            {
                "variant": "NAMPT pathway context",
                "location_group": "Global pattern",
                "summary": "Because NAMPT sits at the center of NAD salvage and stress signaling, cross-population interpretation should emphasize disease context and metabolic state over blanket assumptions about one common allele.",
            },
        ],
    },
    {
        "gene_name": "CDKN2A",
        "cytoband": "9p21.3",
        "chromosome": "9",
        "start": 21967751,
        "end": 21995323,
        "strand": "-",
        "gene_summary": (
            "CDKN2A encodes multiple tumor-suppressor products, including p16INK4A and p14ARF, that regulate the RB and p53 pathways, G1 cell-cycle control, senescence, and oncogenic stress responses."
        ),
        "clinical_context": (
            "The bundled CDKN2A knowledge base is tumor-suppressor focused. Clinically important CDKN2A events are often rare germline variants, deletions, or promoter methylation changes rather than common low-effect polymorphisms."
        ),
        "variant_effect_overview": [
            "For CDKN2A, rare germline and somatic loss-of-function events usually carry more weight than common SNPs.",
            "Promoter methylation and deletion are central recurring mechanisms in cancer biology and senescence research at this locus.",
            "Observed PASS variants in the interval should therefore be treated as locus context unless they match a seeded research SNP or have outside clinical curation.",
        ],
        "condition_research_overview": [
            "CDKN2A is a major tumor-suppressor locus in melanoma, pancreatic cancer predisposition, glioma, and many additional tumors.",
            "The locus is also central to aging and senescence research because p16INK4A expression is a canonical cellular-senescence marker.",
        ],
        "methylation_interpretation": (
            "CDKN2A methylation is often biologically meaningful because promoter hypermethylation is a recurrent mode of tumor-suppressor silencing. "
            "In this workbench it should still be read as regulatory context rather than as a stand-alone diagnosis."
        ),
        "methylation_effects": [
            "Promoter hypermethylation can be consistent with reduced p16INK4A pathway accessibility and is widely studied in cancer cohorts.",
            "At the same time, CDKN2A expression itself is a senescence biomarker, so methylation should be integrated with tissue and disease context.",
        ],
        "methylation_condition_research": [
            "Cancer-associated promoter methylation studies.",
            "Aging and cellular-senescence research.",
            "Tumor-suppressor silencing and biomarker development.",
        ],
        "evidence": [
            _evidence("NCBI Gene 1029: CDKN2A gene summary", "https://www.ncbi.nlm.nih.gov/gene/1029"),
            _evidence("PubMed 26557774: CDKN2A 3'UTR sequence variants in melanoma", "https://pubmed.ncbi.nlm.nih.gov/26557774/"),
            _evidence("PubMed 31270053: CDKN2A methylation and aging", "https://pubmed.ncbi.nlm.nih.gov/31270053/"),
        ],
        "variants": [
            {
                "variant": "rs11515",
                "common_name": "CDKN2A 3'UTR regulatory SNP",
                "region_class": "gene_body",
                "interpretation_scope": "Research association / regulatory marker",
                "clinical_interpretation": (
                    "rs11515 is a 3'UTR CDKN2A polymorphism curated here as a research-level regulatory marker. "
                    "It has been discussed in melanoma and cancer-risk literature but is not treated as a pathogenic stand-alone allele in this local database."
                ),
                "clinical_significance": "Research-level regulatory association marker.",
                "functional_effects": [
                    "Likely reflects post-transcriptional regulatory context at the CDKN2A locus.",
                    "Any disease signal is modest relative to rare truncating, deletion, or methylation-mediated CDKN2A loss.",
                ],
                "associated_conditions": [
                    "Melanoma and cancer-association studies",
                    "Tumor-suppressor regulatory biology",
                ],
                "research_context": [
                    "Common CDKN2A SNP work should be interpreted cautiously because the major clinical burden at this locus is usually rare loss of function or promoter silencing.",
                ],
                "usual_variant_note": "Common research-facing CDKN2A 3'UTR SNP.",
                "methylation_interpretation": (
                    "CDKN2A methylation and rs11515 should be read as different layers of regulatory context, not as interchangeable biomarkers."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PubMed 26557774: CDKN2A 3'UTR sequence variants in melanoma", "https://pubmed.ncbi.nlm.nih.gov/26557774/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Przybyla et al., 2015 (PMID 26557774)",
                        "genotypes": "Common 3'UTR variant carrier states",
                        "phenotype": "Melanoma-focused sequence-variant analysis",
                        "finding": "The study examined CDKN2A 3'UTR variation including rs11515 as part of melanoma-focused regulatory sequence analysis, reinforcing its use as a research-level rather than high-penetrance marker.",
                        "url": "https://pubmed.ncbi.nlm.nih.gov/26557774/",
                    }
                ],
            }
        ],
        "population_intro": "Broader population patterns curated from CDKN2A tumor-suppressor, methylation, and senescence literature.",
        "population_coverage_note": (
            "The bundled CDKN2A population database emphasizes rare-disease and epigenetic patterns because common SNP frequency panels capture only a small part of CDKN2A biology."
        ),
        "population_sources": [
            _evidence("NCBI Gene 1029: CDKN2A gene summary", "https://www.ncbi.nlm.nih.gov/gene/1029"),
            _evidence("PubMed 31270053: CDKN2A methylation and aging", "https://pubmed.ncbi.nlm.nih.gov/31270053/"),
            _evidence("PubMed 26557774: CDKN2A 3'UTR variants in melanoma", "https://pubmed.ncbi.nlm.nih.gov/26557774/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "Rare CDKN2A loss-of-function variants",
                "location_group": "Rare disease families",
                "summary": "Across populations, the most clinically important CDKN2A events are rare germline or somatic loss-of-function alterations rather than common low-effect polymorphisms.",
            },
            {
                "variant": "CDKN2A promoter methylation",
                "location_group": "Cancer cohorts",
                "summary": "Promoter methylation is a recurring mode of CDKN2A silencing across tumor cohorts and is often more informative biologically than common inherited markers at this locus.",
            },
            {
                "variant": "CDKN2A methylation",
                "location_group": "Longevity cohorts",
                "summary": "In healthy-aging studies, CDKN2A methylation and expression are interpreted through senescence biology, which is distinct from the rare-variant cancer-predisposition literature.",
            },
        ],
    },
    {
        "gene_name": "TERT",
        "cytoband": "5p15.33",
        "chromosome": "5",
        "start": 1253282,
        "end": 1295183,
        "strand": "-",
        "gene_summary": (
            "TERT encodes telomerase reverse transcriptase, the catalytic core of telomerase that maintains telomere repeats and influences replicative lifespan, genome stability, stem-cell function, and oncogenesis."
        ),
        "clinical_context": (
            "The local TERT knowledge base spans telomere biology, cancer susceptibility, and promoter regulation. "
            "Common TERT variants are interpreted as research and risk-modifier signals, while the strongest clinical events often involve promoter mutations or rare high-impact pathogenic alleles."
        ),
        "variant_effect_overview": [
            "Common TERT polymorphisms are often interpreted through telomere-length, cancer-risk, or promoter-activity frameworks rather than as stand-alone diagnostic variants.",
            "TERT biology is highly tissue and phenotype specific, so common-variant effects should be read in the context of telomere maintenance and cohort design.",
        ],
        "condition_research_overview": [
            "TERT is central to telomere biology, cellular senescence, stem-cell maintenance, and oncogenesis.",
            "Common germline SNPs and somatic promoter mutations are both heavily studied in cancer biology, but they represent different interpretive layers.",
        ],
        "methylation_interpretation": (
            "TERT methylation is best interpreted as promoter and chromatin context around telomerase regulation. "
            "It can inform whether the locus appears permissive or repressed, but it is not a direct substitute for telomere-length or promoter-mutation testing."
        ),
        "methylation_effects": [
            "Promoter-proximal methylation can contribute to TERT expression modeling, especially in cancer biology.",
            "Interpret methylation alongside telomere biology, promoter mutation status, and tissue context.",
        ],
        "methylation_condition_research": [
            "Cancer and telomerase-reactivation studies.",
            "Aging and cellular-senescence research.",
            "Stem-cell and tissue-renewal biology.",
        ],
        "evidence": [
            _evidence("NCBI Gene 7015: TERT gene summary", "https://www.ncbi.nlm.nih.gov/gene/7015"),
            _evidence("PMCID PMC3329415: rs2736100 and cancer risk", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3329415/"),
            _evidence("PubMed 22994782: rs2736098 meta-analysis", "https://pubmed.ncbi.nlm.nih.gov/22994782/"),
        ],
        "variants": [
            {
                "variant": "rs2736100",
                "common_name": "intronic TERT cancer- and telomere-associated SNP",
                "region_class": "gene_body",
                "interpretation_scope": "Research association / telomere biology marker",
                "clinical_interpretation": (
                    "rs2736100 is one of the best-known common TERT SNPs and is curated here as a telomere-biology and cancer-susceptibility marker. "
                    "It is not treated as a pathogenic diagnostic allele in the local database."
                ),
                "clinical_significance": "Research-level telomere and cancer association marker.",
                "functional_effects": [
                    "Commonly interpreted through telomere-maintenance and cancer-risk frameworks.",
                    "Likely reflects regulatory or linked-locus effects on TERT biology rather than a direct truncating lesion.",
                ],
                "associated_conditions": [
                    "Lung and other cancer susceptibility studies",
                    "Telomere-length and telomerase-biology research",
                ],
                "research_context": [
                    "This is the canonical common inherited TERT association SNP in many cancer and telomere studies.",
                ],
                "usual_variant_note": "Most cited common TERT intronic association SNP.",
                "methylation_interpretation": (
                    "TERT methylation adds promoter and chromatin context but is not a direct biomarker of rs2736100 genotype."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("PMCID PMC3329415: rs2736100 meta-analysis", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3329415/"),
                ],
                "literature_findings": [
                    {
                        "paper": "Zhu et al., 2012 (PMCID PMC3329415)",
                        "genotypes": "Common rs2736100 genotypes",
                        "phenotype": "Cancer risk across pooled case-control studies",
                        "finding": "The pooled analysis supported rs2736100 as a recurrent TERT cancer-risk marker, reinforcing its role as a common inherited telomere-biology signal.",
                        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3329415/",
                    }
                ],
            },
            {
                "variant": "rs2853669",
                "common_name": "promoter regulatory TERT SNP",
                "region_class": "promoter",
                "interpretation_scope": "Research association / promoter modifier",
                "clinical_interpretation": (
                    "rs2853669 is a promoter-region TERT polymorphism frequently discussed as a modifier of promoter activity and cancer outcomes. "
                    "The local database treats it as a research-level promoter modifier rather than as a pathogenic allele."
                ),
                "clinical_significance": "Research-level promoter modifier.",
                "functional_effects": [
                    "Promoter-facing regulatory variant often modeled as a transcription-factor binding modifier.",
                    "Frequently interpreted together with broader TERT promoter biology in oncology studies.",
                ],
                "associated_conditions": [
                    "Cancer-risk and prognosis studies",
                    "Telomerase-promoter regulatory biology",
                ],
                "research_context": [
                    "This SNP is commonly discussed as a modifier of TERT promoter effects, especially in cancer cohorts.",
                ],
                "usual_variant_note": "Common TERT promoter modifier SNP.",
                "methylation_interpretation": (
                    "Promoter methylation and rs2853669 should be interpreted as complementary layers of TERT promoter context."
                ),
                "is_assayable_in_snp_vcf": True,
                "evidence": [
                    _evidence("NCBI Gene 7015: TERT gene summary", "https://www.ncbi.nlm.nih.gov/gene/7015"),
                ],
                "literature_findings": [],
            },
        ],
        "population_intro": "Broader population patterns curated from TERT telomere-biology and cancer-susceptibility literature.",
        "population_coverage_note": (
            "The bundled TERT population database is literature oriented and focuses on cohort patterns in telomere and cancer studies rather than a complete allele-frequency panel."
        ),
        "population_sources": [
            _evidence("NCBI Gene 7015: TERT gene summary", "https://www.ncbi.nlm.nih.gov/gene/7015"),
            _evidence("PMCID PMC3329415: rs2736100 meta-analysis", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3329415/"),
            _evidence("PubMed 22994782: rs2736098 meta-analysis", "https://pubmed.ncbi.nlm.nih.gov/22994782/"),
        ],
        "gene_population_patterns": [
            {
                "variant": "rs2736100",
                "location_group": "Cancer cohorts",
                "summary": "Common TERT inherited variation is most densely studied in cancer cohorts, where rs2736100 is repeatedly used as a telomere-biology susceptibility marker.",
            },
            {
                "variant": "TERT promoter regulation",
                "location_group": "Global pattern",
                "summary": "At the TERT locus, inherited polymorphisms, promoter methylation, and somatic promoter mutations all matter, so population interpretation should distinguish among these different regulatory layers.",
            },
        ],
    },
]


def _build_promoter_region(meta: dict[str, Any]) -> dict[str, int | str]:
    chrom = str(meta["chromosome"])
    start = int(meta["start"])
    end = int(meta["end"])
    if meta["strand"] == "-":
        promoter_start = end + 1
        promoter_end = end + 1000
        definition = (
            "A practical 1 kb upstream window relative to the reverse-strand transcription start, "
            "used to flag promoter-adjacent coverage before the canonical gene interval."
        )
    else:
        promoter_start = max(1, start - 1000)
        promoter_end = start - 1
        definition = (
            "A practical 1 kb upstream window used by this app to flag promoter-adjacent coverage "
            "before the canonical gene interval begins."
        )

    return {
        "label": "Operational promoter review window",
        "start": promoter_start,
        "end": promoter_end,
        "definition": definition,
    }


def _select_relevant_probe_ids(subset_df: pd.DataFrame, meta: dict[str, Any]) -> list[str]:
    tss_coordinate = int(meta["end"] if meta["strand"] == "-" else meta["start"])
    prioritized = subset_df.copy()
    prioritized["distance_to_tss"] = (prioritized["MAPINFO"] - tss_coordinate).abs()

    if "UCSC_RefGene_Group" in prioritized.columns:
        group_text = prioritized["UCSC_RefGene_Group"].fillna("").astype(str)
        tss_mask = group_text.str.contains("TSS|5'UTR|1stExon", case=False, regex=True)
        if tss_mask.any():
            prioritized = prioritized[tss_mask].copy()
            prioritized["distance_to_tss"] = (prioritized["MAPINFO"] - tss_coordinate).abs()

    return (
        prioritized.sort_values(["distance_to_tss", "MAPINFO"])["IlmnID"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .head(10)
        .tolist()
    )


def _build_hotspot_region(subset_df: pd.DataFrame, meta: dict[str, Any], promoter_region: dict[str, Any]) -> dict[str, Any]:
    candidate = subset_df.copy()
    if "UCSC_RefGene_Group" in candidate.columns:
        group_text = candidate["UCSC_RefGene_Group"].fillna("").astype(str)
        tss_mask = group_text.str.contains("TSS|5'UTR|1stExon", case=False, regex=True)
        if tss_mask.any():
            candidate = candidate[tss_mask].copy()

    if "Relation_to_UCSC_CpG_Island" in candidate.columns:
        island_text = candidate["Relation_to_UCSC_CpG_Island"].fillna("").astype(str)
        island_mask = island_text.str.contains("Island|Shore", case=False, regex=True)
        if island_mask.any():
            candidate = candidate[island_mask].copy()

    if candidate.empty:
        hotspot_start = int(promoter_region["start"])
        hotspot_end = int(promoter_region["end"])
        definition = (
            "Promoter-proximal review span derived from the operational promoter window because the bundled EPIC subset does not mark a smaller CpG-focused hotspot."
        )
    else:
        hotspot_start = int(candidate["MAPINFO"].min())
        hotspot_end = int(candidate["MAPINFO"].max())
        definition = (
            "Promoter-proximal CpG-focused span represented by the bundled EPIC subset and used here as the main methylation-review hotspot."
        )

    return {
        "label": "Promoter-associated CpG review window",
        "start": hotspot_start,
        "end": hotspot_end,
        "definition": definition,
    }


def _build_interpretation_database(meta: dict[str, Any], subset_df: pd.DataFrame) -> dict[str, Any]:
    promoter_region = _build_promoter_region(meta)
    hotspot_region = _build_hotspot_region(subset_df, meta, promoter_region)
    probe_ids = _select_relevant_probe_ids(subset_df, meta)

    combined_start = min(
        int(meta["start"]),
        int(meta["end"]),
        int(promoter_region["start"]),
        int(promoter_region["end"]),
    )
    combined_end = max(
        int(meta["start"]),
        int(meta["end"]),
        int(promoter_region["start"]),
        int(promoter_region["end"]),
    )
    recommended_region = f"{meta['chromosome']}:{combined_start}-{combined_end}"
    gene_context = {
        "gene_name": meta["gene_name"],
        "assembly": ASSEMBLY,
        "cytoband": meta["cytoband"],
        "chromosome": meta["chromosome"],
        "gene_region": {
            "label": f"{meta['gene_name']} transcribed interval",
            "start": int(meta["start"]),
            "end": int(meta["end"]),
            "definition": f"Canonical {meta['gene_name']} genomic interval on {ASSEMBLY} from NCBI Gene.",
        },
        "promoter_review_region": promoter_region,
        "promoter_hotspot_region": hotspot_region,
        "recommended_promoter_plus_gene_region": recommended_region,
        "gene_summary": meta["gene_summary"],
        "clinical_context": meta["clinical_context"],
        "variant_effect_overview": meta["variant_effect_overview"],
        "condition_research_overview": meta["condition_research_overview"],
        "relevant_methylation_probe_ids": probe_ids,
        "methylation_interpretation": meta["methylation_interpretation"],
        "methylation_effects": meta["methylation_effects"],
        "methylation_condition_research": meta["methylation_condition_research"],
        "evidence": meta["evidence"],
    }

    variant_records: list[dict[str, Any]] = []
    for raw_variant in meta["variants"]:
        record = deepcopy(raw_variant)
        record["gene_name"] = meta["gene_name"]
        record["chromosome"] = meta["chromosome"]
        record.setdefault("display_name", record["variant"])
        record.setdefault("lookup_keys", [record["variant"], f"{meta['gene_name'].lower()}:{record['variant'].lower()}"])
        record["relevant_methylation_probe_ids"] = probe_ids
        variant_records.append(record)

    return {
        "database_name": f"NophiGene {meta['gene_name']} Interpretation Database",
        "version": VERSION,
        "gene_context": gene_context,
        "variant_records": variant_records,
    }


def _build_population_database(meta: dict[str, Any]) -> dict[str, Any]:
    return {
        "database_name": f"NophiGene {meta['gene_name']} Population Database",
        "version": VERSION,
        "assembly": ASSEMBLY,
        "coverage_note": meta["population_coverage_note"],
        "gene_population_patterns_intro": meta["population_intro"],
        "population_categories": COMMON_POPULATION_CATEGORIES,
        "sources": meta["population_sources"],
        "variant_population_records": [],
        "gene_population_patterns": meta["gene_population_patterns"],
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Required manifest not found: {MANIFEST_PATH}")

    GENE_DATA_DIR.mkdir(parents=True, exist_ok=True)

    generated: list[str] = []
    for meta in GENE_DEFINITIONS:
        region = f"{meta['chromosome']}:{meta['start']}-{meta['end']}"
        selection = save_filtered_manifest(
            gene_name=meta["gene_name"],
            manifest_path=str(MANIFEST_PATH),
            region=region,
            genome_build="hg19",
            output_dir=GENE_DATA_DIR,
        )
        subset_path = Path(selection["output_path"])
        subset_df = pd.read_csv(subset_path)

        interpretation_path = GENE_DATA_DIR / f"{meta['gene_name'].lower()}_interpretation_db.json"
        population_path = GENE_DATA_DIR / f"{meta['gene_name'].lower()}_population_db.json"

        _write_json(interpretation_path, _build_interpretation_database(meta, subset_df))
        _write_json(population_path, _build_population_database(meta))

        generated.extend(
            [
                subset_path.name,
                interpretation_path.name,
                population_path.name,
            ]
        )

    print("Generated knowledge-base artifacts:")
    for name in generated:
        print(f" - {name}")


if __name__ == "__main__":
    main()
