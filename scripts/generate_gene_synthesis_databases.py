#!/usr/bin/env python3
"""Generate gene-level predictive synthesis databases from the bundled interpretation JSON files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENE_DATA_DIR = PROJECT_ROOT / "src" / "gene_data"
VERSION = "2026-04-18"
TARGET_GENES = [
    "HERC2",
    "DRD4",
    "IGF1R",
    "GLP1R",
    "FOXO3",
    "MTOR",
    "RPS6",
    "SIK3",
    "FLCN",
    "SIRT6",
    "PRKAA1",
    "NAMPT",
    "CDKN2A",
    "TERT",
]

METHYLATION_SOURCES = [
    {
        "key": "whitelist",
        "label": "Whitelist mean beta",
        "description": "Uses the curated whitelist probes stored in the interpretation database for the current gene.",
    },
    {
        "key": "gene_name_related",
        "label": "Gene-name-related mean beta",
        "description": "Uses rows whose gene annotation explicitly names the current gene.",
    },
    {
        "key": "all_numeric",
        "label": "All numeric-row mean beta",
        "description": "Uses every numeric beta value that survived preprocessing for the current sample.",
    },
]

BAND_CONTEXT = {
    "high": "This pattern leans toward stronger regulatory restraint or reduced local accessibility around the locus.",
    "medium": "This pattern suggests an intermediate or mixed regulatory state rather than a strongly polarized epigenetic signal.",
    "low": "This pattern leans toward a more permissive or less methylated local regulatory state.",
}

GENE_CONCRETE_VARIANT_PREDICTIONS = {
    "HERC2": (
        "The variant observed in this sample suggests a HERC2/OCA2 pigmentation prediction rather than a HERC2 protein-disease prediction. "
        "If the matched marker is rs12913832, the forward-strand G state used by this workbench supports a lighter or blue-eye tendency because it weakens HERC2 enhancer looping to the OCA2 promoter, lowers OCA2 expression, and reduces iris melanin; the A state supports a darker or brown-eye tendency. "
        "Linked HERC2 markers should be read as light-versus-dark eye-colour haplotype evidence, with ancestry and additional OCA2 variants able to modify the visible result."
    ),
    "DRD4": (
        "The variant observed in this sample suggests a dopamine-D4 regulatory or repeat-background thesis: the individual may carry a research-grade modifier of attentional, novelty-seeking, addiction, social-affective, or neuropsychiatric trait studies because DRD4 variation can alter receptor regulation or dopaminergic signaling context. "
        "This is a behavioral-association prediction, not a diagnosis or deterministic personality call."
    ),
    "IGF1R": (
        "The variant observed in this sample suggests an IGF1 receptor signaling-set-point thesis: the individual may carry a low-penetrance modifier of growth, endocrine, aging, pregnancy-history, cardiovascular, or cancer-cohort biology because IGF1R tunes PI3K-AKT and MAPK survival signaling. "
        "The prediction is pathway-contextual rather than a direct disease call."
    ),
    "GLP1R": (
        "The variant observed in this sample suggests an incretin-receptor response thesis: the individual may carry a research-grade modifier of GLP-1 receptor signaling, glucose-dependent insulin secretion, appetite or gastric-emptying biology, BMI/metabolic traits, or response patterns to GLP-1 receptor agonists and DPP-4 inhibitors. "
        "This is pharmacogenetic and metabolic context, not a medication-selection rule or deterministic diabetes/obesity prediction."
    ),
    "FOXO3": (
        "The variant observed in this sample suggests a FOXO3 healthy-aging and stress-resilience thesis: the individual may carry a longevity-associated regulatory background, especially if the matched marker is rs2802292-like, because FOXO3 controls stress response, autophagy, inflammation, and metabolic adaptation. "
        "This predicts a cohort-level resilience signal, not guaranteed lifespan."
    ),
    "MTOR": (
        "The variant observed in this sample suggests an mTOR pathway-threshold thesis: the individual may carry a regulatory modifier of growth signaling, autophagy restraint, nutrient sensing, or cancer/outcome cohorts because MTOR integrates anabolic and survival signals. "
        "This is best interpreted as pathway tuning, not as a stand-alone diagnosis."
    ),
    "RPS6": (
        "The variant observed in this sample suggests an mTORC1-S6K output-context thesis: any visible RPS6 locus signal should be treated as translational-control context because RPS6 is primarily interpreted through ribosomal protein S6 abundance and phosphorylation rather than through a mature common-variant phenotype catalog."
    ),
    "SIK3": (
        "The variant observed in this sample suggests a sleep-need and AMPK-related kinase thesis: if the matched event is the N783Y-like SIK3 variant, the concrete phenotype prediction is natural short sleep or altered sleep need because SIK3 participates in sleep homeostasis and neuronal signaling. "
        "Other SIK3 variants should be read as emerging sleep, metabolic, and transcriptional signaling context."
    ),
    "FLCN": (
        "The variant observed in this sample suggests a folliculin tumor-suppressor thesis: if the matched event is a pathogenic FLCN frameshift such as c.1285dupC, the concrete clinical prediction is Birt-Hogg-Dube syndrome risk context, including renal tumors, lung cysts, spontaneous pneumothorax, and fibrofolliculoma biology. "
        "Common or unclassified signals should not be upgraded to that clinical prediction without variant-level pathogenic evidence."
    ),
    "SIRT6": (
        "The variant observed in this sample suggests a SIRT6 genome-maintenance and longevity thesis: the individual may carry a modifier of DNA repair, telomere maintenance, inflammation, glucose or lipid metabolism, or centenarian-enriched biology because SIRT6 regulates chromatin and stress resilience. "
        "This is an aging-biology prediction, not a guarantee of exceptional longevity."
    ),
    "PRKAA1": (
        "The variant observed in this sample suggests an AMPK energy-sensing thesis: the individual may carry a low-penetrance modifier of gastric-cancer association studies, metabolic stress adaptation, mitochondrial response, or mTOR restraint because PRKAA1 encodes the AMPK alpha-1 catalytic subunit. "
        "This is a pathway-risk context rather than a deterministic disease prediction."
    ),
    "NAMPT": (
        "The variant observed in this sample suggests an NAD-salvage and inflammatory-metabolic thesis: the individual may carry a regulatory modifier of NAMPT expression, lipid or cardiovascular cohorts, inflammation, cancer metabolism, or cellular resilience because NAMPT controls a major route of NAD regeneration. "
        "This prediction is biomarker and pathway oriented."
    ),
    "CDKN2A": (
        "The variant observed in this sample suggests a cell-cycle checkpoint and tumor-suppressor thesis: the individual may carry a modifier of p16INK4A/p14ARF, RB-p53, senescence, melanoma, pancreatic cancer, glioma, or broader cancer-risk biology because CDKN2A controls G1 arrest and oncogenic-stress response. "
        "Only known loss-of-function or pathogenic CDKN2A variants should be treated as high-risk clinical findings."
    ),
    "TERT": (
        "The variant observed in this sample suggests a telomere-maintenance thesis: the individual may carry a modifier of telomere length, replicative lifespan, stem-cell maintenance, promoter activity, or cancer susceptibility because TERT encodes the catalytic telomerase subunit. "
        "This predicts telomere-biology context, not a direct cancer diagnosis."
    ),
}

VARIANT_CONCRETE_PREDICTION_OVERRIDES = {
    "HERC2": {
        "rs12913832": (
            "Observed rs12913832 is the strongest concrete prediction in this bundle: in the forward-strand A/G representation used here, G supports a lighter or blue-eye tendency and A supports a darker or brown-eye tendency because this enhancer variant changes HERC2-OCA2 looping and OCA2-driven iris melanin production. "
            "A single observed allele without full genotype should be treated as directional evidence rather than a final eye-colour call."
        ),
        "rs1129038": (
            "Observed rs1129038 suggests the sample carries the linked HERC2/OCA2 pigmentation haplotype used in light-versus-dark eye-colour prediction models. "
            "It refines the brown/blue eye-colour thesis mainly through linkage with rs12913832 rather than through a stronger independent mechanism."
        ),
        "rs7170852": (
            "Observed rs7170852 suggests secondary HERC2 eye-colour prediction context: it can support a pigmentation-haplotype thesis, but it should not override rs12913832 because its independent effect is weaker and cohort dependent."
        ),
        "rs916977": (
            "Observed rs916977 suggests a HERC2 pigmentation-background thesis with European iris-colour cline evidence and pigmentation-linked skin-cancer context. "
            "It is a secondary light-versus-dark eye-colour marker rather than the main causal enhancer call."
        ),
        "rs11636232": (
            "Observed rs11636232 suggests a secondary HERC2/OCA2 eye-colour prediction signal that can sharpen multilocus pigmentation models, especially when read with rs12913832 and rs1129038."
        ),
    },
    "SIK3": {
        "SIK3 p.Asn783Tyr": (
            "Observed SIK3 p.Asn783Tyr suggests a natural-short-sleep thesis: the individual may need less sleep than average because this rare functional SIK3 variant has been reported in human short-sleep biology."
        ),
    },
    "FLCN": {
        "FLCN c.1285dupC": (
            "Observed FLCN c.1285dupC suggests a Birt-Hogg-Dube syndrome thesis because this recurrent frameshift disrupts folliculin. "
            "The concrete risk context is renal tumor surveillance, lung cyst or pneumothorax susceptibility, and skin fibrofolliculoma biology."
        ),
    },
    "SIRT6": {
        "centSIRT6 (N308K/A313S)": (
            "Observed centSIRT6-like N308K/A313S suggests a centenarian-enriched SIRT6 thesis: the sample may carry a rare functional background linked to stronger genome-maintenance activity in experimental follow-up."
        ),
    },
    "GLP1R": {
        "rs6923761": (
            "Observed GLP1R rs6923761 / Gly168Ser suggests a pharmacogenetic response-context thesis: the sample carries a GLP-1 receptor missense marker studied in gliptin response and GLP-1 receptor agonist response cohorts. "
            "The concrete prediction is possible heterogeneity in incretin-drug glycemic response, not a universal prediction of semaglutide, liraglutide, dulaglutide, or gliptin benefit."
        ),
        "rs10305420": (
            "Observed GLP1R rs10305420 suggests a liraglutide and metabolic-trait response-context thesis. "
            "In the bundled evidence, T-allele homozygosity was associated with optimal glycemic response to liraglutide in one Iranian type 2 diabetes cohort, but this app cannot infer homozygosity from a single REF -> ALT row."
        ),
        "rs3765467": (
            "Observed GLP1R rs3765467 / p.R131Q suggests a metabolic disease and GLP-1 receptor agonist response-context thesis. "
            "It has been discussed in early-onset type 2 diabetes, dyslipidemia, gestational diabetes, and case-level dulaglutide response literature, so interpretation should stay cohort specific."
        ),
        "rs1042044": (
            "Observed GLP1R rs1042044 suggests a BMI and glucose-insulin trait interaction thesis, especially in developmental or gestational-diabetes exposure cohorts. "
            "This is a metabolic-trait context signal rather than a direct obesity or diabetes prediction."
        ),
    },
}

VARIANT_ALLELE_CHANGE_PREDICTION_OVERRIDES = {
    "HERC2": {
        "rs12913832": [
            {
                "change": "A>G",
                "alt_allele": "G",
                "prediction": (
                    "This sample's actual rs12913832 DNA change is A -> G, so the observed alternate allele is G. "
                    "In the forward-strand hg19 A/G representation used by this workbench, G is the light-eye-associated state: it reduces HERC2 enhancer support for OCA2 expression, lowers iris melanin biology, and therefore supports a blue or lighter-eye tendency. "
                    "Because the current preview uses REF -> ALT and does not infer genotype dosage, treat this as directional allele evidence rather than a complete eye-colour call."
                ),
                "basis": "Observed ALT G at rs12913832, the HERC2/OCA2 enhancer state associated with reduced OCA2 expression and lighter iris pigmentation.",
            },
            {
                "change": "G>A",
                "alt_allele": "A",
                "prediction": (
                    "This sample's actual rs12913832 DNA change is G -> A, so the observed alternate allele is A. "
                    "In the forward-strand hg19 A/G representation used by this workbench, A is the darker-eye-associated state: it is compatible with stronger enhancer-promoter contact, higher OCA2 expression, and greater iris melanin biology, so it supports a brown or darker-eye tendency. "
                    "Because the current preview uses REF -> ALT and does not infer genotype dosage, treat this as directional allele evidence rather than a complete eye-colour call."
                ),
                "basis": "Observed ALT A at rs12913832, the HERC2/OCA2 enhancer state associated with stronger OCA2 expression and darker iris pigmentation.",
            },
        ],
    },
    "DRD4": {
        "rs3758653": [
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at rs3758653. "
                    "The bundled DRD4 literature links T-containing genotypes to higher DRD4 promoter methylation in mQTL work and, in one heroin-use cohort, longer mean latency from first exposure to addiction than CC. "
                    "The resulting prediction is a regulatory/methylation-QTL context signal, not a diagnosis or a deterministic behavioral forecast."
                ),
                "basis": "Observed ALT T at rs3758653, interpreted through DRD4 promoter methylation-QTL and addiction-timing literature.",
            },
            {
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele at rs3758653. "
                    "In the bundled DRD4 studies, C is the contrast allele to the T-associated higher-methylation pattern; CC was also the shortest-latency group in one heroin-use cohort. "
                    "This should be read as promoter-regulatory context rather than as a stand-alone addiction or cognition prediction."
                ),
                "basis": "Observed ALT C at rs3758653, interpreted as the contrast state to the T-linked methylation-QTL pattern.",
            },
        ],
        "rs1800955": [
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at rs1800955 (-521C>T). "
                    "Early reporter assays in the bundled evidence found lower DRD4 promoter transcriptional efficiency for the T allele relative to C, while later brain-expression and association data were mixed. "
                    "The concrete prediction is therefore a possible lower-promoter-efficiency regulatory thesis, not a clinical diagnosis."
                ),
                "basis": "Observed ALT T at rs1800955, a promoter state reported to reduce transcriptional efficiency in early functional assays.",
            },
            {
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele at rs1800955 (-521C>T). "
                    "C is the contrast state to the T allele that reduced transcriptional efficiency in early reporter assays; meta-analysis evidence in the local bundle treats CC as a small-effect susceptibility context in some schizophrenia datasets. "
                    "The prediction remains a modest promoter-regulatory research signal."
                ),
                "basis": "Observed ALT C at rs1800955, interpreted through promoter-activity and small-effect association literature.",
            },
        ],
    },
    "IGF1R": {
        "rs2229765": [
            {
                "alt_allele": "A",
                "prediction": (
                    "This sample row reports A as the observed alternate allele at IGF1R rs2229765. "
                    "In the bundled studies, A-containing states were linked to male longevity and lower free IGF-1 in older AA men, and GA/AA was associated with reduced spontaneous-preterm-birth risk in one cohort. "
                    "The resulting thesis is an IGF1R signaling-set-point modifier, not a deterministic endocrine or pregnancy prediction."
                ),
                "basis": "Observed ALT A at rs2229765, the allele highlighted in local longevity and pregnancy-outcome association evidence.",
            },
            {
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at IGF1R rs2229765. "
                    "G is the contrast state to the A-linked findings in the bundled longevity and pregnancy studies, so the sample should not be described as carrying the A-direction signal from those papers on this row. "
                    "Keep the prediction at IGF1R pathway-context level unless genotype dosage is available."
                ),
                "basis": "Observed ALT G at rs2229765, interpreted as the contrast state to the A-linked cohort findings.",
            },
        ],
        "rs2016347": [
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at IGF1R rs2016347. "
                    "The bundled reproductive-history studies associate T-carrier states with reduced terminal duct lobular unit counts and lower hormone-receptor-positive breast-cancer risk in specific hypertensive-pregnancy or preeclampsia contexts. "
                    "The prediction is a context-dependent IGF-axis regulatory thesis, not an unconditional cancer-risk call."
                ),
                "basis": "Observed ALT T at rs2016347, the allele highlighted in local reproductive-history and breast-involution evidence.",
            },
            {
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at IGF1R rs2016347. "
                    "G is the contrast state to the T-carrier direction emphasized in the bundled breast-involution and hypertensive-pregnancy literature. "
                    "Interpret this as IGF1R regulatory background rather than applying the T-carrier protective-context thesis."
                ),
                "basis": "Observed ALT G at rs2016347, interpreted as the contrast state to the T-carrier cohort findings.",
            },
        ],
    },
    "GLP1R": {
        "rs6923761": [
            {
                "alt_allele": "A",
                "prediction": (
                    "This sample row reports A as the observed alternate allele at GLP1R rs6923761 / Gly168Ser. "
                    "A is the minor allele in the local Ensembl record and represents the Ser-direction state commonly discussed in GLP1R Gly168Ser literature. "
                    "The bundled gliptin study found lower HbA1c reduction in Ser/Ser homozygotes than in Gly-allele carriers, so ALT A supports a possible reduced gliptin-response thesis only when genotype dosage and treatment context are compatible."
                ),
                "basis": "Observed ALT A at rs6923761, interpreted through GLP1R Gly168Ser pharmacogenetic evidence.",
            },
            {
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at GLP1R rs6923761. "
                    "G is the ancestral/reference Gly-direction state in the local Ensembl record, so this row should not be labelled as carrying the Ser/Ser reduced-gliptin-response signal from the bundled study. "
                    "Keep the conclusion at GLP1R incretin-response context unless full genotype dosage is available."
                ),
                "basis": "Observed ALT G at rs6923761, interpreted as the Gly-direction contrast state to the Ser-focused pharmacogenetic evidence.",
            },
        ],
        "rs10305420": [
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at GLP1R rs10305420. "
                    "In the bundled liraglutide pharmacogenetic study, T-allele homozygosity was associated with optimal glycemic response in an Iranian type 2 diabetes cohort. "
                    "Because this app does not infer homozygosity, treat ALT T as directional response-context evidence rather than a final liraglutide-responder call."
                ),
                "basis": "Observed ALT T at rs10305420, the allele highlighted in local liraglutide-response evidence.",
            },
            {
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele at GLP1R rs10305420. "
                    "C is the contrast state to the T-homozygous response signal in the bundled liraglutide study, so the sample should not be described as carrying that T-direction thesis from this row alone."
                ),
                "basis": "Observed ALT C at rs10305420, interpreted as the contrast state to the T-homozygous response evidence.",
            },
        ],
        "rs3765467": [
            {
                "alt_allele": "A",
                "prediction": (
                    "This sample row reports A as the observed alternate allele at GLP1R rs3765467 / p.R131Q. "
                    "The bundled literature discusses rs3765467 in early-onset type 2 diabetes and GLP-1 receptor agonist response contexts, including a dulaglutide case report. "
                    "ALT A therefore supports a cohort-specific metabolic and medication-response context thesis, not a deterministic therapy prediction."
                ),
                "basis": "Observed ALT A at rs3765467, interpreted through p.R131Q diabetes and GLP-1 receptor agonist response-context evidence.",
            },
            {
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at GLP1R rs3765467. "
                    "G is the ancestral/reference state in the local Ensembl record, so avoid applying the alternate p.R131Q direction without full genotype and allele-orientation context."
                ),
                "basis": "Observed ALT G at rs3765467, interpreted as the contrast state to the p.R131Q-focused evidence.",
            },
        ],
        "rs1042044": [
            {
                "alt_allele": "A",
                "prediction": (
                    "This sample row reports A as the observed alternate allele at GLP1R rs1042044. "
                    "The bundled EPOCH studies discuss rs1042044 carrier states as modifiers of BMI growth or metabolic traits in gestational-diabetes exposure and pediatric metabolic contexts. "
                    "ALT A supports a developmental metabolic-trait interaction thesis, not a direct obesity or diabetes prediction."
                ),
                "basis": "Observed ALT A at rs1042044, interpreted through GLP-1R pediatric BMI and metabolic-trait interaction evidence.",
            },
            {
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele at GLP1R rs1042044. "
                    "C is the ancestral/reference state in the local Ensembl record; keep the sample interpretation at general GLP1R metabolic context unless genotype dosage supports a carrier-state thesis."
                ),
                "basis": "Observed ALT C at rs1042044, interpreted as the contrast state to carrier-state EPOCH evidence.",
            },
        ],
    },
    "FOXO3": {
        "rs2802292": [
            {
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at FOXO3 rs2802292. "
                    "The bundled longevity literature treats G as the healthy-aging-associated state and mechanistic follow-up links it to stress-responsive FOXO3 enhancer behavior. "
                    "The prediction is a modest resilience and stress-response pathway thesis, not a guaranteed lifespan prediction."
                ),
                "basis": "Observed ALT G at rs2802292, the FOXO3 allele associated with longevity and stress-responsive enhancer behavior.",
            },
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at FOXO3 rs2802292. "
                    "T is the contrast state to the G-linked longevity and stress-responsive enhancer findings in the local bundle, so the sample should not be labelled as carrying the G-direction FOXO3 longevity signal from this row alone. "
                    "Keep the interpretation as general FOXO3 pathway context."
                ),
                "basis": "Observed ALT T at rs2802292, interpreted as the contrast state to the G-linked longevity evidence.",
            },
        ],
    },
    "MTOR": {
        "rs2295080": [
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at MTOR rs2295080. "
                    "The bundled cancer-association evidence discusses TT as a higher-risk genotype in the Chinese studies pooled by the local database, so ALT T supports a promoter-regulatory risk-context thesis when cohort ancestry and genotype dosage are appropriate. "
                    "It should not be interpreted as a stand-alone cancer diagnosis."
                ),
                "basis": "Observed ALT T at rs2295080, interpreted through promoter-regulatory MTOR cancer-association evidence.",
            },
            {
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at MTOR rs2295080. "
                    "G is the contrast state to the TT-focused higher-risk signal in the bundled cancer-association evidence, so the prediction should stay at pathway-regulatory context unless full genotype and cohort context support a stronger statement."
                ),
                "basis": "Observed ALT G at rs2295080, interpreted as the contrast state to the TT-focused association evidence.",
            },
        ],
    },
    "PRKAA1": {
        "rs13361707": [
            {
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele at PRKAA1 rs13361707. "
                    "The bundled gastric-cancer association study reports the C allele as the higher-risk direction in an eastern Chinese population, so ALT C supports a low-penetrance AMPK-pathway risk-context thesis. "
                    "Do not treat this as diagnostic without ancestry, phenotype, and genotype-dose context."
                ),
                "basis": "Observed ALT C at rs13361707, the PRKAA1 allele highlighted in local gastric-cancer susceptibility evidence.",
            },
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at PRKAA1 rs13361707. "
                    "T is the contrast state to the C-linked gastric-cancer risk direction in the bundled evidence, so the prediction should remain an AMPK-pathway context signal rather than a C-allele risk statement."
                ),
                "basis": "Observed ALT T at rs13361707, interpreted as the contrast state to the C-linked gastric-cancer evidence.",
            },
        ],
    },
    "NAMPT": {
        "rs61330082": [
            {
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at NAMPT rs61330082. "
                    "The bundled cardiovascular follow-up evidence highlights TT as the higher-event group after HCV clearance, so ALT T supports a metabolic-inflammatory risk-context thesis when genotype and clinical context are compatible. "
                    "This remains cohort-specific research evidence."
                ),
                "basis": "Observed ALT T at rs61330082, interpreted through NAMPT cardiovascular-event and inflammatory-metabolic association evidence.",
            },
            {
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele at NAMPT rs61330082. "
                    "C is the contrast state to the TT-focused cardiovascular-event signal in the bundled evidence, so the sample should be framed as NAMPT regulatory context rather than as carrying that T-direction thesis from this row alone."
                ),
                "basis": "Observed ALT C at rs61330082, interpreted as the contrast state to the TT-focused NAMPT evidence.",
            },
        ],
    },
}


def _clean_text(value: Any) -> str:
    """Normalize arbitrary values into compact single-line strings."""
    return " ".join(str(value or "").split())


def _first_nonempty(*values: Any) -> str:
    """Return the first cleaned non-empty value from a list of candidates."""
    for value in values:
        cleaned = _clean_text(value)
        if cleaned:
            return cleaned
    return ""


def _dedupe_text_items(values: list[str]) -> list[str]:
    """Return text items without duplicates while preserving order."""
    deduped: list[str] = []
    seen: set[str] = set()
    for value in values:
        cleaned = _clean_text(value)
        if not cleaned:
            continue
        key = cleaned.casefold()
        if key in seen:
            continue
        seen.add(key)
        deduped.append(cleaned)
    return deduped


def _candidate_interpretation_paths(gene_name: str) -> list[Path]:
    """Return likely interpretation database paths for one gene."""
    return [
        GENE_DATA_DIR / f"{gene_name.lower()}_interpretation_db.json",
        GENE_DATA_DIR / f"{gene_name}_interpretation_db.json",
        GENE_DATA_DIR / f"{gene_name.upper()}_interpretation_db.json",
    ]


def _load_interpretation_database(gene_name: str) -> dict[str, Any]:
    """Load the bundled interpretation database for a gene."""
    for path in _candidate_interpretation_paths(gene_name):
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
    raise FileNotFoundError(f"No interpretation database found for {gene_name}")


def _build_seeded_markers(knowledge_base: dict[str, Any]) -> list[str]:
    """Return a short display list of seeded markers for the synthesis UI."""
    markers: list[str] = []
    for record in knowledge_base.get("variant_records", []):
        display_name = _first_nonempty(record.get("display_name"), record.get("variant"))
        common_name = _clean_text(record.get("common_name"))
        if display_name and common_name:
            markers.append(f"{display_name} ({common_name})")
        elif display_name:
            markers.append(display_name)
    return _dedupe_text_items(markers)[:6]


def _collect_research_focus(knowledge_base: dict[str, Any]) -> list[str]:
    """Collect a concise research-focus list from the interpretation bundle."""
    gene_context = knowledge_base.get("gene_context", {})
    focus_items: list[str] = []
    focus_items.extend(gene_context.get("condition_research_overview", []))
    focus_items.extend(gene_context.get("methylation_condition_research", []))
    for record in knowledge_base.get("variant_records", []):
        focus_items.extend(record.get("associated_conditions", []))
    return _dedupe_text_items(focus_items)[:4]


def _concrete_variant_prediction_for_gene(gene_name: str) -> str:
    """Return a concrete gene-level variant prediction text."""
    return GENE_CONCRETE_VARIANT_PREDICTIONS.get(
        gene_name,
        (
            f"The variant observed in this sample suggests a {gene_name} gene-specific research thesis. "
            "Interpret the result through the bundled gene context and do not treat it as a deterministic clinical prediction."
        ),
    )


def _variant_lookup_values(record: dict[str, Any]) -> list[str]:
    """Return stable lookup values for a curated variant record."""
    return [
        _clean_text(record.get("variant")),
        _clean_text(record.get("display_name")),
        _clean_text(record.get("common_name")),
    ]


def _variant_prediction_override(gene_name: str, record: dict[str, Any]) -> str:
    """Return an override text for a specific variant record when one is curated."""
    overrides = VARIANT_CONCRETE_PREDICTION_OVERRIDES.get(gene_name, {})
    for value in _variant_lookup_values(record):
        if value and value in overrides:
            return overrides[value]
    return ""


def _format_change_display(change: str) -> str:
    """Render a compact allele-change key as the UI-style REF -> ALT text."""
    cleaned = _clean_text(change).replace(" ", "")
    if ">" not in cleaned:
        return cleaned
    ref, alt = cleaned.split(">", 1)
    return f"{ref} -> {alt}"


def _build_sample_change_template(gene_name: str, display_name: str) -> str:
    """Return the generic sample-change template for one variant rule."""
    return _clean_text(
        f"This sample matched {display_name} with the actual DNA change {{change}}. "
        f"Anchor the {gene_name} prediction to that observed REF -> ALT result instead of treating the variant as only a theoretical marker. "
        "The current preview does not infer zygosity or allele dose, so dosage-specific conclusions require genotype fields."
    )


def _build_allele_change_rules(gene_name: str, record: dict[str, Any]) -> list[dict[str, Any]]:
    """Return curated allele or change-specific prediction rules when available."""
    overrides = VARIANT_ALLELE_CHANGE_PREDICTION_OVERRIDES.get(gene_name, {})
    selected_rules: list[dict[str, Any]] = []
    for value in _variant_lookup_values(record):
        if value and value in overrides:
            selected_rules = overrides[value]
            break

    allele_change_rules: list[dict[str, Any]] = []
    for override in selected_rules:
        rule: dict[str, Any] = {
            "prediction": _clean_text(override.get("prediction")),
            "basis": _clean_text(override.get("basis")),
            "evidence": record.get("evidence", []),
        }
        change = _clean_text(override.get("change"))
        alt_allele = _clean_text(override.get("alt_allele")).upper()
        if change:
            rule["change"] = _format_change_display(change)
        if alt_allele:
            rule["alt_allele"] = alt_allele
        allele_change_rules.append(rule)
    return allele_change_rules


def _build_variant_prediction_rules(
    knowledge_base: dict[str, Any],
    *,
    concrete_variant_prediction: str,
) -> list[dict[str, Any]]:
    """Build concrete variant-level prediction rules for matched records."""
    gene_name = _first_nonempty(knowledge_base.get("gene_context", {}).get("gene_name"), "UNKNOWN")
    rules: list[dict[str, Any]] = []
    for record in knowledge_base.get("variant_records", []):
        variant = _first_nonempty(record.get("variant"), record.get("display_name"))
        if not variant:
            continue

        display_name = _first_nonempty(record.get("display_name"), variant)
        common_name = _clean_text(record.get("common_name"))
        record_prediction = _variant_prediction_override(gene_name, record)
        if not record_prediction:
            record_prediction = _clean_text(
                f"Observed {display_name} suggests this sample carries the {gene_name} prediction context described by the gene-level thesis. "
                f"{concrete_variant_prediction} This specific marker is curated as {common_name or 'a gene-specific research marker'}. "
                f"{record.get('clinical_interpretation', '')}"
            )

        rules.append(
            {
                "variant": variant,
                "display_name": display_name,
                "common_name": common_name,
                "lookup_keys": record.get("lookup_keys", []),
                "prediction": _clean_text(record_prediction),
                "sample_change_template": _build_sample_change_template(gene_name, display_name),
                "allele_change_rules": _build_allele_change_rules(gene_name, record),
                "basis": _first_nonempty(
                    *(record.get("functional_effects", []) or []),
                    record.get("clinical_significance"),
                    record.get("interpretation_scope"),
                ),
                "evidence": record.get("evidence", []),
            }
        )
    return rules


def _build_base_case(
    *,
    gene_name: str,
    clinical_context: str,
    variant_focus: str,
    concrete_variant_prediction: str,
    research_focus: list[str],
) -> dict[str, Any]:
    """Build the base variant-only synthesis case."""
    prediction = _clean_text(
        f"{concrete_variant_prediction} "
        f"An observed {gene_name} variant places this sample in the curated {gene_name} research context. "
        f"{clinical_context} {variant_focus}"
    )
    return {
        "case_id": "gene_variant_found",
        "label": "Gene variant found",
        "requires_variant": True,
        "methylation_source": None,
        "methylation_band": None,
        "prediction": prediction,
        "rationale": (
            "This is the base thesis that activates as soon as any promoter or gene-body variant is visible in the current sample."
        ),
        "research_focus": research_focus[:3],
    }


def _build_combined_case(
    *,
    gene_name: str,
    source_key: str,
    source_label: str,
    source_description: str,
    band: str,
    clinical_context: str,
    variant_focus: str,
    methylation_context: str,
    concrete_variant_prediction: str,
    research_focus: list[str],
) -> dict[str, Any]:
    """Build one of the nine variant-plus-methylation synthesis cases."""
    prediction = _clean_text(
        f"When a {gene_name} variant is paired with {band} methylation in the {source_label.lower()}, "
        f"the sample best fits a combined regulatory-context thesis. {BAND_CONTEXT[band]} "
        f"The variant side of the thesis is: {concrete_variant_prediction} "
        f"{methylation_context} {clinical_context} {variant_focus}"
    )
    return {
        "case_id": f"gene_variant_found__{source_key}__{band}",
        "label": f"Gene variant found + {band} {source_label.lower()}",
        "requires_variant": True,
        "methylation_source": source_key,
        "methylation_band": band,
        "prediction": prediction,
        "rationale": source_description,
        "research_focus": research_focus[:3],
    }


def build_synthesis_database(knowledge_base: dict[str, Any]) -> dict[str, Any]:
    """Create the 10-case predictive synthesis matrix for one gene."""
    gene_context = knowledge_base.get("gene_context", {})
    gene_name = _first_nonempty(gene_context.get("gene_name"), "UNKNOWN")
    clinical_context = _first_nonempty(
        gene_context.get("clinical_context"),
        gene_context.get("gene_summary"),
        f"The bundled {gene_name} database is intended for research context.",
    )
    variant_focus = _first_nonempty(
        *(gene_context.get("variant_effect_overview", []) or []),
        "Observed variants should be treated as locus-aware research context unless stronger external evidence exists.",
    )
    methylation_context = _first_nonempty(
        gene_context.get("methylation_interpretation"),
        *(gene_context.get("methylation_effects", []) or []),
        f"{gene_name} methylation is best interpreted as regulatory context rather than as a standalone biomarker.",
    )
    concrete_variant_prediction = _concrete_variant_prediction_for_gene(gene_name)
    variant_prediction_rules = _build_variant_prediction_rules(
        knowledge_base,
        concrete_variant_prediction=concrete_variant_prediction,
    )
    research_focus = _collect_research_focus(knowledge_base)

    cases = [
        _build_base_case(
            gene_name=gene_name,
            clinical_context=clinical_context,
            variant_focus=variant_focus,
            concrete_variant_prediction=concrete_variant_prediction,
            research_focus=research_focus,
        )
    ]

    for source in METHYLATION_SOURCES:
        for band in ("high", "medium", "low"):
            cases.append(
                _build_combined_case(
                    gene_name=gene_name,
                    source_key=source["key"],
                    source_label=source["label"],
                    source_description=source["description"],
                    band=band,
                    clinical_context=clinical_context,
                    variant_focus=variant_focus,
                    methylation_context=methylation_context,
                    concrete_variant_prediction=concrete_variant_prediction,
                    research_focus=research_focus,
                )
            )

    return {
        "database_name": f"NophiGene {gene_name} Predictive Synthesis Database",
        "version": VERSION,
        "gene_name": gene_name,
        "source_interpretation_database": knowledge_base.get(
            "database_name",
            f"NophiGene {gene_name} Interpretation Database",
        ),
        "matching_rule": (
            "One base case matches when a promoter or gene-body variant is visible in the current sample. "
            "Three additional case families match when the whitelist mean beta, the gene-name-related mean beta, "
            "or the all-numeric mean beta resolves to low, medium, or high methylation. "
            "Variant prediction rules also receive the actual sample REF -> ALT change; curated allele-change rules take precedence when the observed ALT allele or exact change has a known direction."
        ),
        "disclaimer": (
            "Predictive theses in this database are literature-guided research summaries derived from the bundled gene interpretation bundle. "
            "They are designed for exploratory synthesis in the UI and should not be treated as diagnostic or therapeutic claims."
        ),
        "seeded_markers": _build_seeded_markers(knowledge_base),
        "concrete_variant_prediction": concrete_variant_prediction,
        "variant_prediction_rules": variant_prediction_rules,
        "case_count": len(cases),
        "methylation_sources": METHYLATION_SOURCES,
        "cases": cases,
    }


def main() -> int:
    """Generate one synthesis JSON file per target gene."""
    GENE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    for gene_name in TARGET_GENES:
        knowledge_base = _load_interpretation_database(gene_name)
        synthesis_database = build_synthesis_database(knowledge_base)
        output_path = GENE_DATA_DIR / f"{gene_name.lower()}_synthesis.json"
        output_path.write_text(
            json.dumps(synthesis_database, indent=2, ensure_ascii=True) + "\n",
            encoding="utf-8",
        )
        print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
