#!/usr/bin/env python3
"""Generate gene-level predictive synthesis databases from the bundled interpretation JSON files."""

from __future__ import annotations

import json
import zipfile
from functools import lru_cache
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENE_DATA_DIR = PROJECT_ROOT / "src" / "gene_data"
GENE_DATA_BUNDLE_PATH = GENE_DATA_DIR / "gene_data_bundle.zip"
VERSION = "2026-06-07"
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
    "CLRN2",
    "ARHGAP10",
    "CCDC66",
    "TYW5",
    "ELOVL7",
    "SH3PXD2B",
    "FRMD3",
    "TMEM218",
    "FAM170A",
    "SYCE3",
    "POTEB3",
    "BLTP3B",
    "CIROP",
    "MT-RNR1",
    "IL33",
    "IL1RL1",
    "ORMDL3",
    "GSDMB",
    "HLA-DQA1",
    "HLA-DQB1",
    "TSLP",
    "IL4R",
    "STAT6",
    "IL13",
    "IL4",
    "FLG",
    "TLR10",
    "TNFRSF8",
    "CD30",
    "MUC5AC",
    "SMAD3",
    "IL18R1",
    "IL18RAP",
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
        "A GT-confirmed HERC2/OCA2 genotype can support a pigmentation prediction rather than a HERC2 protein-disease prediction. "
        "For rs12913832, genotype dosage matters: A/A, A/G, and G/G should be interpreted separately. The forward-strand G allele used by this workbench supports a lighter or blue-eye tendency because it weakens HERC2 enhancer looping to the OCA2 promoter, lowers OCA2 expression, and reduces iris melanin, while A is the darker-eye-compatible contrast state. "
        "Linked HERC2 markers should be read as probabilistic light-versus-dark eye-colour haplotype evidence, with ancestry and additional OCA2/pigmentation variants able to modify the visible result."
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
    "CLRN2": (
        "The variant observed in this sample suggests a CLRN2 stereocilia-maintenance and DFNB117 hearing-loss thesis: the individual may carry a rare clarin-2 marker relevant to auditory hair-bundle maintenance, mechanotransduction, or autosomal recessive nonsyndromic sensorineural hearing-loss review. "
        "Only biallelic pathogenic or clinically confirmed CLRN2 findings should be escalated toward a DFNB117 disease interpretation; single heterozygous or VUS markers remain carrier or research context."
    ),
    "ARHGAP10": (
        "The variant observed in this sample suggests an ARHGAP10 RhoGAP and neuronal-morphology thesis: the individual may carry a research-grade signal relevant to RhoA/Cdc42 signaling, cytoskeletal organization, exonic CNV schizophrenia literature, or cancer-cell migration biology. "
        "This is rare-variant and pathway context, not a schizophrenia diagnosis, cancer prediction, or stand-alone clinical classification."
    ),
    "CCDC66": (
        "The variant observed in this sample suggests a CCDC66 microtubule, ciliary-transition-zone, and retinal-development thesis: the individual may carry a research-grade signal relevant to photoreceptor inner-segment biology, retinal degeneration model systems, high-myopia sequencing literature, centrosome or centriolar-satellite trafficking, CEP290/PCM1 interaction, cilium length and signaling, or mitotic/cytokinetic microtubule organization. "
        "Because ClinGen has not published a CCDC66 clinical curation, CCDC66 findings should stay phenotype-, variant-classification-, zygosity-, and assay-aware rather than being treated as deterministic retinal-disease, myopia, cancer, or ciliopathy diagnoses."
    ),
    "TYW5": (
        "The variant observed in this sample suggests a TYW5 tRNA(Phe) hydroxywybutosine and schizophrenia regulatory-expression thesis: the individual may carry a research-grade marker relevant to Fe(II)/2-oxoglutarate JmjC RNA hydroxylase biology, wybutosine-tRNA modification, translational reading-frame fidelity, brain TYW5 eQTL regulation, neurodevelopment, dendritic spine morphology, or structural-MRI schizophrenia genetics. "
        "Because ClinGen has not published a TYW5 clinical curation, TYW5 findings should stay cohort-, tissue-, expression-, genotype-, and phenotype-aware rather than being treated as deterministic schizophrenia, neurodevelopmental, retinal, cancer, or RNA-modification disorder predictions."
    ),
    "ELOVL7": (
        "The variant observed in this sample suggests an ELOVL7 ER fatty-acid elongase and VLCFA-remodeling thesis: the individual may carry a research-grade signal relevant to saturated or polyunsaturated very-long-chain fatty-acid metabolism, prostate-cancer lipid biosynthesis, HCMV lipidome remodeling, necroptotic membrane biology, liver-fibrosis expression context, or MSA locus research. "
        "This is lipid-metabolism and locus context, not a diagnostic cancer, infection, fibrosis, neurodegenerative-disease, or inherited-disorder prediction."
    ),
    "SH3PXD2B": (
        "The variant observed in this sample suggests a SH3PXD2B/TKS4 podosome-adaptor and Frank-ter Haar syndrome thesis: the individual may carry a research-grade or clinical-review signal relevant to autosomal-recessive skeletal, ocular, cardiac, craniofacial, dermal, collagen-remodeling, or podosome/invadopodia biology. "
        "Only biallelic pathogenic, likely pathogenic, or well-supported loss-of-function SH3PXD2B findings should be escalated toward a Frank-ter Haar or Borrone dermato-cardio-skeletal syndrome interpretation; single heterozygous, benign, or VUS findings remain carrier or research context."
    ),
    "FRMD3": (
        "The variant observed in this sample suggests a FRMD3/protein 4.1O FERM-domain and diabetic-kidney-disease regulatory thesis: the individual may carry a research-grade signal relevant to kidney cytoskeletal architecture, albuminuria or diabetic nephropathy association studies, FRMD3/BMP-pathway promoter regulation, podocyte or tubular expression, or tumor-suppressor and focal-adhesion biology. "
        "This is complex-trait and pathway context, not a diagnostic kidney-disease, diabetes, retinopathy, cancer-risk, or monogenic-disease prediction."
    ),
    "TMEM218": (
        "The variant observed in this sample suggests a TMEM218 ciliary-transition-zone and Joubert-Meckel ciliopathy thesis: the individual may carry a research-grade or clinical-review signal relevant to autosomal-recessive primary-cilium diffusion-barrier biology, Joubert syndrome 39, Meckel syndrome, retinal dystrophy, cystic kidney disease, polydactyly, encephalocele, or NPHP-module interaction studies. "
        "Only biallelic pathogenic, likely pathogenic, or strongly supported TMEM218 loss/reduced-function findings should be escalated toward a ciliopathy interpretation; single heterozygous, benign, or VUS findings remain carrier or research context."
    ),
    "FAM170A": (
        "The variant observed in this sample suggests a FAM170A/ZNFD nuclear zinc-finger and spermiogenesis thesis: the individual may carry a research-grade signal relevant to testis-enriched transcription-factor biology, AP1 or heat-shock-element regulation, sperm chromatin remodeling, or histone-to-protamine exchange. "
        "This is male-fertility and transcriptional-regulation context, not a diagnostic infertility prediction or validated monogenic disease call."
    ),
    "SYCE3": (
        "The variant observed in this sample suggests a SYCE3 synaptonemal-complex central-element thesis: the individual may carry a research-grade signal relevant to meiotic prophase I homolog synapsis, central-element assembly, recombination progression, crossover formation, or spermatocyte/oocyte meiotic-arrest biology. "
        "This is fertility and meiotic-structure context, not a diagnostic infertility, non-obstructive azoospermia, or premature ovarian insufficiency prediction."
    ),
    "POTEB3": (
        "The variant observed in this sample suggests a POTEB3/POTE-family structural-region thesis: the individual may carry a signal in a highly paralogous 15q11.2 cancer-testis gene region where assembly choice, copy number, and read-mapping uniqueness matter more than single-SNV heuristics. "
        "This is a paralog-aware research flag, not a POTEB3-specific diagnosis or methylation prediction."
    ),
    "BLTP3B": (
        "The variant observed in this sample suggests a BLTP3B/UHRF1BP1L membrane-trafficking and MYP3-locus thesis: the individual may carry a research-grade marker tied to bridge-like lipid transfer, early endosome-to-Golgi traffic, or cohort-level high-grade myopia association context. "
        "This is an exploratory locus signal, not a diagnostic myopia prediction or a proven BLTP3B protein-function call."
    ),
    "CIROP": (
        "The variant observed in this sample suggests a rare left-right asymmetry and heterotaxy thesis: the individual may carry a CIROP marker that belongs in developmental laterality-disorder review because CIROP is implicated in ciliated left-right organizer biology and autosomal visceral heterotaxy 12. "
        "This is a rare-disease research flag that depends on zygosity, phase, inheritance, phenotype match, and confirmatory clinical testing, not an app-only diagnosis."
    ),
    "MT-RNR1": (
        "The variant observed in this sample suggests a mitochondrial 12S rRNA pharmacogenetic thesis: the individual may carry an MT-RNR1 genotype that changes aminoglycoside-induced hearing-loss risk classification. "
        "For m.1555A>G, m.1494C>T, and m.1095T>C, the concrete prediction is increased susceptibility to aminoglycoside cochleotoxicity that requires clinical prescribing review, not an app-only medication decision."
    ),
}

VARIANT_CONCRETE_PREDICTION_OVERRIDES = {
    "HERC2": {
        "rs12913832": (
            "GT-decoded rs12913832 is the strongest concrete prediction in this bundle: in the forward-strand A/G representation used here, G dosage supports a lighter or blue-eye tendency and A dosage supports a darker or brown-eye-compatible tendency because this enhancer variant changes HERC2-OCA2 looping and OCA2-driven iris melanin production. "
            "A/G should reduce certainty relative to G/G, and a small SNP subset should never be treated as a final eye-colour call."
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
            "In the bundled evidence, T-allele homozygosity was associated with optimal glycemic response to liraglutide in one Iranian type 2 diabetes cohort, so this app only applies that direction when decoded GT supports the relevant dosage."
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
    "CDKN2A": {
        "rs11515": (
            "Observed CDKN2A rs11515 suggests a common 3'UTR regulatory-context thesis. "
            "This marker is best interpreted as low-effect melanoma/cancer association and haplotype context, with benign germline ClinVar framing, not as high-penetrance CDKN2A cancer predisposition."
        ),
        "rs3088440": (
            "Observed CDKN2A rs3088440 suggests a common 3'UTR melanoma-risk modifier thesis. "
            "Some cohorts report melanoma-risk or haplotype associations, but ClinVar submitter summaries frame the germline variant as benign, so the prediction should remain population- and haplotype-contextual."
        ),
        "rs3731249": (
            "Observed CDKN2A rs3731249 / p.Ala148Thr suggests a benign-polymorphism thesis: this sample carries a historically studied p16INK4A missense marker that current ClinVar curation classifies as benign. "
            "Do not convert this into a CDKN2A cancer-predisposition prediction without other pathogenic evidence."
        ),
        "CDKN2A p.Gly101Trp": (
            "Observed CDKN2A p.Gly101Trp / G101W suggests a rare pathogenic tumor-suppressor thesis. "
            "If the exact GRCh37 C -> A allele is GT-confirmed, the sample carries a high-priority CDKN2A melanoma and melanoma-pancreatic cancer predisposition marker that warrants external clinical confirmation rather than a stand-alone app diagnosis."
        ),
    },
    "BLTP3B": {
        "rs7134216": (
            "Observed BLTP3B/UHRF1BP1L rs7134216 suggests a MYP3 high-grade myopia locus thesis: the sample matched an intronic marker reported in association mapping and replicated for quantitative spherical refractive error in an independent high-grade myopia cohort. "
            "Keep the prediction cohort-level and research-oriented because the bundled evidence does not establish deterministic myopia risk or a direct BLTP3B protein mechanism."
        ),
    },
    "POTEB3": {
        "15q11.1-q11.2 CNV including POTEB3": (
            "Observed 15q11.1-q11.2 copy-number context including POTEB3 suggests a regional structural-variation thesis: the sample matched a broad segmental-duplication CNV context that includes POTEB3 among many genes. "
            "Keep this interpretation regional and assembly-aware; it is not a POTEB3-specific pathogenic or protective call."
        ),
    },
    "CLRN2": {
        "CLRN2 c.494C>A": (
            "Observed CLRN2 c.494C>A / p.Thr165Lys suggests a DFNB117 hearing-loss thesis when genotype dosage and inheritance support a biallelic recessive model. "
            "The practical prediction is rare hearing-loss follow-up context, especially phase, second-allele, phenotype, audiology, and clinical-grade confirmation review, not an app-only diagnosis."
        ),
        "CLRN2 c.236G>T": (
            "Observed CLRN2 c.236G>T / p.Arg79Leu suggests a CLRN2 VUS hearing-loss review thesis. "
            "Keep this as candidate rare-variant context unless newer ClinVar, segregation, functional, or phenotype evidence supports reclassification."
        ),
    },
    "ARHGAP10": {
        "ARHGAP10 exonic CNV": (
            "Observed ARHGAP10 exonic CNV context suggests a rare structural-variant schizophrenia-research thesis: the sample matched a copy-number marker class reported in a Japanese case-control study and biologically tied to RhoGAP/RhoA neuronal-morphology mechanisms. "
            "Keep this interpretation breakpoint-, assay-, and ancestry-aware; it is not a deterministic ARHGAP10 or schizophrenia diagnosis."
        ),
        "ARHGAP10 p.Ser490Pro": (
            "Observed ARHGAP10 p.Ser490Pro / rs483352828 suggests a rare RhoGAP-domain missense thesis. "
            "The strongest bundled evidence is the reported double-hit context with an exonic ARHGAP10 deletion on the other allele, so a single heterozygous VCF row should remain rare-variant research context unless CNV, phase, phenotype, and external review support stronger interpretation."
        ),
    },
    "CCDC66": {
        "CCDC66 c.C172T / p.Q58X": (
            "Observed CCDC66 c.C172T / p.Q58X suggests a high-myopia and retinal-development thesis because this suspected nonsense variant was reported to co-segregate with high myopia in a family, with additional rare CCDC66 variants observed in sporadic high-myopia cases. "
            "Keep this interpretation research-grade and ophthalmology-review oriented because the evidence is emerging, ClinGen has not published a CCDC66 curation, and variant dosage, transcript, segregation, retinal phenotype, and current clinical databases matter."
        ),
        "CCDC66 retinal degeneration/ciliary loss model": (
            "Observed CCDC66 loss, frameshift, deletion, or ciliary functional context suggests a retinal-ciliopathy model thesis: the sample matched a marker class tied to photoreceptor inner-segment expression, canine and mouse retinal degeneration, CEP290/PCM1 centriolar-satellite trafficking, cilium length/signaling, Hedgehog/Wnt pathway response, and mitotic/cytokinetic microtubule organization. "
            "Treat this as model-system and pathway context unless a clinically reviewed human pathogenic CCDC66 allele, zygosity, segregation, and phenotype fit are available."
        ),
    },
    "TYW5": {
        "rs796364 / rs281759": (
            "Observed TYW5-linked rs796364 or rs281759 suggests a 2q33.1 schizophrenia regulatory-expression thesis: the sample matched a functional regulatory marker pair reported to disrupt CTCF/RAD21/FOXP2 binding, physically interact with TYW5, associate with brain TYW5 expression, and affect TYW5 expression after CRISPR deletion of nearby regulatory sequence. "
            "Keep this interpretation cohort-level and non-diagnostic because schizophrenia risk is polygenic, the markers are noncoding and LD-dependent, and individual interpretation depends on ancestry, genotype dosage, brain-expression context, phenotype, and current GWAS/eQTL evidence."
        ),
        "rs203772": (
            "Observed TYW5 rs203772 suggests an integrative schizophrenia eQTL and neuroimaging thesis: the sample matched a marker whose risk allele was associated with higher TYW5 transcription in prefrontal cortex, stronger schizophrenia association in Sherlock/SMR analyses, higher TYW5 expression in schizophrenia brain or induced-neuron data, and gray-matter volume differences in first-episode antipsychotic-naive schizophrenia MRI analyses. "
            "Treat this as research-grade neurogenetic context rather than a diagnostic psychiatric prediction."
        ),
        "TYW5 enzymatic wybutosine-hydroxylase model": (
            "Observed TYW5 enzymatic, loss-of-function, overexpression, or RNA-modification context suggests a JmjC RNA-hydroxylase thesis: the sample matched a model class tied to Fe(II)/2-oxoglutarate-dependent hydroxylation of yW-72 to OHyW* in tRNA(Phe), tRNA-binding Arg residues, homodimer structure, translational reading-frame fidelity, neurodevelopmental expression effects, and dendritic spine morphology. "
            "Keep this as pathway and expression biology unless a clinically reviewed human TYW5 variant, RNA-modification assay, or disease model supports stronger interpretation."
        ),
    },
    "ELOVL7": {
        "rs7715147": (
            "Observed ELOVL7 rs7715147 suggests an MSA GWAS-interest lipid-dyshomeostasis locus thesis: the sample matched an intronic marker retained for neurogenetic research context around ELOVL7 and very-long-chain fatty-acid biology. "
            "Keep this locus-level and non-diagnostic because follow-up sequencing and copy-number work did not support rare ELOVL7 coding variants or CNV as a major MSA risk driver."
        ),
        "ELOVL7 functional lipid-elongation model": (
            "Observed ELOVL7 functional or expression context suggests a VLCFA-remodeling thesis: the sample matched a model class tied to ER fatty-acid elongation, prostate-cancer lipid biosynthesis, HCMV lipidome remodeling, necroptotic membrane disruption, or liver-fibrosis expression studies. "
            "Treat this as pathway and regulatory biology rather than a specific inherited variant or stand-alone disease prediction."
        ),
    },
    "SH3PXD2B": {
        "SH3PXD2B c.76-2A>C": (
            "Observed SH3PXD2B c.76-2A>C / rs775217258 suggests a high-priority Frank-ter Haar syndrome splice-acceptor thesis because ClinVar/OMIM literature curates this germline variant as pathogenic in an affected homozygous patient. "
            "Interpret dosage and phase carefully: a single heterozygous call is carrier context, while biallelic pathogenic SH3PXD2B findings warrant clinical genetics review for the autosomal-recessive FTHS/BDCS spectrum."
        ),
        "SH3PXD2B loss-of-function/deletion model": (
            "Observed SH3PXD2B loss-of-function or exon-deletion context suggests a TKS4 podosome-loss thesis: the sample matched the established autosomal-recessive mechanism for Frank-ter Haar syndrome and overlapping Borrone dermato-cardio-skeletal syndrome, with skeletal, ocular, cardiac, craniofacial, dermal, fibrosis, and collagen-remodeling relevance. "
            "Keep interpretation zygosity-, breakpoint-, transcript-, and phenotype-aware, because heterozygous carrier findings and broad CNV calls should not be treated as diagnostic without confirmatory review."
        ),
    },
    "FRMD3": {
        "rs1888747": (
            "Observed FRMD3 rs1888747 suggests a diabetic-kidney-disease regulatory-locus thesis: the sample matched a promoter-proximal marker reported as the strongest FRMD3-region signal in type 1 diabetes nephropathy GWAS and replicated or re-evaluated in several type 2 diabetes kidney cohorts. "
            "Keep this as probabilistic complex-trait context because the literature is cohort-dependent, C/C has been reported as protective in one T2D DKD cohort, and kidney expression studies did not show a direct rs1888747 genotype effect on FRMD3 mRNA or protein."
        ),
        "FRMD3 tumor-suppressor/cytoskeletal model": (
            "Observed FRMD3 functional or expression context suggests a protein 4.1O cytoskeletal and tumor-suppressor thesis: the sample matched a model class tied to FERM-domain membrane/cytoskeletal organization, NSCLC tumor-suppressor candidate evidence, breast-cancer vimentin degradation, focal-adhesion impairment, and migration or invasion biology. "
            "Treat this as pathway and expression-model context, not a germline cancer-risk or therapy-selection prediction."
        ),
    },
    "TMEM218": {
        "TMEM218 c.111G>T / p.Arg37Ser": (
            "Observed TMEM218 c.111G>T / p.Arg37Ser suggests a high-priority Joubert syndrome 39 and syndromic-ciliopathy thesis because the homozygous missense variant was reported in an affected individual with JBTS/BBS-like ciliopathy features and ClinVar curates the variant as likely pathogenic for JBTS39. "
            "Interpret dosage and phase carefully: a single heterozygous call is carrier context, while biallelic pathogenic or likely pathogenic TMEM218 findings warrant clinical genetics review for the Joubert-Meckel ciliopathy spectrum."
        ),
        "TMEM218 biallelic Joubert-Meckel ciliopathy model": (
            "Observed TMEM218 biallelic loss, reduced-function, or transition-zone variant context suggests a ciliary diffusion-barrier thesis: the sample matched a model class tied to TMEM218/MKS-module biology, TMEM67 interaction, ARL13B/GPR161 ciliary membrane localization, and Joubert-Meckel phenotypes including retinal dystrophy, molar-tooth sign, encephalocele, cystic kidneys, polydactyly, and perinatal lethality at the severe end. "
            "Keep this interpretation transcript-, zygosity-, phase-, phenotype-, and variant-classification-aware because TMEM218 genotype-phenotype severity is variant-specific and heterozygous findings alone do not establish disease."
        ),
    },
    "FAM170A": {
        "FAM170A loss-of-function/deletion model": (
            "Observed FAM170A loss-of-function or deletion context suggests a male-fertility and sperm chromatin-remodeling thesis: the sample matched a model-system marker class tied to Fam170a-deficient mouse infertility, sperm-head abnormality, and impaired histone-to-protamine exchange. "
            "Keep this interpretation research-grade and phenotype-aware because ClinGen has not published a FAM170A clinical gene-disease or dosage curation."
        ),
    },
    "SYCE3": {
        "SYCE3 loss-of-function/deletion model": (
            "Observed SYCE3 loss-of-function or deletion context suggests a meiotic-arrest and synaptonemal-complex central-element thesis: the sample matched a model-system marker class tied to Syce3 knockout infertility, failed synapsis initiation, absent MLH1 crossover foci, and impaired recombination progression. "
            "Keep this interpretation research-grade and phenotype-aware because confirmed human SYCE3 pathogenic mutations have not yet been established in the bundled review evidence."
        ),
    },
    "CIROP": {
        "CIROP c.92C>T": (
            "Observed CIROP c.92C>T / p.Ser31Phe suggests a rare CIROP heterotaxy 12 thesis. "
            "Treat this as developmental laterality-disorder context that needs ancestry-frequency, zygosity, phase, inheritance, and phenotype review before clinical interpretation."
        ),
        "CIROP c.571C>T": (
            "Observed CIROP c.571C>T / p.Arg191Ter suggests the strongest bundled CIROP loss-of-function thesis: a truncating marker classified in ClinVar for autosomal visceral heterotaxy 12. "
            "The practical prediction is rare-disease follow-up context, especially second-allele and inheritance review, not an app-only diagnosis."
        ),
        "CIROP c.1037G>A": (
            "Observed CIROP c.1037G>A / p.Trp346Ter suggests a truncating CIROP heterotaxy 12 thesis. "
            "Prioritize clinical-grade confirmation, phase, and phenotype match because CIROP-associated heterotaxy is interpreted through a recessive developmental-disease model."
        ),
        "CIROP c.1151C>T": (
            "Observed CIROP c.1151C>T / p.Ser384Leu suggests a rare CIROP missense heterotaxy 12 thesis. "
            "This is a laterality-disorder sequencing flag, not a broad adult-trait prediction."
        ),
        "CIROP c.1166G>T": (
            "Observed CIROP c.1166G>T / p.Arg389Ile suggests a rare CIROP missense heterotaxy 12 thesis. "
            "Interpretation should center on phenotype consistency, zygosity, and inheritance rather than methylation or common-variant heuristics."
        ),
        "CIROP c.1364TCT[1]": (
            "Observed CIROP c.1364TCT[1] / p.Phe456del suggests an in-frame deletion CIROP heterotaxy 12 thesis. "
            "Because indel representation can vary, reconcile rsID, transcript HGVS, and normalized coordinates before relying on the match."
        ),
    },
    "MT-RNR1": {
        "m.1555A>G": (
            "Observed MT-RNR1 m.1555A>G suggests the strongest bundled aminoglycoside ototoxicity thesis: CPIC assigns this genotype to increased risk of aminoglycoside-induced hearing loss. "
            "The practical prediction is that systemic aminoglycosides should be avoided unless infection severity and lack of alternatives justify the risk under clinician oversight."
        ),
        "m.1494C>T": (
            "Observed MT-RNR1 m.1494C>T suggests a high-evidence aminoglycoside ototoxicity thesis: CPIC assigns this genotype to increased risk, with ClinVar and hearing-loss literature supporting aminoglycoside-induced deafness relevance."
        ),
        "m.1095T>C": (
            "Observed MT-RNR1 m.1095T>C suggests a CPIC increased-risk aminoglycoside ototoxicity thesis with a more moderate evidence base than m.1555A>G and m.1494C>T. "
            "Interpret the call with mitochondrial heteroplasmy, depth, and test-method caveats."
        ),
        "m.827A>G": (
            "Observed MT-RNR1 m.827A>G suggests a normal-risk contrast thesis for aminoglycoside-induced hearing loss under current CPIC framing. "
            "Do not upgrade this allele to an aminoglycoside contraindication without newer expert guidance or additional variant evidence."
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
                    "This sample's rs12913832 site definition is A -> G, and interpretation now depends on the decoded GT dosage rather than ALT presence alone. "
                    "In the forward-strand hg19 A/G representation used by this workbench, G dosage is the light-eye-associated state: it reduces HERC2 enhancer support for OCA2 expression, lowers iris melanin biology, and therefore directionally supports a blue or lighter-eye tendency. "
                    "If GT is A/G rather than G/G, reduce certainty and keep brown or hazel plausible."
                ),
                "basis": "GT-confirmed G dosage at rs12913832, the HERC2/OCA2 enhancer state associated with reduced OCA2 expression and lighter iris pigmentation.",
            },
            {
                "change": "G>A",
                "alt_allele": "A",
                "prediction": (
                    "This sample's rs12913832 site definition is G -> A, and interpretation now depends on decoded GT dosage rather than ALT presence alone. "
                    "In the forward-strand hg19 A/G representation used by this workbench, A is the darker-eye-compatible state: it is compatible with stronger enhancer-promoter contact, higher OCA2 expression, and greater iris melanin biology, so A dosage supports a brown or darker-eye tendency. "
                    "This remains probabilistic because eye colour is polygenic."
                ),
                "basis": "GT-confirmed A dosage at rs12913832, the HERC2/OCA2 enhancer state associated with stronger OCA2 expression and darker iris pigmentation.",
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
                    "Only treat the T-homozygous liraglutide-response evidence as applicable when decoded GT supports T/T; heterozygous T dosage remains directional response-context evidence rather than a final responder call."
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
                    "G is the ancestral/reference state in the local Ensembl record, so avoid applying the p.R131Q direction unless decoded GT and allele orientation support it."
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
    "CDKN2A": {
        "rs3088440": [
            {
                "change": "G>A",
                "alt_allele": "A",
                "prediction": (
                    "This sample row reports A as the observed alternate allele at CDKN2A rs3088440 in the hg19 forward-strand VCF representation. "
                    "The same marker is often described as 540C>T or c.*69C>T in transcript-oriented literature; interpret A dosage as common 3'UTR regulatory and melanoma-association context, not as a pathogenic CDKN2A finding."
                ),
                "basis": "Observed ALT A at rs3088440, a common CDKN2A 3'UTR marker discussed in 9p21 melanoma-risk and haplotype studies.",
            },
        ],
        "rs3731249": [
            {
                "change": "C>T",
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at CDKN2A rs3731249 / p.Ala148Thr. "
                    "Current ClinVar curation classifies p.Ala148Thr as benign, so the concrete prediction is benign-polymorphism context rather than inherited CDKN2A cancer predisposition."
                ),
                "basis": "Observed ALT T at rs3731249, interpreted through ClinVar benign curation for CDKN2A p.Ala148Thr.",
            },
        ],
        "CDKN2A p.Gly101Trp": [
            {
                "change": "C>A",
                "alt_allele": "A",
                "prediction": (
                    "This sample row reports A as the observed alternate allele at chr9:21,971,057, the hg19 forward-strand representation of CDKN2A c.301G>T / p.Gly101Trp. "
                    "GT-confirmed A dosage is a rare pathogenic G101W signal in the local bundle, so the thesis should escalate to external clinical confirmation and genetics review while avoiding a diagnosis from this app alone."
                ),
                "basis": "Exact GRCh37 C -> A match for CDKN2A p.Gly101Trp / G101W, a ClinVar pathogenic melanoma and melanoma-pancreatic cancer marker.",
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
    "ARHGAP10": {
        "ARHGAP10 p.Ser490Pro": [
            {
                "change": "T>C",
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele for ARHGAP10 c.1468T>C / p.Ser490Pro. "
                    "If GT and QC support non-reference dosage, the sample fits a rare ARHGAP10 RhoGAP-domain missense review thesis. "
                    "The strongest bundled evidence requires checking for a second ARHGAP10 hit such as an exonic deletion, so report this as research context unless CNV, phase, phenotype, and external clinical review support escalation."
                ),
                "basis": "Observed C dosage at NM_024605.4(ARHGAP10):c.1468T>C / p.Ser490Pro, a rare missense marker highlighted in ARHGAP10 schizophrenia double-hit model literature.",
            }
        ],
    },
    "CLRN2": {
        "CLRN2 c.494C>A": [
            {
                "change": "C>A",
                "alt_allele": "A",
                "prediction": (
                    "This sample row reports A as the observed alternate allele at CLRN2 c.494C>A / p.Thr165Lys. "
                    "If GT, phase, and inheritance support biallelic pathogenic CLRN2 dosage, the sample fits a DFNB117 autosomal recessive nonsyndromic hearing-loss review thesis. "
                    "If only one ALT copy is present, report this as carrier or recessive-disease context rather than as a deterministic hearing-loss prediction."
                ),
                "basis": "Observed A dosage at NM_001079827.2(CLRN2):c.494C>A, the pathogenic CLRN2 DFNB117 allele with missense and splicing evidence.",
            }
        ],
        "CLRN2 c.236G>T": [
            {
                "change": "G>T",
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at CLRN2 c.236G>T / p.Arg79Leu. "
                    "Current local curation treats this as a VUS-level CLRN2 hearing-loss candidate, so the prediction should stay in rare-variant review context and should not be upgraded without external evidence."
                ),
                "basis": "Observed T dosage at NM_001079827.2(CLRN2):c.236G>T, a ClinVar VUS in DFNB117-related review.",
            }
        ],
    },
    "MT-RNR1": {
        "m.1555A>G": [
            {
                "change": "A>G",
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at MT-RNR1 m.1555A>G. "
                    "If genotype or heteroplasmy evidence supports the G state, the sample fits the CPIC increased-risk MT-RNR1 phenotype for aminoglycoside-induced hearing loss. "
                    "Use this as a high-priority pharmacogenetic safety flag requiring clinical confirmation and prescribing review."
                ),
                "basis": "Observed G dosage at NC_012920.1:m.1555A>G, a CPIC increased-risk MT-RNR1 aminoglycoside ototoxicity allele.",
            }
        ],
        "m.1494C>T": [
            {
                "change": "C>T",
                "alt_allele": "T",
                "prediction": (
                    "This sample row reports T as the observed alternate allele at MT-RNR1 m.1494C>T. "
                    "If genotype or heteroplasmy evidence supports the T state, the sample fits the CPIC increased-risk MT-RNR1 phenotype for aminoglycoside-induced hearing loss."
                ),
                "basis": "Observed T dosage at NC_012920.1:m.1494C>T, a CPIC increased-risk MT-RNR1 aminoglycoside ototoxicity allele.",
            }
        ],
        "m.1095T>C": [
            {
                "change": "T>C",
                "alt_allele": "C",
                "prediction": (
                    "This sample row reports C as the observed alternate allele at MT-RNR1 m.1095T>C. "
                    "If genotype or heteroplasmy evidence supports the C state, the sample fits CPIC's increased-risk MT-RNR1 phenotype, but with a more moderate evidence caveat than m.1555A>G or m.1494C>T."
                ),
                "basis": "Observed C dosage at NC_012920.1:m.1095T>C, a CPIC increased-risk MT-RNR1 aminoglycoside ototoxicity allele.",
            }
        ],
        "m.827A>G": [
            {
                "change": "A>G",
                "alt_allele": "G",
                "prediction": (
                    "This sample row reports G as the observed alternate allele at MT-RNR1 m.827A>G. "
                    "Current CPIC recommendations use m.827A>G as a normal-risk example for aminoglycoside-induced hearing loss, so this match should dampen rather than escalate aminoglycoside ototoxicity concern."
                ),
                "basis": "Observed G dosage at NC_012920.1:m.827A>G, a CPIC normal-risk example allele.",
            }
        ],
    },
}


ALLERGIC_INFLAMMATION_CONCRETE_PREDICTIONS = {
    "IL33": (
        "The variant observed in this sample suggests an epithelial alarmin thesis: the individual may carry a research-grade modifier of IL-33 release, IL1RL1/ST2 signaling, type 2 inflammation, asthma, allergy, or atopic dermatitis context. "
        "This is allergic-airway susceptibility context, not a deterministic diagnosis."
    ),
    "IL1RL1": (
        "The variant observed in this sample suggests an IL-33/ST2 receptor thesis: the individual may carry a research-grade modifier of soluble or membrane ST2 biology, eosinophilic airway inflammation, asthma, or atopy. "
        "Interpret it with the broader 2q12 IL1 receptor-cluster LD context."
    ),
    "ORMDL3": (
        "The variant observed in this sample suggests a 17q12-q21 asthma-locus thesis: the individual may carry a regulatory haplotype relevant to ORMDL3, GSDMB, sphingolipid homeostasis, childhood asthma, or allergic airway inflammation. "
        "This is locus-level context rather than a single-gene causal assignment."
    ),
    "GSDMB": (
        "The variant observed in this sample suggests a GSDMB/17q12-q21 thesis: the individual may carry a research-grade marker relevant to gasdermin B, epithelial immune biology, asthma, or linked ORMDL3/GSDMB regulatory haplotypes. "
        "Separate direct GSDMB protein hypotheses from linked-locus expression effects."
    ),
    "HLA-DQA1": (
        "The variant observed in this sample suggests an HLA-DQ antigen-presentation thesis: the individual may carry MHC class II haplotype context relevant to immune, autoimmune, allergic, or asthma association studies. "
        "HLA interpretation usually requires haplotype-aware review, not isolated SNP overcalling."
    ),
    "HLA-DQB1": (
        "The variant observed in this sample suggests an HLA-DQ antigen-presentation thesis: the individual may carry MHC class II haplotype context relevant to immune, autoimmune, allergic, or asthma association studies. "
        "HLA interpretation usually requires haplotype-aware review, not isolated SNP overcalling."
    ),
    "TSLP": (
        "The variant observed in this sample suggests an epithelial TSLP alarmin thesis: the individual may carry a research-grade modifier of TSLP expression, dendritic-cell priming, type 2 inflammation, asthma, atopic dermatitis, or allergic disease context. "
        "This does not create a medication-response rule."
    ),
    "IL4R": (
        "The variant observed in this sample suggests an IL-4/IL-13 receptor signaling thesis: the individual may carry a research-grade modifier of type 2 cytokine signaling, IgE biology, asthma, atopy, or eczema context. "
        "Do not use this alone for biologic therapy selection."
    ),
    "STAT6": (
        "The variant observed in this sample suggests a STAT6 type 2 transcriptional signaling thesis: the individual may carry a research-grade modifier of IL-4/IL-13 pathway output, allergic sensitization, eosinophilia, asthma, or epithelial remodeling. "
        "This is pathway tuning context."
    ),
    "IL13": (
        "The variant observed in this sample suggests an IL-13 type 2 cytokine thesis: the individual may carry a research-grade modifier of mucus biology, airway remodeling, IgE, asthma, eczema, or allergic inflammation. "
        "Interpret it with IL4/IL13 5q31 cytokine-cluster context."
    ),
    "IL4": (
        "The variant observed in this sample suggests an IL-4 Th2 cytokine thesis: the individual may carry a research-grade modifier of Th2 differentiation, IgE biology, asthma, atopy, or allergic sensitization. "
        "Interpret it with IL4/IL13 5q31 cytokine-cluster context."
    ),
    "FLG": (
        "The variant observed in this sample suggests a filaggrin epithelial-barrier thesis: the individual may carry a barrier-function signal relevant to atopic dermatitis, ichthyosis vulgaris, allergic sensitization, or atopic asthma context. "
        "Loss-of-function markers need zygosity and clinical-grade confirmation before escalation."
    ),
    "TLR10": (
        "The variant observed in this sample suggests a TLR10 innate-immune thesis: the individual may carry a research-grade modifier of toll-like receptor cluster signaling, microbial exposure response, asthma, allergy, or immune regulation. "
        "Interpret it with nearby TLR1/TLR6/TLR10 LD context."
    ),
    "TNFRSF8": (
        "The variant observed in this sample suggests a TNFRSF8/CD30 immune-activation thesis: the individual may carry research context relevant to activated lymphocytes, CD30 expression, inflammatory immune tone, or CD30-positive lymphoproliferative biology. "
        "Expression and phenotype context matter more than common-SNP determinism."
    ),
    "CD30": (
        "The variant observed in this sample suggests a CD30/TNFRSF8 immune-activation thesis: the individual may carry research context relevant to activated lymphocytes, CD30 expression, inflammatory immune tone, or CD30-positive lymphoproliferative biology. "
        "This alias bundle uses TNFRSF8 coordinates and evidence."
    ),
    "MUC5AC": (
        "The variant observed in this sample suggests a MUC5AC airway mucus thesis: the individual may carry research context relevant to goblet-cell metaplasia, mucus hypersecretion, asthma, chronic airway inflammation, or epithelial regulation. "
        "Methylation and tissue-specific expression context are especially important."
    ),
    "SMAD3": (
        "The variant observed in this sample suggests a SMAD3/TGF-beta remodeling thesis: the individual may carry a research-grade modifier of airway remodeling, fibrosis, immune regulation, asthma, or inflammatory disease context. "
        "Common association markers should not be confused with rare high-impact SMAD3 clinical variants."
    ),
    "IL18R1": (
        "The variant observed in this sample suggests an IL-18 receptor alpha thesis: the individual may carry a research-grade modifier of IL-18 signaling, Th1/NK-cell inflammation, eosinophil traits, asthma, or inflammatory disease context. "
        "Interpret it with the 2q12 IL1 receptor-cluster LD context."
    ),
    "IL18RAP": (
        "The variant observed in this sample suggests an IL-18 receptor accessory-protein thesis: the individual may carry a research-grade modifier of IL-18 receptor-complex signaling, inflammatory disease, asthma, or immune activation context. "
        "Interpret it with IL18R1 and IL1RL1 cluster context."
    ),
}

GENE_CONCRETE_VARIANT_PREDICTIONS.update(ALLERGIC_INFLAMMATION_CONCRETE_PREDICTIONS)


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


@lru_cache(maxsize=1)
def _gene_data_bundle_members() -> tuple[str, ...]:
    """Return sorted gene-data bundle member names, or an empty tuple when absent."""
    if not GENE_DATA_BUNDLE_PATH.exists():
        return ()
    try:
        with zipfile.ZipFile(GENE_DATA_BUNDLE_PATH) as bundle:
            return tuple(sorted(name for name in bundle.namelist() if not name.endswith("/")))
    except zipfile.BadZipFile:
        return ()


def _read_gene_data_bundle_text(member_name: str) -> str | None:
    """Read a text member from the gene-data bundle if present."""
    if member_name not in _gene_data_bundle_members():
        return None
    with zipfile.ZipFile(GENE_DATA_BUNDLE_PATH) as bundle:
        return bundle.read(member_name).decode("utf-8")


def _load_interpretation_database(gene_name: str) -> dict[str, Any]:
    """Load the bundled interpretation database for a gene."""
    for path in _candidate_interpretation_paths(gene_name):
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
    for path in _candidate_interpretation_paths(gene_name):
        payload = _read_gene_data_bundle_text(path.name)
        if payload is not None:
            return json.loads(payload)
    raise FileNotFoundError(f"No interpretation database found for {gene_name}")


def _discover_interpretation_genes() -> list[str]:
    """Return gene names for interpretation databases present in the local bundle directory."""
    discovered: list[str] = []
    for path in sorted(GENE_DATA_DIR.glob("*_interpretation_db.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        gene_name = _clean_text(payload.get("gene_context", {}).get("gene_name"))
        if gene_name:
            discovered.append(gene_name)
    for member_name in _gene_data_bundle_members():
        if not member_name.endswith("_interpretation_db.json"):
            continue
        payload_text = _read_gene_data_bundle_text(member_name)
        if payload_text is None:
            continue
        try:
            payload = json.loads(payload_text)
        except json.JSONDecodeError:
            continue
        gene_name = _clean_text(payload.get("gene_context", {}).get("gene_name"))
        if gene_name:
            discovered.append(gene_name)
    return _dedupe_text_items(discovered)


def _iter_synthesis_gene_names() -> list[str]:
    """Return curated targets plus any available interpretation-only bundles."""
    return _dedupe_text_items([*TARGET_GENES, *_discover_interpretation_genes()])


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


def _concrete_variant_prediction_for_gene(gene_name: str, knowledge_base: dict[str, Any]) -> str:
    """Return a concrete gene-level variant prediction text."""
    curated_prediction = _clean_text(
        knowledge_base.get("gene_context", {}).get("concrete_variant_prediction", "")
    )
    if curated_prediction:
        return curated_prediction

    return GENE_CONCRETE_VARIANT_PREDICTIONS.get(
        gene_name,
        (
            f"A GT-confirmed non-reference genotype in this sample suggests a {gene_name} gene-specific research thesis. "
            "Interpret the result through decoded genotype dosage, call QC, and the bundled gene context; do not treat it as a deterministic clinical prediction."
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
        f"This sample matched {display_name} at site definition {{change}} with GT {{gt_raw}} decoded as {{genotype}} ({{zygosity}}). "
        f"Anchor the {gene_name} prediction to genotype dosage {{allele_dosage}} rather than REF -> ALT presence alone. "
        "Dosage-specific conclusions should be softened when GT, depth, balance, or genotype quality are uncertain."
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
        record_prediction = _first_nonempty(
            record.get("concrete_prediction"),
            _variant_prediction_override(gene_name, record),
        )
        if not record_prediction:
            record_prediction = _clean_text(
                f"GT-confirmed {display_name} dosage suggests this sample carries the {gene_name} prediction context described by the gene-level thesis. "
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
    band_context: str,
    clinical_context: str,
    variant_focus: str,
    methylation_context: str,
    concrete_variant_prediction: str,
    research_focus: list[str],
) -> dict[str, Any]:
    """Build one of the nine variant-plus-methylation synthesis cases."""
    prediction = _clean_text(
        f"When a {gene_name} variant is paired with {band} methylation in the {source_label.lower()}, "
        f"the sample best fits a combined regulatory-context thesis. {band_context} "
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
    methylation_band_interpretations = gene_context.get(
        "methylation_band_interpretations", {}
    )
    concrete_variant_prediction = _concrete_variant_prediction_for_gene(gene_name, knowledge_base)
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
                    band_context=_first_nonempty(
                        methylation_band_interpretations.get(band)
                        if isinstance(methylation_band_interpretations, dict)
                        else "",
                        BAND_CONTEXT[band],
                    ),
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
            "Variant prediction rules receive the decoded sample GT, zygosity, ALT dosage, and site REF -> ALT definition; curated allele rules take precedence only when the sample genotype carries the relevant allele dosage."
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
    skipped: list[str] = []
    for gene_name in _iter_synthesis_gene_names():
        try:
            knowledge_base = _load_interpretation_database(gene_name)
        except FileNotFoundError:
            skipped.append(gene_name)
            continue
        synthesis_database = build_synthesis_database(knowledge_base)
        output_path = GENE_DATA_DIR / f"{gene_name.lower()}_synthesis.json"
        output_path.write_text(
            json.dumps(synthesis_database, indent=2, ensure_ascii=True) + "\n",
            encoding="utf-8",
        )
        print(f"Wrote {output_path}")
    if skipped:
        print("Skipped missing interpretation databases:")
        for gene_name in skipped:
            print(f" - {gene_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
