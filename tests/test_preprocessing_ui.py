"""Regression tests for preprocessing form behavior in the Flask UI."""

from __future__ import annotations

import io
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from src.analysis import load_gene_interpretation_database, load_gene_population_database
from src.webapp import _build_data_sources_payload, _classify_functional_family, app


def _stub_common_discovery(monkeypatch) -> None:
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_bam_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_report_history", lambda: [])


def test_data_sources_payload_combines_curated_and_dynamic_sources() -> None:
    """The Result Viewer Data Sources payload should merge curated and dynamic provenance."""
    knowledge_base = load_gene_interpretation_database("DRD4")
    population_database = load_gene_population_database("DRD4")
    assert knowledge_base is not None
    assert population_database is not None

    payload = _build_data_sources_payload(
        knowledge_base=knowledge_base,
        population_database=population_database,
        population_insights={
            "summary": "Population summary",
            "variant_population_records": [{"variant": "rs1"}],
            "gene_population_patterns": [],
            "sources": [{"label": "Population source", "url": "https://example.com/pop"}],
        },
        methylation_insights={
            "summary": "Methylation summary",
            "clinical_context": "Methylation context",
            "evidence": [{"label": "Methylation paper", "url": "https://example.com/methylation"}],
            "whitelist_probe_reference_rows": [
                {"probe_id": "cg1", "papers": [{"label": "Probe paper", "url": "https://example.com/probe"}]}
            ],
        },
        dynamic_payload={
            "provider_statuses": [
                {
                    "source_key": "clinvar",
                    "name": "ClinVar",
                    "lane": "clinical",
                    "status": "ok",
                    "message": "ClinVar returned one record.",
                    "record_count": 1,
                    "homepage": "https://www.ncbi.nlm.nih.gov/clinvar/",
                },
                {
                    "source_key": "clingen",
                    "name": "ClinGen",
                    "lane": "clinical",
                    "status": "ok",
                    "message": "Queried ClinGen; 1 gene-centered curation record(s) returned for DRD4.",
                    "record_count": 1,
                    "homepage": "https://clinicalgenome.org/",
                },
                {
                    "source_key": "medgen",
                    "name": "MedGen",
                    "lane": "clinical",
                    "status": "ok",
                    "message": "Queried MedGen; 1 condition or phenotype record(s) returned for DRD4.",
                    "record_count": 1,
                    "homepage": "https://www.ncbi.nlm.nih.gov/medgen/",
                },
                {
                    "source_key": "ensembl",
                    "name": "Ensembl",
                    "lane": "clinical",
                    "status": "ok",
                    "message": "Queried Ensembl; 2 record(s).",
                    "record_count": 2,
                    "homepage": "https://rest.ensembl.org/",
                },
                {
                    "source_key": "dbsnp",
                    "name": "dbSNP",
                    "lane": "population",
                    "status": "ok",
                    "message": "Queried NCBI snp; 1 record(s) returned.",
                    "record_count": 1,
                    "homepage": "https://www.ncbi.nlm.nih.gov/snp/",
                },
                {
                    "source_key": "gnomad",
                    "name": "gnomAD",
                    "lane": "population",
                    "status": "ok",
                    "message": "Queried gnomAD; variant frequency and gene constraint records returned for 11-637293-C-T.",
                    "record_count": 2,
                    "homepage": "https://gnomad.broadinstitute.org/",
                },
                {
                    "source_key": "gwas_catalog",
                    "name": "GWAS Catalog",
                    "lane": "population",
                    "status": "ok",
                    "message": (
                        "Queried GWAS Catalog; 1 GWAS association record(s) returned for DRD4; "
                        "no rsID-specific associations were found for rs927984495."
                    ),
                    "record_count": 1,
                    "homepage": "https://www.ebi.ac.uk/gwas/",
                },
                {
                    "source_key": "pgs_catalog",
                    "name": "PGS Catalog",
                    "lane": "population",
                    "status": "ok",
                    "message": "Queried PGS Catalog; 1 linked polygenic score record(s) returned for rs123.",
                    "record_count": 1,
                    "homepage": "https://www.pgscatalog.org/",
                },
                {
                    "source_key": "igsr",
                    "name": "1000 Genomes Project (IGSR)",
                    "lane": "population",
                    "status": "ok",
                    "message": (
                        "Queried IGSR/1000 Genomes FTP release listings; "
                        "1 data-access record returned for rs927984495 at 11:637293 C>T."
                    ),
                    "record_count": 1,
                    "homepage": "https://www.internationalgenome.org/",
                },
                {
                    "source_key": "ucsc",
                    "name": "UCSC Genome Browser",
                    "lane": "regulatory",
                    "status": "ok",
                    "message": (
                        "Queried UCSC Genome Browser API tracks for hg38 chr11:637293-640706; "
                        "1 compact annotation record(s) returned."
                    ),
                    "record_count": 1,
                    "homepage": "https://api.genome.ucsc.edu/",
                },
                {
                    "source_key": "encode",
                    "name": "ENCODE Portal (DCC)",
                    "lane": "regulatory",
                    "status": "ok",
                    "message": "Queried ENCODE Portal DCC; 1 experiment record(s) returned for DRD4.",
                    "record_count": 1,
                    "homepage": "https://www.encodeproject.org/",
                },
                {
                    "source_key": "ewas_catalog",
                    "name": "EWAS Catalog",
                    "lane": "regulatory",
                    "status": "ok",
                    "message": "Queried EWAS Catalog gene results for DRD4; 1 compact association record(s) returned.",
                    "record_count": 1,
                    "homepage": "https://www.ewascatalog.org/",
                },
                {
                    "source_key": "ewas_atlas",
                    "name": "EWAS Atlas / EWAS Open Platform",
                    "lane": "regulatory",
                    "status": "ok",
                    "message": (
                        "Queried EWAS Atlas REST position query hg19 chr11:637293-640706; "
                        "1 compact association record(s) returned for DRD4."
                    ),
                    "record_count": 1,
                    "homepage": "https://ngdc.cncb.ac.cn/ewas/",
                },
                {
                    "source_key": "screen",
                    "name": "SCREEN",
                    "lane": "regulatory",
                    "status": "ok",
                    "message": (
                        "Queried SCREEN cCRE Registry for chr11:637293-640706; "
                        "1 overlapping cCRE record(s) returned for DRD4."
                    ),
                    "record_count": 1,
                    "homepage": "https://screen.encodeproject.org/",
                },
                {
                    "source_key": "civic",
                    "name": "CIViC",
                    "lane": "clinical",
                    "status": "ok",
                    "message": "Queried CIViC; 1 variant evidence record(s) returned for DRD4.",
                    "record_count": 1,
                    "homepage": "https://civicdb.org/",
                },
                {
                    "source_key": "panelapp",
                    "name": "PanelApp",
                    "lane": "clinical",
                    "status": "ok",
                    "message": "Queried PanelApp; 1 exact gene panel record(s) returned for BRCA1.",
                    "record_count": 1,
                    "homepage": "https://panelapp.genomicsengland.co.uk/",
                },
                {
                    "source_key": "mavedb",
                    "name": "MaveDB",
                    "lane": "clinical",
                    "status": "ok",
                    "message": "Queried MaveDB; 1 published score set record(s) returned for BRCA1.",
                    "record_count": 1,
                    "homepage": "https://www.mavedb.org/",
                },
                {
                    "source_key": "hgmd",
                    "name": "HGMD",
                    "lane": "licensed",
                    "status": "needs_export",
                    "message": "Upload a permitted export.",
                    "record_count": 0,
                    "warnings": ["HGMD needs a licensed export."],
                    "license_note": "License-gated source.",
                },
            ],
            "source_records": [
                {
                    "source_key": "clinvar",
                    "category": "clinical",
                    "label": "DRD4 c.1A>G",
                    "summary": "Clinical significance: Likely benign; Variant type: single nucleotide variant",
                    "source_id": "123",
                    "variant": "rs1800955",
                    "rsid": "rs1800955",
                    "url": "https://example.com/clinvar-drd4",
                },
                {
                    "source_key": "clingen",
                    "category": "gene_disease_validity",
                    "label": "ClinGen validity: DRD4 - Example syndrome",
                    "summary": "Definitive gene-disease validity for Example syndrome (AD)",
                    "classification": "Definitive",
                    "disease": "Example syndrome",
                    "url": "https://search.clinicalgenome.org/kb/gene-validity/DRD4",
                },
                {
                    "source_key": "medgen",
                    "category": "clinical_condition",
                    "label": "DRD4-related phenotype",
                    "summary": "MedGen clinical condition associated with DRD4.",
                    "concept_id": "C123456",
                    "url": "https://www.ncbi.nlm.nih.gov/medgen/C123456",
                    "research_links": [
                        {
                            "label": "MedGen GTR clinical-test records for DRD4",
                            "url": "https://www.ncbi.nlm.nih.gov/medgen/?term=%22medgen%20gtr%20tests%20clinical%22%5BFilter%5D%20AND%20DRD4%5Bgene%5D",
                        }
                    ],
                },
                {
                    "source_key": "ensembl",
                    "category": "variant_annotation",
                    "label": "rs927984495",
                    "summary": (
                        "rs927984495 (C/T) at GRCh37 11:637293-637293: "
                        "5_prime_UTR_variant in DRD4 canonical transcript ENST00000176183, "
                        "exon 1/4, cDNA position 1; variant class SNP; evidence: Frequency, TOPMed, gnomAD."
                    ),
                    "source_id": "rs927984495",
                    "variant": "11:637293:C>T",
                    "rsid": "rs927984495",
                    "url": "https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs927984495",
                },
                {
                    "source_key": "dbsnp",
                    "category": "population",
                    "label": "rs2533154733",
                    "summary": (
                        "rs2533154733: deletion at 11:637373 (NC_000011.10:637372:T:); "
                        "DRD4 coding_sequence_variant / frameshift_variant; "
                        "HGVS NC_000011.10:g.637373del, NM_000797.4:c.69del; "
                        "frequency GnomAD_exomes 1/998848; validated by-frequency; submitted by GNOMAD."
                    ),
                    "source_id": "2533154733",
                    "variant": "rs2533154733",
                    "rsid": "rs2533154733",
                    "url": "https://www.ncbi.nlm.nih.gov/snp/2533154733",
                },
                {
                    "source_key": "gnomad",
                    "category": "population_frequency",
                    "label": "rs927984495",
                    "summary": (
                        "gnomAD v4 rs927984495 (11-637293-C-T) at GRCh38 11:637293 C>T: "
                        "joint AF 1.67e-05 (AC 20/AN 1,195,848; hom 0; FAF95 nfe 1.35e-05); "
                        "genomes AF 3.97e-05 (AC 6/AN 151,146; hom 0; FAF95 nfe 3.78e-05); "
                        "exomes AF 1.34e-05 (AC 14/AN 1,044,702; hom 0; FAF95 nfe 9.05e-06); "
                        "highest observed population genome nfe AF 8.86e-05 (AC 6/AN 67,718); "
                        "consequence 5_prime_UTR_variant in DRD4 transcript ENST00000176183; "
                        "HGVS ENST00000176183.5:c.-11C>T; predictors: CADD 3.2; "
                        "filters: discrepant_frequencies."
                    ),
                    "source_id": "11-637293-C-T",
                    "variant": "rs927984495",
                    "rsid": "rs927984495",
                    "url": "https://gnomad.broadinstitute.org/variant/11-637293-C-T?dataset=gnomad_r4",
                },
                {
                    "source_key": "gnomad",
                    "category": "gene_constraint",
                    "label": "DRD4 gnomAD constraint",
                    "summary": (
                        "gnomAD v4 DRD4 gene constraint (ENSG00000069696): pLI 1.03e-10; "
                        "LoF O/E 1.175 (0.512-2.211), Z -0.6971; "
                        "missense O/E 1.444 (1.18-1.76), Z -4.759; "
                        "LoF observed/expected 4/3.405; missense observed/expected 114/79.21; "
                        "constraint flags: outlier_mis, outlier_syn; gene location GRCh38 11:637269-640706."
                    ),
                    "source_id": "ENSG00000069696",
                    "gene": "DRD4",
                    "url": "https://gnomad.broadinstitute.org/gene/DRD4?dataset=gnomad_r4",
                },
                {
                    "source_key": "gwas_catalog",
                    "category": "population_association",
                    "label": "rs1870723-A - Hypothyroidism",
                    "summary": (
                        "GWAS Catalog association rs1870723-A with Hypothyroidism / hypothyroidism "
                        "mapped to DRD4 at 11:640349: p=5e-18; "
                        "beta 0.0352 unit decrease; CI 0.027-0.043; risk frequency 0.2688; "
                        "study GCST90627750, PMID 41644669, first author White SL; "
                        "initial sample 257,365 cases, 2,186,763 controls; "
                        "full summary statistics available; strongest risk allele rs1870723-A; "
                        "data release 2026-06-22."
                    ),
                    "source_id": "216548355",
                    "variant": "rs927984495",
                    "rsids": ["rs1870723"],
                    "accession_id": "GCST90627750",
                    "url": "https://www.ebi.ac.uk/gwas/rest/api/v2/associations/216548355",
                },
                {
                    "source_key": "pgs_catalog",
                    "category": "polygenic_score",
                    "label": "PGS000001 - PRS77_BC - Breast cancer",
                    "summary": (
                        "PGS Catalog score PGS000001 (PRS77_BC) includes variant rs123 and predicts "
                        "Breast cancer: mapped traits: breast carcinoma; 77 variants; "
                        "method SNPs passing genome-wide significance (P<5x10-8); weight type beta; "
                        "variant-source sample 22,627 European (Finland, Sweden, U.S., Australia, "
                        "Netherlands, Germany, U.K.) source GCST001937; "
                        "GWAS ancestry n=22,627 (EUR 100%); evaluation ancestry n=10 (EUR 80%, NR 20%); "
                        "performance All breast cancer: OR 1.55 (1.52-1.58), "
                        "C-index 0.622 (0.619-0.627); evaluated sample 67,054 European; "
                        "publication Mavaddat N, 2015-04-08, PMID 25855707, DOI 10.1093/jnci/djv036; "
                        "harmonized scoring files: GRCh37, GRCh38; released 2019-10-14."
                    ),
                    "source_id": "PGS000001",
                    "variant": "rs123",
                    "pgs_id": "PGS000001",
                    "url": "https://www.pgscatalog.org/score/PGS000001/",
                },
                {
                    "source_key": "igsr",
                    "category": "population_reference_panel",
                    "label": "1000 Genomes 2504 high-coverage phased callset chr11",
                    "summary": (
                        "IGSR/1000 Genomes high-coverage GRCh38 data-access context for "
                        "rs927984495 at 11:637293 C>T: chromosome VCF "
                        "CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz "
                        "(1.6G, updated 2020-10-29 16:36) and tabix index are available "
                        "for the 30x NYGC 2504-sample Phase 3 panel; use the indexed VCF "
                        "to extract exact genotypes, AC/AN, or allele counts for this coordinate; "
                        "legacy Phase 3 GRCh37 integrated data are also available via "
                        "ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "
                        "(701M) with global and superpopulation AF tags; rsIDs were removed "
                        "from the Phase 3 v5b VCF, so coordinate lookup or Ensembl rsID mapping is needed."
                    ),
                    "source_id": "1000G_2504_high_coverage_GRCh38_chr11",
                    "variant": "rs927984495",
                    "rsid": "rs927984495",
                    "url": (
                        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
                        "1000G_2504_high_coverage/working/20201028_3202_phased/"
                        "CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz"
                    ),
                },
                {
                    "source_key": "ucsc",
                    "category": "gene_model",
                    "label": "UCSC ncbiRefSeq NM_000797.4",
                    "summary": (
                        "UCSC ncbiRefSeq NM_000797.4 for DRD4 on hg38 chr11:637269-640706 (+): "
                        "4 exons; CDS chr11:637305-640603; complete CDS start/end; "
                        "query window hg38 chr11:637293-640706."
                    ),
                    "source_id": "NM_000797.4",
                    "gene": "DRD4",
                    "transcript": "NM_000797.4",
                    "url": "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr11:637293-640706",
                },
                {
                    "source_key": "encode",
                    "category": "regulatory_experiment",
                    "label": "ENCSRDRD4TF - TF ChIP-seq",
                    "summary": (
                        "ENCODE DCC ENCSRDRD4TF TF ChIP-seq matched DRD4 through target gene metadata: "
                        "target DRD4 (DRD4); biosample Homo sapiens K562 cell line; "
                        "Homo sapiens; 2 biological replicate(s); "
                        "2/2 released file(s); outputs optimal IDR thresholded peaks, signal p-value; "
                        "assemblies GRCh38; lab ENCODE Processing Pipeline; project ENCODE; "
                        "status released; released 2026-01-15."
                    ),
                    "source_id": "ENCSRDRD4TF",
                    "gene": "DRD4",
                    "assay": "TF ChIP-seq",
                    "target": "DRD4",
                    "url": "https://www.encodeproject.org/experiments/ENCSRDRD4TF/",
                },
                {
                    "source_key": "ewas_catalog",
                    "category": "ewas_association",
                    "label": "cg02762115 - age",
                    "summary": (
                        "EWAS Catalog DRD4 CpG cg02762115 at chr11:640446: "
                        "outcome DNA methylation; exposure/trait age; tissue Whole blood; "
                        "N=2,338; beta 0.0021; p=0E+00; PMID 33450751, Mulder RH; "
                        "analysis Model 1 with age as a fixed effect."
                    ),
                    "source_id": "EWAS0001",
                    "gene": "DRD4",
                    "cpg": "cg02762115",
                    "location": "chr11:640446",
                    "url": "https://www.ewascatalog.org/?cpg=cg02762115",
                },
                {
                    "source_key": "ewas_atlas",
                    "category": "ewas_association",
                    "label": "cg02762115 - cognitive function",
                    "summary": (
                        "EWAS Atlas DRD4 probe cg02762115 at hg19 chr11:640446: "
                        "Shelf; DRD4 transcript ENST00000176183.5, 3153 bp from TSS; "
                        "trait cognitive function; negative methylation-trait correlation; "
                        "rank 4; study ES00743; PMID 29311653; "
                        "position query hg19 chr11:637293-640706."
                    ),
                    "source_id": "ES00743:cg02762115",
                    "gene": "DRD4",
                    "probe_id": "cg02762115",
                    "trait": "cognitive function",
                    "rank": "4",
                    "url": "https://ngdc.cncb.ac.cn/ewas/search?item=cg02762115&term=Probe+Id",
                },
                {
                    "source_key": "screen",
                    "category": "candidate_regulatory_element",
                    "label": "EH38E2937824 - PLS",
                    "summary": (
                        "SCREEN cCRE EH38E2937824 overlaps DRD4 query window at "
                        "GRCh38 chr11:637026-637375: promoter-like signature (PLS); "
                        "CTCF-bound; intersecting genes: DRD4 (protein_coding); "
                        "max assay Z-scores DNase 3.52, H3K4me3 4.00, H3K27ac 2.03, CTCF 1.88; "
                        "rDHS EH38D4573084; query window chr11:637293-640706."
                    ),
                    "source_id": "EH38E2937824",
                    "gene": "DRD4",
                    "ccre_group": "PLS",
                    "location": "GRCh38 chr11:637026-637375",
                    "url": "https://screen.encodeproject.org/search/?q=EH38E2937824&assembly=GRCh38",
                },
                {
                    "source_key": "civic",
                    "category": "cancer_variant",
                    "label": "DRD4 V194G",
                    "summary": (
                        "CIViC DRD4 variant DRD4 V194G (variant type missense_variant; aliases: rs1800955; "
                        "HGVS NM_000797.4:c.581T>G): PREDICTIVE evidence, level B, rating 4, SUPPORTS; "
                        "significance SENSITIVITYRESPONSE; disease Example cancer; therapies Example therapy; "
                        "source PMID 123456; status ACCEPTED."
                    ),
                    "source_id": "101",
                    "variant": "DRD4 V194G",
                    "evidence_type": "PREDICTIVE",
                    "evidence_level": "B",
                    "disease": "Example cancer",
                    "therapies": ["Example therapy"],
                    "url": "https://civicdb.org/variants/101",
                },
                {
                    "source_key": "panelapp",
                    "category": "gene_panel",
                    "label": "BRCA1 in Inherited ovarian cancer (without breast cancer)",
                    "summary": (
                        "BRCA1 is listed on PanelApp panel Inherited ovarian cancer (without breast cancer) "
                        "(public, v5.1); confidence 3 (green/high evidence); "
                        "phenotypes: {Breast-ovarian cancer, familial, 1}, OMIM:604370; "
                        "inheritance: BOTH monoallelic and biallelic, autosomal or pseudoautosomal; "
                        "evidence: NHS GMS, Expert Review Green, Expert list."
                    ),
                    "source_id": "143:BRCA1",
                    "gene": "BRCA1",
                    "panel_name": "Inherited ovarian cancer (without breast cancer)",
                    "confidence_level": "3",
                    "confidence_label": "green/high evidence",
                    "url": "https://panelapp.genomicsengland.co.uk/panels/143/gene/BRCA1/",
                },
                {
                    "source_key": "mavedb",
                    "category": "functional_assay_score_set",
                    "label": "Scores from multiplexed functional assay of BRCA1 variants",
                    "summary": (
                        "MaveDB score set urn:mavedb:00001222-b-2 for BRCA1 "
                        "(BRCA1 DNA repair associated): "
                        "Scores from multiplexed functional assay of BRCA1 variants; "
                        "2,271 scored variants; published 2025-10-22; "
                        "experiment: Multiplexed functional assay of BRCA1 variants; "
                        "assay summary: Multiplexed assay of BRCA1 variants measuring "
                        "homology directed repair activity; target genes: BRCA1; "
                        "publications: PMID 39999999; license CC BY 4.0."
                    ),
                    "source_id": "urn:mavedb:00001222-b-2",
                    "gene": "BRCA1",
                    "score_set_urn": "urn:mavedb:00001222-b-2",
                    "num_variants": 2271,
                    "url": "https://www.mavedb.org/score-sets/urn:mavedb:00001222-b-2",
                },
            ],
            "literature_records": [
                {
                    "source_key": "europe_pmc",
                    "category": "literature",
                    "title": "DRD4 literature",
                    "summary": "Literature summary",
                    "url": "https://example.com/lit",
                }
            ],
            "local_article_evidence": {
                "status": "ok",
                "message": "Extracted one local snippet.",
                "records": [
                    {
                        "source_key": "local_pdf_articles",
                        "title": "Local DRD4 PDF",
                        "snippet": "DRD4 local finding",
                        "url": "https://example.com/local",
                    }
                ],
                "provenance": {"warnings": [], "errors": []},
            },
        },
        selected_source_keys=["clinvar", "hgmd"],
    )

    cards = [card for group in payload["groups"] for card in group["cards"]]
    by_key = {card["source_key"]: card for card in cards}

    assert payload["dynamic_status"] == "available"
    assert by_key["curated_gene_bundle"]["status"] == "ok"
    assert any("NCBI Gene 1815" in link["label"] for link in by_key["curated_gene_bundle"]["links"])
    assert by_key["clinvar"]["status"] == "ok"
    assert by_key["clinvar"]["findings"] == [
        "DRD4 c.1A>G (Variation ID 123; rs1800955): "
        "Clinical significance: Likely benign; Variant type: single nucleotide variant"
    ]
    assert by_key["clingen"]["status"] == "ok"
    assert by_key["clingen"]["findings"] == ["Definitive gene-disease validity for Example syndrome (AD)"]
    assert "Expert assessment of gene" not in by_key["clingen"]["summary"]
    assert by_key["medgen"]["status"] == "ok"
    assert by_key["medgen"]["findings"] == ["MedGen clinical condition associated with DRD4."]
    assert any("GTR clinical-test" in link["label"] for link in by_key["medgen"]["links"])
    assert by_key["ensembl"]["status"] == "ok"
    assert by_key["ensembl"]["findings"] == [
        "rs927984495 (C/T) at GRCh37 11:637293-637293: "
        "5_prime_UTR_variant in DRD4 canonical transcript ENST00000176183, "
        "exon 1/4, cDNA position 1; variant class SNP; evidence: Frequency, TOPMed, gnomAD."
    ]
    assert by_key["dbsnp"]["status"] == "ok"
    assert by_key["dbsnp"]["findings"] == [
        "rs2533154733: deletion at 11:637373 (NC_000011.10:637372:T:); "
        "DRD4 coding_sequence_variant / frameshift_variant; "
        "HGVS NC_000011.10:g.637373del, NM_000797.4:c.69del; "
        "frequency GnomAD_exomes 1/998848; validated by-frequency; submitted by GNOMAD."
    ]
    assert by_key["gnomad"]["status"] == "ok"
    assert by_key["gnomad"]["findings"] == [
        "gnomAD v4 rs927984495 (11-637293-C-T) at GRCh38 11:637293 C>T: "
        "joint AF 1.67e-05 (AC 20/AN 1,195,848; hom 0; FAF95 nfe 1.35e-05); "
        "genomes AF 3.97e-05 (AC 6/AN 151,146; hom 0; FAF95 nfe 3.78e-05); "
        "exomes AF 1.34e-05 (AC 14/AN 1,044,702; hom 0; FAF95 nfe 9.05e-06); "
        "highest observed population genome nfe AF 8.86e-05 (AC 6/AN 67,718); "
        "consequence 5_prime_UTR_variant in DRD4 transcript ENST00000176183; "
        "HGVS ENST00000176183.5:c.-11C>T; predictors: CADD 3.2; "
        "filters: discrepant_frequencies.",
        "gnomAD v4 DRD4 gene constraint (ENSG00000069696): pLI 1.03e-10; "
        "LoF O/E 1.175 (0.512-2.211), Z -0.6971; "
        "missense O/E 1.444 (1.18-1.76), Z -4.759; "
        "LoF observed/expected 4/3.405; missense observed/expected 114/79.21; "
        "constraint flags: outlier_mis, outlier_syn; gene location GRCh38 11:637269-640706.",
    ]
    assert by_key["gwas_catalog"]["status"] == "ok"
    assert by_key["gwas_catalog"]["findings"] == [
        "GWAS Catalog association rs1870723-A with Hypothyroidism / hypothyroidism "
        "mapped to DRD4 at 11:640349: p=5e-18; "
        "beta 0.0352 unit decrease; CI 0.027-0.043; risk frequency 0.2688; "
        "study GCST90627750, PMID 41644669, first author White SL; "
        "initial sample 257,365 cases, 2,186,763 controls; "
        "full summary statistics available; strongest risk allele rs1870723-A; "
        "data release 2026-06-22."
    ]
    assert by_key["pgs_catalog"]["status"] == "ok"
    assert by_key["pgs_catalog"]["findings"] == [
        "PGS Catalog score PGS000001 (PRS77_BC) includes variant rs123 and predicts "
        "Breast cancer: mapped traits: breast carcinoma; 77 variants; "
        "method SNPs passing genome-wide significance (P<5x10-8); weight type beta; "
        "variant-source sample 22,627 European (Finland, Sweden, U.S., Australia, "
        "Netherlands, Germany, U.K.) source GCST001937; "
        "GWAS ancestry n=22,627 (EUR 100%); evaluation ancestry n=10 (EUR 80%, NR 20%); "
        "performance All breast cancer: OR 1.55 (1.52-1.58), "
        "C-index 0.622 (0.619-0.627); evaluated sample 67,054 European; "
        "publication Mavaddat N, 2015-04-08, PMID 25855707, DOI 10.1093/jnci/djv036; "
        "harmonized scoring files: GRCh37, GRCh38; released 2019-10-14."
    ]
    assert by_key["igsr"]["status"] == "ok"
    assert by_key["igsr"]["findings"] == [
        "IGSR/1000 Genomes high-coverage GRCh38 data-access context for "
        "rs927984495 at 11:637293 C>T: chromosome VCF "
        "CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz "
        "(1.6G, updated 2020-10-29 16:36) and tabix index are available "
        "for the 30x NYGC 2504-sample Phase 3 panel; use the indexed VCF "
        "to extract exact genotypes, AC/AN, or allele counts for this coordinate; "
        "legacy Phase 3 GRCh37 integrated data are also available via "
        "ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "
        "(701M) with global and superpopulation AF tags; rsIDs were removed "
        "from the Phase 3 v5b VCF, so coordinate lookup or Ensembl rsID mapping is needed."
    ]
    assert by_key["ucsc"]["status"] == "ok"
    assert by_key["ucsc"]["findings"] == [
        "UCSC ncbiRefSeq NM_000797.4 for DRD4 on hg38 chr11:637269-640706 (+): "
        "4 exons; CDS chr11:637305-640603; complete CDS start/end; "
        "query window hg38 chr11:637293-640706."
    ]
    assert by_key["encode"]["status"] == "ok"
    assert by_key["encode"]["findings"] == [
        "ENCODE DCC ENCSRDRD4TF TF ChIP-seq matched DRD4 through target gene metadata: "
        "target DRD4 (DRD4); biosample Homo sapiens K562 cell line; "
        "Homo sapiens; 2 biological replicate(s); "
        "2/2 released file(s); outputs optimal IDR thresholded peaks, signal p-value; "
        "assemblies GRCh38; lab ENCODE Processing Pipeline; project ENCODE; "
        "status released; released 2026-01-15."
    ]
    assert by_key["ewas_catalog"]["status"] == "ok"
    assert by_key["ewas_catalog"]["findings"] == [
        "EWAS Catalog DRD4 CpG cg02762115 at chr11:640446: "
        "outcome DNA methylation; exposure/trait age; tissue Whole blood; "
        "N=2,338; beta 0.0021; p=0E+00; PMID 33450751, Mulder RH; "
        "analysis Model 1 with age as a fixed effect."
    ]
    assert by_key["ewas_atlas"]["status"] == "ok"
    assert by_key["ewas_atlas"]["findings"] == [
        "EWAS Atlas DRD4 probe cg02762115 at hg19 chr11:640446: "
        "Shelf; DRD4 transcript ENST00000176183.5, 3153 bp from TSS; "
        "trait cognitive function; negative methylation-trait correlation; "
        "rank 4; study ES00743; PMID 29311653; "
        "position query hg19 chr11:637293-640706."
    ]
    assert by_key["screen"]["status"] == "ok"
    assert by_key["screen"]["findings"] == [
        "SCREEN cCRE EH38E2937824 overlaps DRD4 query window at "
        "GRCh38 chr11:637026-637375: promoter-like signature (PLS); "
        "CTCF-bound; intersecting genes: DRD4 (protein_coding); "
        "max assay Z-scores DNase 3.52, H3K4me3 4.00, H3K27ac 2.03, CTCF 1.88; "
        "rDHS EH38D4573084; query window chr11:637293-640706."
    ]
    assert by_key["civic"]["status"] == "ok"
    assert by_key["civic"]["findings"] == [
        "CIViC DRD4 variant DRD4 V194G (variant type missense_variant; aliases: rs1800955; "
        "HGVS NM_000797.4:c.581T>G): PREDICTIVE evidence, level B, rating 4, SUPPORTS; "
        "significance SENSITIVITYRESPONSE; disease Example cancer; therapies Example therapy; "
        "source PMID 123456; status ACCEPTED."
    ]
    assert by_key["panelapp"]["status"] == "ok"
    assert by_key["panelapp"]["findings"] == [
        "BRCA1 is listed on PanelApp panel Inherited ovarian cancer (without breast cancer) "
        "(public, v5.1); confidence 3 (green/high evidence); "
        "phenotypes: {Breast-ovarian cancer, familial, 1}, OMIM:604370; "
        "inheritance: BOTH monoallelic and biallelic, autosomal or pseudoautosomal; "
        "evidence: NHS GMS, Expert Review Green, Expert list."
    ]
    assert by_key["mavedb"]["status"] == "ok"
    assert by_key["mavedb"]["findings"] == [
        "MaveDB score set urn:mavedb:00001222-b-2 for BRCA1 "
        "(BRCA1 DNA repair associated): "
        "Scores from multiplexed functional assay of BRCA1 variants; "
        "2,271 scored variants; published 2025-10-22; "
        "experiment: Multiplexed functional assay of BRCA1 variants; "
        "assay summary: Multiplexed assay of BRCA1 variants measuring "
        "homology directed repair activity; target genes: BRCA1; "
        "publications: PMID 39999999; license CC BY 4.0."
    ]
    assert by_key["hgmd"]["status"] == "needs_export"
    assert by_key["hgmd"]["license_note"] == "License-gated source."
    assert "europe_pmc_literature" in by_key
    assert by_key["local_pdf_articles"]["record_count"] == 1
    assert any(link["url"] == "https://example.com/probe" for link in by_key["methylation_evidence"]["links"])


def test_data_sources_payload_marks_selected_dynamic_sources_not_run() -> None:
    """When no dynamic KB is present, selected providers should appear as not_run cards."""
    payload = _build_data_sources_payload(
        knowledge_base={"database_name": "Mock KB", "gene_context": {}, "variant_records": []},
        population_database={"database_name": "Mock population DB"},
        population_insights={},
        methylation_insights={},
        dynamic_payload=None,
        selected_source_keys=["clinvar", "hgmd"],
    )

    cards = [card for group in payload["groups"] for card in group["cards"]]
    statuses = {card["source_key"]: card["status"] for card in cards}

    assert payload["dynamic_status"] == "not_run"
    assert statuses["clinvar"] == "not_run"
    assert statuses["hgmd"] == "not_run"


def test_preprocess_find_region_submission_updates_session(monkeypatch) -> None:
    """A valid preprocessing action should resolve the gene interval and persist it."""

    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr(
        "src.webapp.find_gene_region",
        lambda gene_symbol: {
            "gene_name": gene_symbol.upper(),
            "selected_region": "chr11:639677-643057",
            "selected_sources": ["NCBI RefSeq"],
            "candidate_regions": [{"source": "NCBI RefSeq", "region": "chr11:639677-643057"}],
        },
    )

    client = app.test_client()
    response = client.post(
        "/",
        data={
            "workflow": "preprocess",
            "gene_name": "drd4",
            "preprocess_action": "find_region",
        },
    )

    assert response.status_code == 200
    assert "Resolved DRD4 to standard promoter+gene region 11:636269-640706." in response.get_data(as_text=True)

    with client.session_transaction() as session_state:
        preprocess_state = session_state["preprocess_state"]

    assert preprocess_state["gene_name"] == "DRD4"
    assert preprocess_state["region"] == "11:636269-640706"
    assert preprocess_state["analysis_scope"] == "promoter_plus_gene"
    assert preprocess_state["scope_regions"]["promoter_plus_gene"] == "11:636269-640706"
    assert preprocess_state["scope_regions"]["promoter_only"] == "11:636269-637268"
    assert preprocess_state["scope_regions"]["gene_only"] == "11:637269-640706"
    assert preprocess_state["region_ready"] is True
    assert preprocess_state["manifest_ready"] is False
    assert preprocess_state["analysis_ready"] is False
    assert preprocess_state["selected_sources"] == ["NCBI RefSeq", "Local curated promoter/gene intervals"]


def test_mt_rnr1_zero_probe_preprocessing_unlocks_analysis(monkeypatch, tmp_path: Path) -> None:
    """Curated mitochondrial genes should unlock analysis even with a zero-row EPIC subset."""
    manifest_path = tmp_path / "manifest.csv"
    manifest_path.write_text("IlmnID,CHR,MAPINFO,UCSC_RefGene_Name\n", encoding="utf-8")
    output_path = tmp_path / "MT-RNR1_epigenetics_hg19.csv"
    captured_call: dict[str, object] = {}

    def fake_save_filtered_manifest(**kwargs: object) -> dict[str, object]:
        captured_call.update(kwargs)
        output_path.write_text("IlmnID,CHR,MAPINFO,UCSC_RefGene_Name\n", encoding="utf-8")
        return {"output_path": output_path, "probe_count": 0}

    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.save_filtered_manifest", fake_save_filtered_manifest)

    client = app.test_client()
    with client.session_transaction() as session_state:
        session_state["preprocess_state"] = {
            "gene_name": "MT-RNR1",
            "region": "MT:1-1601",
            "analysis_scope": "promoter_plus_gene",
            "scope_regions": {
                "promoter_plus_gene": "MT:1-1601",
                "promoter_only": "MT:1-647",
                "gene_only": "MT:648-1601",
            },
            "scope_region_source": "Local curated promoter/gene intervals",
            "manifest_source": str(manifest_path),
            "filtered_manifest": "",
            "region_candidates": [],
            "selected_sources": [],
            "region_ready": True,
            "manifest_ready": False,
            "analysis_ready": False,
            "probe_count": 0,
            "build": "hg19",
            "logs": [],
            "region_recently_updated": False,
            "overwrite_filtered_manifest": True,
        }

    response = client.post(
        "/",
        data={
            "workflow": "preprocess",
            "gene_name": "MT-RNR1",
            "preprocess_region": "MT:1-1601",
            "manifest_source": str(manifest_path),
            "overwrite_filtered_manifest": "1",
            "preprocess_action": "select_methylation",
        },
    )

    assert response.status_code == 200
    assert captured_call["allow_empty"] is True

    with client.session_transaction() as session_state:
        preprocess_state = session_state["preprocess_state"]

    assert preprocess_state["manifest_ready"] is True
    assert preprocess_state["analysis_ready"] is True
    assert preprocess_state["probe_count"] == 0


def test_functional_map_groups_genes_and_links_processed_reports(monkeypatch, tmp_path: Path) -> None:
    """The functional map should group knowledge-base genes and link completed reports."""
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    report_path = results_dir / "sirt6_report.html"
    report_path.write_text("<html><body>SIRT6 report</body></html>", encoding="utf-8")

    monkeypatch.setattr("src.webapp.RESULTS_DIR", results_dir)
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_bam_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])

    client = app.test_client()
    response = client.get("/functional-map")
    page = response.get_data(as_text=True)

    assert response.status_code == 200
    assert "NophiGene Functional Map" in page
    assert "Longevity &amp; Healthy Aging" in page
    assert "Senses &amp; Sensory Signaling" in page
    assert "Asthma, Allergy &amp; Airways" in page
    assert 'value="SIRT6"' in page
    assert "/results/sirt6_report.html" in page
    assert "Open latest report" in page


def test_functional_map_prioritizes_requested_families() -> None:
    """High-signal descriptions should land in the intended functional families."""
    assert _classify_functional_family("human longevity and centenarian healthy aging") == "longevity"
    assert _classify_functional_family("human sensory biology and visual phototransduction") == "senses"
    assert _classify_functional_family("asthma allergy airway eosinophil biology") == "asthma_allergy"


def test_v2_navigation_replaces_monolithic_workspace(monkeypatch) -> None:
    monkeypatch.setattr("src.webapp.discover_vcf_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_idat_prefixes", lambda: [])
    monkeypatch.setattr("src.webapp.discover_population_stats_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_manifest_files", lambda: [])
    monkeypatch.setattr("src.webapp.discover_report_history", lambda: [])

    page = app.test_client().get("/").get_data(as_text=True)
    for label in ("Run", "Results", "History", "Data Explorer", "Settings"):
        assert label in page
    assert "App Structure" not in page
    assert "Predictive Theses" not in page
