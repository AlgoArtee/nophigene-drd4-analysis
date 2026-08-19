"""Database connectors used by dynamic variant knowledge-base generation."""

from __future__ import annotations

import csv
import html
import io
import os
import re
import time
from html.parser import HTMLParser
from typing import Any
from urllib.parse import quote, urljoin
from xml.etree import ElementTree

from .client import KnowledgeRequestError, RequestClient
from .credentials import ResolvedCredential
from .models import KnowledgeQuery, SourceResult, SourceSpec


def _elapsed_ms(started: float) -> int:
    return int((time.monotonic() - started) * 1000)


def _first_variant(query: KnowledgeQuery):
    return query.variants[0] if query.variants else None


def _first_rsid(query: KnowledgeQuery) -> str:
    return query.rsids[0] if query.rsids else ""


def _literature_query(query: KnowledgeQuery) -> str:
    rsid = _first_rsid(query)
    if rsid:
        return f"{query.gene} {rsid}"
    return query.gene


NCBI_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_REFSNP_BASE_URL = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp"
NCBI_REFSNP_TIMEOUT_SECONDS = 8
NCBI_TOOL_NAME = "NophiGeneDynamicKB"
NLM_CLINICAL_TABLES_VARIANTS_URL = "https://clinicaltables.nlm.nih.gov/api/variants/v4/search"
CLINGEN_VALIDITY_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
CLINGEN_DOSAGE_URL = "https://search.clinicalgenome.org/kb/gene-dosage/download"
CLINGEN_SUMMARY_URL = "https://search.clinicalgenome.org/kb/reports/curation-activity-summary-report"
CLINGEN_ACTIONABILITY_URLS = (
    ("Adult", "https://actionability.clinicalgenome.org/ac/Adult/api/summ?flavor=flat"),
    ("Pediatric", "https://actionability.clinicalgenome.org/ac/Pediatric/api/summ?flavor=flat"),
)
CLINGEN_PRIMARY_TIMEOUT_SECONDS = 20
CLINGEN_SUMMARY_TIMEOUT_SECONDS = 20
CLINGEN_OPTIONAL_TIMEOUT_SECONDS = 8
MEDGEN_RETMAX = 20
MEDGEN_GTR_CLINICAL_FILTER = '"medgen gtr tests clinical"[Filter]'
GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api/"
GNOMAD_DATASET = "gnomad_r4"
GWAS_CATALOG_API_V2_BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"
GWAS_CATALOG_MAX_ASSOCIATIONS = 5
PGS_CATALOG_REST_BASE = "https://www.pgscatalog.org/rest"
PGS_CATALOG_MAX_LINKED_SCORES = 5
PGS_CATALOG_MAX_PERFORMANCE_RECORDS = 3
IGSR_PORTAL_URL = "https://www.internationalgenome.org/data-portal/"
IGSR_FTP_BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
IGSR_PHASE3_RELEASE_URL = f"{IGSR_FTP_BASE}/release/20130502/"
IGSR_PHASE3_README_URL = f"{IGSR_PHASE3_RELEASE_URL}README_phase3_callset_20150220"
IGSR_PHASE3_ANNOTATION_README_URL = f"{IGSR_PHASE3_RELEASE_URL}README_vcf_info_annotation.20141104"
IGSR_HIGH_COVERAGE_URL = f"{IGSR_FTP_BASE}/data_collections/1000G_2504_high_coverage/"
IGSR_HIGH_COVERAGE_PHASED_URL = f"{IGSR_HIGH_COVERAGE_URL}working/20201028_3202_phased/"
IGSR_HIGH_COVERAGE_README_URL = f"{IGSR_HIGH_COVERAGE_URL}20190405_1000G_2504_high_cov_README.md"
UCSC_API_BASE = "https://api.genome.ucsc.edu"
UCSC_TRACK_MAX_ITEMS = 5
ENCODE_PORTAL_BASE_URL = "https://www.encodeproject.org"
ENCODE_SEARCH_URL = f"{ENCODE_PORTAL_BASE_URL}/search/"
ENCODE_REGION_SEARCH_URL = f"{ENCODE_PORTAL_BASE_URL}/region-search/"
ENCODE_MAX_RECORDS = 5
SCREEN_PORTAL_BASE_URL = "https://screen.encodeproject.org"
SCREEN_GRAPHQL_URL = f"{SCREEN_PORTAL_BASE_URL}/api/screen-graphql"
SCREEN_MAX_RECORDS = 5
SCREEN_GRCH38_CCRE_BED_URL = "https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.bed"
EWAS_CATALOG_BASE_URL = "https://www.ewascatalog.org"
EWAS_CATALOG_MAX_RECORDS = 5
EWAS_ATLAS_REST_BASE = "https://ngdc.cncb.ac.cn/ewas/rest"
EWAS_ATLAS_PORTAL_BASE = "https://ngdc.cncb.ac.cn/ewas"
EWAS_ATLAS_MAX_RECORDS = 5
STRING_V12_API_BASE = "https://version-12-0.string-db.org/api/json"
STRING_V12_NETWORK_URL = f"{STRING_V12_API_BASE}/network"
STRING_DIRECT_PARTNER_LIMIT = 10
STRING_REQUIRED_SCORE = 400


def _ncbi_request_params(params: dict[str, Any]) -> dict[str, Any]:
    out = dict(params)
    tool = os.environ.get("NOPHIGENE_NCBI_TOOL", NCBI_TOOL_NAME).strip() or NCBI_TOOL_NAME
    out["tool"] = tool.replace(" ", "_")
    email = os.environ.get("NOPHIGENE_NCBI_EMAIL", "").strip()
    if email:
        out["email"] = email
    api_key = os.environ.get("NOPHIGENE_NCBI_API_KEY", "").strip() or os.environ.get("NCBI_API_KEY", "").strip()
    if api_key:
        out["api_key"] = api_key
    return out


def _request_failure_result(
    spec: SourceSpec,
    exc: KnowledgeRequestError,
    *,
    queried_urls: list[str],
    started: float,
) -> SourceResult:
    remediation = getattr(exc, "remediation", "")
    code = getattr(exc, "code", "request_failed")
    message = str(exc)
    return SourceResult(
        source_key=spec.key,
        status="failed",
        message=message,
        errors=[message],
        queried_urls=queried_urls,
        elapsed_ms=_elapsed_ms(started),
        error_code=code,
        remediation=remediation,
    )


def _clean_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(_clean_cell(item) for item in value if _clean_cell(item))
    if isinstance(value, dict):
        for key in ("description", "name", "label", "title", "trait_name", "variation_name"):
            cleaned = _clean_cell(value.get(key))
            if cleaned:
                return cleaned
        return ""
    return str(value).strip()


class BaseConnector:
    """Base connector contract."""

    def __init__(self, spec: SourceSpec, client: RequestClient, credential: ResolvedCredential) -> None:
        self.spec = spec
        self.client = client
        self.credential = credential

    def query(self, context: KnowledgeQuery) -> SourceResult:
        return SourceResult(
            source_key=self.spec.key,
            status="metadata_only",
            message=self.spec.license_note,
            records=[self._metadata_record()],
        )

    def _metadata_record(self) -> dict[str, Any]:
        return {
            "category": "source_metadata",
            "source": self.spec.name,
            "label": self.spec.name,
            "summary": self.spec.description or self.spec.license_note,
            "url": self.spec.homepage,
            "license_note": self.spec.license_note,
        }


class AuthMetadataConnector(BaseConnector):
    """Credential-aware connector for sources whose official API is not implemented in v1."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        if not self.credential.present:
            return SourceResult(
                source_key=self.spec.key,
                status="needs_credentials",
                message=f"Set {self.spec.env_var} or enter a session token before querying {self.spec.name}.",
                records=[self._metadata_record()],
            )
        return SourceResult(
            source_key=self.spec.key,
            status="metadata_only",
            message=(
                f"Credentials are available for {self.spec.name}, but this v1 connector only records "
                "metadata until the licensed endpoint contract is configured."
            ),
            records=[self._metadata_record()],
        )


class LicensedMetadataConnector(BaseConnector):
    """Explicitly non-scraping connector for licensed databases/search portals."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        return SourceResult(
            source_key=self.spec.key,
            status="metadata_only",
            message="Licensed or non-open source: no scraping is performed.",
            records=[self._metadata_record()],
        )


class _EwasCatalogDownloadParser(HTMLParser):
    """Extract the query-specific TSV download generated by EWAS Catalog."""

    def __init__(self) -> None:
        super().__init__()
        self.tsv_href = ""

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if self.tsv_href or tag.lower() != "a":
            return
        attr_map = {key.lower(): value or "" for key, value in attrs}
        href = attr_map.get("href", "")
        if href.endswith(".tsv") and ("download" in attr_map or "/tmp/" in href):
            self.tsv_href = href


class NcbiEutilsConnector(BaseConnector):
    """NCBI E-Utilities connector for ClinVar, dbSNP, Gene, PubMed, PMC, GEO, and LitVar-like evidence."""

    DB_BY_KIND = {
        "clinvar": "clinvar",
        "dbsnp": "snp",
        "ncbi_gene": "gene",
        "pubmed": "pubmed",
        "pmc": "pmc",
        "geo": "gds",
        "litvar": "pubmed",
    }

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        db = self.DB_BY_KIND.get(self.spec.connector_kind, "pubmed")
        term = self._term(context)
        records: list[dict[str, Any]] = []
        urls: list[str] = []
        warnings: list[str] = []
        try:
            search_url = f"{NCBI_EUTILS_BASE_URL}/esearch.fcgi"
            search = self.client.get_json(
                search_url,
                params=_ncbi_request_params({"db": db, "term": term, "retmode": "json", "retmax": 5}),
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(search_url)
            ids = list(search.get("esearchresult", {}).get("idlist", []))
            if ids:
                summary_url = f"{NCBI_EUTILS_BASE_URL}/esummary.fcgi"
                summary = self.client.get_json(
                    summary_url,
                    params=_ncbi_request_params({"db": db, "id": ",".join(ids[:5]), "retmode": "json"}),
                    rate_limit_per_second=self.spec.rate_limit_per_second,
                )
                urls.append(summary_url)
                result = summary.get("result", {})
                for item_id in ids[:5]:
                    item = result.get(str(item_id), {})
                    if self.spec.connector_kind == "dbsnp":
                        records.append(self._dbsnp_record(str(item_id), item, context, urls, warnings))
                        continue
                    title = str(item.get("title") or item.get("name") or item.get("uid") or item_id)
                    records.append(
                        {
                            "category": self._category(),
                            "source": self.spec.name,
                            "label": title,
                            "summary": str(item.get("description") or item.get("fulljournalname") or title),
                            "source_id": str(item_id),
                            "url": self._record_url(db, str(item_id)),
                            "variant": _first_rsid(context),
                        }
                    )
            return SourceResult(
                source_key=self.spec.key,
                status="ok",
                message=f"Queried NCBI {db}; {len(records)} record(s) returned.",
                records=records,
                warnings=warnings,
                queried_urls=urls,
                elapsed_ms=_elapsed_ms(started),
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=urls, started=started)

    def _dbsnp_record(
        self,
        item_id: str,
        item: dict[str, Any],
        context: KnowledgeQuery,
        urls: list[str],
        warnings: list[str],
    ) -> dict[str, Any]:
        rsid = self._dbsnp_rsid(item_id)
        esummary_fields = self._dbsnp_esummary_fields(item, item_id)
        detail_fields = self._dbsnp_refsnp_fields(item_id, context, urls, warnings)
        fields = self._merge_dbsnp_fields(esummary_fields, detail_fields)
        fields["rsid"] = fields.get("rsid") or rsid
        fields["source_id"] = fields.get("source_id") or str(item_id)
        label = fields["rsid"]
        record = {
            "category": self._category(),
            "source": self.spec.name,
            "label": label,
            "summary": self._dbsnp_summary(label, fields),
            "source_id": str(item_id),
            "url": self._record_url("snp", str(item_id)),
            "variant": _first_rsid(context) or label,
        }
        record.update(fields)
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _dbsnp_esummary_fields(self, item: dict[str, Any], item_id: str) -> dict[str, Any]:
        fields: dict[str, Any] = {
            "rsid": self._dbsnp_rsid(item.get("snp_id") or item.get("uid") or item_id),
            "source_id": _clean_cell(item.get("snp_id") or item.get("uid") or item_id),
            "snp_class": _clean_cell(item.get("snp_class")),
            "chromosome": _clean_cell(item.get("chr")),
            "chrpos": _clean_cell(item.get("chrpos")),
            "spdi": _clean_cell(item.get("spdi")),
            "validated": self._clean_string_list(item.get("validated")),
            "submitter_handle": _clean_cell(item.get("handle")),
            "createdate": _clean_cell(item.get("createdate")),
            "updatedate": _clean_cell(item.get("updatedate")),
        }
        if fields["chrpos"] and ":" in fields["chrpos"]:
            chromosome, position = fields["chrpos"].split(":", 1)
            fields["chromosome"] = fields["chromosome"] or chromosome
            fields["position"] = position
        fields["functional_classes"] = self._split_dbsnp_terms(item.get("fxn_class"))
        fields["genes"] = self._dbsnp_esummary_genes(item)
        fields["global_mafs"] = self._dbsnp_global_mafs(item.get("global_mafs"))
        if fields["global_mafs"]:
            fields["frequencies"] = list(fields["global_mafs"])

        docsum = self._parse_dbsnp_docsum(_clean_cell(item.get("docsum")))
        if docsum.get("HGVS"):
            fields["hgvs"] = self._dedupe_clean(docsum["HGVS"].split(","))
        if docsum.get("SEQ"):
            fields["alleles"] = docsum["SEQ"].strip("[]")
        if docsum.get("GENE"):
            gene_name, _, gene_id = docsum["GENE"].partition(":")
            gene_record = {
                "name": _clean_cell(gene_name),
                "gene_id": _clean_cell(gene_id),
            }
            fields["genes"] = self._dedupe_gene_records([*fields["genes"], gene_record])
        return {key: value for key, value in fields.items() if value not in ("", [], {}, None)}

    def _dbsnp_refsnp_fields(
        self,
        item_id: str,
        context: KnowledgeQuery,
        urls: list[str],
        warnings: list[str],
    ) -> dict[str, Any]:
        rsid = self._dbsnp_rsid(item_id)
        detail_url = f"{NCBI_REFSNP_BASE_URL}/{quote(str(item_id))}"
        fields: dict[str, Any] = {"refsnp_detail_url": detail_url}
        try:
            payload = self.client.get_json(
                detail_url,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
                timeout=NCBI_REFSNP_TIMEOUT_SECONDS,
            )
            urls.append(detail_url)
        except KnowledgeRequestError:
            warnings.append(f"Optional dbSNP RefSNP detail lookup failed for {rsid}; ESummary details were still used.")
            return fields

        if not isinstance(payload, dict):
            return fields
        source_id = _clean_cell(payload.get("refsnp_id") or item_id)
        if source_id:
            fields["source_id"] = source_id
            fields["rsid"] = self._dbsnp_rsid(source_id)
        fields["createdate"] = _clean_cell(payload.get("create_date")) or fields.get("createdate", "")
        fields["updatedate"] = _clean_cell(payload.get("last_update_date")) or fields.get("updatedate", "")

        snapshot = payload.get("primary_snapshot_data") if isinstance(payload.get("primary_snapshot_data"), dict) else {}
        fields["variant_type"] = _clean_cell(snapshot.get("variant_type"))
        fields["mane_select_ids"] = self._clean_string_list(snapshot.get("mane_select_ids"))

        placement = self._dbsnp_primary_placement(snapshot, context)
        fields.update(self._dbsnp_placement_fields(placement))
        fields.update(self._dbsnp_annotation_fields(snapshot.get("allele_annotations")))

        support = self._dbsnp_support(snapshot.get("support"))
        if support:
            fields["support"] = support
            fields.setdefault("submitter_handle", support[0])

        merged_rsids = self._dbsnp_merged_rsids(payload.get("dbsnp1_merges"))
        if merged_rsids:
            fields["merged_rsids"] = merged_rsids
        citations = payload.get("citations")
        if isinstance(citations, list):
            fields["citation_count"] = len(citations)
        return {key: value for key, value in fields.items() if value not in ("", [], {}, None)}

    def _dbsnp_primary_placement(self, snapshot: dict[str, Any], context: KnowledgeQuery) -> dict[str, Any]:
        placements = snapshot.get("placements_with_allele")
        if not isinstance(placements, list):
            return {}
        target_assembly = "GRCh38" if context.genome_build == "hg38" else "GRCh37"
        first: dict[str, Any] = {}
        ptlp: dict[str, Any] = {}
        for placement in placements:
            if not isinstance(placement, dict):
                continue
            first = first or placement
            assemblies = self._dbsnp_placement_assemblies(placement)
            if any(assembly.startswith(target_assembly) for assembly in assemblies):
                return placement
            annot = placement.get("placement_annot") if isinstance(placement.get("placement_annot"), dict) else {}
            if annot.get("is_ptlp") or placement.get("is_ptlp"):
                ptlp = ptlp or placement
        return ptlp or first

    def _dbsnp_placement_assemblies(self, placement: dict[str, Any]) -> list[str]:
        annot = placement.get("placement_annot") if isinstance(placement.get("placement_annot"), dict) else {}
        traits = annot.get("seq_id_traits_by_assembly")
        assemblies: list[str] = []
        if isinstance(traits, list):
            for trait in traits:
                if isinstance(trait, dict):
                    assemblies.extend(self._clean_string_list(trait.get("assembly_name")))
        assemblies.extend(self._clean_string_list(placement.get("assembly_name")))
        return self._dedupe_clean(assemblies)

    def _dbsnp_placement_fields(self, placement: dict[str, Any]) -> dict[str, Any]:
        if not placement:
            return {}
        fields: dict[str, Any] = {}
        assemblies = self._dbsnp_placement_assemblies(placement)
        if assemblies:
            fields["assembly"] = assemblies[0]
        hgvs: list[str] = []
        spdis: list[str] = []
        alleles: list[str] = []
        placement_alleles = placement.get("alleles")
        if isinstance(placement_alleles, list):
            for item in placement_alleles:
                if not isinstance(item, dict):
                    continue
                hgvs.extend(self._clean_string_list(item.get("hgvs")))
                spdi = self._dbsnp_spdi_dict(item)
                spdi_text = self._dbsnp_spdi_text(spdi)
                if spdi_text:
                    spdis.append(spdi_text)
                allele_text = self._dbsnp_allele_text(spdi)
                if allele_text:
                    alleles.append(allele_text)
        if hgvs:
            fields["hgvs"] = self._dedupe_clean(hgvs)[:6]
        if spdis:
            fields["spdi"] = self._dedupe_clean(spdis)[0]
            fields["spdis"] = self._dedupe_clean(spdis)[:6]
        if alleles:
            fields["alleles"] = "/".join(self._dedupe_clean(alleles)[:4])
        return fields

    def _dbsnp_annotation_fields(self, annotations: Any) -> dict[str, Any]:
        if not isinstance(annotations, list):
            return {}
        genes: list[dict[str, str]] = []
        transcripts: list[dict[str, Any]] = []
        proteins: list[dict[str, Any]] = []
        sequence_ontology: list[str] = []
        hgvs: list[str] = []
        frequencies: list[dict[str, str]] = []
        clinical: list[str] = []

        for annotation in annotations:
            if not isinstance(annotation, dict):
                continue
            frequencies.extend(self._dbsnp_frequency_rows(annotation.get("frequency")))
            clinical.extend(self._dbsnp_clinical_terms(annotation.get("clinical")))
            assembly_annotations = annotation.get("assembly_annotation")
            if not isinstance(assembly_annotations, list):
                continue
            for assembly_annotation in assembly_annotations:
                if not isinstance(assembly_annotation, dict):
                    continue
                sequence_ontology.extend(self._dbsnp_sequence_ontology_terms(assembly_annotation))
                for gene in assembly_annotation.get("genes") or []:
                    if not isinstance(gene, dict):
                        continue
                    gene_record = {
                        "name": _clean_cell(gene.get("name") or gene.get("symbol")),
                        "gene_id": _clean_cell(gene.get("id") or gene.get("gene_id")),
                    }
                    genes.append(gene_record)
                    for rna in gene.get("rnas") or []:
                        if not isinstance(rna, dict):
                            continue
                        rna_hgvs = self._dbsnp_hgvs_values(rna)
                        rna_so = self._dbsnp_sequence_ontology_terms(rna)
                        hgvs.extend(rna_hgvs)
                        sequence_ontology.extend(rna_so)
                        transcript = {
                            "accession": _clean_cell(rna.get("accession_version") or rna.get("id")),
                            "hgvs": rna_hgvs[:3],
                            "sequence_ontology": rna_so,
                        }
                        transcripts.append({key: value for key, value in transcript.items() if value not in ("", [], None)})
                        protein = rna.get("protein") or rna.get("product")
                        if isinstance(protein, dict):
                            protein_hgvs = self._dbsnp_hgvs_values(protein)
                            protein_so = self._dbsnp_sequence_ontology_terms(protein)
                            hgvs.extend(protein_hgvs)
                            sequence_ontology.extend(protein_so)
                            protein_record = {
                                "accession": _clean_cell(protein.get("accession_version") or protein.get("id")),
                                "hgvs": protein_hgvs[:3],
                                "sequence_ontology": protein_so,
                            }
                            proteins.append(
                                {key: value for key, value in protein_record.items() if value not in ("", [], None)}
                            )

        fields: dict[str, Any] = {
            "genes": self._dedupe_gene_records(genes),
            "transcripts": self._dedupe_record_list(transcripts, ("accession", "hgvs")),
            "proteins": self._dedupe_record_list(proteins, ("accession", "hgvs")),
            "sequence_ontology": self._dedupe_clean(sequence_ontology),
            "hgvs": self._dedupe_clean(hgvs)[:8],
            "frequencies": self._dedupe_record_list(frequencies, ("study", "allele_count", "total_count", "frequency")),
            "clinical_significance": self._dedupe_clean(clinical),
        }
        return {key: value for key, value in fields.items() if value not in ("", [], {}, None)}

    def _dbsnp_summary(self, label: str, fields: dict[str, Any]) -> str:
        coordinate = _clean_cell(fields.get("chrpos"))
        if not coordinate and fields.get("chromosome") and fields.get("position"):
            coordinate = f"{fields['chromosome']}:{fields['position']}"
        variant_type = self._dbsnp_variant_type_text(fields.get("variant_type") or fields.get("snp_class"))
        main = variant_type or "dbSNP variant"
        if coordinate:
            main = f"{main} at {coordinate}"
        spdi = _clean_cell(fields.get("spdi"))
        if spdi:
            main = f"{main} ({spdi})"

        gene_text = ", ".join(self._dbsnp_gene_names(fields.get("genes"))[:3])
        consequences = self._dedupe_clean(
            [
                *self._clean_string_list(fields.get("sequence_ontology")),
                *self._clean_string_list(fields.get("functional_classes")),
            ]
        )
        if gene_text and consequences:
            main = f"{main}; {gene_text} {' / '.join(consequences[:4])}"
        elif gene_text:
            main = f"{main}; gene {gene_text}"
        elif consequences:
            main = f"{main}; {' / '.join(consequences[:4])}"

        details: list[str] = []
        hgvs = self._clean_string_list(fields.get("hgvs"))
        if hgvs:
            details.append(f"HGVS {', '.join(hgvs[:3])}")
        frequency_text = self._dbsnp_frequency_summary(fields.get("frequencies") or fields.get("global_mafs"))
        if frequency_text:
            details.append(f"frequency {frequency_text}")
        clinical = self._clean_string_list(fields.get("clinical_significance"))
        if clinical:
            details.append(f"clinical significance: {', '.join(clinical[:3])}")
        validated = self._clean_string_list(fields.get("validated"))
        if validated:
            details.append(f"validated {', '.join(validated[:3])}")
        submitter = _clean_cell(fields.get("submitter_handle"))
        if submitter:
            details.append(f"submitted by {submitter}")
        suffix = f"; {'; '.join(details)}" if details else ""
        return f"{label}: {main}{suffix}."

    def _merge_dbsnp_fields(self, base: dict[str, Any], enrichment: dict[str, Any]) -> dict[str, Any]:
        merged = dict(base)
        for key, value in enrichment.items():
            if value in ("", [], {}, None):
                continue
            if key == "frequencies" and isinstance(value, list) and isinstance(merged.get(key), list):
                merged[key] = self._dedupe_record_list(
                    [*value, *merged[key]],
                    ("study", "allele", "allele_count", "total_count", "frequency"),
                )
            elif isinstance(value, list) and isinstance(merged.get(key), list):
                if value and all(isinstance(item, dict) for item in value):
                    merged[key] = self._dedupe_record_list([*merged[key], *value], ("name", "accession", "study", "hgvs"))
                else:
                    merged[key] = self._dedupe_clean([*merged[key], *value])
            elif key == "genes" and isinstance(value, list):
                merged[key] = self._dedupe_gene_records([*merged.get(key, []), *value])
            elif key == "hgvs" and isinstance(value, list) and isinstance(merged.get(key), list):
                merged[key] = self._dedupe_clean([*value, *merged[key]])[:8]
            else:
                merged[key] = value
        return merged

    def _dbsnp_esummary_genes(self, item: dict[str, Any]) -> list[dict[str, str]]:
        genes: list[dict[str, str]] = []
        raw_genes = item.get("genes")
        if isinstance(raw_genes, list):
            for gene in raw_genes:
                if not isinstance(gene, dict):
                    continue
                genes.append(
                    {
                        "name": _clean_cell(gene.get("name") or gene.get("symbol")),
                        "gene_id": _clean_cell(gene.get("gene_id") or gene.get("id")),
                    }
                )
        return self._dedupe_gene_records(genes)

    def _dbsnp_global_mafs(self, values: Any) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "study": _clean_cell(value.get("study") or value.get("study_name")),
                "frequency": _clean_cell(value.get("freq") or value.get("frequency")),
            }
            if row["study"] or row["frequency"]:
                rows.append({key: item for key, item in row.items() if item})
        return self._dedupe_record_list(rows, ("study", "frequency"))

    def _dbsnp_frequency_rows(self, values: Any) -> list[dict[str, str]]:
        if isinstance(values, dict):
            values = [values]
        if not isinstance(values, list):
            return []
        rows: list[dict[str, str]] = []
        for value in values:
            if not isinstance(value, dict):
                continue
            study_value = value.get("study")
            if isinstance(study_value, dict):
                study = _clean_cell(study_value.get("name") or study_value.get("study_name"))
            else:
                study = _clean_cell(value.get("study_name") or study_value or value.get("source"))
            row = {
                "study": study,
                "allele": _clean_cell(value.get("allele") or value.get("observation")),
                "allele_count": _clean_cell(value.get("allele_count") or value.get("local_obs_count") or value.get("count")),
                "total_count": _clean_cell(
                    value.get("total_count") or value.get("local_allele_count") or value.get("allele_total")
                ),
                "frequency": _clean_cell(value.get("frequency") or value.get("freq") or value.get("allele_frequency")),
            }
            cleaned = {key: item for key, item in row.items() if item}
            if cleaned:
                rows.append(cleaned)
        return self._dedupe_record_list(rows, ("study", "allele", "allele_count", "total_count", "frequency"))

    def _dbsnp_frequency_summary(self, values: Any) -> str:
        rows = self._dbsnp_frequency_rows(values)
        if not rows and isinstance(values, list):
            rows = [row for row in values if isinstance(row, dict)]
        if not rows:
            return ""
        row = rows[0]
        study = _clean_cell(row.get("study"))
        allele_count = _clean_cell(row.get("allele_count"))
        total_count = _clean_cell(row.get("total_count"))
        frequency = _clean_cell(row.get("frequency"))
        allele = _clean_cell(row.get("allele"))
        prefix = f"{study} " if study else ""
        if allele_count and total_count:
            text = f"{prefix}{allele_count}/{total_count}"
        elif frequency:
            text = f"{prefix}{frequency}"
        else:
            text = study
        if allele and text:
            text = f"{text} allele {allele}"
        return text.strip()

    def _dbsnp_clinical_terms(self, values: Any) -> list[str]:
        if isinstance(values, dict):
            values = [values]
        if not isinstance(values, list):
            return []
        terms: list[str] = []
        for value in values:
            if isinstance(value, dict):
                for key in ("clinical_significance", "significance", "review_status", "trait", "disease_names"):
                    terms.extend(self._clean_string_list(value.get(key)))
            else:
                terms.extend(self._clean_string_list(value))
        return self._dedupe_clean(terms)

    def _dbsnp_sequence_ontology_terms(self, value: Any) -> list[str]:
        terms: list[str] = []
        if isinstance(value, list):
            for item in value:
                terms.extend(self._dbsnp_sequence_ontology_terms(item))
        elif isinstance(value, dict):
            name = _clean_cell(value.get("name") or value.get("term") or value.get("so_term"))
            accession = _clean_cell(value.get("accession") or value.get("id"))
            if name and (name.endswith("_variant") or accession.upper().startswith("SO:")):
                terms.append(name)
            for key, child in value.items():
                if key in {"sequence_ontology", "so_terms", "genes", "rnas", "protein", "product", "assembly_annotation"}:
                    terms.extend(self._dbsnp_sequence_ontology_terms(child))
        return self._dedupe_clean(terms)

    def _dbsnp_hgvs_values(self, value: Any) -> list[str]:
        hgvs: list[str] = []
        if isinstance(value, list):
            for item in value:
                hgvs.extend(self._dbsnp_hgvs_values(item))
        elif isinstance(value, dict):
            for key, item in value.items():
                if key.lower().startswith("hgvs"):
                    hgvs.extend(self._clean_string_list(item))
                elif key in {"protein", "product", "rnas", "genes", "assembly_annotation"}:
                    hgvs.extend(self._dbsnp_hgvs_values(item))
        return self._dedupe_clean(hgvs)

    def _dbsnp_spdi_dict(self, item: dict[str, Any]) -> dict[str, Any]:
        allele = item.get("allele") if isinstance(item.get("allele"), dict) else {}
        spdi = allele.get("spdi") if isinstance(allele.get("spdi"), dict) else item.get("spdi")
        return spdi if isinstance(spdi, dict) else {}

    def _dbsnp_spdi_text(self, spdi: dict[str, Any]) -> str:
        if not spdi:
            return ""
        parts = [
            _clean_cell(spdi.get("seq_id")),
            _clean_cell(spdi.get("position")),
            _clean_cell(spdi.get("deleted_sequence")),
            _clean_cell(spdi.get("inserted_sequence")),
        ]
        if not parts[0] or not parts[1]:
            return ""
        return ":".join(parts)

    def _dbsnp_allele_text(self, spdi: dict[str, Any]) -> str:
        deleted = _clean_cell(spdi.get("deleted_sequence"))
        inserted = _clean_cell(spdi.get("inserted_sequence"))
        if deleted and inserted:
            return f"{deleted}>{inserted}"
        if deleted and inserted == "":
            return f"{deleted}>deletion"
        if inserted:
            return f"insertion {inserted}"
        return ""

    def _dbsnp_support(self, value: Any) -> list[str]:
        handles: list[str] = []
        if isinstance(value, list):
            for item in value:
                handles.extend(self._dbsnp_support(item))
        elif isinstance(value, dict):
            for key, item in value.items():
                if "handle" in str(key).lower() or "submitter" in str(key).lower():
                    handles.extend(self._clean_string_list(item))
                elif isinstance(item, (dict, list)):
                    handles.extend(self._dbsnp_support(item))
        return self._dedupe_clean(handles)

    def _dbsnp_merged_rsids(self, values: Any) -> list[str]:
        if not isinstance(values, list):
            return []
        rsids: list[str] = []
        for value in values:
            if isinstance(value, dict):
                rsids.extend(
                    self._dbsnp_rsid(item)
                    for item in (
                        value.get("merged_rsid"),
                        value.get("merged_refsnp_id"),
                        value.get("rsid"),
                        value.get("refsnp_id"),
                    )
                    if _clean_cell(item)
                )
        return self._dedupe_clean(rsids)

    def _parse_dbsnp_docsum(self, docsum: str) -> dict[str, str]:
        if not docsum:
            return {}
        fields: dict[str, str] = {}
        for match in re.finditer(r"(HGVS|SEQ|LEN|GENE)=([^|]+)", docsum):
            fields[match.group(1)] = match.group(2).strip()
        return fields

    def _split_dbsnp_terms(self, value: Any) -> list[str]:
        terms: list[str] = []
        for item in self._clean_string_list(value):
            terms.extend(part.strip() for part in re.split(r"[,;]", item) if part.strip())
        return self._dedupe_clean(terms)

    def _dbsnp_gene_names(self, values: Any) -> list[str]:
        names: list[str] = []
        if isinstance(values, list):
            for value in values:
                if isinstance(value, dict):
                    names.append(_clean_cell(value.get("name") or value.get("symbol")))
                else:
                    names.append(_clean_cell(value))
        return self._dedupe_clean(names)

    def _dbsnp_variant_type_text(self, value: Any) -> str:
        text = _clean_cell(value)
        lookup = {
            "del": "deletion",
            "ins": "insertion",
            "mnv": "multi-nucleotide variant",
            "snv": "SNV",
            "snp": "SNP",
        }
        return lookup.get(text.lower(), text)

    def _dbsnp_rsid(self, value: Any) -> str:
        text = _clean_cell(value)
        if not text:
            return ""
        return text if text.lower().startswith("rs") else f"rs{text}"

    def _clean_string_list(self, values: Any) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str):
            values = re.split(r"[,;]", values)
        elif not isinstance(values, list):
            values = [values]
        return self._dedupe_clean(_clean_cell(value) for value in values)

    def _dedupe_clean(self, values: Any) -> list[str]:
        seen: set[str] = set()
        cleaned: list[str] = []
        for value in values or []:
            text = _clean_cell(value)
            if not text:
                continue
            key = text.casefold()
            if key in seen:
                continue
            seen.add(key)
            cleaned.append(text)
        return cleaned

    def _dedupe_gene_records(self, records: list[dict[str, str]]) -> list[dict[str, str]]:
        return self._dedupe_record_list(
            [{key: value for key, value in record.items() if value} for record in records],
            ("name", "gene_id"),
        )

    def _dedupe_record_list(self, records: list[dict[str, Any]], keys: tuple[str, ...]) -> list[dict[str, Any]]:
        seen: set[tuple[str, ...]] = set()
        deduped: list[dict[str, Any]] = []
        for record in records:
            if not isinstance(record, dict):
                continue
            cleaned = {key: value for key, value in record.items() if value not in ("", [], {}, None)}
            if not cleaned:
                continue
            marker = tuple(_clean_cell(cleaned.get(key)).casefold() for key in keys)
            if marker in seen:
                continue
            seen.add(marker)
            deduped.append(cleaned)
        return deduped

    def _term(self, context: KnowledgeQuery) -> str:
        rsid = _first_rsid(context)
        if self.spec.connector_kind in {"clinvar", "dbsnp"} and rsid:
            return rsid
        if self.spec.connector_kind == "ncbi_gene":
            return f"{context.gene}[sym] AND Homo sapiens[orgn]"
        if self.spec.connector_kind == "litvar":
            return f'("{context.gene}"[Title/Abstract]) AND ({rsid or context.region})'
        return _literature_query(context)

    def _category(self) -> str:
        if self.spec.connector_kind in {"pubmed", "pmc", "litvar", "geo"}:
            return "literature"
        if self.spec.connector_kind == "dbsnp":
            return "population"
        return "clinical_variant"

    def _record_url(self, db: str, item_id: str) -> str:
        if db == "pubmed":
            return f"https://pubmed.ncbi.nlm.nih.gov/{item_id}/"
        if db == "pmc":
            return f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{item_id}/"
        if db == "snp":
            return f"https://www.ncbi.nlm.nih.gov/snp/{item_id}"
        if db == "gene":
            return f"https://www.ncbi.nlm.nih.gov/gene/{item_id}"
        if db == "clinvar":
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{item_id}/"
        return f"https://www.ncbi.nlm.nih.gov/{db}/{item_id}"


class ClinVarConnector(BaseConnector):
    """ClinVar connector using official NCBI query patterns."""

    RETMAX = 10

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        queried_urls: list[str] = []
        warnings: list[str] = []
        seen_ids: set[str] = set()
        records: list[dict[str, Any]] = []
        try:
            for rsid in context.rsids:
                ids = self._esearch_ids(rsid, queried_urls)
                records.extend(self._esummary_records(ids, context, seen_ids, queried_urls))

            gene_ids = self._esearch_ids(f"{context.gene}[gene] AND single_gene[prop]", queried_urls)
            if not gene_ids:
                gene_ids = self._esearch_ids(f"{context.gene}[gene]", queried_urls)
            records.extend(self._esummary_records(gene_ids, context, seen_ids, queried_urls))

            if not records:
                try:
                    records.extend(self._clinical_tables_records(context, seen_ids, queried_urls))
                except KnowledgeRequestError as exc:
                    warnings.append(str(exc))

            return SourceResult(
                source_key=self.spec.key,
                status="ok",
                message=f"Queried ClinVar; {len(records)} clinical variant record(s) returned.",
                records=records,
                warnings=warnings,
                queried_urls=queried_urls,
                elapsed_ms=_elapsed_ms(started),
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=queried_urls, started=started)

    def _esearch_ids(self, term: str, queried_urls: list[str]) -> list[str]:
        url = f"{NCBI_EUTILS_BASE_URL}/esearch.fcgi"
        payload = self.client.get_json(
            url,
            params=_ncbi_request_params(
                {
                    "db": "clinvar",
                    "term": term,
                    "retmode": "json",
                    "retmax": self.RETMAX,
                }
            ),
            rate_limit_per_second=self.spec.rate_limit_per_second,
        )
        queried_urls.append(url)
        idlist = payload.get("esearchresult", {}).get("idlist", [])
        return [str(item_id) for item_id in idlist[: self.RETMAX]]

    def _esummary_records(
        self,
        ids: list[str],
        context: KnowledgeQuery,
        seen_ids: set[str],
        queried_urls: list[str],
    ) -> list[dict[str, Any]]:
        if not ids:
            return []
        url = f"{NCBI_EUTILS_BASE_URL}/esummary.fcgi"
        payload = self.client.get_json(
            url,
            params=_ncbi_request_params(
                {
                    "db": "clinvar",
                    "id": ",".join(ids[: self.RETMAX]),
                    "retmode": "json",
                }
            ),
            rate_limit_per_second=self.spec.rate_limit_per_second,
        )
        queried_urls.append(url)
        result = payload.get("result", {})
        records: list[dict[str, Any]] = []
        for item_id in ids[: self.RETMAX]:
            item = result.get(str(item_id), {})
            if not isinstance(item, dict):
                continue
            source_id = self._variation_id(item, item_id)
            if source_id in seen_ids:
                continue
            seen_ids.add(source_id)
            records.append(self._record_from_summary(item, source_id, context))
        return records

    def _clinical_tables_records(
        self,
        context: KnowledgeQuery,
        seen_ids: set[str],
        queried_urls: list[str],
    ) -> list[dict[str, Any]]:
        terms = context.rsids[0] if context.rsids else context.gene
        payload = self.client.get_json(
            NLM_CLINICAL_TABLES_VARIANTS_URL,
            params={
                "terms": terms,
                "maxList": 5,
                "df": "VariationID,GeneSymbol,dbSNP,Name",
            },
            rate_limit_per_second=self.spec.rate_limit_per_second,
        )
        queried_urls.append(NLM_CLINICAL_TABLES_VARIANTS_URL)
        if not isinstance(payload, list) or len(payload) < 4 or not isinstance(payload[3], list):
            return []
        records: list[dict[str, Any]] = []
        for row in payload[3][:5]:
            if not isinstance(row, list):
                continue
            cells = [_clean_cell(value) for value in row]
            source_id = cells[0] if cells else ""
            if not source_id or source_id in seen_ids:
                continue
            seen_ids.add(source_id)
            gene_symbol = cells[1] if len(cells) > 1 else context.gene
            rsid = cells[2] if len(cells) > 2 else ""
            label = cells[3] if len(cells) > 3 else f"ClinVar variation {source_id}"
            records.append(
                {
                    "category": "clinical_variant",
                    "source": self.spec.name,
                    "label": label or f"ClinVar variation {source_id}",
                    "summary": f"ClinVar Clinical Tables record for {gene_symbol or context.gene}.",
                    "source_id": source_id,
                    "url": self._record_url(source_id),
                    "variant": rsid or source_id,
                    "rsid": rsid,
                    "gene": gene_symbol or context.gene,
                }
            )
        return records

    def _record_from_summary(
        self,
        item: dict[str, Any],
        source_id: str,
        context: KnowledgeQuery,
    ) -> dict[str, Any]:
        title = (
            _clean_cell(item.get("title"))
            or self._first_variation_name(item)
            or _clean_cell(item.get("accession"))
            or f"ClinVar variation {source_id}"
        )
        significance = self._clinical_significance(item)
        traits = self._traits(item)
        rsid = self._rsid(item) or _first_rsid(context)
        clinical_metadata = self._clinical_metadata(item)
        summary_parts = []
        if significance:
            summary_parts.append(f"Clinical significance: {significance}")
        if traits:
            summary_parts.append(f"Phenotype/trait: {'; '.join(traits[:3])}")
        variant_type = _clean_cell(item.get("variant_type") or item.get("obj_type"))
        if variant_type:
            summary_parts.append(f"Variant type: {variant_type}")
        summary = "; ".join(summary_parts) or _clean_cell(item.get("description")) or title
        record = {
            "category": "clinical_variant",
            "source": self.spec.name,
            "label": title,
            "summary": summary,
            "source_id": source_id,
            "url": self._record_url(source_id),
            "variant": rsid or source_id,
            "rsid": rsid,
            "gene": context.gene,
        }
        if significance:
            record["clinical_significance"] = significance
            record["assertion"] = significance
        if traits:
            record["phenotype"] = "; ".join(traits[:5])
        record.update(clinical_metadata)
        return record

    def _clinical_metadata(self, item: dict[str, Any]) -> dict[str, Any]:
        """Preserve ClinVar assertion metadata instead of flattening its label."""
        metadata: dict[str, Any] = {}
        for key in (
            "clinical_significance",
            "germline_classification",
            "somatic_clinical_impact",
            "oncogenicity_classification",
        ):
            value = item.get(key)
            if not isinstance(value, dict):
                continue
            review_status = _clean_cell(value.get("review_status"))
            if review_status:
                metadata["review_status"] = review_status
            assertion_criteria = _clean_cell(
                value.get("assertion_criteria")
                or value.get("assertion_criteria_url")
                or value.get("criteria")
            )
            if assertion_criteria:
                metadata["assertion_criteria"] = assertion_criteria
            last_evaluated = _clean_cell(
                value.get("last_evaluated") or value.get("last_evaluated_date") or value.get("date_last_evaluated")
            )
            if last_evaluated:
                metadata["last_evaluated"] = last_evaluated
            conflicts = _clean_cell(
                value.get("conflicting_interpretations") or value.get("conflict_status")
            )
            if conflicts:
                metadata["conflict_status"] = conflicts
            submitter_count = value.get("number_of_submitters") or value.get("submission_count")
            if submitter_count not in (None, ""):
                metadata["submitter_count"] = submitter_count
            break
        return metadata

    def _variation_id(self, item: dict[str, Any], fallback: str) -> str:
        for key in ("variation_id", "uid", "id"):
            value = _clean_cell(item.get(key))
            if value:
                return value
        return str(fallback)

    def _first_variation_name(self, item: dict[str, Any]) -> str:
        variation_set = item.get("variation_set")
        if not isinstance(variation_set, list):
            return ""
        for variation in variation_set:
            if isinstance(variation, dict):
                name = _clean_cell(variation.get("variation_name"))
                if name:
                    return name
        return ""

    def _clinical_significance(self, item: dict[str, Any]) -> str:
        for key in (
            "clinical_significance",
            "germline_classification",
            "somatic_clinical_impact",
            "oncogenicity_classification",
        ):
            value = item.get(key)
            if isinstance(value, dict):
                cleaned = _clean_cell(value.get("description") or value.get("review_status"))
            else:
                cleaned = _clean_cell(value)
            if cleaned:
                return cleaned
        return ""

    def _traits(self, item: dict[str, Any]) -> list[str]:
        traits: list[str] = []
        for key in ("trait_set", "traits"):
            value = item.get(key)
            if isinstance(value, list):
                for trait in value:
                    cleaned = _clean_cell(trait)
                    if cleaned:
                        traits.append(cleaned)
            else:
                cleaned = _clean_cell(value)
                if cleaned:
                    traits.append(cleaned)
        seen: set[str] = set()
        ordered: list[str] = []
        for trait in traits:
            normalized = trait.lower()
            if normalized in seen:
                continue
            seen.add(normalized)
            ordered.append(trait)
        return ordered

    def _rsid(self, item: dict[str, Any]) -> str:
        for key in ("rsid", "dbsnp", "refsnp_id"):
            value = _clean_cell(item.get(key))
            if value:
                return value if value.lower().startswith("rs") else f"rs{value}"
        variation_set = item.get("variation_set")
        if not isinstance(variation_set, list):
            return ""
        for variation in variation_set:
            if not isinstance(variation, dict):
                continue
            xrefs = variation.get("variation_xrefs")
            if not isinstance(xrefs, list):
                continue
            for xref in xrefs:
                if not isinstance(xref, dict):
                    continue
                source = _clean_cell(xref.get("db_source") or xref.get("db"))
                if source.lower() != "dbsnp":
                    continue
                value = _clean_cell(xref.get("db_id") or xref.get("id"))
                if value:
                    return value if value.lower().startswith("rs") else f"rs{value}"
        return ""

    def _record_url(self, source_id: str) -> str:
        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{source_id}/"


class MedGenConnector(BaseConnector):
    """MedGen connector using official NCBI ESearch and ESummary endpoints."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        queried_urls: list[str] = []
        warnings: list[str] = []
        ordered_ids: list[str] = []
        query_contexts_by_id: dict[str, list[str]] = {}

        try:
            for query_context, term in self._search_terms(context):
                ids = self._esearch_ids(term, queried_urls)
                for item_id in ids:
                    contexts = query_contexts_by_id.setdefault(item_id, [])
                    if query_context not in contexts:
                        contexts.append(query_context)
                    if item_id not in ordered_ids:
                        ordered_ids.append(item_id)

            records = self._esummary_records(
                ordered_ids,
                context,
                query_contexts_by_id,
                queried_urls,
                warnings,
            )
            if records:
                message = f"Queried MedGen; {len(records)} condition or phenotype record(s) returned for {context.gene}."
            else:
                message = f"Queried MedGen; no MedGen records found for {context.gene}."
            return SourceResult(
                source_key=self.spec.key,
                status="ok",
                message=message,
                records=records,
                warnings=warnings,
                queried_urls=queried_urls,
                elapsed_ms=_elapsed_ms(started),
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=queried_urls, started=started)

    def _search_terms(self, context: KnowledgeQuery) -> tuple[tuple[str, str], ...]:
        gene_term = f"{context.gene}[gene]"
        return (
            ("gene", gene_term),
            ("gtr_clinical_tests", f"{MEDGEN_GTR_CLINICAL_FILTER} AND {gene_term}"),
        )

    def _esearch_ids(self, term: str, queried_urls: list[str]) -> list[str]:
        url = f"{NCBI_EUTILS_BASE_URL}/esearch.fcgi"
        payload = self.client.get_json(
            url,
            params=_ncbi_request_params(
                {
                    "db": "medgen",
                    "term": term,
                    "retmode": "json",
                    "retmax": MEDGEN_RETMAX,
                }
            ),
            rate_limit_per_second=self.spec.rate_limit_per_second,
        )
        queried_urls.append(url)
        idlist = payload.get("esearchresult", {}).get("idlist", [])
        return [str(item_id) for item_id in idlist[:MEDGEN_RETMAX]]

    def _esummary_records(
        self,
        ids: list[str],
        context: KnowledgeQuery,
        query_contexts_by_id: dict[str, list[str]],
        queried_urls: list[str],
        warnings: list[str],
    ) -> list[dict[str, Any]]:
        if not ids:
            return []
        url = f"{NCBI_EUTILS_BASE_URL}/esummary.fcgi"
        payload = self.client.get_json(
            url,
            params=_ncbi_request_params(
                {
                    "db": "medgen",
                    "id": ",".join(ids[:MEDGEN_RETMAX]),
                    "retmode": "json",
                }
            ),
            rate_limit_per_second=self.spec.rate_limit_per_second,
        )
        queried_urls.append(url)
        result = payload.get("result", {})
        records: list[dict[str, Any]] = []
        for item_id in ids[:MEDGEN_RETMAX]:
            item = result.get(str(item_id), {})
            if not isinstance(item, dict):
                continue
            records.append(self._record_from_summary(item, item_id, context, query_contexts_by_id, warnings))
        return records

    def _record_from_summary(
        self,
        item: dict[str, Any],
        item_id: str,
        context: KnowledgeQuery,
        query_contexts_by_id: dict[str, list[str]],
        warnings: list[str],
    ) -> dict[str, Any]:
        title = self._summary_value(item, "title", "Title") or f"MedGen record {item_id}"
        concept_id = self._summary_value(item, "conceptid", "concept_id", "ConceptId", "cui")
        definition = self._summary_value(item, "definition", "Definition")
        semantic_id = self._summary_value(item, "semanticid", "SemanticId")
        semantic_type = self._summary_value(item, "semantictype", "SemanticType")
        modification_date = self._summary_value(item, "modificationdate", "ModificationDate")
        record_url = self._record_url(concept_id or item_id)
        meta = self._parse_concept_meta(self._summary_value(item, "conceptmeta", "ConceptMeta"), item_id, warnings)
        summary = definition or f"MedGen {semantic_type or 'medical genetics'} concept associated with {context.gene}: {title}."
        return {
            "category": "clinical_condition",
            "source": self.spec.name,
            "source_key": self.spec.key,
            "source_id": item_id,
            "gene": context.gene,
            "medgen_uid": item_id,
            "concept_id": concept_id,
            "title": title,
            "definition": definition,
            "semantic_id": semantic_id,
            "semantic_type": semantic_type,
            "modification_date": modification_date,
            "summary": summary,
            "label": title,
            "url": record_url,
            "query_contexts": query_contexts_by_id.get(item_id, []),
            "research_links": self._research_links(context, record_url),
            **meta,
        }

    def _parse_concept_meta(self, raw_meta: Any, item_id: str, warnings: list[str]) -> dict[str, Any]:
        text = html.unescape(_clean_cell(raw_meta))
        if not text:
            return {}
        parsed: dict[str, Any] = {}
        parsed["hpo_ids"] = sorted({match.upper() for match in re.findall(r"\bHP:\d{7}\b", text, flags=re.IGNORECASE)})
        parsed["orphanet_ids"] = sorted(
            {
                match.lower().replace(":", "_")
                for match in re.findall(r"\borphanet[:_]\d+\b", text, flags=re.IGNORECASE)
            }
        )
        parsed["omim_ids"] = sorted(
            {
                match
                for match in re.findall(r"\b(?:OMIM|MIM)[:_\s-]*(\d{3,6})\b", text, flags=re.IGNORECASE)
            }
        )
        parsed["related_genes"] = sorted(
            {
                match.strip()
                for match in re.findall(r"<(?:GeneSymbol|Gene|Symbol)>([^<]+)</(?:GeneSymbol|Gene|Symbol)>", text)
                if match.strip()
            }
        )
        lower = text.lower()
        if "clinvar" in lower:
            parsed["has_clinvar"] = True
        if "gtr" in lower or "genetic testing registry" in lower:
            parsed["has_gtr"] = True
        if "<" in text and ">" in text:
            try:
                ElementTree.fromstring(f"<root>{text}</root>")
            except ElementTree.ParseError as exc:
                warnings.append(f"MedGen ConceptMeta for UID {item_id} could not be fully parsed: {exc}")
        return {key: value for key, value in parsed.items() if value not in ("", [], {}, False)}

    def _summary_value(self, item: dict[str, Any], *keys: str) -> str:
        for key in keys:
            if key in item:
                value = _clean_cell(item.get(key))
                if value:
                    return value
        requested = {key.lower() for key in keys}
        for key, value in item.items():
            if str(key).lower() in requested:
                cleaned = _clean_cell(value)
                if cleaned:
                    return cleaned
        return ""

    def _research_links(self, context: KnowledgeQuery, record_url: str) -> list[dict[str, str]]:
        gene_term = f"{context.gene}[gene]"
        gtr_term = f"{MEDGEN_GTR_CLINICAL_FILTER} AND {gene_term}"
        return [
            {
                "label": f"MedGen records for {context.gene}",
                "url": self._search_url(gene_term),
            },
            {
                "label": f"MedGen GTR clinical-test records for {context.gene}",
                "url": self._search_url(gtr_term),
            },
            {
                "label": "MedGen additional descriptions",
                "url": f"{record_url}#Additional_description",
            },
        ]

    def _search_url(self, term: str) -> str:
        return f"https://www.ncbi.nlm.nih.gov/medgen/?term={quote(term)}"

    def _record_url(self, identifier: str) -> str:
        clean_identifier = identifier.strip()
        if not clean_identifier:
            return self.spec.homepage
        return f"https://www.ncbi.nlm.nih.gov/medgen/{quote(clean_identifier)}"


class EnsemblConnector(BaseConnector):
    """Ensembl REST connector for gene and regional variant context."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        if context.genome_build not in {"hg19", "hg38"}:
            return SourceResult(
                self.spec.key,
                "skipped",
                "Ensembl/VEP allele annotation requires a declared hg19 or hg38 assembly; source was not queried.",
                elapsed_ms=_elapsed_ms(started),
            )
        server = "https://rest.ensembl.org" if context.genome_build == "hg38" else "https://grch37.rest.ensembl.org"
        assembly_name = "GRCh38" if context.genome_build == "hg38" else "GRCh37"
        headers = {"Accept": "application/json", "Content-Type": "application/json"}
        urls: list[str] = []
        records: list[dict[str, Any]] = []
        warnings: list[str] = []
        try:
            lookup_url = f"{server}/lookup/symbol/homo_sapiens/{context.gene}"
            gene_payload = self.client.get_json(lookup_url, params={"expand": 1}, headers=headers)
            urls.append(lookup_url)
            records.append(
                {
                    "category": "gene_annotation",
                    "source": self.spec.name,
                    "label": f"Ensembl {context.gene}",
                    "summary": self._gene_summary(context, gene_payload, assembly_name),
                    "source_id": gene_payload.get("id"),
                    "url": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_payload.get('id')}",
                }
            )
            # ``KnowledgeQuery.variants`` is already split to one record per
            # observed ALT. Keep external annotation aligned with that exact
            # identity rather than silently annotating only the first site.
            for variant in context.variants[:25]:
                region = f"{variant.chrom.removeprefix('chr')}:{variant.pos}-{variant.pos}"
                overlap_url = f"{server}/overlap/region/homo_sapiens/{region}"
                overlap = self.client.get_json(overlap_url, params={"feature": "variation"}, headers=headers)
                urls.append(overlap_url)
                vep_payload = self._vep_payload(context, variant, server, headers, urls, warnings)
                for item in self._dedupe_overlap_items(overlap)[:5]:
                    records.append(
                        self._variant_record(
                            context,
                            variant,
                            item,
                            server=server,
                            headers=headers,
                            urls=urls,
                            warnings=warnings,
                            assembly_name=assembly_name,
                            vep_payload=vep_payload,
                        )
                    )
            return SourceResult(
                self.spec.key,
                "ok",
                f"Queried Ensembl; {len(records)} record(s).",
                records,
                warnings=warnings,
                queried_urls=urls,
                elapsed_ms=_elapsed_ms(started),
            )
        except KnowledgeRequestError as exc:
            return SourceResult(self.spec.key, "failed", str(exc), errors=[str(exc)], queried_urls=urls, elapsed_ms=_elapsed_ms(started))

    def _gene_summary(self, context: KnowledgeQuery, gene_payload: dict[str, Any], assembly_name: str) -> str:
        region = f"{gene_payload.get('seq_region_name')}:{gene_payload.get('start')}-{gene_payload.get('end')}"
        details = [
            _clean_cell(gene_payload.get("id")),
            _clean_cell(gene_payload.get("canonical_transcript")),
            _clean_cell(gene_payload.get("biotype")),
        ]
        detail_text = "; ".join(detail for detail in details if detail)
        if detail_text:
            return f"{context.gene} resolved on {assembly_name} to {region} ({detail_text})."
        return f"{context.gene} resolved on {assembly_name} to {region}."

    def _vep_payload(
        self,
        context: KnowledgeQuery,
        variant: Any,
        server: str,
        headers: dict[str, str],
        urls: list[str],
        warnings: list[str],
    ) -> list[dict[str, Any]]:
        if not self._is_simple_snv(variant):
            return []
        chrom = variant.chrom.removeprefix("chr")
        vep_url = f"{server}/vep/homo_sapiens/region/{chrom}:{variant.pos}:{variant.pos}/{quote(variant.alt)}"
        try:
            payload = self.client.get_json(
                vep_url,
                params={"canonical": 1, "mane": 1, "numbers": 1, "variant_class": 1},
                headers=headers,
            )
            urls.append(vep_url)
            return payload if isinstance(payload, list) else []
        except KnowledgeRequestError as exc:
            warnings.append(f"Optional Ensembl VEP annotation failed for {variant.label}; overlap details were still used.")
            return []

    def _variant_record(
        self,
        context: KnowledgeQuery,
        variant: Any,
        overlap_item: dict[str, Any],
        *,
        server: str,
        headers: dict[str, str],
        urls: list[str],
        warnings: list[str],
        assembly_name: str,
        vep_payload: list[dict[str, Any]],
    ) -> dict[str, Any]:
        source_id = _clean_cell(overlap_item.get("id"))
        variation_payload = self._variation_payload(source_id, variant, server, headers, urls, warnings)
        vep_record = self._matching_vep_record(source_id, vep_payload)
        transcript_consequence = self._canonical_transcript_consequence(context, vep_record)
        colocated_variant = self._matching_colocated_variant(source_id, vep_record)

        location = self._variant_location(variation_payload, overlap_item, assembly_name)
        alleles = self._variant_alleles(variation_payload, overlap_item, assembly_name)
        consequence = _clean_cell(
            (vep_record or {}).get("most_severe_consequence")
            or variation_payload.get("most_severe_consequence")
            or overlap_item.get("consequence_type")
        )
        variant_class = _clean_cell(
            variation_payload.get("var_class")
            or (vep_record or {}).get("variant_class")
            or overlap_item.get("feature_type")
        )
        evidence = self._clean_string_list(variation_payload.get("evidence"))
        clinical_significance = self._clean_string_list(
            variation_payload.get("clinical_significance") or overlap_item.get("clinical_significance")
        )
        synonyms = self._clean_string_list(variation_payload.get("synonyms"))
        phenotypes = self._phenotype_summaries(variation_payload.get("phenotypes"))
        label = source_id or variant.label
        record = {
            "category": "variant_annotation",
            "source": self.spec.name,
            "label": label,
            "summary": self._variant_summary(
                context=context,
                variant=variant,
                label=label,
                location=location,
                alleles=alleles,
                consequence=consequence,
                variant_class=variant_class,
                clinical_significance=clinical_significance,
                evidence=evidence,
                source=_clean_cell(overlap_item.get("source") or variation_payload.get("source")),
                transcript_consequence=transcript_consequence,
                phenotypes=phenotypes,
                assembly_name=assembly_name,
            ),
            "source_id": source_id,
            "url": f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={quote(label)}",
            "variant": variant.label,
            "rsid": source_id if source_id.lower().startswith("rs") else "",
            "chromosome": variant.chrom.removeprefix("chr"),
            "position": variant.pos,
            "ref": _clean_cell(getattr(variant, "ref", "")).upper(),
            "alt": _clean_cell(getattr(variant, "alt", "")).upper(),
            "genome_build": context.genome_build,
            "reference_allele_validation": self._reference_allele_validation(
                variant,
                variation_payload,
                overlap_item,
                assembly_name,
            ),
            "location": location,
            "alleles": alleles,
            "consequence": consequence,
            "variant_class": variant_class,
            "clinical_significance": clinical_significance,
            "evidence": evidence,
            "synonyms": synonyms,
            "phenotypes": phenotypes,
            "minor_allele": _clean_cell(variation_payload.get("minor_allele")),
            "maf": variation_payload.get("MAF"),
            "transcript_consequence": transcript_consequence,
            "colocated_variant": colocated_variant,
        }
        return {key: value for key, value in record.items() if value not in ("", [], None)}

    def _variation_payload(
        self,
        source_id: str,
        variant: Any,
        server: str,
        headers: dict[str, str],
        urls: list[str],
        warnings: list[str],
    ) -> dict[str, Any]:
        if not source_id:
            return {}
        variation_url = f"{server}/variation/homo_sapiens/{quote(source_id)}"
        try:
            payload = self.client.get_json(variation_url, params={"phenotypes": 1}, headers=headers)
            urls.append(variation_url)
            return payload if isinstance(payload, dict) else {}
        except KnowledgeRequestError:
            warnings.append(f"Optional Ensembl variation detail failed for {source_id}; overlap details were still used.")
            return {}

    def _dedupe_overlap_items(self, overlap: Any) -> list[dict[str, Any]]:
        if not isinstance(overlap, list):
            return []
        seen: set[tuple[str, str, str, str]] = set()
        deduped: list[dict[str, Any]] = []
        for item in overlap:
            if not isinstance(item, dict):
                continue
            key = (
                _clean_cell(item.get("id")),
                _clean_cell(item.get("seq_region_name")),
                _clean_cell(item.get("start")),
                _clean_cell(item.get("end")),
            )
            if key in seen:
                continue
            seen.add(key)
            deduped.append(item)
        return deduped

    def _is_simple_snv(self, variant: Any) -> bool:
        return (
            len(str(getattr(variant, "ref", "") or "")) == 1
            and len(str(getattr(variant, "alt", "") or "")) == 1
            and str(getattr(variant, "ref", "")).upper() in {"A", "C", "G", "T"}
            and str(getattr(variant, "alt", "")).upper() in {"A", "C", "G", "T"}
        )

    def _variant_location(self, variation_payload: dict[str, Any], overlap_item: dict[str, Any], assembly_name: str) -> str:
        for mapping in variation_payload.get("mappings") or []:
            if not isinstance(mapping, dict):
                continue
            if _clean_cell(mapping.get("assembly_name")) == assembly_name and _clean_cell(mapping.get("location")):
                return f"{assembly_name} {_clean_cell(mapping.get('location'))}"
        seq_region = _clean_cell(overlap_item.get("seq_region_name"))
        start = _clean_cell(overlap_item.get("start"))
        end = _clean_cell(overlap_item.get("end"))
        if seq_region and start and end:
            return f"{assembly_name} {seq_region}:{start}-{end}"
        return ""

    def _variant_alleles(self, variation_payload: dict[str, Any], overlap_item: dict[str, Any], assembly_name: str) -> str:
        for mapping in variation_payload.get("mappings") or []:
            if not isinstance(mapping, dict):
                continue
            if _clean_cell(mapping.get("assembly_name")) == assembly_name and _clean_cell(mapping.get("allele_string")):
                return _clean_cell(mapping.get("allele_string"))
        alleles = self._clean_string_list(overlap_item.get("alleles"))
        return "/".join(alleles)

    def _matching_vep_record(self, source_id: str, vep_payload: list[dict[str, Any]]) -> dict[str, Any]:
        for record in vep_payload:
            if not isinstance(record, dict):
                continue
            colocated = record.get("colocated_variants") or []
            if source_id and any(_clean_cell(item.get("id")) == source_id for item in colocated if isinstance(item, dict)):
                return record
        for record in vep_payload:
            if isinstance(record, dict):
                return record
        return {}

    def _matching_colocated_variant(self, source_id: str, vep_record: dict[str, Any]) -> dict[str, Any]:
        for item in vep_record.get("colocated_variants") or []:
            if not isinstance(item, dict):
                continue
            if source_id and _clean_cell(item.get("id")) == source_id:
                return item
        colocated = vep_record.get("colocated_variants") or []
        return colocated[0] if colocated and isinstance(colocated[0], dict) else {}

    def _canonical_transcript_consequence(self, context: KnowledgeQuery, vep_record: dict[str, Any]) -> dict[str, Any]:
        consequences = [item for item in vep_record.get("transcript_consequences") or [] if isinstance(item, dict)]
        if not consequences:
            return {}
        gene = context.gene.casefold()
        for item in consequences:
            if (
                _clean_cell(item.get("gene_symbol")).casefold() == gene
                and (_clean_cell(item.get("mane_select")) or _clean_cell(item.get("mane_plus_clinical")))
            ):
                return self._annotate_transcript_selection(item)
        for item in consequences:
            if item.get("canonical") == 1 and _clean_cell(item.get("gene_symbol")).casefold() == gene:
                return self._annotate_transcript_selection(item)
        for item in consequences:
            if _clean_cell(item.get("gene_symbol")).casefold() == gene:
                return self._annotate_transcript_selection(item)
        return self._annotate_transcript_selection(consequences[0])

    def _annotate_transcript_selection(self, item: dict[str, Any]) -> dict[str, Any]:
        """Expose VEP's MANE/canonical flags in report-neutral fields."""
        selected = dict(item)
        selected["is_mane_select"] = bool(_clean_cell(item.get("mane_select")))
        selected["is_mane_plus_clinical"] = bool(_clean_cell(item.get("mane_plus_clinical")))
        selected["is_canonical"] = item.get("canonical") == 1 or bool(item.get("canonical"))
        return selected

    def _reference_allele_validation(
        self,
        variant: Any,
        variation_payload: dict[str, Any],
        overlap_item: dict[str, Any],
        assembly_name: str,
    ) -> str:
        """Validate REF/ALT against the assembly-specific Ensembl mapping."""
        ref = _clean_cell(getattr(variant, "ref", "")).upper()
        alt = _clean_cell(getattr(variant, "alt", "")).upper()
        if not ref or not alt:
            return "not_assessed"
        allele_strings: list[str] = []
        for mapping in variation_payload.get("mappings") or []:
            if not isinstance(mapping, dict):
                continue
            if _clean_cell(mapping.get("assembly_name")) == assembly_name:
                allele_strings.append(_clean_cell(mapping.get("allele_string")))
        if not allele_strings:
            alleles = self._clean_string_list(overlap_item.get("alleles"))
            if alleles:
                allele_strings.append("/".join(alleles))
        expected = {ref, alt}
        for allele_string in allele_strings:
            observed = {part.strip().upper() for part in allele_string.replace("|", "/").split("/") if part.strip()}
            if expected <= observed:
                return "validated"
        return "mismatch_or_not_returned"

    def _phenotype_summaries(self, phenotypes: Any) -> list[dict[str, str]]:
        summaries: list[dict[str, str]] = []
        if not isinstance(phenotypes, list):
            return summaries
        for item in phenotypes[:3]:
            if not isinstance(item, dict):
                continue
            summary = {
                "trait": _clean_cell(item.get("trait")),
                "source": _clean_cell(item.get("source")),
                "pvalue": _clean_cell(item.get("pvalue")),
                "risk_allele": _clean_cell(item.get("risk_allele")),
                "study": _clean_cell(item.get("study")),
            }
            summaries.append({key: value for key, value in summary.items() if value})
        return summaries

    def _clean_string_list(self, values: Any) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str):
            return [values] if values.strip() else []
        if not isinstance(values, list):
            return [_clean_cell(values)] if _clean_cell(values) else []
        cleaned: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in cleaned:
                cleaned.append(text)
        return cleaned

    def _variant_summary(
        self,
        *,
        context: KnowledgeQuery,
        variant: Any,
        label: str,
        location: str,
        alleles: str,
        consequence: str,
        variant_class: str,
        clinical_significance: list[str],
        evidence: list[str],
        source: str,
        transcript_consequence: dict[str, Any],
        phenotypes: list[dict[str, str]],
        assembly_name: str,
    ) -> str:
        identity = label
        if alleles and alleles not in identity:
            identity = f"{identity} ({alleles})"
        if location:
            identity = f"{identity} at {location}"

        main = consequence or "Ensembl variation overlap"
        transcript_text = self._transcript_summary(context, transcript_consequence)
        if transcript_text:
            main = f"{main} {transcript_text}"

        details: list[str] = []
        if variant_class:
            details.append(f"variant class {variant_class}")
        if clinical_significance:
            details.append(f"clinical significance: {', '.join(clinical_significance[:3])}")
        if evidence:
            details.append(f"evidence: {', '.join(evidence[:4])}")
        if phenotypes:
            phenotype_text = self._phenotype_summary_text(phenotypes)
            if phenotype_text:
                details.append(f"phenotype annotations: {phenotype_text}")
        if source:
            details.append(f"source {source}")
        details.append(f"input variant {variant.label}")

        suffix = f"; {'; '.join(details)}" if details else ""
        return f"{identity}: {main}{suffix}."

    def _transcript_summary(self, context: KnowledgeQuery, transcript_consequence: dict[str, Any]) -> str:
        if not transcript_consequence:
            return ""
        gene_symbol = _clean_cell(transcript_consequence.get("gene_symbol")) or context.gene
        transcript_id = _clean_cell(transcript_consequence.get("transcript_id"))
        canonical = " canonical" if transcript_consequence.get("canonical") == 1 else ""
        parts = [f"in {gene_symbol}{canonical} transcript {transcript_id}".strip()]
        exon = _clean_cell(transcript_consequence.get("exon"))
        if exon:
            parts.append(f"exon {exon}")
        cdna_start = _clean_cell(transcript_consequence.get("cdna_start"))
        cdna_end = _clean_cell(transcript_consequence.get("cdna_end"))
        if cdna_start and cdna_end:
            cdna = cdna_start if cdna_start == cdna_end else f"{cdna_start}-{cdna_end}"
            parts.append(f"cDNA position {cdna}")
        impact = _clean_cell(transcript_consequence.get("impact"))
        if impact:
            parts.append(f"impact {impact}")
        return ", ".join(part for part in parts if part)

    def _phenotype_summary_text(self, phenotypes: list[dict[str, str]]) -> str:
        entries: list[str] = []
        for item in phenotypes[:2]:
            text = item.get("trait", "")
            if item.get("pvalue"):
                text = f"{text} p={item['pvalue']}"
            if item.get("source"):
                text = f"{text} ({item['source']})"
            if text:
                entries.append(text)
        return "; ".join(entries)


class UcscConnector(BaseConnector):
    """UCSC Genome Browser API connector for compact track context."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        genome = self._genome(context)
        headers = {"Accept": "application/json"}
        urls: list[str] = []
        records: list[dict[str, Any]] = []
        warnings: list[str] = []

        variant = _first_variant(context)
        region = self._parse_region(context.region)
        search_payload: dict[str, Any] = {}
        if not region:
            search_payload = self._get_json(
                f"{UCSC_API_BASE}/search",
                {"search": context.gene, "genome": genome},
                headers,
                urls,
                warnings,
                "search",
            )
            region = self._region_from_search(search_payload, context.gene)
        if not region:
            region = self._variant_window(variant)
        if not region:
            return SourceResult(
                self.spec.key,
                "skipped",
                "No genome region or variant coordinate was available for UCSC Genome Browser lookup.",
                warnings=warnings,
                queried_urls=urls,
                elapsed_ms=_elapsed_ms(started),
            )

        variant_interval = self._variant_interval(variant, region)
        track_url = f"{UCSC_API_BASE}/getData/track"
        sequence_url = f"{UCSC_API_BASE}/getData/sequence"

        refseq_payload = self._track_payload(track_url, genome, "ncbiRefSeq", region, headers, urls, warnings)
        gene_models = self._track_items(refseq_payload, "ncbiRefSeq")
        if gene_models:
            records.append(self._gene_model_record(context, genome, region, gene_models[0]))
        elif search_payload:
            records.extend(self._search_records(context, genome, search_payload))

        if variant_interval:
            sequence_payload = self._get_json(
                sequence_url,
                {
                    "genome": genome,
                    "chrom": variant_interval["chrom"],
                    "start": variant_interval["start0"],
                    "end": variant_interval["end1"],
                },
                headers,
                urls,
                warnings,
                "reference sequence",
            )
            sequence_record = self._sequence_record(context, genome, variant, variant_interval, sequence_payload)
            if sequence_record:
                records.append(sequence_record)

        cpg_payload = self._track_payload(track_url, genome, "cpgIslandExt", region, headers, urls, warnings)
        cpg_items = self._track_items(cpg_payload, "cpgIslandExt")
        if cpg_items:
            records.append(self._cpg_record(context, genome, region, cpg_items[0]))

        ccre_payload = self._track_payload(track_url, genome, "encodeCcreCombined", region, headers, urls, warnings)
        ccre_items = self._track_items(ccre_payload, "encodeCcreCombined")
        if ccre_items:
            records.append(self._ccre_record(context, genome, region, ccre_items, bool(ccre_payload.get("maxItemsLimit"))))

        tfbs_payload = self._track_payload(track_url, genome, "encRegTfbsClustered", region, headers, urls, warnings)
        tfbs_items = self._track_items(tfbs_payload, "encRegTfbsClustered")
        if tfbs_items:
            records.append(self._tfbs_record(context, genome, region, tfbs_items, bool(tfbs_payload.get("maxItemsLimit"))))

        repeat_payload = self._track_payload(track_url, genome, "rmsk", region, headers, urls, warnings)
        repeat_items = self._track_items(repeat_payload, "rmsk")
        if repeat_items:
            records.append(self._repeat_record(context, genome, region, repeat_items, bool(repeat_payload.get("maxItemsLimit"))))

        snp_items: list[dict[str, Any]] = []
        snp_scope = "query window"
        if variant_interval:
            snp_exact_payload = self._track_payload(
                track_url,
                genome,
                "snp151Common",
                variant_interval,
                headers,
                urls,
                warnings,
                warn_label="snp151Common exact-variant track",
            )
            snp_items = self._track_items(snp_exact_payload, "snp151Common")
            snp_scope = "exact variant coordinate"
        if not snp_items:
            snp_region_payload = self._track_payload(track_url, genome, "snp151Common", region, headers, urls, warnings)
            snp_items = self._track_items(snp_region_payload, "snp151Common")
            snp_scope = "query window"
        if snp_items:
            records.append(self._snp_record(context, genome, region, variant, snp_items, snp_scope))

        if not records and search_payload:
            records.extend(self._search_records(context, genome, search_payload))

        status = "failed" if not records and warnings else "ok"
        if not records and not warnings:
            message = f"Queried UCSC Genome Browser API tracks for {self._region_text(genome, region)}, but no compact annotations were returned."
        elif status == "failed":
            message = "UCSC Genome Browser lookup failed; no track annotations were returned."
        else:
            message = (
                f"Queried UCSC Genome Browser API tracks for {self._region_text(genome, region)}; "
                f"{len(records)} compact annotation record(s) returned."
            )
        return SourceResult(
            self.spec.key,
            status,
            message,
            records,
            warnings=warnings,
            errors=warnings if status == "failed" else [],
            queried_urls=urls,
            elapsed_ms=_elapsed_ms(started),
        )

    def _genome(self, context: KnowledgeQuery) -> str:
        build = (context.genome_build or "").casefold()
        return "hg38" if "38" in build else "hg19"

    def _get_json(
        self,
        url: str,
        params: dict[str, Any],
        headers: dict[str, str],
        urls: list[str],
        warnings: list[str],
        label: str,
    ) -> dict[str, Any]:
        urls.append(url)
        try:
            payload = self.client.get_json(
                url,
                params=params,
                headers=headers,
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError:
            warnings.append(f"Optional UCSC {label} lookup failed; other UCSC annotations were still used.")
            return {}
        if not isinstance(payload, dict):
            warnings.append(f"Optional UCSC {label} lookup returned no usable data; other UCSC annotations were still used.")
            return {}
        if payload.get("error"):
            warnings.append(f"Optional UCSC {label} lookup failed; other UCSC annotations were still used.")
            return {}
        return payload

    def _track_payload(
        self,
        url: str,
        genome: str,
        track: str,
        interval: dict[str, Any],
        headers: dict[str, str],
        urls: list[str],
        warnings: list[str],
        *,
        warn_label: str = "",
    ) -> dict[str, Any]:
        return self._get_json(
            url,
            {
                "genome": genome,
                "track": track,
                "chrom": interval["chrom"],
                "start": interval["start0"],
                "end": interval["end1"],
                "maxItemsOutput": UCSC_TRACK_MAX_ITEMS,
            },
            headers,
            urls,
            warnings,
            warn_label or f"{track} track",
        )

    def _track_items(self, payload: dict[str, Any], track: str) -> list[dict[str, Any]]:
        items = payload.get(track) or payload.get("items") or []
        if not isinstance(items, list):
            return []
        deduped: list[dict[str, Any]] = []
        seen: set[tuple[str, str, str]] = set()
        for item in items:
            if not isinstance(item, dict):
                continue
            key = (
                _clean_cell(item.get("name") or item.get("chrom")),
                _clean_cell(item.get("chromStart") or item.get("txStart") or item.get("start")),
                _clean_cell(item.get("chromEnd") or item.get("txEnd") or item.get("end")),
            )
            if key in seen:
                continue
            seen.add(key)
            deduped.append(item)
        return deduped

    def _parse_region(self, value: Any) -> dict[str, Any]:
        text = _clean_cell(value).replace(",", "")
        match = re.search(r"(?P<chrom>(?:chr)?[A-Za-z0-9_.]+)\s*:\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)", text)
        if not match:
            return {}
        start1 = self._as_int(match.group("start"))
        end1 = self._as_int(match.group("end"))
        if start1 is None or end1 is None:
            return {}
        if end1 < start1:
            start1, end1 = end1, start1
        return self._interval(match.group("chrom"), start1, end1)

    def _region_from_search(self, payload: dict[str, Any], gene: str) -> dict[str, Any]:
        gene_lower = gene.casefold()
        fallback: dict[str, Any] = {}
        for group in payload.get("positionMatches") or []:
            if not isinstance(group, dict):
                continue
            for match in group.get("matches") or []:
                if not isinstance(match, dict):
                    continue
                interval = self._parse_region(match.get("position"))
                if not interval:
                    continue
                if not fallback:
                    fallback = interval
                name = _clean_cell(match.get("posName") or match.get("name")).casefold()
                if name == gene_lower or gene_lower in name:
                    return interval
        return fallback

    def _variant_window(self, variant: Any) -> dict[str, Any]:
        if variant is None:
            return {}
        chrom = _clean_cell(getattr(variant, "chrom", ""))
        pos = self._as_int(getattr(variant, "pos", None))
        if not (chrom and pos):
            return {}
        return self._interval(chrom, max(1, pos - 250), pos + 250)

    def _variant_interval(self, variant: Any, fallback_region: dict[str, Any]) -> dict[str, Any]:
        if variant is None:
            return {}
        chrom = _clean_cell(getattr(variant, "chrom", "")) or _clean_cell(fallback_region.get("chrom"))
        pos = self._as_int(getattr(variant, "pos", None))
        if not (chrom and pos):
            return {}
        return self._interval(chrom, pos, pos)

    def _interval(self, chrom: str, start1: int, end1: int) -> dict[str, Any]:
        return {
            "chrom": self._ucsc_chrom(chrom),
            "start1": max(1, start1),
            "end1": max(1, end1),
            "start0": max(0, start1 - 1),
        }

    def _ucsc_chrom(self, chrom: str) -> str:
        text = _clean_cell(chrom)
        if not text:
            return ""
        if text.lower().startswith("chr"):
            return f"chr{text[3:]}"
        if text in {"M", "MT"}:
            return "chrM"
        return f"chr{text}"

    def _gene_model_record(
        self,
        context: KnowledgeQuery,
        genome: str,
        region: dict[str, Any],
        item: dict[str, Any],
    ) -> dict[str, Any]:
        transcript = _clean_cell(item.get("name"))
        gene = _clean_cell(item.get("name2")) or context.gene
        tx_interval = self._item_interval(item, "txStart", "txEnd")
        cds_interval = self._item_interval(item, "cdsStart", "cdsEnd")
        exons = self._exon_intervals(item)
        record = {
            "category": "gene_model",
            "source": self.spec.name,
            "label": f"UCSC ncbiRefSeq {transcript or gene}",
            "summary": self._gene_model_summary(context, genome, region, item, tx_interval, cds_interval, exons),
            "source_id": transcript,
            "url": self._browser_url(genome, region),
            "track": "ncbiRefSeq",
            "genome": genome,
            "gene": gene,
            "transcript": transcript,
            "strand": _clean_cell(item.get("strand")),
            "chromosome": _clean_cell(item.get("chrom")),
            "transcript_interval": tx_interval,
            "cds_interval": cds_interval,
            "exon_count": self._as_int(item.get("exonCount")),
            "exons": exons,
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _gene_model_summary(
        self,
        context: KnowledgeQuery,
        genome: str,
        region: dict[str, Any],
        item: dict[str, Any],
        tx_interval: str,
        cds_interval: str,
        exons: list[str],
    ) -> str:
        transcript = _clean_cell(item.get("name")) or "RefSeq transcript"
        gene = _clean_cell(item.get("name2")) or context.gene
        strand = _clean_cell(item.get("strand"))
        exon_count = self._as_int(item.get("exonCount")) or len(exons)
        parts = [f"{exon_count} exons" if exon_count else ""]
        if cds_interval:
            parts.append(f"CDS {cds_interval}")
        cds_status = self._cds_status(item)
        if cds_status:
            parts.append(cds_status)
        parts.append(f"query window {self._region_text(genome, region)}")
        detail = "; ".join(part for part in parts if part)
        strand_text = f" ({strand})" if strand else ""
        return f"UCSC ncbiRefSeq {transcript} for {gene} on {genome} {tx_interval}{strand_text}: {detail}."

    def _sequence_record(
        self,
        context: KnowledgeQuery,
        genome: str,
        variant: Any,
        interval: dict[str, Any],
        payload: dict[str, Any],
    ) -> dict[str, Any]:
        dna = _clean_cell(payload.get("dna")).upper()
        if not dna:
            return {}
        label = _clean_cell(getattr(variant, "label", "")) or self._region_text(genome, interval)
        allele_text = self._variant_allele_text(variant)
        allele_clause = f"; sample allele {allele_text} can be compared against this UCSC assembly coordinate" if allele_text else ""
        position_text = f"{interval['chrom']}:{interval['start1']}"
        return {
            "category": "reference_sequence",
            "source": self.spec.name,
            "label": f"UCSC reference base {position_text}",
            "summary": f"UCSC {genome} reference base at {label} / {position_text} is {dna}{allele_clause}.",
            "source_id": f"{genome}:{position_text}",
            "url": self._browser_url(genome, interval),
            "genome": genome,
            "chromosome": interval["chrom"],
            "position": interval["start1"],
            "reference_base": dna,
            "variant": label,
            "gene": context.gene,
        }

    def _cpg_record(self, context: KnowledgeQuery, genome: str, region: dict[str, Any], item: dict[str, Any]) -> dict[str, Any]:
        interval = self._item_interval(item, "chromStart", "chromEnd")
        label = _clean_cell(item.get("name")) or "CpG island"
        length = self._format_int(item.get("length"))
        cpg_num = self._format_int(item.get("cpgNum"))
        gc_percent = self._format_float(item.get("perGc"))
        obs_exp = self._format_float(item.get("obsExp"))
        details = []
        if length:
            details.append(f"length {length} bp")
        if cpg_num:
            details.append(f"CpG count {cpg_num}")
        if gc_percent:
            details.append(f"GC {gc_percent}%")
        if obs_exp:
            details.append(f"observed/expected CpG {obs_exp}")
        summary = f"UCSC CpG island {label} overlaps the {context.gene} query window at {genome} {interval}"
        if details:
            summary = f"{summary}: {', '.join(details)}"
        record = {
            "category": "cpg_island",
            "source": self.spec.name,
            "label": label,
            "summary": f"{summary}.",
            "source_id": label,
            "url": self._browser_url(genome, region),
            "track": "cpgIslandExt",
            "genome": genome,
            "location": f"{genome} {interval}",
            "length": item.get("length"),
            "cpg_count": item.get("cpgNum"),
            "gc_percent": item.get("perGc"),
            "observed_expected": item.get("obsExp"),
            "gene": context.gene,
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _ccre_record(
        self,
        context: KnowledgeQuery,
        genome: str,
        region: dict[str, Any],
        items: list[dict[str, Any]],
        capped: bool,
    ) -> dict[str, Any]:
        first = items[0]
        first_id = _clean_cell(first.get("name"))
        first_desc = _clean_cell(first.get("description")) or _clean_cell(first.get("ccre"))
        first_interval = self._item_interval(first, "chromStart", "chromEnd")
        labels = self._dedupe_clean(
            _clean_cell(item.get("encodeLabel") or item.get("ucscLabel") or item.get("ccre"))
            for item in items
        )
        z_score = self._format_float(first.get("zScore"))
        first_parts = [first_id, first_desc]
        first_text = " ".join(part for part in first_parts if part)
        label_text = _clean_cell(first.get("ccre"))
        if label_text:
            first_text = f"{first_text} ({label_text}"
            if z_score:
                first_text = f"{first_text}, z={z_score}"
            first_text = f"{first_text})"
        count_text = self._count_phrase(len(items), "regulatory element")
        cap_text = ", capped by maxItemsOutput" if capped else ""
        summary = f"UCSC ENCODE cCRE returned {count_text} in the {context.gene} window{cap_text}; first {first_text} at {genome} {first_interval}"
        if labels:
            summary = f"{summary}; labels include {', '.join(labels[:4])}"
        record = {
            "category": "regulatory_element",
            "source": self.spec.name,
            "label": f"UCSC ENCODE cCRE {first_id or context.gene}",
            "summary": f"{summary}.",
            "source_id": first_id,
            "url": self._browser_url(genome, region),
            "track": "encodeCcreCombined",
            "genome": genome,
            "gene": context.gene,
            "count": len(items),
            "capped": capped,
            "first_element": self._compact_items(items[:1])[0] if items else {},
            "labels": labels,
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _tfbs_record(
        self,
        context: KnowledgeQuery,
        genome: str,
        region: dict[str, Any],
        items: list[dict[str, Any]],
        capped: bool,
    ) -> dict[str, Any]:
        ranked = sorted(items, key=lambda item: self._as_int(item.get("score")) or 0, reverse=True)
        signals = []
        for item in ranked[:4]:
            name = _clean_cell(item.get("name"))
            score = _clean_cell(item.get("score"))
            sources = self._format_int(item.get("sourceCount"))
            if not name:
                continue
            text = f"{name} score {score}" if score else name
            if sources:
                text = f"{text} ({sources} sources)"
            signals.append(text)
        count_text = self._count_phrase(len(items), "transcription-factor binding cluster")
        cap_text = ", capped by maxItemsOutput" if capped else ""
        summary = f"UCSC ENCODE TFBS clustered track returned {count_text} in the {context.gene} window{cap_text}"
        if signals:
            summary = f"{summary}; top signals: {', '.join(signals)}"
        record = {
            "category": "transcription_factor_binding",
            "source": self.spec.name,
            "label": f"UCSC ENCODE TFBS clusters for {context.gene}",
            "summary": f"{summary}.",
            "source_id": "encRegTfbsClustered",
            "url": self._browser_url(genome, region),
            "track": "encRegTfbsClustered",
            "genome": genome,
            "gene": context.gene,
            "count": len(items),
            "capped": capped,
            "top_signals": self._compact_items(ranked[:4]),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _repeat_record(
        self,
        context: KnowledgeQuery,
        genome: str,
        region: dict[str, Any],
        items: list[dict[str, Any]],
        capped: bool,
    ) -> dict[str, Any]:
        examples = []
        for item in items[:3]:
            name = _clean_cell(item.get("repName") or item.get("name"))
            repeat_class = _clean_cell(item.get("repClass"))
            family = _clean_cell(item.get("repFamily"))
            interval = self._item_interval(item, "chromStart", "chromEnd")
            label = name
            class_text = "/".join(part for part in (repeat_class, family) if part)
            if class_text:
                label = f"{label} {class_text}" if label else class_text
            if label and interval:
                examples.append(f"{label} at {genome} {interval}")
        count_text = self._count_phrase(len(items), "repeat/low-complexity annotation")
        cap_text = ", capped by maxItemsOutput" if capped else ""
        summary = f"UCSC RepeatMasker returned {count_text} in the {context.gene} window{cap_text}"
        if examples:
            summary = f"{summary}; examples: {'; '.join(examples)}"
        record = {
            "category": "repeat_annotation",
            "source": self.spec.name,
            "label": f"UCSC RepeatMasker annotations for {context.gene}",
            "summary": f"{summary}.",
            "source_id": "rmsk",
            "url": self._browser_url(genome, region),
            "track": "rmsk",
            "genome": genome,
            "gene": context.gene,
            "count": len(items),
            "capped": capped,
            "examples": self._compact_items(items[:3]),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _snp_record(
        self,
        context: KnowledgeQuery,
        genome: str,
        region: dict[str, Any],
        variant: Any,
        items: list[dict[str, Any]],
        scope: str,
    ) -> dict[str, Any]:
        examples = []
        for item in items[:3]:
            name = _clean_cell(item.get("name"))
            observed = _clean_cell(item.get("observed") or item.get("alleles"))
            snp_class = _clean_cell(item.get("class"))
            valid = _clean_cell(item.get("valid")).replace(",", "/")
            interval = self._item_interval(item, "chromStart", "chromEnd")
            parts = [name, f"at {genome} {interval}" if interval else "", observed, snp_class]
            if valid:
                parts.append(f"validated {valid}")
            examples.append(" ".join(part for part in parts if part))
        variant_label = _clean_cell(getattr(variant, "label", "")) if variant is not None else ""
        nearby = f" near {variant_label}" if variant_label and scope != "exact variant coordinate" else ""
        count_text = self._count_phrase(len(items), "common dbSNP variant")
        summary = f"UCSC dbSNP Common (snp151Common) found {count_text} at the {scope}{nearby}"
        if examples:
            summary = f"{summary}; examples: {'; '.join(examples)}"
        record = {
            "category": "common_variant",
            "source": self.spec.name,
            "label": f"UCSC dbSNP Common variants for {context.gene}",
            "summary": f"{summary}.",
            "source_id": "snp151Common",
            "url": self._browser_url(genome, region),
            "track": "snp151Common",
            "genome": genome,
            "gene": context.gene,
            "variant": variant_label,
            "count": len(items),
            "scope": scope,
            "examples": self._compact_items(items[:3]),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _search_records(self, context: KnowledgeQuery, genome: str, payload: dict[str, Any]) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for group in payload.get("positionMatches") or []:
            if not isinstance(group, dict):
                continue
            for match in group.get("matches") or []:
                if not isinstance(match, dict):
                    continue
                position = _clean_cell(match.get("position"))
                if not position:
                    continue
                interval = self._parse_region(position)
                records.append(
                    {
                        "category": "gene_annotation",
                        "source": self.spec.name,
                        "label": _clean_cell(match.get("posName")) or context.gene,
                        "summary": f"UCSC search resolved {context.gene} to {genome} {position}.",
                        "url": self._browser_url(genome, interval) if interval else self.spec.homepage,
                    }
                )
                if len(records) >= 3:
                    return records
        return records

    def _item_interval(self, item: dict[str, Any], start_key: str, end_key: str) -> str:
        chrom = _clean_cell(item.get("chrom"))
        start0 = self._as_int(item.get(start_key))
        end1 = self._as_int(item.get(end_key))
        if not (chrom and start0 is not None and end1 is not None):
            return ""
        return f"{chrom}:{start0 + 1}-{end1}"

    def _exon_intervals(self, item: dict[str, Any]) -> list[str]:
        chrom = _clean_cell(item.get("chrom"))
        starts = self._split_ints(item.get("exonStarts"))
        ends = self._split_ints(item.get("exonEnds"))
        intervals = []
        for start0, end1 in zip(starts, ends):
            if chrom:
                intervals.append(f"{chrom}:{start0 + 1}-{end1}")
        return intervals

    def _cds_status(self, item: dict[str, Any]) -> str:
        start_status = _clean_cell(item.get("cdsStartStat"))
        end_status = _clean_cell(item.get("cdsEndStat"))
        if start_status == "cmpl" and end_status == "cmpl":
            return "complete CDS start/end"
        statuses = [status for status in (start_status, end_status) if status]
        if statuses:
            return f"CDS status {'/'.join(statuses)}"
        return ""

    def _variant_allele_text(self, variant: Any) -> str:
        if variant is None:
            return ""
        ref = _clean_cell(getattr(variant, "ref", ""))
        alt = _clean_cell(getattr(variant, "alt", ""))
        if ref and alt:
            return f"{ref}>{alt}"
        return ""

    def _region_text(self, genome: str, interval: dict[str, Any]) -> str:
        return f"{genome} {interval['chrom']}:{interval['start1']}-{interval['end1']}"

    def _browser_url(self, genome: str, interval: dict[str, Any]) -> str:
        if not interval:
            return self.spec.homepage
        position = f"{interval['chrom']}:{interval['start1']}-{interval['end1']}"
        return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={quote(genome)}&position={quote(position, safe=':-,')}"

    def _compact_items(self, items: list[dict[str, Any]]) -> list[dict[str, Any]]:
        compact: list[dict[str, Any]] = []
        keys = (
            "name",
            "chrom",
            "chromStart",
            "chromEnd",
            "score",
            "sourceCount",
            "ccre",
            "encodeLabel",
            "ucscLabel",
            "repName",
            "repClass",
            "repFamily",
            "observed",
            "class",
            "valid",
        )
        for item in items:
            row = {key: item.get(key) for key in keys if item.get(key) not in ("", [], {}, None)}
            if row:
                compact.append(row)
        return compact

    def _split_ints(self, value: Any) -> list[int]:
        ints: list[int] = []
        for part in _clean_cell(value).split(","):
            number = self._as_int(part)
            if number is not None:
                ints.append(number)
        return ints

    def _dedupe_clean(self, values: Any) -> list[str]:
        cleaned: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in cleaned:
                cleaned.append(text)
        return cleaned

    def _count_phrase(self, count: int, singular: str) -> str:
        noun = singular if count == 1 else f"{singular}(s)"
        return f"{count} {noun}"

    def _as_int(self, value: Any) -> int | None:
        try:
            return int(str(value).replace(",", "").strip())
        except (TypeError, ValueError):
            return None

    def _format_int(self, value: Any) -> str:
        number = self._as_int(value)
        if number is None:
            return _clean_cell(value)
        return f"{number:,}"

    def _format_float(self, value: Any) -> str:
        try:
            return f"{float(value):.3g}"
        except (TypeError, ValueError):
            return _clean_cell(value)


class LiteratureSearchConnector(BaseConnector):
    """Connector for open literature APIs with straightforward search endpoints."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        kind = self.spec.connector_kind
        term = _literature_query(context)
        urls: list[str] = []
        records: list[dict[str, Any]] = []
        try:
            if kind == "europe_pmc":
                url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
                payload = self.client.get_json(url, params={"query": term, "format": "json", "pageSize": 5})
                urls.append(url)
                for item in payload.get("resultList", {}).get("result", [])[:5]:
                    records.append(
                        {
                            "category": "literature",
                            "source": self.spec.name,
                            "label": str(item.get("title") or item.get("id") or term),
                            "summary": str(item.get("abstractText") or item.get("journalTitle") or ""),
                            "source_id": item.get("id"),
                            "url": str(item.get("fullTextUrlList", {}).get("fullTextUrl", [{}])[0].get("url") or f"https://europepmc.org/article/{item.get('source', 'MED')}/{item.get('id', '')}"),
                            "variant": _first_rsid(context),
                        }
                    )
            elif kind == "openalex":
                url = "https://api.openalex.org/works"
                payload = self.client.get_json(url, params={"search": term, "per-page": 5})
                urls.append(url)
                for item in payload.get("results", [])[:5]:
                    records.append(
                        {
                            "category": "literature",
                            "source": self.spec.name,
                            "label": str(item.get("title") or item.get("id") or term),
                            "summary": str((item.get("primary_location") or {}).get("source", {}).get("display_name") or ""),
                            "source_id": item.get("id"),
                            "url": str(item.get("doi") or item.get("id") or ""),
                            "variant": _first_rsid(context),
                        }
                    )
            elif kind == "crossref":
                url = "https://api.crossref.org/works"
                payload = self.client.get_json(url, params={"query": term, "rows": 5})
                urls.append(url)
                for item in payload.get("message", {}).get("items", [])[:5]:
                    title = " ".join(item.get("title", [])[:1])
                    records.append(
                        {
                            "category": "literature",
                            "source": self.spec.name,
                            "label": title or str(item.get("DOI") or term),
                            "summary": str(item.get("container-title", [""])[0] if item.get("container-title") else ""),
                            "source_id": item.get("DOI"),
                            "url": str(item.get("URL") or ""),
                            "variant": _first_rsid(context),
                        }
                    )
            elif kind == "semantic_scholar":
                url = "https://api.semanticscholar.org/graph/v1/paper/search"
                payload = self.client.get_json(url, params={"query": term, "limit": 5, "fields": "title,year,url,abstract"})
                urls.append(url)
                for item in payload.get("data", [])[:5]:
                    records.append(
                        {
                            "category": "literature",
                            "source": self.spec.name,
                            "label": str(item.get("title") or term),
                            "summary": str(item.get("abstract") or ""),
                            "source_id": item.get("paperId"),
                            "url": str(item.get("url") or ""),
                            "variant": _first_rsid(context),
                        }
                    )
            elif kind in {"biorxiv", "medrxiv"}:
                # The public details endpoint is date-based, so v1 records a precise search link.
                records.append(
                    {
                        "category": "literature",
                        "source": self.spec.name,
                        "label": f"{self.spec.name} search for {term}",
                        "summary": "Preprint search linkout; article API is date-window based.",
                        "url": f"https://www.{kind}.org/search/{term.replace(' ', '%2520')}",
                        "variant": _first_rsid(context),
                    }
                )
            return SourceResult(self.spec.key, "ok", f"Queried {self.spec.name}; {len(records)} record(s).", records, queried_urls=urls, elapsed_ms=_elapsed_ms(started))
        except KnowledgeRequestError as exc:
            return SourceResult(self.spec.key, "failed", str(exc), errors=[str(exc)], queried_urls=urls, elapsed_ms=_elapsed_ms(started))


class GwasCatalogConnector(BaseConnector):
    """NHGRI-EBI GWAS Catalog connector."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        rsid = _first_rsid(context)
        records: list[dict[str, Any]] = []
        urls: list[str] = []
        warnings: list[str] = []
        try:
            metadata = self._gwas_metadata(urls, warnings)
            query_context = "gene"
            payload: dict[str, Any] = {}
            rsid_total = 0
            if rsid:
                query_context = "rsid"
                payload = self._gwas_associations({"rs_id": rsid, "size": GWAS_CATALOG_MAX_ASSOCIATIONS}, urls)
                rsid_total = self._gwas_total(payload)
            if not self._gwas_association_rows(payload) and context.gene:
                query_context = "gene_fallback" if rsid else "gene"
                payload = self._gwas_associations(
                    {"mapped_gene": context.gene, "size": GWAS_CATALOG_MAX_ASSOCIATIONS},
                    urls,
                )

            for item in self._gwas_association_rows(payload)[:GWAS_CATALOG_MAX_ASSOCIATIONS]:
                records.append(self._gwas_record(item, context, metadata, query_context, urls, warnings))

            if not records:
                label = f"{rsid} / {context.gene}" if rsid else context.gene
                message = f"Queried GWAS Catalog; no association records found for {label}."
            else:
                total = self._gwas_total(payload)
                message = f"Queried GWAS Catalog; {len(records)} GWAS association record(s) returned"
                if total and total != len(records):
                    message = f"{message} ({total} total available)"
                if query_context == "gene_fallback":
                    message = f"{message} for {context.gene}; no rsID-specific associations were found for {rsid}."
                elif query_context == "rsid":
                    message = f"{message} for {rsid}."
                else:
                    message = f"{message} for {context.gene}."
                if query_context == "gene_fallback" and rsid_total:
                    message = f"{message} rsID query returned {rsid_total} record(s)."
            return SourceResult(
                self.spec.key,
                "ok",
                message,
                records,
                warnings=warnings,
                queried_urls=urls,
                elapsed_ms=_elapsed_ms(started),
            )
        except KnowledgeRequestError as exc:
            return SourceResult(self.spec.key, "failed", str(exc), errors=[str(exc)], queried_urls=urls, elapsed_ms=_elapsed_ms(started))

    def _gwas_metadata(self, urls: list[str], warnings: list[str]) -> dict[str, Any]:
        url = f"{GWAS_CATALOG_API_V2_BASE}/metadata"
        try:
            payload = self.client.get_json(
                url,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(url)
            return payload if isinstance(payload, dict) else {}
        except KnowledgeRequestError:
            warnings.append("Optional GWAS Catalog metadata lookup failed; association details were still used.")
            return {}

    def _gwas_associations(self, params: dict[str, Any], urls: list[str]) -> dict[str, Any]:
        url = f"{GWAS_CATALOG_API_V2_BASE}/associations"
        payload = self.client.get_json(
            url,
            params=params,
            headers={"Accept": "application/json"},
            rate_limit_per_second=self.spec.rate_limit_per_second,
        )
        urls.append(url)
        return payload if isinstance(payload, dict) else {}

    def _gwas_record(
        self,
        item: dict[str, Any],
        context: KnowledgeQuery,
        metadata: dict[str, Any],
        query_context: str,
        urls: list[str],
        warnings: list[str],
    ) -> dict[str, Any]:
        association_id = _clean_cell(item.get("association_id"))
        accession_id = _clean_cell(item.get("accession_id"))
        study = self._gwas_study(accession_id, urls, warnings) if accession_id else {}
        loci = self._gwas_loci(item, association_id, urls, warnings)
        snp_alleles = self._gwas_snp_alleles(item.get("snp_allele"))
        effect_alleles = self._gwas_clean_list(item.get("snp_effect_allele"))
        rsids = self._gwas_rsids(snp_alleles, effect_alleles, item)
        traits = self._gwas_clean_list(item.get("reported_trait"))
        efo_traits = self._gwas_efo_traits(item.get("efo_traits"))
        label_variant = effect_alleles[0] if effect_alleles else (rsids[0] if rsids else context.gene)
        label_trait = traits[0] if traits else (efo_traits[0]["trait"] if efo_traits else "GWAS association")
        record = {
            "category": "population_association",
            "source": self.spec.name,
            "label": f"{label_variant} - {label_trait}",
            "summary": self._gwas_summary(item, study, loci, metadata, label_variant),
            "source_id": association_id,
            "url": self._gwas_link(item, "self", f"{GWAS_CATALOG_API_V2_BASE}/associations/{association_id}"),
            "browser_url": f"https://www.ebi.ac.uk/gwas/associations/{quote(association_id)}" if association_id else "",
            "variant": _first_rsid(context) or (rsids[0] if rsids else ""),
            "query_context": query_context,
            "association_id": association_id,
            "accession_id": accession_id,
            "study_url": f"https://www.ebi.ac.uk/gwas/studies/{quote(accession_id)}" if accession_id else "",
            "snp_url": self._gwas_link(item, "snp", ""),
            "rsids": rsids,
            "snp_effect_alleles": effect_alleles,
            "snp_alleles": snp_alleles,
            "risk_frequency": _clean_cell(item.get("risk_frequency")),
            "p_value": item.get("p_value"),
            "pvalue_mantissa": item.get("pvalue_mantissa"),
            "pvalue_exponent": item.get("pvalue_exponent"),
            "pvalue_description": _clean_cell(item.get("pvalue_description")),
            "beta": _clean_cell(item.get("beta")),
            "odds_ratio": _clean_cell(item.get("or_per_copy_number") or item.get("odds_ratio")),
            "range": _clean_cell(item.get("range")),
            "ci_lower": item.get("ci_lower"),
            "ci_upper": item.get("ci_upper"),
            "reported_traits": traits,
            "efo_traits": efo_traits,
            "background_efo_traits": self._gwas_efo_traits(item.get("bg_efo_traits")),
            "mapped_genes": self._gwas_clean_list(item.get("mapped_genes")),
            "locations": self._gwas_clean_list(item.get("locations")),
            "pubmed_id": _clean_cell(item.get("pubmed_id")),
            "first_author": _clean_cell(item.get("first_author")),
            "multi_snp_haplotype": item.get("multi_snp_haplotype"),
            "snp_interaction": item.get("snp_interaction"),
            "study": study,
            "loci": loci,
            "data_release_date": _clean_cell(metadata.get("data_release_date")),
            "api_release_date": _clean_cell(metadata.get("api_release_date")),
            "dbsnp_build": _clean_cell(metadata.get("dbsnp_build")),
            "gene_build": _clean_cell(metadata.get("gene_build")),
            "efo_version": _clean_cell(metadata.get("efo_version")),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _gwas_study(self, accession_id: str, urls: list[str], warnings: list[str]) -> dict[str, Any]:
        url = f"{GWAS_CATALOG_API_V2_BASE}/studies/{quote(accession_id)}"
        try:
            payload = self.client.get_json(
                url,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(url)
        except KnowledgeRequestError:
            warnings.append(f"Optional GWAS Catalog study detail lookup failed for {accession_id}; association details were still used.")
            return {}
        if not isinstance(payload, dict):
            return {}
        row = {
            "accession_id": _clean_cell(payload.get("accession_id")),
            "disease_trait": _clean_cell(payload.get("disease_trait")),
            "initial_sample_size": _clean_cell(payload.get("initial_sample_size")),
            "replication_sample_size": _clean_cell(payload.get("replication_sample_size")),
            "discovery_ancestry": self._gwas_clean_list(payload.get("discovery_ancestry")),
            "replication_ancestry": self._gwas_clean_list(payload.get("replication_ancestry")),
            "full_summary_stats_available": payload.get("full_summary_stats_available"),
            "full_summary_stats": _clean_cell(payload.get("full_summary_stats")),
            "terms_of_license": _clean_cell(payload.get("terms_of_license")),
            "snp_count": payload.get("snp_count"),
            "imputed": payload.get("imputed"),
            "pooled": payload.get("pooled"),
            "platforms": _clean_cell(payload.get("platforms")),
            "genotyping_technologies": self._gwas_clean_list(payload.get("genotyping_technologies")),
            "cohort": self._gwas_clean_list(payload.get("cohort")),
        }
        return {key: value for key, value in row.items() if value not in ("", [], None)}

    def _gwas_loci(
        self,
        item: dict[str, Any],
        association_id: str,
        urls: list[str],
        warnings: list[str],
    ) -> list[dict[str, Any]]:
        url = self._gwas_link(item, "loci", "")
        if not url:
            return []
        try:
            payload = self.client.get_json(
                url,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(url)
        except KnowledgeRequestError:
            warnings.append(f"Optional GWAS Catalog loci lookup failed for association {association_id}; association details were still used.")
            return []
        rows = ((payload.get("_embedded") or {}).get("loci") or []) if isinstance(payload, dict) else []
        loci: list[dict[str, Any]] = []
        for row in rows:
            if not isinstance(row, dict):
                continue
            locus = {
                "description": _clean_cell(row.get("description")),
                "strongest_risk_alleles": self._gwas_risk_alleles(row.get("strongest_risk_alleles")),
                "author_reported_genes": self._gwas_clean_list(row.get("author_reported_genes")),
                "url": self._gwas_link(row, "self", ""),
            }
            cleaned = {key: value for key, value in locus.items() if value not in ("", [], None)}
            if cleaned:
                loci.append(cleaned)
        return loci

    def _gwas_summary(
        self,
        item: dict[str, Any],
        study: dict[str, Any],
        loci: list[dict[str, Any]],
        metadata: dict[str, Any],
        label_variant: str,
    ) -> str:
        traits = self._gwas_clean_list(item.get("reported_trait"))
        efo_traits = [row["trait"] for row in self._gwas_efo_traits(item.get("efo_traits")) if row.get("trait")]
        trait_text = self._gwas_join([*traits[:2], *[trait for trait in efo_traits[:2] if trait not in traits]])
        if not trait_text:
            trait_text = _clean_cell(study.get("disease_trait")) or "reported trait"
        mapped_genes = self._gwas_clean_list(item.get("mapped_genes"))
        locations = self._gwas_clean_list(item.get("locations"))
        lead = f"GWAS Catalog association {label_variant} with {trait_text}"
        if mapped_genes:
            lead = f"{lead} mapped to {', '.join(mapped_genes[:3])}"
        if locations:
            lead = f"{lead} at {', '.join(locations[:2])}"

        parts: list[str] = []
        p_value = self._gwas_pvalue_text(item)
        if p_value:
            parts.append(f"p={p_value}")
        effect = self._gwas_effect_text(item)
        if effect:
            parts.append(effect)
        risk_frequency = _clean_cell(item.get("risk_frequency"))
        if risk_frequency:
            parts.append(f"risk frequency {risk_frequency}")
        accession_id = _clean_cell(item.get("accession_id"))
        pubmed_id = _clean_cell(item.get("pubmed_id"))
        first_author = _clean_cell(item.get("first_author"))
        study_parts = []
        if accession_id:
            study_parts.append(accession_id)
        if pubmed_id:
            study_parts.append(f"PMID {pubmed_id}")
        if first_author:
            study_parts.append(f"first author {first_author}")
        if study_parts:
            parts.append(f"study {', '.join(study_parts)}")
        sample_size = _clean_cell(study.get("initial_sample_size"))
        if sample_size:
            parts.append(f"initial sample {self._gwas_clip(sample_size, 180)}")
        discovery_ancestry = self._gwas_clean_list(study.get("discovery_ancestry"))
        if discovery_ancestry:
            parts.append(f"discovery ancestry {self._gwas_clip(discovery_ancestry[0], 180)}")
        if study.get("full_summary_stats_available") is True:
            parts.append("full summary statistics available")
        elif study.get("full_summary_stats_available") is False:
            parts.append("full summary statistics not marked available")
        strongest = self._gwas_loci_text(loci)
        if strongest:
            parts.append(strongest)
        release = _clean_cell(metadata.get("data_release_date"))
        if release:
            parts.append(f"data release {release}")
        return f"{lead}: {'; '.join(parts)}." if parts else f"{lead}."

    def _gwas_association_rows(self, payload: dict[str, Any]) -> list[dict[str, Any]]:
        rows = ((payload.get("_embedded") or {}).get("associations") or []) if isinstance(payload, dict) else []
        return [row for row in rows if isinstance(row, dict)]

    def _gwas_total(self, payload: dict[str, Any]) -> int:
        page = payload.get("page") if isinstance(payload, dict) else {}
        try:
            return int(page.get("totalElements") or 0)
        except (TypeError, ValueError):
            return 0

    def _gwas_snp_alleles(self, values: Any) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "rsid": _clean_cell(value.get("rs_id")),
                "effect_allele": _clean_cell(value.get("effect_allele")),
            }
            cleaned = {key: item for key, item in row.items() if item}
            if cleaned and cleaned not in rows:
                rows.append(cleaned)
        return rows

    def _gwas_rsids(self, snp_alleles: list[dict[str, str]], effect_alleles: list[str], item: dict[str, Any]) -> list[str]:
        rsids = [row["rsid"] for row in snp_alleles if row.get("rsid")]
        for effect_allele in effect_alleles:
            rsid, _, _allele = effect_allele.partition("-")
            if rsid.lower().startswith("rs"):
                rsids.append(rsid)
        snp_href = self._gwas_link(item, "snp", "")
        if "/single-nucleotide-polymorphisms/" in snp_href:
            rsids.append(snp_href.rsplit("/", 1)[-1])
        return self._gwas_clean_list(rsids)

    def _gwas_efo_traits(self, values: Any) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "id": _clean_cell(value.get("efo_id")),
                "trait": _clean_cell(value.get("efo_trait")),
            }
            cleaned = {key: item for key, item in row.items() if item}
            if cleaned and cleaned not in rows:
                rows.append(cleaned)
        return rows

    def _gwas_risk_alleles(self, values: Any) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "risk_allele_name": _clean_cell(value.get("risk_allele_name")),
                "genome_wide": value.get("genome_wide"),
                "limited_list": value.get("limited_list"),
            }
            cleaned = {key: item for key, item in row.items() if item not in ("", None)}
            if cleaned:
                rows.append(cleaned)
        return rows

    def _gwas_pvalue_text(self, item: dict[str, Any]) -> str:
        value = item.get("p_value")
        try:
            number = float(value)
        except (TypeError, ValueError):
            mantissa = _clean_cell(item.get("pvalue_mantissa"))
            exponent = _clean_cell(item.get("pvalue_exponent"))
            return f"{mantissa}e{exponent}" if mantissa and exponent else _clean_cell(value)
        return f"{number:.3g}"

    def _gwas_effect_text(self, item: dict[str, Any]) -> str:
        parts: list[str] = []
        beta = _clean_cell(item.get("beta"))
        odds_ratio = _clean_cell(item.get("or_per_copy_number") or item.get("odds_ratio"))
        if beta:
            parts.append(f"beta {beta}")
        if odds_ratio:
            parts.append(f"OR {odds_ratio}")
        ci_lower = _clean_cell(item.get("ci_lower"))
        ci_upper = _clean_cell(item.get("ci_upper"))
        range_text = _clean_cell(item.get("range"))
        if ci_lower and ci_upper:
            parts.append(f"CI {ci_lower}-{ci_upper}")
        elif range_text and range_text != "-":
            parts.append(f"range {range_text}")
        return "; ".join(parts)

    def _gwas_loci_text(self, loci: list[dict[str, Any]]) -> str:
        if not loci:
            return ""
        first = loci[0]
        alleles = first.get("strongest_risk_alleles") if isinstance(first, dict) else []
        allele_names = [
            _clean_cell(allele.get("risk_allele_name"))
            for allele in alleles
            if isinstance(allele, dict) and _clean_cell(allele.get("risk_allele_name"))
        ]
        author_genes = self._gwas_clean_list(first.get("author_reported_genes"))
        parts = []
        if allele_names:
            parts.append(f"strongest risk allele {', '.join(allele_names[:3])}")
        if author_genes:
            parts.append(f"author-reported genes {', '.join(author_genes[:3])}")
        return "; ".join(parts)

    def _gwas_link(self, item: dict[str, Any], rel: str, fallback: str) -> str:
        links = item.get("_links") if isinstance(item.get("_links"), dict) else {}
        link = links.get(rel) if isinstance(links.get(rel), dict) else {}
        return _clean_cell(link.get("href")) or fallback

    def _gwas_clean_list(self, values: Any) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, list):
            values = [values]
        cleaned: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in cleaned:
                cleaned.append(text)
        return cleaned

    def _gwas_join(self, values: list[str]) -> str:
        cleaned = self._gwas_clean_list(values)
        return " / ".join(cleaned[:4])

    def _gwas_clip(self, value: str, limit: int) -> str:
        if len(value) <= limit:
            return value
        return value[: limit - 1].rstrip() + "..."


class PgsCatalogConnector(BaseConnector):
    """PGS Catalog connector."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        rsid = _first_rsid(context)
        urls: list[str] = []
        warnings: list[str] = []
        if not rsid:
            return SourceResult(self.spec.key, "skipped", "No rsID was available for PGS Catalog lookup.")
        url = f"{PGS_CATALOG_REST_BASE}/variant/{quote(rsid)}"
        try:
            payload = self.client.get_json(
                url,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(url)
        except KnowledgeRequestError as exc:
            urls.append(url)
            if self._pgs_is_not_found(exc):
                return SourceResult(
                    self.spec.key,
                    "ok",
                    f"Queried PGS Catalog; no deposited score variant record found for {rsid}.",
                    [],
                    queried_urls=urls,
                    elapsed_ms=_elapsed_ms(started),
                )
            return SourceResult(self.spec.key, "failed", str(exc), errors=[str(exc)], queried_urls=urls, elapsed_ms=_elapsed_ms(started))

        pgs_ids = self._pgs_associated_ids(payload)
        if not pgs_ids:
            return SourceResult(
                self.spec.key,
                "ok",
                f"Queried PGS Catalog; variant {rsid} is present but has no linked polygenic scores.",
                [],
                queried_urls=urls,
                elapsed_ms=_elapsed_ms(started),
            )

        records: list[dict[str, Any]] = []
        for pgs_id in pgs_ids[:PGS_CATALOG_MAX_LINKED_SCORES]:
            score = self._pgs_score(pgs_id, urls, warnings)
            performance_records = self._pgs_performance(pgs_id, urls, warnings) if score else []
            records.append(self._pgs_record(rsid, pgs_id, payload, score, performance_records))

        message = f"Queried PGS Catalog; {len(records)} linked polygenic score record(s) returned for {rsid}."
        if len(pgs_ids) > len(records):
            message = f"{message} {len(pgs_ids)} total linked score(s) are available."
        return SourceResult(
            self.spec.key,
            "ok",
            message,
            records,
            warnings=warnings,
            queried_urls=urls,
            elapsed_ms=_elapsed_ms(started),
        )

    def _pgs_score(self, pgs_id: str, urls: list[str], warnings: list[str]) -> dict[str, Any]:
        url = f"{PGS_CATALOG_REST_BASE}/score/{quote(pgs_id)}"
        try:
            payload = self.client.get_json(
                url,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(url)
            return payload if isinstance(payload, dict) else {}
        except KnowledgeRequestError:
            urls.append(url)
            warnings.append(f"Optional PGS Catalog score detail lookup failed for {pgs_id}; variant link was still used.")
            return {}

    def _pgs_performance(self, pgs_id: str, urls: list[str], warnings: list[str]) -> list[dict[str, Any]]:
        url = f"{PGS_CATALOG_REST_BASE}/performance/search"
        try:
            payload = self.client.get_json(
                url,
                params={"pgs_id": pgs_id},
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(url)
        except KnowledgeRequestError:
            urls.append(url)
            warnings.append(f"Optional PGS Catalog performance lookup failed for {pgs_id}; score details were still used.")
            return []
        rows = payload.get("results") if isinstance(payload, dict) else []
        performance: list[dict[str, Any]] = []
        for row in (rows or [])[:PGS_CATALOG_MAX_PERFORMANCE_RECORDS]:
            if isinstance(row, dict):
                performance.append(self._pgs_performance_record(row))
        return performance

    def _pgs_record(
        self,
        rsid: str,
        pgs_id: str,
        variant_payload: dict[str, Any],
        score: dict[str, Any],
        performance_records: list[dict[str, Any]],
    ) -> dict[str, Any]:
        score_id = _clean_cell(score.get("id")) or pgs_id
        name = _clean_cell(score.get("name"))
        trait = _clean_cell(score.get("trait_reported"))
        label_parts = [score_id]
        if name:
            label_parts.append(name)
        if trait:
            label_parts.append(trait)
        record = {
            "category": "polygenic_score",
            "source": self.spec.name,
            "label": " - ".join(label_parts),
            "summary": self._pgs_summary(rsid, pgs_id, score, performance_records),
            "source_id": score_id,
            "url": f"https://www.pgscatalog.org/score/{quote(score_id)}/",
            "variant_url": f"https://www.pgscatalog.org/variant/{quote(rsid)}/",
            "variant": rsid,
            "linked_variant": rsid,
            "pgs_id": score_id,
            "name": name,
            "trait_reported": trait,
            "trait_additional": _clean_cell(score.get("trait_additional")),
            "trait_efo": self._pgs_traits(score.get("trait_efo")),
            "method_name": _clean_cell(score.get("method_name")),
            "method_params": _clean_cell(score.get("method_params")),
            "variants_number": score.get("variants_number"),
            "variants_interactions": score.get("variants_interactions"),
            "variants_genomebuild": _clean_cell(score.get("variants_genomebuild")),
            "weight_type": _clean_cell(score.get("weight_type")),
            "publication": self._pgs_publication(score.get("publication")),
            "matches_publication": score.get("matches_publication"),
            "samples_variants": self._pgs_samples(score.get("samples_variants")),
            "samples_training": self._pgs_samples(score.get("samples_training")),
            "ancestry_distribution": score.get("ancestry_distribution"),
            "performance": performance_records,
            "ftp_scoring_file": _clean_cell(score.get("ftp_scoring_file")),
            "ftp_harmonized_scoring_files": self._pgs_harmonized_files(score.get("ftp_harmonized_scoring_files")),
            "date_release": _clean_cell(score.get("date_release")),
            "license": _clean_cell(score.get("license")),
            "associated_pgs_ids": self._pgs_associated_ids(variant_payload),
        }
        if not score:
            record["summary"] = (
                f"PGS Catalog variant {rsid} is linked to polygenic score {pgs_id}, "
                "but score details were not returned by the optional lookup."
            )
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _pgs_summary(
        self,
        rsid: str,
        pgs_id: str,
        score: dict[str, Any],
        performance_records: list[dict[str, Any]],
    ) -> str:
        score_id = _clean_cell(score.get("id")) or pgs_id
        name = _clean_cell(score.get("name"))
        trait = _clean_cell(score.get("trait_reported")) or "reported trait"
        lead = f"PGS Catalog score {score_id}"
        if name:
            lead = f"{lead} ({name})"
        lead = f"{lead} includes variant {rsid} and predicts {trait}"

        parts: list[str] = []
        mapped_traits = [trait_row["label"] for trait_row in self._pgs_traits(score.get("trait_efo")) if trait_row.get("label")]
        if mapped_traits:
            parts.append(f"mapped traits: {', '.join(mapped_traits[:3])}")
        variants_number = self._pgs_count(score.get("variants_number"))
        if variants_number:
            parts.append(f"{variants_number} variants")
        interactions = self._pgs_count(score.get("variants_interactions"))
        if interactions and interactions != "0":
            parts.append(f"{interactions} interaction terms")
        method = _clean_cell(score.get("method_name"))
        method_params = _clean_cell(score.get("method_params"))
        if method:
            method_text = f"method {method}"
            if method_params:
                method_text = f"{method_text} ({method_params})"
            parts.append(method_text)
        weight_type = _clean_cell(score.get("weight_type"))
        if weight_type and weight_type != "NR":
            parts.append(f"weight type {weight_type}")
        genome_build = _clean_cell(score.get("variants_genomebuild"))
        if genome_build and genome_build != "NR":
            parts.append(f"variant genome build {genome_build}")
        sample_text = self._pgs_sample_text(score.get("samples_variants"), "variant-source sample")
        if sample_text:
            parts.append(sample_text)
        training_text = self._pgs_sample_text(score.get("samples_training"), "training sample")
        if training_text:
            parts.append(training_text)
        ancestry_text = self._pgs_ancestry_distribution_text(score.get("ancestry_distribution"))
        if ancestry_text:
            parts.append(ancestry_text)
        performance_text = self._pgs_performance_text(performance_records)
        if performance_text:
            parts.append(performance_text)
        publication_text = self._pgs_publication_text(score.get("publication"))
        if publication_text:
            parts.append(publication_text)
        harmonized = self._pgs_harmonized_files(score.get("ftp_harmonized_scoring_files"))
        if harmonized:
            parts.append(f"harmonized scoring files: {', '.join(harmonized.keys())}")
        release = _clean_cell(score.get("date_release"))
        if release:
            parts.append(f"released {release}")
        return f"{lead}: {'; '.join(parts)}." if parts else f"{lead}."

    def _pgs_associated_ids(self, payload: Any) -> list[str]:
        if not isinstance(payload, dict):
            return []
        values = payload.get("associated_pgs_ids") or payload.get("associated_PGS_ids") or payload.get("pgs_ids")
        ids: list[str] = []
        if isinstance(values, dict):
            for nested in values.values():
                ids.extend(self._pgs_clean_list(nested))
        else:
            ids.extend(self._pgs_clean_list(values))
        return [value for value in ids if value.upper().startswith("PGS")]

    def _pgs_traits(self, values: Any) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "id": _clean_cell(value.get("id")),
                "label": _clean_cell(value.get("label")),
                "description": _clean_cell(value.get("description")),
                "url": _clean_cell(value.get("url")),
            }
            cleaned = {key: item for key, item in row.items() if item}
            if cleaned and cleaned not in rows:
                rows.append(cleaned)
        return rows

    def _pgs_publication(self, value: Any) -> dict[str, Any]:
        if not isinstance(value, dict):
            return {}
        row = {
            "id": _clean_cell(value.get("id")),
            "title": _clean_cell(value.get("title")),
            "doi": _clean_cell(value.get("doi")),
            "pmid": _clean_cell(value.get("PMID") or value.get("pmid")),
            "journal": _clean_cell(value.get("journal")),
            "first_author": _clean_cell(value.get("firstauthor") or value.get("first_author")),
            "date_publication": _clean_cell(value.get("date_publication")),
        }
        return {key: item for key, item in row.items() if item}

    def _pgs_samples(self, values: Any) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        if not isinstance(values, list):
            return rows
        for value in values[:5]:
            if not isinstance(value, dict):
                continue
            row = {
                "sample_number": value.get("sample_number"),
                "sample_cases": value.get("sample_cases"),
                "sample_controls": value.get("sample_controls"),
                "sample_percent_male": value.get("sample_percent_male"),
                "sample_age": _clean_cell(value.get("sample_age")),
                "phenotyping": _clean_cell(value.get("phenotyping_free")),
                "followup_time": _clean_cell(value.get("followup_time")),
                "ancestry_broad": _clean_cell(value.get("ancestry_broad")),
                "ancestry_free": _clean_cell(value.get("ancestry_free")),
                "ancestry_country": _clean_cell(value.get("ancestry_country")),
                "ancestry_additional": _clean_cell(value.get("ancestry_additional")),
                "source_gwas_catalog": _clean_cell(value.get("source_GWAS_catalog")),
                "source_pmid": _clean_cell(value.get("source_PMID")),
                "source_doi": _clean_cell(value.get("source_DOI")),
                "cohorts": self._pgs_cohorts(value.get("cohorts")),
                "cohorts_additional": _clean_cell(value.get("cohorts_additional")),
            }
            rows.append({key: item for key, item in row.items() if item not in ("", [], None)})
        return rows

    def _pgs_cohorts(self, values: Any) -> list[str]:
        if not isinstance(values, list):
            return []
        cohorts: list[str] = []
        for value in values:
            if isinstance(value, dict):
                text = _clean_cell(value.get("name_short") or value.get("name_full") or value.get("name_others"))
            else:
                text = _clean_cell(value)
            if text and text not in cohorts:
                cohorts.append(text)
        return cohorts

    def _pgs_performance_record(self, value: dict[str, Any]) -> dict[str, Any]:
        sampleset = value.get("sampleset") if isinstance(value.get("sampleset"), dict) else {}
        row = {
            "id": _clean_cell(value.get("id")),
            "associated_pgs_id": _clean_cell(value.get("associated_pgs_id")),
            "phenotyping_reported": _clean_cell(value.get("phenotyping_reported")),
            "publication": self._pgs_publication(value.get("publication")),
            "sampleset_id": _clean_cell(sampleset.get("id")),
            "samples": self._pgs_samples(sampleset.get("samples")),
            "metrics": self._pgs_performance_metrics(value.get("performance_metrics")),
            "covariates": _clean_cell(value.get("covariates")),
            "performance_comments": _clean_cell(value.get("performance_comments")),
        }
        return {key: item for key, item in row.items() if item not in ("", [], {}, None)}

    def _pgs_performance_metrics(self, value: Any) -> list[dict[str, Any]]:
        if not isinstance(value, dict):
            return []
        rows: list[dict[str, Any]] = []
        for category, label in (
            ("effect_sizes", "effect_size"),
            ("class_acc", "classification_accuracy"),
            ("othermetrics", "other"),
        ):
            for metric in value.get(category) or []:
                if not isinstance(metric, dict):
                    continue
                row = {
                    "category": label,
                    "name_long": _clean_cell(metric.get("name_long")),
                    "name_short": _clean_cell(metric.get("name_short")),
                    "estimate": metric.get("estimate"),
                    "ci_lower": metric.get("ci_lower"),
                    "ci_upper": metric.get("ci_upper"),
                }
                rows.append({key: item for key, item in row.items() if item not in ("", None)})
        return rows

    def _pgs_harmonized_files(self, value: Any) -> dict[str, str]:
        if not isinstance(value, dict):
            return {}
        rows: dict[str, str] = {}
        for build, payload in value.items():
            if not isinstance(payload, dict):
                continue
            url = _clean_cell(payload.get("positions") or payload.get("effect_allele") or payload.get("url"))
            if url:
                rows[_clean_cell(build)] = url
        return rows

    def _pgs_sample_text(self, values: Any, label: str) -> str:
        samples = self._pgs_samples(values)
        if not samples:
            return ""
        first = samples[0]
        sample_number = self._pgs_count(first.get("sample_number"))
        ancestry = _clean_cell(first.get("ancestry_broad"))
        countries = _clean_cell(first.get("ancestry_country"))
        parts = [label]
        if sample_number:
            parts.append(sample_number)
        if ancestry:
            parts.append(ancestry)
        if countries:
            parts.append(f"({self._pgs_clip(countries, 120)})")
        source = _clean_cell(first.get("source_gwas_catalog") or first.get("source_pmid") or first.get("source_doi"))
        if source:
            parts.append(f"source {source}")
        return " ".join(parts)

    def _pgs_ancestry_distribution_text(self, value: Any) -> str:
        if not isinstance(value, dict):
            return ""
        parts: list[str] = []
        for key, label in (("gwas", "GWAS"), ("eval", "evaluation")):
            payload = value.get(key)
            if not isinstance(payload, dict):
                continue
            count = self._pgs_count(payload.get("count"))
            dist = payload.get("dist") if isinstance(payload.get("dist"), dict) else {}
            dist_text = ", ".join(f"{group} {amount}%" for group, amount in dist.items())
            text = f"{label} ancestry"
            if count:
                text = f"{text} n={count}"
            if dist_text:
                text = f"{text} ({dist_text})"
            parts.append(text)
        return "; ".join(parts)

    def _pgs_performance_text(self, performance_records: list[dict[str, Any]]) -> str:
        if not performance_records:
            return ""
        first = performance_records[0]
        metrics = first.get("metrics") if isinstance(first.get("metrics"), list) else []
        metric_texts = []
        for metric in metrics[:3]:
            name = _clean_cell(metric.get("name_short") or metric.get("name_long"))
            estimate = _clean_cell(metric.get("estimate"))
            ci_lower = _clean_cell(metric.get("ci_lower"))
            ci_upper = _clean_cell(metric.get("ci_upper"))
            if not (name and estimate):
                continue
            text = f"{name} {estimate}"
            if ci_lower and ci_upper:
                text = f"{text} ({ci_lower}-{ci_upper})"
            metric_texts.append(text)
        if not metric_texts:
            return ""
        phenotype = _clean_cell(first.get("phenotyping_reported"))
        sample_text = self._pgs_sample_text(first.get("samples"), "evaluated sample")
        details = ", ".join(metric_texts)
        if phenotype:
            details = f"{phenotype}: {details}"
        if sample_text:
            details = f"{details}; {sample_text}"
        return f"performance {details}"

    def _pgs_publication_text(self, value: Any) -> str:
        publication = self._pgs_publication(value)
        if not publication:
            return ""
        parts = []
        if publication.get("first_author"):
            parts.append(publication["first_author"])
        if publication.get("date_publication"):
            parts.append(publication["date_publication"])
        if publication.get("pmid"):
            parts.append(f"PMID {publication['pmid']}")
        if publication.get("doi"):
            parts.append(f"DOI {publication['doi']}")
        return f"publication {', '.join(parts)}" if parts else ""

    def _pgs_clean_list(self, values: Any) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, list):
            values = [values]
        cleaned: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in cleaned:
                cleaned.append(text)
        return cleaned

    def _pgs_count(self, value: Any) -> str:
        try:
            return f"{int(value):,}"
        except (TypeError, ValueError):
            return _clean_cell(value)

    def _pgs_clip(self, value: str, limit: int) -> str:
        if len(value) <= limit:
            return value
        return value[: limit - 1].rstrip() + "..."

    def _pgs_is_not_found(self, exc: KnowledgeRequestError) -> bool:
        text = str(exc).lower()
        return "404" in text or "not found" in text


class IgsrConnector(BaseConnector):
    """1000 Genomes Project / IGSR data-access connector."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        variant = _first_variant(context)
        chrom = self._chromosome(context, variant)
        if not chrom:
            return SourceResult(
                self.spec.key,
                "skipped",
                "No chromosome or variant coordinate was available for 1000 Genomes/IGSR lookup.",
                elapsed_ms=_elapsed_ms(started),
            )

        urls: list[str] = []
        warnings: list[str] = []
        high_listing = self._get_listing(
            IGSR_HIGH_COVERAGE_PHASED_URL,
            urls,
            warnings,
            "high-coverage GRCh38 phased VCF listing",
        )
        phase3_listing = self._get_listing(
            IGSR_PHASE3_RELEASE_URL,
            urls,
            warnings,
            "Phase 3 GRCh37 release listing",
        )
        high_file = self._file_meta(high_listing, IGSR_HIGH_COVERAGE_PHASED_URL, self._high_coverage_filename(chrom))
        phase3_file = self._file_meta(phase3_listing, IGSR_PHASE3_RELEASE_URL, self._phase3_filename(chrom))
        phase3_site_file = self._file_meta(
            phase3_listing,
            IGSR_PHASE3_RELEASE_URL,
            "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz",
        )
        sample_panel = self._file_meta(
            phase3_listing,
            IGSR_PHASE3_RELEASE_URL,
            "integrated_call_samples_v3.20130502.ALL.panel",
        )

        primary_kind = self._primary_kind(context, high_file, phase3_file)
        primary_file = high_file if primary_kind == "high_coverage" else phase3_file
        if not primary_file:
            primary_file = high_file or phase3_file
            primary_kind = "high_coverage" if primary_file == high_file else "phase3"

        if not primary_file:
            message = "Queried IGSR/1000 Genomes release listings, but no chromosome VCF was identified for chromosome "
            message = f"{message}{chrom}."
            status = "failed" if warnings and len(warnings) >= 2 else "ok"
            return SourceResult(
                self.spec.key,
                status,
                message,
                [],
                warnings=warnings,
                errors=warnings if status == "failed" else [],
                queried_urls=urls,
                elapsed_ms=_elapsed_ms(started),
            )

        record = self._record(
            context,
            variant,
            chrom,
            primary_kind,
            primary_file,
            high_file,
            phase3_file,
            phase3_site_file,
            sample_panel,
        )
        label = self._variant_label(context, variant, chrom)
        return SourceResult(
            self.spec.key,
            "ok",
            f"Queried IGSR/1000 Genomes FTP release listings; 1 data-access record returned for {label}.",
            [record],
            warnings=warnings,
            queried_urls=urls,
            elapsed_ms=_elapsed_ms(started),
        )

    def _get_listing(self, url: str, urls: list[str], warnings: list[str], label: str) -> str:
        try:
            text = self.client.get_text(
                url,
                headers={"Accept": "text/html,text/plain,*/*"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            urls.append(url)
            return text
        except KnowledgeRequestError:
            urls.append(url)
            warnings.append(f"Optional IGSR {label} lookup failed; available release details were still used.")
            return ""

    def _record(
        self,
        context: KnowledgeQuery,
        variant: Any,
        chrom: str,
        primary_kind: str,
        primary_file: dict[str, str],
        high_file: dict[str, str],
        phase3_file: dict[str, str],
        phase3_site_file: dict[str, str],
        sample_panel: dict[str, str],
    ) -> dict[str, Any]:
        is_high = primary_kind == "high_coverage"
        assembly = "GRCh38" if is_high else "GRCh37"
        dataset = (
            "1000 Genomes 2504 high-coverage phased callset"
            if is_high
            else "1000 Genomes Phase 3 integrated callset"
        )
        record = {
            "category": "population_reference_panel",
            "source": self.spec.name,
            "label": f"{dataset} chr{chrom}",
            "summary": self._summary(context, variant, chrom, primary_kind, primary_file, high_file, phase3_file),
            "source_id": self._source_id(primary_kind, chrom),
            "url": primary_file.get("url"),
            "data_portal_url": IGSR_PORTAL_URL,
            "variant": _clean_cell(getattr(variant, "label", "")) or self._variant_label(context, variant, chrom),
            "rsid": _clean_cell(getattr(variant, "rsid", "")),
            "gene": context.gene,
            "chromosome": chrom,
            "position": getattr(variant, "pos", None),
            "query_build": context.genome_build,
            "assembly": assembly,
            "dataset": dataset,
            "samples": 2504,
            "populations": 26,
            "population_groups": ["AFR", "AMR", "EAS", "EUR", "SAS"],
            "file_url": primary_file.get("url"),
            "index_url": primary_file.get("index_url"),
            "file_name": primary_file.get("name"),
            "file_size": primary_file.get("size"),
            "last_modified": primary_file.get("last_modified"),
            "data_fields": self._data_fields(primary_kind),
            "readme_urls": [IGSR_HIGH_COVERAGE_README_URL, IGSR_PHASE3_README_URL, IGSR_PHASE3_ANNOTATION_README_URL],
            "related_files": self._related_files(primary_kind, high_file, phase3_file, phase3_site_file, sample_panel),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _summary(
        self,
        context: KnowledgeQuery,
        variant: Any,
        chrom: str,
        primary_kind: str,
        primary_file: dict[str, str],
        high_file: dict[str, str],
        phase3_file: dict[str, str],
    ) -> str:
        variant_label = self._variant_label(context, variant, chrom)
        file_text = self._file_text(primary_file)
        if primary_kind == "high_coverage":
            lead = (
                f"IGSR/1000 Genomes high-coverage GRCh38 data-access context for {variant_label}: "
                f"chromosome VCF {file_text} and tabix index are available for the 30x NYGC 2504-sample "
                "Phase 3 panel"
            )
            parts = [
                "use the indexed VCF to extract exact genotypes, AC/AN, or allele counts for this coordinate",
            ]
            if phase3_file:
                parts.append(
                    "legacy Phase 3 GRCh37 integrated data are also available via "
                    f"{phase3_file.get('name')} ({phase3_file.get('size')}) with global and superpopulation AF tags"
                )
            parts.append("rsIDs were removed from the Phase 3 v5b VCF, so coordinate lookup or Ensembl rsID mapping is needed")
            return f"{lead}; {'; '.join(parts)}."

        lead = (
            f"IGSR/1000 Genomes Phase 3 GRCh37 data-access context for {variant_label}: "
            f"chromosome genotype VCF {file_text} and tabix index are available"
        )
        parts = [
            "final Phase 3 release has phased genotypes for 2,504 individuals from 26 populations",
            "VCF INFO contains AN, AC, global AF, and EAS/EUR/AFR/AMR/SAS allele-frequency tags",
            "rsIDs were removed from the v5b VCF, so coordinate lookup or Ensembl rsID mapping is needed",
        ]
        if self._is_hg38(context):
            parts.append("input build is GRCh38; lift over coordinates or prefer the high-coverage GRCh38 file before direct extraction")
        elif high_file:
            parts.append(f"a GRCh38 high-coverage phased file is also available as {high_file.get('name')} ({high_file.get('size')})")
        return f"{lead}; {'; '.join(parts)}."

    def _data_fields(self, primary_kind: str) -> list[str]:
        if primary_kind == "phase3":
            return ["phased genotypes", "AN", "AC", "AF", "EAS_AF", "EUR_AF", "AFR_AF", "AMR_AF", "SAS_AF"]
        return ["phased genotypes", "indexed chromosome VCF", "sample-level genotype extraction"]

    def _related_files(
        self,
        primary_kind: str,
        high_file: dict[str, str],
        phase3_file: dict[str, str],
        phase3_site_file: dict[str, str],
        sample_panel: dict[str, str],
    ) -> list[dict[str, str]]:
        files: list[dict[str, str]] = []
        for label, row in (
            ("High-coverage GRCh38 phased chromosome VCF", high_file),
            ("Phase 3 GRCh37 integrated chromosome VCF", phase3_file),
            ("Phase 3 WGS sites VCF", phase3_site_file),
            ("Phase 3 sample population panel", sample_panel),
        ):
            if not row:
                continue
            files.append(
                {
                    "label": label,
                    "name": row.get("name", ""),
                    "url": row.get("url", ""),
                    "index_url": row.get("index_url", ""),
                    "size": row.get("size", ""),
                    "last_modified": row.get("last_modified", ""),
                    "primary": label.startswith("High-coverage") if primary_kind == "high_coverage" else label.startswith("Phase 3 GRCh37"),
                }
            )
        return [{key: value for key, value in row.items() if value not in ("", None)} for row in files]

    def _file_meta(self, listing: str, base_url: str, filename: str) -> dict[str, str]:
        if not listing or not filename:
            return {}
        entries = self._listing_entries(listing, base_url)
        row = entries.get(filename)
        if not row:
            return {}
        index_name = f"{filename}.tbi"
        index_row = entries.get(index_name, {})
        if index_row:
            row["index_url"] = index_row.get("url", "")
            row["index_size"] = index_row.get("size", "")
            row["index_last_modified"] = index_row.get("last_modified", "")
        return row

    def _listing_entries(self, listing: str, base_url: str) -> dict[str, dict[str, str]]:
        entries: dict[str, dict[str, str]] = {}
        pattern = re.compile(
            r'<tr><td[^>]*>.*?</td><td><a href="(?P<href>[^"]+)">(?P<name>.*?)</a></td>'
            r'<td[^>]*>(?P<modified>.*?)</td><td[^>]*>(?P<size>.*?)</td>',
            flags=re.IGNORECASE | re.DOTALL,
        )
        for match in pattern.finditer(listing):
            name = html.unescape(re.sub(r"<.*?>", "", match.group("name"))).strip()
            href = html.unescape(match.group("href")).strip()
            if not name or name == "Parent Directory" or href.endswith("/"):
                continue
            url = href if href.startswith(("http://", "https://")) else f"{base_url.rstrip('/')}/{href.lstrip('/')}"
            entries[name] = {
                "name": name,
                "url": url,
                "last_modified": html.unescape(re.sub(r"<.*?>", "", match.group("modified"))).strip(),
                "size": html.unescape(re.sub(r"<.*?>", "", match.group("size"))).strip(),
            }
        return entries

    def _primary_kind(self, context: KnowledgeQuery, high_file: dict[str, str], phase3_file: dict[str, str]) -> str:
        if self._is_hg38(context) and high_file:
            return "high_coverage"
        if phase3_file:
            return "phase3"
        return "high_coverage" if high_file else "phase3"

    def _source_id(self, primary_kind: str, chrom: str) -> str:
        if primary_kind == "high_coverage":
            return f"1000G_2504_high_coverage_GRCh38_chr{chrom}"
        return f"1000G_phase3_20130502_GRCh37_chr{chrom}"

    def _file_text(self, row: dict[str, str]) -> str:
        name = row.get("name", "the chromosome VCF")
        size = row.get("size")
        modified = row.get("last_modified")
        details = []
        if size:
            details.append(size)
        if modified:
            details.append(f"updated {modified}")
        return f"{name} ({', '.join(details)})" if details else name

    def _variant_label(self, context: KnowledgeQuery, variant: Any, chrom: str) -> str:
        label = _clean_cell(getattr(variant, "rsid", ""))
        position = getattr(variant, "pos", None)
        coordinate = f"{chrom}:{position}" if position else context.region or f"chr{chrom}"
        change = ""
        ref = _clean_cell(getattr(variant, "ref", ""))
        alt = _clean_cell(getattr(variant, "alt", ""))
        if ref and alt:
            change = f" {ref}>{alt}"
        return f"{label} at {coordinate}{change}" if label else f"{coordinate}{change}"

    def _chromosome(self, context: KnowledgeQuery, variant: Any) -> str:
        chrom = _clean_cell(getattr(variant, "chrom", ""))
        if not chrom and context.region:
            chrom = context.region.split(":", 1)[0]
        chrom = chrom.removeprefix("chr").removeprefix("Chr").strip()
        upper = chrom.upper()
        if upper in {"M", "MT", "MITO", "MITOCHONDRIAL"}:
            return "MT"
        return upper if upper in {"X", "Y"} else chrom

    def _is_hg38(self, context: KnowledgeQuery) -> bool:
        build = context.genome_build.lower()
        return "38" in build or "grch38" in build or "b38" in build

    def _high_coverage_filename(self, chrom: str) -> str:
        if chrom == "X":
            return "CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
        if chrom in {"Y", "MT"}:
            return ""
        return f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"

    def _phase3_filename(self, chrom: str) -> str:
        if chrom == "X":
            return "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
        if chrom == "Y":
            return "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
        if chrom == "MT":
            return "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
        return f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"


class EncodeConnector(BaseConnector):
    """ENCODE Portal search connector."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        if self.spec.connector_kind == "screen":
            return ScreenConnector(self.spec, self.client, self.credential).query(context)
        return self._query_encode_portal(context, started)

    def _query_encode_portal(self, context: KnowledgeQuery, started: float) -> SourceResult:
        headers = {"Accept": "application/json"}
        urls: list[str] = []
        warnings: list[str] = []
        records: list[dict[str, Any]] = []
        failures = 0
        seen: set[str] = set()

        target_payload, failed = self._encode_get(
            ENCODE_SEARCH_URL,
            {
                "type": "Experiment",
                "target.genes.symbol": context.gene,
                "format": "json",
                "frame": "object",
                "limit": ENCODE_MAX_RECORDS,
            },
            headers,
            urls,
            warnings,
            "target-gene experiment search",
        )
        failures += 1 if failed else 0
        records.extend(self._experiment_records(context, target_payload, seen, "target_gene"))

        text_payload, failed = self._encode_get(
            ENCODE_SEARCH_URL,
            {
                "type": "Experiment",
                "searchTerm": context.gene,
                "format": "json",
                "frame": "object",
                "limit": ENCODE_MAX_RECORDS,
            },
            headers,
            urls,
            warnings,
            "gene-text experiment search",
        )
        failures += 1 if failed else 0
        records.extend(self._experiment_records(context, text_payload, seen, "gene_text"))

        region_payload: dict[str, Any] = {}
        region = self._encode_region(context)
        if region:
            region_payload, _ = self._encode_get(
                ENCODE_REGION_SEARCH_URL,
                {
                    "region": region,
                    "format": "json",
                    "limit": ENCODE_MAX_RECORDS,
                },
                headers,
                urls,
                warnings,
                "region search",
                optional=True,
            )
            records.extend(self._region_records(context, region_payload, seen, region))
            notification = _clean_cell(region_payload.get("notification"))
            if notification and notification.lower().startswith("error"):
                warnings.append("Optional ENCODE region search did not return compact DCC rows; experiment searches were still used.")

        if not records and failures >= 2:
            message = "ENCODE Portal DCC lookup failed; no experiment or region-search records were returned."
            return SourceResult(
                self.spec.key,
                "failed",
                message,
                [],
                warnings=warnings,
                errors=warnings or [message],
                queried_urls=urls,
                elapsed_ms=_elapsed_ms(started),
            )

        if not records:
            records.append(self._query_context_record(context, region, target_payload, text_payload, region_payload))

        direct_count = sum(1 for record in records if record.get("category") == "regulatory_experiment")
        region_count = sum(1 for record in records if record.get("category") == "regulatory_region_hit")
        context_count = sum(1 for record in records if record.get("category") == "encode_query_context")
        parts: list[str] = []
        if direct_count:
            parts.append(f"{direct_count} experiment record(s)")
        if region_count:
            parts.append(f"{region_count} region-search record(s)")
        if context_count:
            parts.append("query context recorded")
        detail = ", ".join(parts) if parts else f"{len(records)} record(s)"
        return SourceResult(
            self.spec.key,
            "ok",
            f"Queried ENCODE Portal DCC; {detail} returned for {context.gene}.",
            records,
            warnings=warnings,
            queried_urls=urls,
            elapsed_ms=_elapsed_ms(started),
        )

    def _query_screen_legacy(self, context: KnowledgeQuery, started: float) -> SourceResult:
        try:
            payload = self.client.get_json(
                ENCODE_SEARCH_URL,
                params={"type": "Experiment", "searchTerm": context.gene, "format": "json", "limit": ENCODE_MAX_RECORDS},
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            records = [
                {
                    "category": "regulatory_experiment",
                    "source": self.spec.name,
                    "label": str(item.get("accession") or item.get("@id") or context.gene),
                    "summary": str(item.get("assay_title") or item.get("description") or "ENCODE experiment"),
                    "source_id": item.get("accession"),
                    "url": self._portal_url(item.get("@id")),
                }
                for item in self._graph_items(payload)[:ENCODE_MAX_RECORDS]
            ]
            return SourceResult(
                self.spec.key,
                "ok",
                f"Queried SCREEN/ENCODE search; {len(records)} record(s).",
                records,
                queried_urls=[ENCODE_SEARCH_URL],
                elapsed_ms=_elapsed_ms(started),
            )
        except KnowledgeRequestError as exc:
            return SourceResult(
                self.spec.key,
                "failed",
                str(exc),
                errors=[str(exc)],
                queried_urls=[ENCODE_SEARCH_URL],
                elapsed_ms=_elapsed_ms(started),
            )

    def _encode_get(
        self,
        url: str,
        params: dict[str, Any],
        headers: dict[str, str],
        urls: list[str],
        warnings: list[str],
        label: str,
        *,
        optional: bool = False,
    ) -> tuple[dict[str, Any], bool]:
        urls.append(url)
        try:
            payload = self.client.get_json(
                url,
                params=params,
                headers=headers,
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError:
            warnings.append(f"Optional ENCODE {label} failed; other ENCODE Portal results were still used.")
            return {}, True
        if not isinstance(payload, dict):
            warnings.append(f"Optional ENCODE {label} returned no usable data; other ENCODE Portal results were still used.")
            return {}, not optional
        return payload, False

    def _experiment_records(
        self,
        context: KnowledgeQuery,
        payload: dict[str, Any],
        seen: set[str],
        match_kind: str,
    ) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for item in self._graph_items(payload):
            accession = _clean_cell(item.get("accession") or item.get("@id"))
            dedupe_key = accession or _clean_cell(item.get("@id"))
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            records.append(self._experiment_record(context, item, match_kind))
            if len(records) >= ENCODE_MAX_RECORDS:
                break
        return records

    def _region_records(
        self,
        context: KnowledgeQuery,
        payload: dict[str, Any],
        seen: set[str],
        region: str,
    ) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for item in self._graph_items(payload):
            accession = _clean_cell(item.get("accession") or item.get("@id") or item.get("uuid"))
            dedupe_key = f"region:{accession}"
            if not accession or dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            records.append(self._region_record(context, item, region))
            if len(records) >= ENCODE_MAX_RECORDS:
                break
        return records

    def _experiment_record(self, context: KnowledgeQuery, item: dict[str, Any], match_kind: str) -> dict[str, Any]:
        accession = _clean_cell(item.get("accession") or item.get("@id") or context.gene)
        assay = _clean_cell(item.get("assay_title") or item.get("assay_term_name"))
        target = self._target_label(item)
        target_genes = self._target_gene_symbols(item)
        biosample = self._biosample_text(item)
        files = self._file_summary(item.get("files"))
        record = {
            "category": "regulatory_experiment",
            "source": self.spec.name,
            "label": f"{accession} - {assay}" if assay else accession,
            "summary": self._experiment_summary(context, item, match_kind, accession, assay, target, target_genes, biosample, files),
            "source_id": accession,
            "url": self._portal_url(item.get("@id") or f"/experiments/{accession}/"),
            "gene": context.gene,
            "match_type": match_kind,
            "assay": assay,
            "target": target,
            "target_genes": target_genes,
            "biosample": biosample,
            "organism": self._organism_text(item),
            "status": _clean_cell(item.get("status")),
            "lab": self._lab_text(item),
            "project": self._project_text(item),
            "date_released": _clean_cell(item.get("date_released")),
            "replication_type": _clean_cell(item.get("replication_type")),
            "biological_replicates": self._biological_replicates(item.get("replicates")),
            "files": files,
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _region_record(self, context: KnowledgeQuery, item: dict[str, Any], region: str) -> dict[str, Any]:
        accession = _clean_cell(item.get("accession") or item.get("@id") or region)
        assay = _clean_cell(item.get("assay_title") or item.get("assay_term_name") or item.get("annotation_type"))
        biosample = self._biosample_text(item)
        status = _clean_cell(item.get("status"))
        details = [f"assay {assay}" if assay else "", biosample, f"status {status}" if status else ""]
        summary = (
            f"ENCODE DCC region search returned {accession} overlapping {region} for {context.gene}: "
            f"{'; '.join(part for part in details if part)}."
        )
        return {
            "category": "regulatory_region_hit",
            "source": self.spec.name,
            "label": f"{accession} region hit",
            "summary": summary,
            "source_id": accession,
            "url": self._portal_url(item.get("@id")),
            "gene": context.gene,
            "region": region,
            "assay": assay,
            "biosample": biosample,
            "status": status,
        }

    def _query_context_record(
        self,
        context: KnowledgeQuery,
        region: str,
        target_payload: dict[str, Any],
        text_payload: dict[str, Any],
        region_payload: dict[str, Any],
    ) -> dict[str, Any]:
        target_total = self._payload_total(target_payload)
        text_total = self._payload_total(text_payload)
        notification = _clean_cell(region_payload.get("notification"))
        region_text = f" and region {region}" if region else ""
        detail = f"target-gene search total {target_total}; gene-text search total {text_total}"
        if notification:
            detail = f"{detail}; region search note: {notification}"
        return {
            "category": "encode_query_context",
            "source": self.spec.name,
            "label": f"ENCODE Portal query context for {context.gene}",
            "summary": (
                f"ENCODE Portal DCC searched released experiment metadata for {context.gene}{region_text}; "
                f"no direct experiment records were returned ({detail})."
            ),
            "source_id": f"encode:{context.gene}",
            "url": f"{ENCODE_SEARCH_URL}?type=Experiment&searchTerm={quote(context.gene)}",
            "gene": context.gene,
            "region": region,
            "target_search_total": target_total,
            "text_search_total": text_total,
            "region_search_notification": notification,
        }

    def _experiment_summary(
        self,
        context: KnowledgeQuery,
        item: dict[str, Any],
        match_kind: str,
        accession: str,
        assay: str,
        target: str,
        target_genes: list[str],
        biosample: str,
        files: dict[str, Any],
    ) -> str:
        lead = f"ENCODE DCC {accession}"
        if assay:
            lead = f"{lead} {assay}"
        if match_kind == "target_gene":
            lead = f"{lead} matched {context.gene} through target gene metadata"
        else:
            lead = f"{lead} matched {context.gene} through portal text search"

        details: list[str] = []
        target_text = target
        if target_genes:
            gene_text = ", ".join(target_genes[:4])
            target_text = f"{target_text} ({gene_text})" if target_text else gene_text
        if target_text:
            details.append(f"target {target_text}")
        if biosample:
            details.append(f"biosample {biosample}")
        organism = self._organism_text(item)
        if organism:
            details.append(organism)
        replicate_count = self._biological_replicates(item.get("replicates"))
        if replicate_count:
            details.append(f"{replicate_count} biological replicate(s)")
        file_text = self._file_text(files)
        if file_text:
            details.append(file_text)
        lab = self._lab_text(item)
        if lab:
            details.append(f"lab {lab}")
        project = self._project_text(item)
        if project:
            details.append(f"project {project}")
        status = _clean_cell(item.get("status"))
        if status:
            details.append(f"status {status}")
        released = _clean_cell(item.get("date_released"))
        if released:
            details.append(f"released {released}")
        return f"{lead}: {'; '.join(details)}." if details else f"{lead}."

    def _file_summary(self, values: Any) -> dict[str, Any]:
        files = values if isinstance(values, list) else []
        released = 0
        assemblies: list[str] = []
        output_types: list[str] = []
        formats: list[str] = []
        for file_item in files:
            if not isinstance(file_item, dict):
                continue
            if _clean_cell(file_item.get("status")) == "released":
                released += 1
            assemblies.extend(self._clean_list(file_item.get("assembly")))
            output_types.extend(self._clean_list(file_item.get("output_type") or file_item.get("file_type")))
            formats.extend(self._clean_list(file_item.get("file_format")))
        return {
            "total": len(files),
            "released": released or None,
            "assemblies": self._dedupe_text(assemblies),
            "output_types": self._dedupe_text(output_types),
            "formats": self._dedupe_text(formats),
        }

    def _file_text(self, files: dict[str, Any]) -> str:
        total = files.get("total")
        if not total:
            return ""
        parts = [f"{files.get('released') or total}/{total} released file(s)"]
        output_types = files.get("output_types") or []
        if output_types:
            parts.append(f"outputs {', '.join(output_types[:4])}")
        assemblies = files.get("assemblies") or []
        if assemblies:
            parts.append(f"assemblies {', '.join(assemblies[:3])}")
        return "; ".join(parts)

    def _target_label(self, item: dict[str, Any]) -> str:
        target = item.get("target")
        if isinstance(target, dict):
            return _clean_cell(target.get("label") or target.get("name") or target.get("title"))
        return _clean_cell(target)

    def _target_gene_symbols(self, item: dict[str, Any]) -> list[str]:
        targets = []
        target = item.get("target")
        if isinstance(target, dict):
            targets.extend(self._target_genes_from_value(target.get("genes")))
        targets.extend(self._target_genes_from_value(item.get("targets")))
        return self._dedupe_text(targets)

    def _target_genes_from_value(self, value: Any) -> list[str]:
        genes: list[str] = []
        values = value if isinstance(value, list) else [value]
        for item in values:
            if isinstance(item, dict):
                symbol = _clean_cell(item.get("symbol") or item.get("gene_symbol") or item.get("name") or item.get("label"))
                if symbol:
                    genes.append(symbol)
                genes.extend(self._target_genes_from_value(item.get("genes")))
            else:
                text = _clean_cell(item)
                if text:
                    genes.append(text)
        return genes

    def _biosample_text(self, item: dict[str, Any]) -> str:
        summary = _clean_cell(item.get("biosample_summary") or item.get("simple_biosample_summary"))
        ontology = item.get("biosample_ontology") if isinstance(item.get("biosample_ontology"), dict) else {}
        term = _clean_cell(ontology.get("term_name"))
        classification = _clean_cell(ontology.get("classification"))
        if summary:
            return summary
        if term and classification:
            return f"{term} {classification}"
        return term or classification

    def _organism_text(self, item: dict[str, Any]) -> str:
        organisms: list[str] = []
        replicates = item.get("replicates") if isinstance(item.get("replicates"), list) else []
        for replicate in replicates:
            if not isinstance(replicate, dict):
                continue
            biosample = ((replicate.get("library") or {}).get("biosample") or {}) if isinstance(replicate.get("library"), dict) else {}
            organism = biosample.get("organism") if isinstance(biosample.get("organism"), dict) else {}
            text = _clean_cell(organism.get("scientific_name"))
            if text:
                organisms.append(text)
        return ", ".join(self._dedupe_text(organisms)[:3])

    def _lab_text(self, item: dict[str, Any]) -> str:
        lab = item.get("lab")
        if isinstance(lab, dict):
            return _clean_cell(lab.get("title") or lab.get("name"))
        return _clean_cell(lab)

    def _project_text(self, item: dict[str, Any]) -> str:
        award = item.get("award")
        if isinstance(award, dict):
            return _clean_cell(award.get("project") or award.get("name"))
        return _clean_cell(award)

    def _biological_replicates(self, values: Any) -> int | None:
        if not isinstance(values, list):
            return None
        reps: set[str] = set()
        for item in values:
            if not isinstance(item, dict):
                continue
            number = _clean_cell(item.get("biological_replicate_number"))
            if number:
                reps.add(number)
        return len(reps) if reps else None

    def _payload_total(self, payload: dict[str, Any]) -> int:
        try:
            return int(payload.get("total") or 0)
        except (TypeError, ValueError):
            return 0

    def _graph_items(self, payload: Any) -> list[dict[str, Any]]:
        if not isinstance(payload, dict):
            return []
        graph = payload.get("@graph")
        if not isinstance(graph, list):
            return []
        return [item for item in graph if isinstance(item, dict)]

    def _encode_region(self, context: KnowledgeQuery) -> str:
        text = _clean_cell(context.region).replace(",", "")
        match = re.search(r"(?P<chrom>(?:chr)?[A-Za-z0-9_.]+)\s*:\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)", text)
        if not match:
            return ""
        chrom = match.group("chrom")
        if not chrom.lower().startswith("chr"):
            chrom = f"chr{chrom}"
        return f"{chrom}:{match.group('start')}-{match.group('end')}"

    def _portal_url(self, path: Any) -> str:
        text = _clean_cell(path)
        if not text:
            return self.spec.homepage or ENCODE_PORTAL_BASE_URL
        if text.startswith(("http://", "https://")):
            return text
        return f"{ENCODE_PORTAL_BASE_URL}{text if text.startswith('/') else f'/{text}'}"

    def _clean_list(self, value: Any) -> list[str]:
        if value is None:
            return []
        values = value if isinstance(value, list) else [value]
        return [text for text in (_clean_cell(item) for item in values) if text]

    def _dedupe_text(self, values: list[str]) -> list[str]:
        deduped: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in deduped:
                deduped.append(text)
        return deduped


class EwasCatalogConnector(BaseConnector):
    """MRC-IEU EWAS Catalog query connector."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        warnings: list[str] = []
        urls: list[str] = [EWAS_CATALOG_BASE_URL]
        if not self._is_hg19(context):
            warnings.append(
                f"EWAS Catalog reports CpG locations in hg19 coordinates; input build {context.genome_build or 'unknown'} was used as provided."
            )

        try:
            page = self.client.get_text(
                f"{EWAS_CATALOG_BASE_URL}/",
                params={"gene": context.gene},
                headers={"Accept": "text/html"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=urls, started=started)

        tsv_url = self._download_url(page)
        records: list[dict[str, Any]] = []
        if tsv_url:
            urls.append(tsv_url)
            try:
                tsv_text = self.client.get_text(
                    tsv_url,
                    headers={"Accept": "text/tab-separated-values,text/plain,*/*"},
                    rate_limit_per_second=self.spec.rate_limit_per_second,
                )
                records = self._records_from_tsv(context, tsv_text)
            except KnowledgeRequestError:
                warnings.append("Optional EWAS Catalog query-specific TSV download failed; query context was still recorded.")

        if not records:
            records.append(self._context_record(context, bool(tsv_url)))

        count = sum(1 for record in records if record.get("category") == "ewas_association")
        if count:
            message = f"Queried EWAS Catalog gene results for {context.gene}; {count} compact association record(s) returned."
        else:
            message = f"EWAS Catalog query context recorded for {context.gene}; no compact association rows were parsed."
        return SourceResult(
            self.spec.key,
            "ok",
            message,
            records,
            warnings=warnings,
            queried_urls=urls,
            elapsed_ms=_elapsed_ms(started),
        )

    def _download_url(self, page: str) -> str:
        parser = _EwasCatalogDownloadParser()
        parser.feed(page or "")
        return urljoin(EWAS_CATALOG_BASE_URL, parser.tsv_href) if parser.tsv_href else ""

    def _records_from_tsv(self, context: KnowledgeQuery, tsv_text: str) -> list[dict[str, Any]]:
        reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
        rows = [row for row in reader if any(_clean_cell(value) for value in row.values())]
        rows.sort(key=self._catalog_sort_key)
        records: list[dict[str, Any]] = []
        seen: set[str] = set()
        for row in rows:
            cpg = self._row_get(row, "CpG", "cpg")
            pmid = self._row_get(row, "PMID", "pmid")
            exposure = self._row_get(row, "Exposure", "exposure", "Trait", "trait")
            key = f"{cpg}:{pmid}:{exposure}:{self._row_get(row, 'P', 'p')}"
            if not cpg or key in seen:
                continue
            seen.add(key)
            records.append(self._association_record(context, row))
            if len(records) >= EWAS_CATALOG_MAX_RECORDS:
                break
        return records

    def _association_record(self, context: KnowledgeQuery, row: dict[str, Any]) -> dict[str, Any]:
        cpg = self._row_get(row, "CpG", "cpg")
        exposure = self._row_get(row, "Exposure", "exposure", "Trait", "trait")
        outcome = self._row_get(row, "Outcome", "outcome")
        tissue = self._row_get(row, "Tissue", "tissue")
        p_value = self._row_get(row, "P", "p")
        beta = self._row_get(row, "Beta", "beta")
        location = self._row_get(row, "Location", "location")
        pmid = self._row_get(row, "PMID", "pmid")
        author = self._row_get(row, "Author", "author")
        study_id = self._row_get(row, "StudyID", "studyid", "StudyId", "study_id")
        label_bits = [cpg]
        if exposure:
            label_bits.append(exposure)
        details: list[str] = []
        if outcome:
            details.append(f"outcome {outcome}")
        if exposure:
            details.append(f"exposure/trait {exposure}")
        if tissue:
            details.append(f"tissue {tissue}")
        n_value = self._row_get(row, "N", "n")
        if n_value:
            details.append(f"N={self._format_int(n_value)}")
        if beta:
            details.append(f"beta {beta}")
        if p_value:
            details.append(f"p={p_value}")
        citation = ", ".join(part for part in (f"PMID {pmid}" if pmid else "", author) if part)
        if citation:
            details.append(citation)
        analysis = self._row_get(row, "Analysis", "analysis")
        if analysis:
            details.append(f"analysis {analysis}")
        location_text = f" at {location}" if location else ""
        summary = f"EWAS Catalog {context.gene} CpG {cpg}{location_text}: {'; '.join(details)}."
        record = {
            "category": "ewas_association",
            "source": self.spec.name,
            "label": " - ".join(label_bits),
            "summary": summary,
            "source_id": study_id or f"{pmid}:{cpg}",
            "url": f"{EWAS_CATALOG_BASE_URL}/?cpg={quote(cpg)}" if cpg else f"{EWAS_CATALOG_BASE_URL}/?gene={quote(context.gene)}",
            "gene": context.gene,
            "cpg": cpg,
            "probe_id": cpg,
            "location": location,
            "tissue": tissue,
            "outcome": outcome,
            "exposure": exposure,
            "beta": beta,
            "p_value": p_value,
            "n": n_value,
            "pmid": pmid,
            "author": author,
            "study_id": study_id,
            "analysis": analysis,
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _context_record(self, context: KnowledgeQuery, has_download: bool) -> dict[str, Any]:
        availability = "generated a query-specific TSV link" if has_download else "did not expose a query-specific TSV link"
        return {
            "category": "ewas_catalog_context",
            "source": self.spec.name,
            "label": f"EWAS Catalog query context for {context.gene}",
            "summary": (
                f"EWAS Catalog was queried for {context.gene}; the website {availability}. "
                "Catalog rows include CpG, hg19 location, gene, beta, standard error, p-value, study metadata, tissue, "
                "sample size, exposure/outcome, PMID, and analysis fields when available."
            ),
            "source_id": f"ewas_catalog:{context.gene}",
            "url": f"{EWAS_CATALOG_BASE_URL}/?gene={quote(context.gene)}",
            "gene": context.gene,
        }

    def _catalog_sort_key(self, row: dict[str, Any]) -> tuple[float, str]:
        return (self._float_sort_value(self._row_get(row, "P", "p")), self._row_get(row, "CpG", "cpg"))

    def _row_get(self, row: dict[str, Any], *keys: str) -> str:
        normalized = {str(key).strip().lower(): value for key, value in row.items()}
        for key in keys:
            value = _clean_cell(normalized.get(key.lower()))
            if value and value.upper() != "NA":
                return value
        return ""

    def _format_int(self, value: Any) -> str:
        try:
            return f"{int(str(value).replace(',', '').strip()):,}"
        except (TypeError, ValueError):
            return _clean_cell(value)

    def _float_sort_value(self, value: Any) -> float:
        try:
            return float(str(value).replace(",", "").strip())
        except (TypeError, ValueError):
            return float("inf")

    def _is_hg19(self, context: KnowledgeQuery) -> bool:
        build = _clean_cell(context.genome_build).lower()
        return not build or "19" in build or "37" in build or "grch37" in build or "hg19" in build


class EwasAtlasConnector(BaseConnector):
    """EWAS Atlas / EWAS Open Platform REST connector."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        warnings: list[str] = []
        urls: list[str] = []
        if not self._is_hg19(context):
            warnings.append(
                f"EWAS Atlas REST reports CpG locations as chrHg19/posHg19; input build {context.genome_build or 'unknown'} was used as provided."
            )

        endpoint, params, scope = self._query_endpoint(context)
        if not endpoint:
            record = self._context_record(context, "No gene, region, or probe query target was available.")
            return SourceResult(
                self.spec.key,
                "metadata_only",
                "EWAS Atlas context recorded; no compact query target was available.",
                [record],
                warnings=warnings,
                queried_urls=[EWAS_ATLAS_PORTAL_BASE],
                elapsed_ms=_elapsed_ms(started),
            )

        urls.append(endpoint)
        try:
            payload = self.client.get_json(
                endpoint,
                params=params,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=urls, started=started)

        records = self._records_from_payload(context, payload, scope)
        if not records:
            records.append(self._context_record(context, f"No association rows were returned for {scope}."))

        count = sum(1 for record in records if record.get("category") == "ewas_association")
        if count:
            message = f"Queried EWAS Atlas REST {scope}; {count} compact association record(s) returned for {context.gene}."
        else:
            message = f"EWAS Atlas REST context recorded for {context.gene}; no compact association records were parsed."
        return SourceResult(
            self.spec.key,
            "ok",
            message,
            records,
            warnings=warnings,
            queried_urls=urls,
            elapsed_ms=_elapsed_ms(started),
        )

    def _query_endpoint(self, context: KnowledgeQuery) -> tuple[str, dict[str, Any], str]:
        region = self._parse_region(context.region)
        if region:
            endpoint = f"{EWAS_ATLAS_REST_BASE}/pos"
            params = {"chr": region["chrom"], "start": region["start"], "end": region["end"]}
            return endpoint, params, f"position query hg19 chr{region['chrom']}:{region['start']}-{region['end']}"
        probe_id = self._first_probe(context)
        if probe_id:
            return f"{EWAS_ATLAS_REST_BASE}/probe", {"probeId": probe_id}, f"probe query {probe_id}"
        if context.gene:
            return f"{EWAS_ATLAS_REST_BASE}/gene", {"geneSymbol": context.gene}, f"gene query {context.gene}"
        return "", {}, ""

    def _records_from_payload(self, context: KnowledgeQuery, payload: Any, scope: str) -> list[dict[str, Any]]:
        probes = self._probe_items(payload)
        associations: list[tuple[dict[str, Any], dict[str, Any]]] = []
        for probe in probes:
            assoc_list = probe.get("associationList") if isinstance(probe.get("associationList"), list) else []
            for assoc in assoc_list:
                if isinstance(assoc, dict):
                    associations.append((probe, assoc))
        associations.sort(key=lambda item: (self._rank_sort_value(item[1].get("rank")), _clean_cell(item[0].get("probeId"))))

        records: list[dict[str, Any]] = []
        seen: set[str] = set()
        for probe, association in associations:
            key = ":".join(
                _clean_cell(value)
                for value in (
                    probe.get("probeId"),
                    association.get("studyId"),
                    association.get("trait"),
                    association.get("pmid"),
                )
            )
            if key in seen:
                continue
            seen.add(key)
            records.append(self._association_record(context, probe, association, scope))
            if len(records) >= EWAS_ATLAS_MAX_RECORDS:
                break
        return records

    def _probe_items(self, payload: Any) -> list[dict[str, Any]]:
        if not isinstance(payload, dict):
            return []
        data = payload.get("data")
        if isinstance(data, list):
            return [item for item in data if isinstance(item, dict)]
        if isinstance(data, dict):
            probe_list = data.get("probeList")
            if isinstance(probe_list, list):
                return [item for item in probe_list if isinstance(item, dict)]
            if data.get("probeId"):
                return [data]
        return []

    def _association_record(
        self,
        context: KnowledgeQuery,
        probe: dict[str, Any],
        association: dict[str, Any],
        scope: str,
    ) -> dict[str, Any]:
        probe_id = _clean_cell(probe.get("probeId"))
        chrom = _clean_cell(probe.get("chrHg19"))
        pos = _clean_cell(probe.get("posHg19"))
        location = f"hg19 chr{chrom}:{pos}" if chrom and pos else ""
        cpg_island = _clean_cell(probe.get("cpgIsland"))
        trait = _clean_cell(association.get("trait"))
        correlation = _clean_cell(association.get("correlation"))
        rank = _clean_cell(association.get("rank"))
        study_id = _clean_cell(association.get("studyId"))
        pmid = _clean_cell(association.get("pmid"))
        transcripts = self._transcripts(probe.get("relatedTranscription"))
        details: list[str] = []
        if cpg_island:
            details.append(cpg_island)
        transcript_text = self._transcript_text(transcripts)
        if transcript_text:
            details.append(transcript_text)
        if trait:
            details.append(f"trait {trait}")
        direction = self._correlation_text(correlation)
        if direction:
            details.append(direction)
        if rank:
            details.append(f"rank {rank}")
        if study_id:
            details.append(f"study {study_id}")
        if pmid:
            details.append(f"PMID {pmid}")
        details.append(scope)
        summary = f"EWAS Atlas {context.gene} probe {probe_id}"
        if location:
            summary = f"{summary} at {location}"
        if details:
            summary = f"{summary}: {'; '.join(details)}."
        else:
            summary = f"{summary}."
        record = {
            "category": "ewas_association",
            "source": self.spec.name,
            "label": f"{probe_id} - {trait}" if trait else probe_id,
            "summary": summary,
            "source_id": f"{study_id}:{probe_id}" if study_id else probe_id,
            "url": f"{EWAS_ATLAS_PORTAL_BASE}/search?item={quote(probe_id)}&term=Probe+Id",
            "gene": context.gene,
            "probe_id": probe_id,
            "cpg": probe_id,
            "chr_hg19": chrom,
            "pos_hg19": pos,
            "location": location,
            "cpg_island": cpg_island,
            "trait": trait,
            "correlation": correlation,
            "rank": rank,
            "study_id": study_id,
            "pmid": pmid,
            "transcripts": transcripts,
            "query_scope": scope,
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _context_record(self, context: KnowledgeQuery, note: str) -> dict[str, Any]:
        return {
            "category": "ewas_atlas_context",
            "source": self.spec.name,
            "label": f"EWAS Atlas query context for {context.gene}",
            "summary": (
                f"EWAS Atlas / EWAS Open Platform was queried for {context.gene}. "
                "Its REST API can return probe, gene, position, study, and publication objects with curated EWAS "
                f"associations and causality records. {note}"
            ),
            "source_id": f"ewas_atlas:{context.gene}",
            "url": f"{EWAS_ATLAS_PORTAL_BASE}/search?item={quote(context.gene)}&term=Gene",
            "gene": context.gene,
        }

    def _parse_region(self, text: str) -> dict[str, int | str]:
        match = re.search(r"(?P<chrom>(?:chr)?[A-Za-z0-9_.]+)\s*:\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)", _clean_cell(text).replace(",", ""))
        if not match:
            return {}
        chrom = match.group("chrom").removeprefix("chr").removeprefix("Chr")
        try:
            start = int(match.group("start"))
            end = int(match.group("end"))
        except ValueError:
            return {}
        if start > end:
            start, end = end, start
        return {"chrom": chrom, "start": start, "end": end}

    def _first_probe(self, context: KnowledgeQuery) -> str:
        for locus in context.epigenetic_loci:
            probe_id = _clean_cell(getattr(locus, "probe_id", ""))
            if probe_id:
                return probe_id
        return ""

    def _transcripts(self, value: Any) -> list[dict[str, Any]]:
        if not isinstance(value, list):
            return []
        transcripts: list[dict[str, Any]] = []
        seen: set[str] = set()
        for item in value:
            if not isinstance(item, dict):
                continue
            gene_name = _clean_cell(item.get("geneName"))
            transcript_id = _clean_cell(item.get("ensemblTranscriptId"))
            key = f"{gene_name}:{transcript_id}:{_clean_cell(item.get('posToTss'))}"
            if not gene_name or key in seen:
                continue
            seen.add(key)
            transcripts.append(
                {
                    "gene": gene_name,
                    "transcript": transcript_id,
                    "pos_to_tss": item.get("posToTss"),
                }
            )
        return transcripts

    def _transcript_text(self, transcripts: list[dict[str, Any]]) -> str:
        for transcript in transcripts:
            transcript_id = _clean_cell(transcript.get("transcript"))
            pos_to_tss = _clean_cell(transcript.get("pos_to_tss"))
            gene = _clean_cell(transcript.get("gene"))
            if transcript_id and pos_to_tss:
                return f"{gene} transcript {transcript_id}, {pos_to_tss} bp from TSS"
            if transcript_id:
                return f"{gene} transcript {transcript_id}"
        return ""

    def _correlation_text(self, value: str) -> str:
        lower = value.lower()
        if lower == "pos":
            return "positive methylation-trait correlation"
        if lower == "neg":
            return "negative methylation-trait correlation"
        if value and lower != "na":
            return f"correlation {value}"
        return ""

    def _rank_sort_value(self, value: Any) -> int:
        try:
            return int(value)
        except (TypeError, ValueError):
            return 10**12

    def _is_hg19(self, context: KnowledgeQuery) -> bool:
        build = _clean_cell(context.genome_build).lower()
        return not build or "19" in build or "37" in build or "grch37" in build or "hg19" in build


class ScreenConnector(BaseConnector):
    """SCREEN cCRE Registry connector."""

    CCRE_QUERY = """
    query ScreenCcreByRegion($assembly: String!, $coordinates: [GenomicRangeInput!]) {
      cCREQuery(assembly: $assembly, coordinates: $coordinates) {
        accession
        assembly
        rDHS
        group
        ctcf_bound
        dnaseMax: maxZ(assay: "DNase")
        h3k4me3Max: maxZ(assay: "H3K4me3")
        h3k27acMax: maxZ(assay: "H3K27ac")
        ctcfMax: maxZ(assay: "CTCF")
        coordinates {
          chromosome
          start
          end
        }
        nearby_genes {
          intersecting_genes {
            id
            name
            gene_type
            strand
          }
        }
      }
    }
    """

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        warnings: list[str] = []
        region = self._query_region(context)
        if not region:
            record = self._context_record(
                context,
                "",
                "No compact coordinate window was available for a SCREEN cCRE overlap query.",
            )
            return SourceResult(
                self.spec.key,
                "metadata_only",
                "SCREEN cCRE Registry context recorded; no region was available for cCRE overlap lookup.",
                [record],
                queried_urls=[SCREEN_PORTAL_BASE_URL],
                elapsed_ms=_elapsed_ms(started),
            )

        assembly = "GRCh38"
        if not self._is_hg38(context):
            warnings.append(
                f"SCREEN cCRE Registry is queried in GRCh38; provided build {context.genome_build or 'unknown'} may require liftover."
            )

        records: list[dict[str, Any]] = []
        payload: dict[str, Any] = {}
        try:
            payload = self.client.post_json(
                SCREEN_GRAPHQL_URL,
                json_payload={
                    "query": self.CCRE_QUERY,
                    "variables": {
                        "assembly": assembly,
                        "coordinates": [
                            {
                                "chromosome": region["chromosome"],
                                "start": region["query_start"],
                                "end": region["query_end"],
                            }
                        ],
                    },
                },
                headers={"Accept": "application/json", "Content-Type": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            records = self._ccre_records(context, payload, region)
            if self._graphql_errors(payload):
                warnings.append("Optional SCREEN cCRE GraphQL query returned errors; registry context was still recorded.")
        except KnowledgeRequestError:
            warnings.append("Optional SCREEN cCRE GraphQL query failed; registry context was still recorded.")

        total_hits = self._ccre_total(payload)
        if not records:
            records.append(
                self._context_record(
                    context,
                    region["display_region"],
                    "No overlapping cCRE records were returned by the compact SCREEN query.",
                )
            )

        ccre_count = sum(1 for record in records if record.get("category") == "candidate_regulatory_element")
        if ccre_count:
            if total_hits > ccre_count:
                detail = f"{ccre_count} of {total_hits} overlapping cCRE record(s)"
            else:
                detail = f"{ccre_count} overlapping cCRE record(s)"
            message = f"Queried SCREEN cCRE Registry for {region['display_region']}; {detail} returned for {context.gene}."
            status = "ok"
        else:
            message = f"SCREEN cCRE Registry context recorded for {context.gene} at {region['display_region']}."
            status = "ok"

        return SourceResult(
            self.spec.key,
            status,
            message,
            records,
            warnings=warnings,
            queried_urls=[SCREEN_GRAPHQL_URL],
            elapsed_ms=_elapsed_ms(started),
        )

    def _ccre_records(self, context: KnowledgeQuery, payload: Any, region: dict[str, Any]) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        if not isinstance(payload, dict):
            return records
        data = payload.get("data")
        if not isinstance(data, dict):
            return records
        items = data.get("cCREQuery")
        if not isinstance(items, list):
            return records

        seen: set[str] = set()
        for item in items:
            if not isinstance(item, dict):
                continue
            accession = _clean_cell(item.get("accession") or item.get("id"))
            if not accession or accession in seen:
                continue
            seen.add(accession)
            records.append(self._ccre_record(context, item, region))
            if len(records) >= SCREEN_MAX_RECORDS:
                break
        return records

    def _ccre_record(self, context: KnowledgeQuery, item: dict[str, Any], region: dict[str, Any]) -> dict[str, Any]:
        accession = _clean_cell(item.get("accession") or item.get("id"))
        group = _clean_cell(item.get("group"))
        group_label = self._group_label(group)
        location = self._display_location(item.get("coordinates"))
        genes = self._nearby_genes(item.get("nearby_genes"))
        max_scores = self._max_scores(item)
        ctcf_bound = item.get("ctcf_bound") if isinstance(item.get("ctcf_bound"), bool) else None
        label_suffix = group if group else "cCRE"
        record = {
            "category": "candidate_regulatory_element",
            "source": self.spec.name,
            "label": f"{accession} - {label_suffix}",
            "summary": self._ccre_summary(context, item, region, accession, group_label, location, genes, max_scores, ctcf_bound),
            "source_id": accession,
            "url": self._screen_url(accession),
            "gene": context.gene,
            "region": region["display_region"],
            "assembly": _clean_cell(item.get("assembly")) or "GRCh38",
            "accession": accession,
            "rdhs": _clean_cell(item.get("rDHS")),
            "ccre_group": group,
            "ccre_group_label": group_label,
            "ctcf_bound": ctcf_bound,
            "coordinates": item.get("coordinates") if isinstance(item.get("coordinates"), dict) else {},
            "location": location,
            "nearby_genes": genes,
            "max_z_scores": max_scores,
            "registry_download": SCREEN_GRCH38_CCRE_BED_URL,
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _ccre_summary(
        self,
        context: KnowledgeQuery,
        item: dict[str, Any],
        region: dict[str, Any],
        accession: str,
        group_label: str,
        location: str,
        genes: list[dict[str, str]],
        max_scores: dict[str, float],
        ctcf_bound: bool | None,
    ) -> str:
        lead = f"SCREEN cCRE {accession} overlaps {context.gene} query window"
        if location:
            lead = f"{lead} at {location}"
        details: list[str] = []
        group = _clean_cell(item.get("group"))
        if group_label and group:
            details.append(f"{group_label} ({group})")
        elif group:
            details.append(group)
        if ctcf_bound is not None:
            details.append("CTCF-bound" if ctcf_bound else "not CTCF-bound")
        gene_text = self._gene_text(genes)
        if gene_text:
            details.append(f"intersecting genes: {gene_text}")
        score_text = self._score_text(max_scores)
        if score_text:
            details.append(f"max assay Z-scores {score_text}")
        rdhs = _clean_cell(item.get("rDHS"))
        if rdhs:
            details.append(f"rDHS {rdhs}")
        details.append(f"query window {region['display_region']}")
        return f"{lead}: {'; '.join(details)}."

    def _context_record(self, context: KnowledgeQuery, region: str, note: str) -> dict[str, Any]:
        region_text = f" at {region}" if region else ""
        return {
            "category": "screen_registry_context",
            "source": self.spec.name,
            "label": f"SCREEN cCRE Registry context for {context.gene}",
            "summary": (
                f"SCREEN provides the ENCODE candidate cis-regulatory element registry for human GRCh38; "
                f"{context.gene}{region_text} can be queried against cCRE classes such as PLS, pELS, dELS, "
                f"DNase-H3K4me3, and CTCF-only. {note}"
            ),
            "source_id": f"screen:{context.gene}",
            "url": self._screen_region_url(region) if region else SCREEN_PORTAL_BASE_URL,
            "gene": context.gene,
            "region": region,
            "assembly": "GRCh38",
            "registry_download": SCREEN_GRCH38_CCRE_BED_URL,
        }

    def _query_region(self, context: KnowledgeQuery) -> dict[str, Any]:
        text = _clean_cell(context.region).replace(",", "")
        match = re.search(r"(?P<chrom>(?:chr)?[A-Za-z0-9_.]+)\s*:\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)", text)
        if match:
            chrom = match.group("chrom")
            if not chrom.lower().startswith("chr"):
                chrom = f"chr{chrom}"
            start = int(match.group("start"))
            end = int(match.group("end"))
            if start > end:
                start, end = end, start
            return {
                "chromosome": chrom,
                "query_start": max(0, start - 1),
                "query_end": end,
                "display_region": f"{chrom}:{start}-{end}",
            }

        variant = _first_variant(context)
        chrom = _clean_cell(getattr(variant, "chrom", "") if variant else "")
        pos = getattr(variant, "pos", None) if variant else None
        if chrom and pos:
            if not chrom.lower().startswith("chr"):
                chrom = f"chr{chrom}"
            try:
                position = int(pos)
            except (TypeError, ValueError):
                return {}
            return {
                "chromosome": chrom,
                "query_start": max(0, position - 1),
                "query_end": position,
                "display_region": f"{chrom}:{position}-{position}",
            }
        return {}

    def _display_location(self, coordinates: Any) -> str:
        if not isinstance(coordinates, dict):
            return ""
        chrom = _clean_cell(coordinates.get("chromosome"))
        start = self._int_value(coordinates.get("start"))
        end = self._int_value(coordinates.get("end"))
        if not chrom or start is None or end is None:
            return ""
        return f"GRCh38 {chrom}:{start + 1}-{end}"

    def _nearby_genes(self, value: Any) -> list[dict[str, str]]:
        if not isinstance(value, dict):
            return []
        genes = value.get("intersecting_genes")
        if not isinstance(genes, list):
            return []
        deduped: list[dict[str, str]] = []
        seen: set[str] = set()
        for gene in genes:
            if not isinstance(gene, dict):
                continue
            name = _clean_cell(gene.get("name"))
            gene_id = _clean_cell(gene.get("id"))
            key = name or gene_id
            if not key or key in seen:
                continue
            seen.add(key)
            deduped.append(
                {
                    "name": name,
                    "id": gene_id,
                    "gene_type": _clean_cell(gene.get("gene_type")),
                    "strand": _clean_cell(gene.get("strand")),
                }
            )
        return deduped

    def _max_scores(self, item: dict[str, Any]) -> dict[str, float]:
        scores: dict[str, float] = {}
        for key, label in (
            ("dnaseMax", "DNase"),
            ("h3k4me3Max", "H3K4me3"),
            ("h3k27acMax", "H3K27ac"),
            ("ctcfMax", "CTCF"),
        ):
            value = self._float_value(item.get(key))
            if value is not None:
                scores[label] = value
        return scores

    def _group_label(self, group: str) -> str:
        labels = {
            "PLS": "promoter-like signature",
            "pELS": "proximal enhancer-like signature",
            "dELS": "distal enhancer-like signature",
            "DNase-H3K4me3": "DNase-H3K4me3 signature",
            "CTCF-only": "CTCF-only element",
        }
        return labels.get(group, group)

    def _gene_text(self, genes: list[dict[str, str]]) -> str:
        parts: list[str] = []
        for gene in genes[:4]:
            name = gene.get("name") or gene.get("id") or ""
            gene_type = gene.get("gene_type") or ""
            if name and gene_type:
                parts.append(f"{name} ({gene_type})")
            elif name:
                parts.append(name)
        return ", ".join(parts)

    def _score_text(self, scores: dict[str, float]) -> str:
        return ", ".join(f"{label} {value:.2f}" for label, value in scores.items())

    def _ccre_total(self, payload: Any) -> int:
        if not isinstance(payload, dict):
            return 0
        data = payload.get("data")
        if not isinstance(data, dict):
            return 0
        items = data.get("cCREQuery")
        return len(items) if isinstance(items, list) else 0

    def _graphql_errors(self, payload: Any) -> bool:
        return isinstance(payload, dict) and bool(payload.get("errors"))

    def _screen_url(self, accession: str) -> str:
        return f"{SCREEN_PORTAL_BASE_URL}/search/?q={quote(accession)}&assembly=GRCh38"

    def _screen_region_url(self, region: str) -> str:
        return f"{SCREEN_PORTAL_BASE_URL}/search/?q={quote(region)}&assembly=GRCh38"

    def _is_hg38(self, context: KnowledgeQuery) -> bool:
        build = _clean_cell(context.genome_build).lower()
        return not build or "38" in build or "grch38" in build or "hg38" in build or "b38" in build

    def _int_value(self, value: Any) -> int | None:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _float_value(self, value: Any) -> float | None:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None


class BioStudiesConnector(BaseConnector):
    """BioStudies search connector."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        url = "https://www.ebi.ac.uk/biostudies/api/v1/search"
        try:
            payload = self.client.get_json(url, params={"query": context.gene, "pageSize": 5})
            records = [
                {
                    "category": "public_dataset",
                    "source": self.spec.name,
                    "label": str(item.get("title") or item.get("accession") or context.gene),
                    "summary": str(item.get("description") or ""),
                    "source_id": item.get("accession"),
                    "url": f"https://www.ebi.ac.uk/biostudies/studies/{item.get('accession', '')}",
                }
                for item in payload.get("hits", [])[:5]
            ]
            return SourceResult(self.spec.key, "ok", f"Queried BioStudies; {len(records)} record(s).", records, queried_urls=[url], elapsed_ms=_elapsed_ms(started))
        except KnowledgeRequestError as exc:
            return SourceResult(self.spec.key, "failed", str(exc), errors=[str(exc)], queried_urls=[url], elapsed_ms=_elapsed_ms(started))


class GraphqlConnector(BaseConnector):
    """Best-effort GraphQL connector for gnomAD and CIViC."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        if self.spec.connector_kind == "gnomad":
            return self._query_gnomad(context, started)

        return self._query_civic(context, started)

    def _query_gnomad(self, context: KnowledgeQuery, started: float) -> SourceResult:
        variant = _first_variant(context)
        variant_id = self._gnomad_variant_id(variant)
        rsid = _clean_cell(variant.rsid if variant else _first_rsid(context))
        has_variant_lookup = bool(variant_id or rsid)
        variables = {
            "dataset": GNOMAD_DATASET,
            "geneSymbol": context.gene,
            "referenceGenome": self._gnomad_reference_genome(context),
        }
        if has_variant_lookup:
            variables.update({"variantId": variant_id or None, "rsid": None if variant_id else rsid})
            query = self._gnomad_variant_gene_query()
        else:
            query = self._gnomad_gene_query()

        try:
            payload = self.client.post_json(
                GNOMAD_API_URL,
                json_payload={"query": query, "variables": variables},
                headers={"Accept": "application/json", "Content-Type": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError as exc:
            return SourceResult(
                self.spec.key,
                "failed",
                str(exc),
                errors=[str(exc)],
                queried_urls=[GNOMAD_API_URL],
                elapsed_ms=_elapsed_ms(started),
            )

        errors = payload.get("errors") if isinstance(payload, dict) else None
        if errors:
            message = self._graphql_error_message(errors, "gnomAD GraphQL query failed.", "gnomAD")
            return SourceResult(
                self.spec.key,
                "failed",
                message,
                errors=[message],
                queried_urls=[GNOMAD_API_URL],
                elapsed_ms=_elapsed_ms(started),
            )

        data = payload.get("data", {}) if isinstance(payload, dict) else {}
        if not isinstance(data, dict):
            data = {}
        meta = data.get("meta") if isinstance(data.get("meta"), dict) else {}
        variant_payload = data.get("variant") if isinstance(data.get("variant"), dict) else {}
        gene_payload = data.get("gene") if isinstance(data.get("gene"), dict) else {}

        records: list[dict[str, Any]] = []
        if variant_payload:
            records.append(self._gnomad_variant_record(context, variant, variant_payload, gene_payload, meta))
        if gene_payload:
            records.append(self._gnomad_gene_record(context, gene_payload, meta))

        lookup_label = variant_id or rsid
        if not records:
            if has_variant_lookup:
                message = f"Queried gnomAD; no variant or gene records found for {lookup_label} / {context.gene}."
            else:
                message = f"Queried gnomAD; no gene constraint record found for {context.gene}."
        elif variant_payload and gene_payload:
            message = f"Queried gnomAD; variant frequency and gene constraint records returned for {lookup_label}."
        elif variant_payload:
            message = f"Queried gnomAD; variant frequency record returned for {lookup_label}."
        else:
            message = f"Queried gnomAD; gene constraint record returned for {context.gene}."
        return SourceResult(
            self.spec.key,
            "ok",
            message,
            records,
            queried_urls=[GNOMAD_API_URL],
            elapsed_ms=_elapsed_ms(started),
        )

    def _gnomad_variant_gene_query(self) -> str:
        return """
        query GnomadVariantAndGene(
          $variantId: String,
          $rsid: String,
          $dataset: DatasetId!,
          $geneSymbol: String!,
          $referenceGenome: ReferenceGenomeId!
        ) {
          meta { clinvar_release_date }
          variant(variantId: $variantId, rsid: $rsid, dataset: $dataset) {
            variantId
            variant_id
            reference_genome
            chrom
            pos
            ref
            alt
            caid
            rsid
            rsids
            flags
            exome {
              ac an af homozygote_count hemizygote_count ac_hom ac_hemi filters flags
              faf95 { popmax popmax_population }
              faf99 { popmax popmax_population }
              populations { id ac an homozygote_count hemizygote_count ac_hom ac_hemi }
            }
            genome {
              ac an af homozygote_count hemizygote_count ac_hom ac_hemi filters flags
              faf95 { popmax popmax_population }
              faf99 { popmax popmax_population }
              populations { id ac an homozygote_count hemizygote_count ac_hom ac_hemi }
            }
            joint {
              ac an homozygote_count hemizygote_count filters
              faf95 { popmax popmax_population }
              faf99 { popmax popmax_population }
              populations { id ac an homozygote_count hemizygote_count ac_hom ac_hemi }
            }
            transcript_consequences {
              gene_id
              gene_symbol
              transcript_id
              transcript_version
              is_canonical
              is_mane_select
              major_consequence
              consequence_terms
              hgvsc
              hgvsp
              hgvs
              lof
              lof_filter
              lof_flags
              polyphen_prediction
              sift_prediction
              domains
            }
            in_silico_predictors { id value flags }
            lof_curations { gene_id gene_symbol verdict flags project }
            non_coding_constraint { chrom start stop element_id possible observed expected oe z }
          }
          gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
            gene_id
            gene_version
            symbol
            gencode_symbol
            hgnc_id
            ncbi_id
            omim_id
            name
            reference_genome
            chrom
            start
            stop
            strand
            canonical_transcript_id
            mane_select_transcript { ensembl_id ensembl_version refseq_id refseq_version }
            flags
            gnomad_constraint {
              exp_lof exp_mis exp_syn
              obs_lof obs_mis obs_syn
              oe_lof oe_lof_lower oe_lof_upper oe_lof_percentile
              oe_mis oe_mis_lower oe_mis_upper
              oe_syn oe_syn_lower oe_syn_upper
              lof_z mis_z syn_z pli pLI flags
            }
          }
        }
        """

    def _gnomad_gene_query(self) -> str:
        return """
        query GnomadGene($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
          meta { clinvar_release_date }
          gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
            gene_id
            gene_version
            symbol
            gencode_symbol
            hgnc_id
            ncbi_id
            omim_id
            name
            reference_genome
            chrom
            start
            stop
            strand
            canonical_transcript_id
            mane_select_transcript { ensembl_id ensembl_version refseq_id refseq_version }
            flags
            gnomad_constraint {
              exp_lof exp_mis exp_syn
              obs_lof obs_mis obs_syn
              oe_lof oe_lof_lower oe_lof_upper oe_lof_percentile
              oe_mis oe_mis_lower oe_mis_upper
              oe_syn oe_syn_lower oe_syn_upper
              lof_z mis_z syn_z pli pLI flags
            }
          }
        }
        """

    def _gnomad_variant_record(
        self,
        context: KnowledgeQuery,
        input_variant: Any,
        item: dict[str, Any],
        gene: dict[str, Any],
        meta: dict[str, Any],
    ) -> dict[str, Any]:
        variant_id = _clean_cell(item.get("variantId") or item.get("variant_id"))
        rsids = self._gnomad_clean_list(item.get("rsids"))
        rsid = _clean_cell(item.get("rsid")) or (rsids[0] if rsids else "")
        transcript_consequences = self._gnomad_transcript_consequences(item.get("transcript_consequences"), context)
        primary_transcript = self._gnomad_primary_transcript(transcript_consequences, context)
        frequencies = {
            "joint": self._gnomad_frequency_block(item.get("joint"), compute_af=True),
            "genome": self._gnomad_frequency_block(item.get("genome")),
            "exome": self._gnomad_frequency_block(item.get("exome")),
        }
        frequencies = {key: value for key, value in frequencies.items() if value}
        record = {
            "category": "population_frequency",
            "source": self.spec.name,
            "label": rsid or variant_id or (input_variant.label if input_variant else context.gene),
            "summary": self._gnomad_variant_summary(item, frequencies, primary_transcript),
            "source_id": variant_id,
            "url": self._gnomad_variant_url(variant_id or (input_variant.label if input_variant else "")),
            "variant": input_variant.label if input_variant else (rsid or variant_id),
            "variant_id": variant_id,
            "rsid": rsid,
            "rsids": rsids,
            "caid": _clean_cell(item.get("caid")),
            "dataset": GNOMAD_DATASET,
            "reference_genome": _clean_cell(item.get("reference_genome")),
            "chromosome": _clean_cell(item.get("chrom")),
            "position": item.get("pos"),
            "ref": _clean_cell(item.get("ref")),
            "alt": _clean_cell(item.get("alt")),
            "filters": self._gnomad_variant_filters(frequencies),
            "flags": self._gnomad_clean_list(item.get("flags")),
            "frequencies": frequencies,
            "top_populations": self._gnomad_top_populations(frequencies),
            "faf": self._gnomad_faf_summary(frequencies),
            "transcript_consequence": primary_transcript,
            "transcript_consequences": transcript_consequences[:5],
            "in_silico_predictors": self._gnomad_predictors(item.get("in_silico_predictors")),
            "lof_curations": self._gnomad_lof_curations(item.get("lof_curations")),
            "non_coding_constraint": self._gnomad_non_coding_constraint(item.get("non_coding_constraint")),
            "gene_constraint": self._gnomad_constraint(gene.get("gnomad_constraint") if gene else {}),
            "clinvar_release_date": _clean_cell(meta.get("clinvar_release_date")),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _gnomad_gene_record(self, context: KnowledgeQuery, gene: dict[str, Any], meta: dict[str, Any]) -> dict[str, Any]:
        symbol = _clean_cell(gene.get("symbol")) or context.gene
        record = {
            "category": "gene_constraint",
            "source": self.spec.name,
            "label": f"{symbol} gnomAD constraint",
            "summary": self._gnomad_gene_summary(gene),
            "source_id": _clean_cell(gene.get("gene_id")),
            "url": self._gnomad_gene_url(symbol),
            "gene": symbol,
            "gene_id": _clean_cell(gene.get("gene_id")),
            "gene_version": _clean_cell(gene.get("gene_version")),
            "gene_name": _clean_cell(gene.get("name")),
            "gencode_symbol": _clean_cell(gene.get("gencode_symbol")),
            "hgnc_id": _clean_cell(gene.get("hgnc_id")),
            "ncbi_id": _clean_cell(gene.get("ncbi_id")),
            "omim_id": _clean_cell(gene.get("omim_id")),
            "reference_genome": _clean_cell(gene.get("reference_genome")),
            "location": self._gnomad_gene_location(gene),
            "canonical_transcript_id": _clean_cell(gene.get("canonical_transcript_id")),
            "mane_select_transcript": self._gnomad_mane_transcript(gene.get("mane_select_transcript")),
            "flags": self._gnomad_clean_list(gene.get("flags")),
            "constraint": self._gnomad_constraint(gene.get("gnomad_constraint")),
            "clinvar_release_date": _clean_cell(meta.get("clinvar_release_date")),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _gnomad_variant_id(self, variant: Any) -> str:
        if not variant or not variant.chrom or not variant.pos or not variant.ref or not variant.alt:
            return ""
        alt = str(variant.alt).split(",", 1)[0].strip()
        if not alt:
            return ""
        return f"{variant.chrom.removeprefix('chr')}-{variant.pos}-{variant.ref}-{alt}"

    def _gnomad_reference_genome(self, context: KnowledgeQuery) -> str:
        return "GRCh38" if str(context.genome_build).lower() in {"hg38", "grch38"} else "GRCh37"

    def _gnomad_variant_summary(
        self,
        item: dict[str, Any],
        frequencies: dict[str, dict[str, Any]],
        transcript: dict[str, Any],
    ) -> str:
        variant_id = _clean_cell(item.get("variantId") or item.get("variant_id"))
        rsid = _clean_cell(item.get("rsid"))
        identity = rsid or variant_id or "variant"
        coordinate = self._gnomad_variant_coordinate(item)
        prefix = f"gnomAD v4 {identity}"
        if variant_id and variant_id != identity:
            prefix = f"{prefix} ({variant_id})"
        if coordinate:
            prefix = f"{prefix} at {coordinate}"

        parts: list[str] = []
        for key, label in (("joint", "joint"), ("genome", "genomes"), ("exome", "exomes")):
            text = self._gnomad_frequency_text(label, frequencies.get(key))
            if text:
                parts.append(text)
        pop_text = self._gnomad_top_population_text(frequencies)
        if pop_text:
            parts.append(pop_text)
        if transcript:
            transcript_text = self._gnomad_transcript_text(transcript)
            if transcript_text:
                parts.append(transcript_text)
        predictors_text = self._gnomad_predictor_text(item.get("in_silico_predictors"))
        if predictors_text:
            parts.append(predictors_text)
        filters = self._gnomad_variant_filters(frequencies)
        if filters:
            parts.append(f"filters: {', '.join(filters[:4])}")
        flags = self._gnomad_clean_list(item.get("flags"))
        if flags:
            parts.append(f"flags: {', '.join(flags[:4])}")
        return f"{prefix}: {'; '.join(parts)}." if parts else f"{prefix}: no population frequency counts returned."

    def _gnomad_gene_summary(self, gene: dict[str, Any]) -> str:
        symbol = _clean_cell(gene.get("symbol")) or "gene"
        gene_id = _clean_cell(gene.get("gene_id"))
        constraint = self._gnomad_constraint(gene.get("gnomad_constraint"))
        prefix = f"gnomAD v4 {symbol} gene constraint"
        if gene_id:
            prefix = f"{prefix} ({gene_id})"
        parts: list[str] = []
        pli = self._gnomad_float_text(constraint.get("pLI") if constraint else None)
        if pli:
            parts.append(f"pLI {pli}")
        oe_lof = self._gnomad_float_text(constraint.get("oe_lof") if constraint else None)
        if oe_lof:
            ci = self._gnomad_ci_text(constraint, "oe_lof")
            lof_z = self._gnomad_float_text(constraint.get("lof_z"))
            text = f"LoF O/E {oe_lof}{ci}"
            if lof_z:
                text = f"{text}, Z {lof_z}"
            parts.append(text)
        oe_mis = self._gnomad_float_text(constraint.get("oe_mis") if constraint else None)
        if oe_mis:
            ci = self._gnomad_ci_text(constraint, "oe_mis")
            mis_z = self._gnomad_float_text(constraint.get("mis_z"))
            text = f"missense O/E {oe_mis}{ci}"
            if mis_z:
                text = f"{text}, Z {mis_z}"
            parts.append(text)
        observed_expected = self._gnomad_observed_expected_text(constraint)
        if observed_expected:
            parts.append(observed_expected)
        flags = self._gnomad_clean_list(constraint.get("flags") if constraint else [])
        if flags:
            parts.append(f"constraint flags: {', '.join(flags[:4])}")
        location = self._gnomad_gene_location(gene)
        if location:
            parts.append(f"gene location {location}")
        return f"{prefix}: {'; '.join(parts)}." if parts else f"{prefix}: no constraint metrics returned."

    def _gnomad_frequency_block(self, value: Any, *, compute_af: bool = False) -> dict[str, Any]:
        if not isinstance(value, dict):
            return {}
        ac = self._gnomad_int(value.get("ac"))
        an = self._gnomad_int(value.get("an"))
        af = value.get("af")
        if af is None and compute_af and ac is not None and an:
            af = ac / an
        homozygote_count = self._gnomad_int(value.get("homozygote_count"))
        if homozygote_count is None:
            homozygote_count = self._gnomad_int(value.get("ac_hom"))
        hemizygote_count = self._gnomad_int(value.get("hemizygote_count"))
        if hemizygote_count is None:
            hemizygote_count = self._gnomad_int(value.get("ac_hemi"))
        block = {
            "ac": ac,
            "an": an,
            "af": af,
            "homozygote_count": homozygote_count,
            "hemizygote_count": hemizygote_count,
            "filters": self._gnomad_clean_list(value.get("filters")),
            "flags": self._gnomad_clean_list(value.get("flags")),
            "faf95": self._gnomad_faf(value.get("faf95")),
            "faf99": self._gnomad_faf(value.get("faf99")),
            "populations": self._gnomad_population_rows(value.get("populations")),
        }
        return {key: item for key, item in block.items() if item not in ("", [], {}, None)}

    def _gnomad_population_rows(self, values: Any) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            ac = self._gnomad_int(value.get("ac"))
            an = self._gnomad_int(value.get("an"))
            homozygote_count = self._gnomad_int(value.get("homozygote_count"))
            if homozygote_count is None:
                homozygote_count = self._gnomad_int(value.get("ac_hom"))
            hemizygote_count = self._gnomad_int(value.get("hemizygote_count"))
            if hemizygote_count is None:
                hemizygote_count = self._gnomad_int(value.get("ac_hemi"))
            row = {
                "id": _clean_cell(value.get("id")),
                "ac": ac,
                "an": an,
                "af": ac / an if ac is not None and an else None,
                "homozygote_count": homozygote_count,
                "hemizygote_count": hemizygote_count,
            }
            rows.append({key: item for key, item in row.items() if item not in ("", None)})
        return rows

    def _gnomad_faf(self, value: Any) -> dict[str, Any]:
        if not isinstance(value, dict):
            return {}
        row = {
            "popmax": value.get("popmax"),
            "popmax_population": _clean_cell(value.get("popmax_population")),
        }
        return {key: item for key, item in row.items() if item not in ("", None)}

    def _gnomad_faf_summary(self, frequencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
        summary: dict[str, Any] = {}
        for key, block in frequencies.items():
            for faf_key in ("faf95", "faf99"):
                if block.get(faf_key):
                    summary[f"{key}_{faf_key}"] = block[faf_key]
        return summary

    def _gnomad_top_populations(self, frequencies: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
        top_rows: list[dict[str, Any]] = []
        for dataset, block in frequencies.items():
            populations = block.get("populations") if isinstance(block, dict) else []
            if not isinstance(populations, list):
                continue
            eligible = [
                row
                for row in populations
                if isinstance(row, dict)
                and self._gnomad_int(row.get("ac"))
                and not self._gnomad_population_is_sex_specific(row.get("id"))
            ]
            if not eligible:
                continue
            top = sorted(
                eligible,
                key=lambda row: (float(row.get("af") or 0), self._gnomad_int(row.get("ac")) or 0),
                reverse=True,
            )[0]
            top_rows.append({"dataset": dataset, **top})
        return sorted(
            top_rows,
            key=lambda row: (float(row.get("af") or 0), self._gnomad_int(row.get("ac")) or 0),
            reverse=True,
        )

    def _gnomad_variant_filters(self, frequencies: dict[str, dict[str, Any]]) -> list[str]:
        filters: list[str] = []
        for block in frequencies.values():
            if not isinstance(block, dict):
                continue
            filters.extend(self._gnomad_clean_list(block.get("filters")))
        return self._gnomad_clean_list(filters)

    def _gnomad_frequency_text(self, label: str, block: dict[str, Any] | None) -> str:
        if not block:
            return ""
        ac = self._gnomad_int(block.get("ac"))
        an = self._gnomad_int(block.get("an"))
        af = block.get("af")
        if af is None and ac is not None and an:
            af = ac / an
        if ac is None and an is None and af is None:
            return ""
        parts = [f"{label} AF {self._gnomad_float_text(af)}" if af is not None else f"{label} counts"]
        if ac is not None and an is not None:
            parts.append(f"AC {self._gnomad_count(ac)}/AN {self._gnomad_count(an)}")
        hom = self._gnomad_int(block.get("homozygote_count"))
        hemi = self._gnomad_int(block.get("hemizygote_count"))
        if hom is not None:
            parts.append(f"hom {self._gnomad_count(hom)}")
        if hemi:
            parts.append(f"hemi {self._gnomad_count(hemi)}")
        faf95 = block.get("faf95") if isinstance(block.get("faf95"), dict) else {}
        faf_pop = _clean_cell(faf95.get("popmax_population"))
        faf_value = self._gnomad_float_text(faf95.get("popmax"))
        if faf_pop and faf_value:
            parts.append(f"FAF95 {faf_pop} {faf_value}")
        return " ".join([parts[0], f"({'; '.join(parts[1:])})" if len(parts) > 1 else ""]).strip()

    def _gnomad_top_population_text(self, frequencies: dict[str, dict[str, Any]]) -> str:
        top_rows = self._gnomad_top_populations(frequencies)
        if not top_rows:
            return ""
        best = sorted(
            top_rows,
            key=lambda row: (float(row.get("af") or 0), self._gnomad_int(row.get("ac")) or 0),
            reverse=True,
        )[0]
        return (
            f"highest observed population {best.get('dataset')} {best.get('id')} "
            f"AF {self._gnomad_float_text(best.get('af'))} "
            f"(AC {self._gnomad_count(best.get('ac'))}/AN {self._gnomad_count(best.get('an'))})"
        )

    def _gnomad_transcript_consequences(self, values: Any, context: KnowledgeQuery) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "gene_id": _clean_cell(value.get("gene_id")),
                "gene_symbol": _clean_cell(value.get("gene_symbol")),
                "transcript_id": _clean_cell(value.get("transcript_id")),
                "transcript_version": _clean_cell(value.get("transcript_version")),
                "is_canonical": value.get("is_canonical"),
                "is_mane_select": value.get("is_mane_select"),
                "major_consequence": _clean_cell(value.get("major_consequence")),
                "consequence_terms": self._gnomad_clean_list(value.get("consequence_terms")),
                "hgvsc": _clean_cell(value.get("hgvsc")),
                "hgvsp": _clean_cell(value.get("hgvsp")),
                "hgvs": _clean_cell(value.get("hgvs")),
                "lof": _clean_cell(value.get("lof")),
                "lof_filter": _clean_cell(value.get("lof_filter")),
                "lof_flags": _clean_cell(value.get("lof_flags")),
                "polyphen_prediction": _clean_cell(value.get("polyphen_prediction")),
                "sift_prediction": _clean_cell(value.get("sift_prediction")),
                "domains": self._gnomad_clean_list(value.get("domains")),
            }
            rows.append({key: item for key, item in row.items() if item not in ("", [], None)})
        return sorted(
            rows,
            key=lambda row: (
                not bool(row.get("is_mane_select")),
                not bool(row.get("is_canonical")),
                _clean_cell(row.get("gene_symbol")).upper() != context.gene.upper(),
            ),
        )

    def _gnomad_primary_transcript(self, rows: list[dict[str, Any]], context: KnowledgeQuery) -> dict[str, Any]:
        if not rows:
            return {}
        for row in rows:
            if row.get("is_mane_select"):
                return row
        for row in rows:
            if row.get("is_canonical"):
                return row
        for row in rows:
            if _clean_cell(row.get("gene_symbol")).upper() == context.gene.upper():
                return row
        return rows[0]

    def _gnomad_transcript_text(self, transcript: dict[str, Any]) -> str:
        consequences = self._gnomad_clean_list(
            transcript.get("consequence_terms") or transcript.get("major_consequence")
        )
        text = ", ".join(consequences[:3]) if consequences else _clean_cell(transcript.get("major_consequence"))
        gene = _clean_cell(transcript.get("gene_symbol"))
        transcript_id = _clean_cell(transcript.get("transcript_id"))
        if text and gene and transcript_id:
            text = f"consequence {text} in {gene} transcript {transcript_id}"
        elif text and transcript_id:
            text = f"consequence {text} in transcript {transcript_id}"
        elif text:
            text = f"consequence {text}"
        hgvs = _clean_cell(transcript.get("hgvsc") or transcript.get("hgvsp") or transcript.get("hgvs"))
        if hgvs:
            text = f"{text}; HGVS {hgvs}" if text else f"HGVS {hgvs}"
        lof = _clean_cell(transcript.get("lof"))
        if lof:
            text = f"{text}; LoF {lof}" if text else f"LoF {lof}"
        return text

    def _gnomad_predictors(self, values: Any) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "id": _clean_cell(value.get("id")),
                "value": _clean_cell(value.get("value")),
                "flags": self._gnomad_clean_list(value.get("flags")),
            }
            cleaned = {key: item for key, item in row.items() if item not in ("", [], None)}
            if cleaned:
                rows.append(cleaned)
        return rows

    def _gnomad_predictor_text(self, values: Any) -> str:
        predictors = self._gnomad_predictors(values)
        if not predictors:
            return ""
        parts = []
        for predictor in predictors[:3]:
            if predictor.get("id") and predictor.get("value"):
                parts.append(f"{predictor['id']} {predictor['value']}")
        return f"predictors: {', '.join(parts)}" if parts else ""

    def _gnomad_lof_curations(self, values: Any) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        if not isinstance(values, list):
            return rows
        for value in values:
            if not isinstance(value, dict):
                continue
            row = {
                "gene_id": _clean_cell(value.get("gene_id")),
                "gene_symbol": _clean_cell(value.get("gene_symbol")),
                "verdict": _clean_cell(value.get("verdict")),
                "flags": self._gnomad_clean_list(value.get("flags")),
                "project": _clean_cell(value.get("project")),
            }
            cleaned = {key: item for key, item in row.items() if item not in ("", [], None)}
            if cleaned:
                rows.append(cleaned)
        return rows

    def _gnomad_non_coding_constraint(self, value: Any) -> dict[str, Any]:
        if not isinstance(value, dict):
            return {}
        keys = ("chrom", "start", "stop", "element_id", "possible", "observed", "expected", "oe", "z")
        row = {key: (_clean_cell(value.get(key)) if key in {"chrom", "element_id"} else value.get(key)) for key in keys}
        return {key: item for key, item in row.items() if item not in ("", None)}

    def _gnomad_constraint(self, value: Any) -> dict[str, Any]:
        if not isinstance(value, dict):
            return {}
        keys = (
            "exp_lof",
            "exp_mis",
            "exp_syn",
            "obs_lof",
            "obs_mis",
            "obs_syn",
            "oe_lof",
            "oe_lof_lower",
            "oe_lof_upper",
            "oe_lof_percentile",
            "oe_mis",
            "oe_mis_lower",
            "oe_mis_upper",
            "oe_syn",
            "oe_syn_lower",
            "oe_syn_upper",
            "lof_z",
            "mis_z",
            "syn_z",
            "pli",
            "pLI",
        )
        row = {key: value.get(key) for key in keys}
        row["flags"] = self._gnomad_clean_list(value.get("flags"))
        return {key: item for key, item in row.items() if item not in ("", [], None)}

    def _gnomad_mane_transcript(self, value: Any) -> dict[str, str]:
        if not isinstance(value, dict):
            return {}
        row = {
            "ensembl_id": _clean_cell(value.get("ensembl_id")),
            "ensembl_version": _clean_cell(value.get("ensembl_version")),
            "refseq_id": _clean_cell(value.get("refseq_id")),
            "refseq_version": _clean_cell(value.get("refseq_version")),
        }
        return {key: item for key, item in row.items() if item}

    def _gnomad_gene_location(self, gene: dict[str, Any]) -> str:
        chrom = _clean_cell(gene.get("chrom"))
        start = _clean_cell(gene.get("start"))
        stop = _clean_cell(gene.get("stop"))
        reference_genome = _clean_cell(gene.get("reference_genome"))
        if chrom and start and stop:
            prefix = f"{reference_genome} " if reference_genome else ""
            return f"{prefix}{chrom}:{start}-{stop}"
        return ""

    def _gnomad_variant_coordinate(self, item: dict[str, Any]) -> str:
        reference = _clean_cell(item.get("reference_genome"))
        chrom = _clean_cell(item.get("chrom"))
        pos = _clean_cell(item.get("pos"))
        ref = _clean_cell(item.get("ref"))
        alt = _clean_cell(item.get("alt"))
        if not (chrom and pos):
            return ""
        allele = f" {ref}>{alt}" if ref or alt else ""
        prefix = f"{reference} " if reference else ""
        return f"{prefix}{chrom}:{pos}{allele}"

    def _gnomad_observed_expected_text(self, constraint: dict[str, Any]) -> str:
        if not constraint:
            return ""
        parts: list[str] = []
        obs_lof = self._gnomad_count(constraint.get("obs_lof"))
        exp_lof = self._gnomad_float_text(constraint.get("exp_lof"))
        if obs_lof and exp_lof:
            parts.append(f"LoF observed/expected {obs_lof}/{exp_lof}")
        obs_mis = self._gnomad_count(constraint.get("obs_mis"))
        exp_mis = self._gnomad_float_text(constraint.get("exp_mis"))
        if obs_mis and exp_mis:
            parts.append(f"missense observed/expected {obs_mis}/{exp_mis}")
        return "; ".join(parts)

    def _gnomad_ci_text(self, constraint: dict[str, Any], key: str) -> str:
        lower = self._gnomad_float_text(constraint.get(f"{key}_lower"))
        upper = self._gnomad_float_text(constraint.get(f"{key}_upper"))
        return f" ({lower}-{upper})" if lower and upper else ""

    def _gnomad_clean_list(self, values: Any) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, list):
            values = [values]
        cleaned: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in cleaned:
                cleaned.append(text)
        return cleaned

    def _gnomad_int(self, value: Any) -> int | None:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _gnomad_count(self, value: Any) -> str:
        integer = self._gnomad_int(value)
        if integer is None:
            return _clean_cell(value)
        return f"{integer:,}"

    def _gnomad_float_text(self, value: Any) -> str:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return _clean_cell(value)
        if number == 0:
            return "0"
        if abs(number) < 0.001 or abs(number) >= 1000:
            return f"{number:.3g}"
        return f"{number:.4g}"

    def _gnomad_population_is_sex_specific(self, value: Any) -> bool:
        text = _clean_cell(value)
        return text in {"XX", "XY"} or text.endswith("_XX") or text.endswith("_XY")

    def _gnomad_variant_url(self, variant_id: str) -> str:
        if variant_id:
            return f"https://gnomad.broadinstitute.org/variant/{quote(variant_id)}?dataset={GNOMAD_DATASET}"
        return f"https://gnomad.broadinstitute.org/?dataset={GNOMAD_DATASET}"

    def _gnomad_gene_url(self, symbol: str) -> str:
        return f"https://gnomad.broadinstitute.org/gene/{quote(symbol)}?dataset={GNOMAD_DATASET}"

    def _query_civic(self, context: KnowledgeQuery, started: float) -> SourceResult:
        url = "https://civicdb.org/api/graphql"
        variables = {"entrezSymbol": context.gene, "variantFirst": 5, "evidenceFirst": 3}
        try:
            payload = self.client.post_json(
                url,
                json_payload={"query": self._civic_gene_evidence_query(), "variables": variables},
                headers={"Accept": "application/json", "Content-Type": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
            errors = payload.get("errors") if isinstance(payload, dict) else None
            if errors:
                message = self._civic_graphql_error_message(errors, "CIViC GraphQL query failed.")
                return SourceResult(
                    self.spec.key,
                    "failed",
                    message,
                    errors=[message],
                    queried_urls=[url],
                    elapsed_ms=_elapsed_ms(started),
                )

            gene = (payload.get("data", {}) if isinstance(payload, dict) else {}).get("gene") or {}
            if not gene:
                return SourceResult(
                    self.spec.key,
                    "ok",
                    f"Queried CIViC; no gene record found for {context.gene}.",
                    [],
                    queried_urls=[url],
                    elapsed_ms=_elapsed_ms(started),
                )

            records = self._civic_records(context, gene)
            return SourceResult(
                self.spec.key,
                "ok",
                f"Queried CIViC; {len(records)} variant evidence record(s) returned for {context.gene}.",
                records,
                queried_urls=[url],
                elapsed_ms=_elapsed_ms(started),
            )
        except KnowledgeRequestError as exc:
            return SourceResult(self.spec.key, "failed", str(exc), errors=[str(exc)], queried_urls=[url], elapsed_ms=_elapsed_ms(started))

    def _civic_gene_evidence_query(self) -> str:
        return """
        query CivicGeneEvidence($entrezSymbol: String!, $variantFirst: Int!, $evidenceFirst: Int!) {
          gene(entrezSymbol: $entrezSymbol) {
            id
            name
            link
            description
            featureAliases
            variants(first: $variantFirst) {
              totalCount
              nodes {
                id
                name
                link
                variantAliases
                variantTypes { id name }
                hgvsDescriptions
                clinvarIds
                alleleRegistryId
                maneSelectTranscript
                singleVariantMolecularProfile {
                  id
                  name
                  link
                  description
                  evidenceItems(first: $evidenceFirst, includeRejected: false) {
                    totalCount
                    nodes {
                      id
                      name
                      link
                      description
                      descriptionWithNames
                      status
                      evidenceType
                      evidenceLevel
                      evidenceRating
                      evidenceDirection
                      significance
                      variantOrigin
                      variantHgvs
                      therapyInteractionType
                      disease { id name link }
                      therapies { id name link }
                      source { id citation citationId sourceType link }
                      phenotypes { id name link }
                    }
                  }
                }
              }
            }
          }
        }
        """

    def _civic_records(self, context: KnowledgeQuery, gene: dict[str, Any]) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        variants = ((gene.get("variants") or {}).get("nodes") or []) if isinstance(gene, dict) else []
        for variant in variants[:5]:
            if not isinstance(variant, dict):
                continue
            molecular_profile = variant.get("singleVariantMolecularProfile") or {}
            if not isinstance(molecular_profile, dict):
                molecular_profile = {}
            evidence_connection = molecular_profile.get("evidenceItems") or {}
            if not isinstance(evidence_connection, dict):
                evidence_connection = {}
            evidence_nodes = [item for item in (evidence_connection.get("nodes") or [])[:3] if isinstance(item, dict)]
            evidence_items = [self._civic_evidence_item(item) for item in evidence_nodes]
            evidence_count = self._civic_int(evidence_connection.get("totalCount")) or len(evidence_items)
            first_evidence = evidence_items[0] if evidence_items else {}

            variant_name = _clean_cell(variant.get("name")) or context.gene
            molecular_profile_name = _clean_cell(molecular_profile.get("name"))
            variant_types = self._civic_names(variant.get("variantTypes"))
            record = {
                "category": "cancer_variant",
                "source": self.spec.name,
                "label": molecular_profile_name or variant_name,
                "summary": self._civic_summary(
                    context=context,
                    gene=gene,
                    variant=variant,
                    molecular_profile=molecular_profile,
                    evidence_items=evidence_items,
                    evidence_count=evidence_count,
                ),
                "source_id": _clean_cell(variant.get("id")),
                "url": self._civic_url(variant.get("link"), f"https://civicdb.org/variants/{_clean_cell(variant.get('id'))}"),
                "gene": _clean_cell(gene.get("name")) or context.gene,
                "gene_id": _clean_cell(gene.get("id")),
                "variant": variant_name,
                "variant_aliases": self._civic_clean_list(variant.get("variantAliases")),
                "variant_types": variant_types,
                "hgvs_descriptions": self._civic_clean_list(variant.get("hgvsDescriptions")),
                "clinvar_ids": self._civic_clean_list(variant.get("clinvarIds")),
                "allele_registry_id": _clean_cell(variant.get("alleleRegistryId")),
                "mane_select_transcript": _clean_cell(variant.get("maneSelectTranscript")),
                "molecular_profile": molecular_profile_name,
                "molecular_profile_id": _clean_cell(molecular_profile.get("id")),
                "molecular_profile_description": _clean_cell(molecular_profile.get("description")),
                "evidence_count": evidence_count,
                "evidence_items": evidence_items,
            }
            record.update(
                {
                    "evidence_type": first_evidence.get("evidence_type", ""),
                    "evidence_level": first_evidence.get("evidence_level", ""),
                    "evidence_rating": first_evidence.get("evidence_rating", ""),
                    "evidence_direction": first_evidence.get("evidence_direction", ""),
                    "significance": first_evidence.get("significance", ""),
                    "variant_origin": first_evidence.get("variant_origin", ""),
                    "variant_hgvs": first_evidence.get("variant_hgvs", ""),
                    "disease": first_evidence.get("disease", ""),
                    "therapies": first_evidence.get("therapies", []),
                    "therapy_interaction_type": first_evidence.get("therapy_interaction_type", ""),
                    "citation": first_evidence.get("citation", ""),
                    "citation_id": first_evidence.get("citation_id", ""),
                    "source_type": first_evidence.get("source_type", ""),
                }
            )
            records.append({key: value for key, value in record.items() if value not in ("", [], None)})
        return records

    def _civic_summary(
        self,
        *,
        context: KnowledgeQuery,
        gene: dict[str, Any],
        variant: dict[str, Any],
        molecular_profile: dict[str, Any],
        evidence_items: list[dict[str, Any]],
        evidence_count: int,
    ) -> str:
        gene_name = _clean_cell(gene.get("name")) or context.gene
        variant_name = _clean_cell(variant.get("name")) or gene_name
        molecular_profile_name = _clean_cell(molecular_profile.get("name"))
        identity = variant_name
        if molecular_profile_name and molecular_profile_name != variant_name:
            identity = f"{identity} / {molecular_profile_name}"

        annotations: list[str] = []
        variant_types = self._civic_names(variant.get("variantTypes"))
        if variant_types:
            annotations.append(f"variant type {', '.join(variant_types[:3])}")
        aliases = self._civic_clean_list(variant.get("variantAliases"))
        if aliases:
            annotations.append(f"aliases: {', '.join(aliases[:3])}")
        hgvs = self._civic_clean_list(variant.get("hgvsDescriptions"))
        if hgvs:
            annotations.append(f"HGVS {', '.join(hgvs[:3])}")
        clinvar_ids = self._civic_clean_list(variant.get("clinvarIds"))
        if clinvar_ids:
            annotations.append(f"ClinVar {', '.join(clinvar_ids[:3])}")
        mane_transcript = _clean_cell(variant.get("maneSelectTranscript"))
        if mane_transcript:
            annotations.append(f"MANE {mane_transcript}")
        allele_registry_id = _clean_cell(variant.get("alleleRegistryId"))
        if allele_registry_id:
            annotations.append(f"Allele Registry {allele_registry_id}")

        prefix = f"CIViC {gene_name} variant {identity}"
        if annotations:
            prefix = f"{prefix} ({'; '.join(annotations)})"

        if not evidence_items:
            if evidence_count:
                return f"{prefix}: {evidence_count} accepted evidence item(s) available in CIViC."
            return f"{prefix}: no accepted single-variant evidence items returned by this query."

        evidence = evidence_items[0]
        clauses: list[str] = []
        evidence_label_parts = []
        evidence_type = evidence.get("evidence_type", "")
        if evidence_type:
            evidence_label_parts.append(f"{evidence_type} evidence")
        evidence_level = evidence.get("evidence_level", "")
        if evidence_level:
            evidence_label_parts.append(f"level {evidence_level}")
        evidence_rating = evidence.get("evidence_rating", "")
        if evidence_rating:
            evidence_label_parts.append(f"rating {evidence_rating}")
        evidence_direction = evidence.get("evidence_direction", "")
        if evidence_direction:
            evidence_label_parts.append(evidence_direction)
        if evidence_label_parts:
            clauses.append(", ".join(evidence_label_parts))

        for label, key in (
            ("significance", "significance"),
            ("disease", "disease"),
            ("variant origin", "variant_origin"),
            ("variant HGVS", "variant_hgvs"),
            ("therapy interaction", "therapy_interaction_type"),
            ("status", "status"),
        ):
            value = evidence.get(key, "")
            if value:
                clauses.append(f"{label} {value}")
        therapies = evidence.get("therapies") or []
        if therapies:
            clauses.append(f"therapies {', '.join(therapies[:3])}")
        source_text = self._civic_source_text(evidence)
        if source_text:
            clauses.append(source_text)
        if evidence_count > 1:
            clauses.append(f"{evidence_count} total accepted evidence item(s)")

        summary = f"{prefix}: {'; '.join(clauses)}."
        description = self._civic_clip(evidence.get("description_with_names") or evidence.get("description"), 220)
        if description:
            summary = f"{summary} {description}"
        return summary

    def _civic_evidence_item(self, item: dict[str, Any]) -> dict[str, Any]:
        source = item.get("source") or {}
        if not isinstance(source, dict):
            source = {}
        record = {
            "id": _clean_cell(item.get("id")),
            "name": _clean_cell(item.get("name")),
            "status": _clean_cell(item.get("status")),
            "evidence_type": _clean_cell(item.get("evidenceType")),
            "evidence_level": _clean_cell(item.get("evidenceLevel")),
            "evidence_rating": _clean_cell(item.get("evidenceRating")),
            "evidence_direction": _clean_cell(item.get("evidenceDirection")),
            "significance": _clean_cell(item.get("significance")),
            "variant_origin": _clean_cell(item.get("variantOrigin")),
            "variant_hgvs": _clean_cell(item.get("variantHgvs")),
            "therapy_interaction_type": _clean_cell(item.get("therapyInteractionType")),
            "disease": self._civic_name(item.get("disease")),
            "therapies": self._civic_names(item.get("therapies")),
            "phenotypes": self._civic_names(item.get("phenotypes")),
            "citation": _clean_cell(source.get("citation")),
            "citation_id": _clean_cell(source.get("citationId")),
            "source_type": _clean_cell(source.get("sourceType")),
            "source_id": _clean_cell(source.get("id")),
            "source_url": self._civic_url(source.get("link"), ""),
            "url": self._civic_url(item.get("link"), f"https://civicdb.org/evidence/{_clean_cell(item.get('id'))}"),
            "description": _clean_cell(item.get("description")),
            "description_with_names": _clean_cell(item.get("descriptionWithNames")),
        }
        return {key: value for key, value in record.items() if value not in ("", [], None)}

    def _civic_source_text(self, evidence: dict[str, Any]) -> str:
        citation_id = _clean_cell(evidence.get("citation_id"))
        source_type = _clean_cell(evidence.get("source_type"))
        citation = _clean_cell(evidence.get("citation"))
        if citation_id and source_type.upper() == "PUBMED":
            return f"source PMID {citation_id}"
        if citation_id and source_type:
            return f"source {source_type} {citation_id}"
        if citation:
            return f"source {citation}"
        return ""

    def _graphql_error_message(self, errors: Any, fallback: str, source_name: str) -> str:
        if isinstance(errors, list) and errors:
            first_error = errors[0]
            if isinstance(first_error, dict):
                message = _clean_cell(first_error.get("message"))
                if message:
                    return f"{source_name} GraphQL query failed: {message}"
            message = _clean_cell(first_error)
            if message:
                return f"{source_name} GraphQL query failed: {message}"
        return fallback

    def _civic_graphql_error_message(self, errors: Any, fallback: str) -> str:
        return self._graphql_error_message(errors, fallback, "CIViC")

    def _civic_url(self, value: Any, fallback: str) -> str:
        text = _clean_cell(value)
        if text.startswith("http://") or text.startswith("https://"):
            return text
        if text.startswith("/"):
            return f"https://civicdb.org{text}"
        return fallback

    def _civic_name(self, value: Any) -> str:
        if isinstance(value, dict):
            return _clean_cell(value.get("name")) or _clean_cell(value.get("id"))
        return _clean_cell(value)

    def _civic_names(self, values: Any) -> list[str]:
        if values is None:
            return []
        if not isinstance(values, list):
            values = [values]
        names: list[str] = []
        for value in values:
            text = self._civic_name(value)
            if text and text not in names:
                names.append(text)
        return names

    def _civic_clean_list(self, values: Any) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str):
            return [values.strip()] if values.strip() else []
        if not isinstance(values, list):
            values = [values]
        cleaned: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in cleaned:
                cleaned.append(text)
        return cleaned

    def _civic_int(self, value: Any) -> int:
        try:
            return int(value)
        except (TypeError, ValueError):
            return 0

    def _civic_clip(self, value: Any, limit: int) -> str:
        text = _clean_cell(value)
        if len(text) <= limit:
            return text
        return f"{text[: limit - 3].rstrip()}..."


class ClinGenConnector(BaseConnector):
    """ClinGen connector backed by official gene-centered downloads and APIs."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        urls: list[str] = []
        records: list[dict[str, Any]] = []
        warnings: list[str] = []
        request_errors: list[KnowledgeRequestError] = []

        download_steps = (
            ("Gene-Disease Validity", CLINGEN_VALIDITY_URL, self._validity_records, CLINGEN_PRIMARY_TIMEOUT_SECONDS),
            ("Dosage Sensitivity", CLINGEN_DOSAGE_URL, self._dosage_records, CLINGEN_PRIMARY_TIMEOUT_SECONDS),
        )
        for label, url, parser, timeout_seconds in download_steps:
            try:
                text = self.client.get_text(
                    url,
                    rate_limit_per_second=self.spec.rate_limit_per_second,
                    timeout=timeout_seconds,
                )
                urls.append(url)
                records.extend(parser(text, context))
            except KnowledgeRequestError as exc:
                request_errors.append(exc)
                warnings.append(f"{label} request failed: {exc}")
            except (ValueError, csv.Error) as exc:
                warnings.append(f"{label} response could not be parsed: {exc}")

        try:
            text = self.client.get_text(
                CLINGEN_SUMMARY_URL,
                rate_limit_per_second=self.spec.rate_limit_per_second,
                timeout=CLINGEN_SUMMARY_TIMEOUT_SECONDS,
            )
            urls.append(CLINGEN_SUMMARY_URL)
            records.extend(self._summary_records(text, context))
        except KnowledgeRequestError as exc:
            request_errors.append(exc)
            if "timed out" in str(exc).lower() or "timeout" in str(exc).lower():
                warnings.append(
                    "Optional ClinGen curation activity summary timed out; primary ClinGen feeds were still used."
                )
            else:
                warnings.append(
                    "Optional ClinGen curation activity summary request failed; primary ClinGen feeds were still used."
                )
        except (ValueError, csv.Error) as exc:
            warnings.append(f"Optional ClinGen curation activity summary response could not be parsed: {exc}")

        for context_label, url in CLINGEN_ACTIONABILITY_URLS:
            try:
                payload = self.client.get_json(
                    url,
                    rate_limit_per_second=self.spec.rate_limit_per_second,
                    timeout=CLINGEN_OPTIONAL_TIMEOUT_SECONDS,
                )
                urls.append(url)
                records.extend(self._actionability_records(payload, context, context_label))
            except KnowledgeRequestError as exc:
                request_errors.append(exc)
                warnings.append(f"Clinical Actionability {context_label} request failed: {exc}")
            except (TypeError, ValueError) as exc:
                warnings.append(f"Clinical Actionability {context_label} response could not be parsed: {exc}")

        records = self._dedupe_records(records)
        expected_request_count = len(download_steps) + 1 + len(CLINGEN_ACTIONABILITY_URLS)
        if request_errors and not records and len(request_errors) == expected_request_count:
            first_error = request_errors[0]
            return _request_failure_result(self.spec, first_error, queried_urls=urls, started=started)

        if records:
            message = f"Queried ClinGen; {len(records)} gene-centered curation record(s) returned for {context.gene}."
        else:
            message = f"Queried ClinGen; no ClinGen curation records found for {context.gene}."
        return SourceResult(
            source_key=self.spec.key,
            status="ok",
            message=message,
            records=records,
            warnings=warnings,
            queried_urls=urls,
            elapsed_ms=_elapsed_ms(started),
        )

    def _validity_records(self, text: str, context: KnowledgeQuery) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for row in self._csv_dict_rows(text, "GENE SYMBOL"):
            if not self._gene_matches(row.get("GENE SYMBOL"), context):
                continue
            disease = _clean_cell(row.get("DISEASE LABEL"))
            classification = _clean_cell(row.get("CLASSIFICATION"))
            moi = _clean_cell(row.get("MOI"))
            gcep = _clean_cell(row.get("GCEP"))
            date = _clean_cell(row.get("CLASSIFICATION DATE"))
            report = _clean_cell(row.get("ONLINE REPORT"))
            summary_parts = [
                f"{classification} gene-disease validity" if classification else "Gene-disease validity curation",
                f"for {disease}" if disease else "",
                f"({moi})" if moi else "",
            ]
            record = {
                "category": "gene_disease_validity",
                "source": self.spec.name,
                "label": f"ClinGen validity: {context.gene}{' - ' + disease if disease else ''}",
                "summary": " ".join(part for part in summary_parts if part).strip(),
                "gene": context.gene,
                "hgnc_id": _clean_cell(row.get("GENE ID (HGNC)")),
                "disease": disease,
                "mondo_id": _clean_cell(row.get("DISEASE ID (MONDO)")),
                "mode_of_inheritance": moi,
                "classification": classification,
                "assertion": classification,
                "sop": _clean_cell(row.get("SOP")),
                "expert_panel": gcep,
                "date": date,
                "url": report or self.spec.homepage,
                "provenance_url": CLINGEN_VALIDITY_URL,
            }
            records.append(record)
        return records

    def _dosage_records(self, text: str, context: KnowledgeQuery) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for row in self._csv_dict_rows(text, "GENE SYMBOL"):
            if not self._gene_matches(row.get("GENE SYMBOL"), context):
                continue
            haplo = _clean_cell(row.get("HAPLOINSUFFICIENCY"))
            triplo = _clean_cell(row.get("TRIPLOSENSITIVITY"))
            report = _clean_cell(row.get("ONLINE REPORT"))
            date = _clean_cell(row.get("DATE"))
            records.append(
                {
                    "category": "dosage_sensitivity",
                    "source": self.spec.name,
                    "label": f"ClinGen dosage sensitivity: {context.gene}",
                    "summary": f"Haploinsufficiency: {haplo or 'not reported'}; triplosensitivity: {triplo or 'not reported'}.",
                    "gene": context.gene,
                    "hgnc_id": _clean_cell(row.get("HGNC ID")),
                    "haploinsufficiency": haplo,
                    "triplosensitivity": triplo,
                    "date": date,
                    "url": report or self.spec.homepage,
                    "provenance_url": CLINGEN_DOSAGE_URL,
                }
            )
        return records

    def _summary_records(self, text: str, context: KnowledgeQuery) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for row in self._csv_dict_rows(text, "gene_symbol"):
            if not self._gene_matches(row.get("gene_symbol"), context):
                continue
            disease = _clean_cell(row.get("disease_label"))
            validity = _clean_cell(row.get("gene_disease_validity_assertion_classifications"))
            actionability = _clean_cell(row.get("actionability_assertion_classifications"))
            dosage_haplo = _clean_cell(row.get("dosage_haploinsufficiency_assertion"))
            dosage_triplo = _clean_cell(row.get("dosage_triplosensitivity_assertion"))
            parts = []
            if validity:
                parts.append(f"validity: {validity}")
            if actionability:
                parts.append(f"actionability: {actionability}")
            if dosage_haplo or dosage_triplo:
                parts.append(f"dosage: HI {dosage_haplo or 'n/a'}, TS {dosage_triplo or 'n/a'}")
            records.append(
                {
                    "category": "clinical_gene_curation_summary",
                    "source": self.spec.name,
                    "label": f"ClinGen curation summary: {context.gene}{' - ' + disease if disease else ''}",
                    "summary": "; ".join(parts) or "ClinGen curation summary row.",
                    "gene": context.gene,
                    "hgnc_id": _clean_cell(row.get("hgnc_id")),
                    "disease": disease,
                    "mondo_id": _clean_cell(row.get("mondo_id")),
                    "mode_of_inheritance": _clean_cell(row.get("mode_of_inheritance")),
                    "classification": validity,
                    "actionability": actionability,
                    "haploinsufficiency": dosage_haplo,
                    "triplosensitivity": dosage_triplo,
                    "expert_panel": _clean_cell(row.get("gene_disease_validity_gceps")),
                    "actionability_group": _clean_cell(row.get("actionability_groups")),
                    "url": _clean_cell(row.get("gene_url")) or self.spec.homepage,
                    "report_url": _clean_cell(row.get("gene_disease_validity_assertion_reports")),
                    "actionability_report_url": _clean_cell(row.get("actionability_assertion_reports")),
                    "dosage_report_url": _clean_cell(row.get("dosage_report")),
                    "provenance_url": CLINGEN_SUMMARY_URL,
                }
            )
        return records

    def _actionability_records(
        self,
        payload: Any,
        context: KnowledgeQuery,
        actionability_context: str,
    ) -> list[dict[str, Any]]:
        rows = self._flat_json_rows(payload)
        records: list[dict[str, Any]] = []
        for row in rows:
            if not self._gene_matches(row.get("geneOrVariant"), context):
                continue
            disease = _clean_cell(row.get("disease"))
            overall = _clean_cell(row.get("overall"))
            intervention = _clean_cell(row.get("intervention"))
            doc_id = _clean_cell(row.get("docId"))
            records.append(
                {
                    "category": "clinical_actionability",
                    "source": self.spec.name,
                    "label": f"ClinGen actionability ({actionability_context}): {context.gene}{' - ' + disease if disease else ''}",
                    "summary": (
                        f"Overall actionability score {overall or 'not reported'}"
                        f"{'; intervention: ' + intervention if intervention else ''}."
                    ),
                    "gene": context.gene,
                    "disease": disease,
                    "classification": _clean_cell(row.get("status-overall")),
                    "actionability_score": overall,
                    "outcome": _clean_cell(row.get("outcome")),
                    "intervention": intervention,
                    "severity": _clean_cell(row.get("severity")),
                    "likelihood": _clean_cell(row.get("likelihood")),
                    "nature_of_intervention": _clean_cell(row.get("natureOfIntervention")),
                    "effectiveness": _clean_cell(row.get("effectiveness")),
                    "context": actionability_context,
                    "date": _clean_cell(row.get("releaseDate") or row.get("lastUpdated")),
                    "url": _clean_cell(row.get("contextIri")) or self.spec.homepage,
                    "source_id": doc_id,
                    "provenance_url": dict(CLINGEN_ACTIONABILITY_URLS).get(actionability_context, ""),
                }
            )
        return records

    def _csv_dict_rows(self, text: str, first_column: str) -> list[dict[str, str]]:
        reader = csv.reader(io.StringIO(text))
        header: list[str] | None = None
        rows: list[dict[str, str]] = []
        for raw_row in reader:
            row = [cell.strip() for cell in raw_row]
            if not any(row):
                continue
            if header is None:
                if row and row[0] == first_column:
                    header = row
                continue
            if row[0].startswith("+"):
                continue
            normalized = row + [""] * max(0, len(header) - len(row))
            rows.append(dict(zip(header, normalized[: len(header)])))
        if header is None:
            raise ValueError(f"Could not find ClinGen CSV header starting with {first_column!r}.")
        return rows

    def _flat_json_rows(self, payload: Any) -> list[dict[str, Any]]:
        if not isinstance(payload, dict):
            raise ValueError("ClinGen actionability response was not a JSON object.")
        columns = payload.get("columns")
        rows = payload.get("rows")
        if not isinstance(columns, list) or not isinstance(rows, list):
            raise ValueError("ClinGen actionability response did not include columns and rows.")
        normalized_rows: list[dict[str, Any]] = []
        for row in rows:
            if not isinstance(row, list):
                continue
            padded = row + [""] * max(0, len(columns) - len(row))
            normalized_rows.append(dict(zip([str(column) for column in columns], padded[: len(columns)])))
        return normalized_rows

    def _gene_matches(self, value: Any, context: KnowledgeQuery) -> bool:
        return _clean_cell(value).upper() == context.gene.upper()

    def _dedupe_records(self, records: list[dict[str, Any]]) -> list[dict[str, Any]]:
        seen: set[tuple[str, str, str, str]] = set()
        deduped: list[dict[str, Any]] = []
        for record in records:
            key = (
                _clean_cell(record.get("category")),
                _clean_cell(record.get("label")).lower(),
                _clean_cell(record.get("url")),
                _clean_cell(record.get("provenance_url")),
            )
            if key in seen:
                continue
            seen.add(key)
            deduped.append(record)
        return deduped


class PanelAppConnector(BaseConnector):
    """Genomics England PanelApp connector for exact gene panel memberships."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        url = "https://panelapp.genomicsengland.co.uk/api/v1/genes/"
        params = {"entity_name": context.gene}
        try:
            payload = self.client.get_json(
                url,
                params=params,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=[url], started=started)

        raw_results = payload.get("results") if isinstance(payload, dict) else []
        exact_results = [item for item in raw_results or [] if self._matches_gene(item, context.gene)]
        records = [self._record_from_entry(item, context) for item in exact_results[:5]]
        if not records:
            return SourceResult(
                self.spec.key,
                "ok",
                f"Queried PanelApp; no exact gene panel entries found for {context.gene}.",
                [],
                queried_urls=[url],
                elapsed_ms=_elapsed_ms(started),
            )
        return SourceResult(
            self.spec.key,
            "ok",
            f"Queried PanelApp; {len(records)} exact gene panel record(s) returned for {context.gene}.",
            records,
            queried_urls=[url],
            elapsed_ms=_elapsed_ms(started),
        )

    def _matches_gene(self, item: Any, gene: str) -> bool:
        if not isinstance(item, dict):
            return False
        gene_data = item.get("gene_data") if isinstance(item.get("gene_data"), dict) else {}
        candidates = (
            item.get("entity_name"),
            gene_data.get("gene_symbol"),
            gene_data.get("hgnc_symbol"),
        )
        return any(_clean_cell(value).upper() == gene.upper() for value in candidates)

    def _record_from_entry(self, item: dict[str, Any], context: KnowledgeQuery) -> dict[str, Any]:
        gene_data = item.get("gene_data") if isinstance(item.get("gene_data"), dict) else {}
        panel = item.get("panel") if isinstance(item.get("panel"), dict) else {}
        panel_id = _clean_cell(panel.get("id"))
        gene_symbol = _clean_cell(item.get("entity_name")) or _clean_cell(gene_data.get("gene_symbol")) or context.gene
        panel_name = _clean_cell(panel.get("name")) or "PanelApp panel"
        confidence_level = _clean_cell(item.get("confidence_level"))
        confidence_label = self._confidence_label(confidence_level)
        record = {
            "category": "gene_panel",
            "source": self.spec.name,
            "label": f"{gene_symbol} in {panel_name}",
            "summary": self._summary(item, gene_data, panel, gene_symbol, panel_name, confidence_level, confidence_label),
            "source_id": f"{panel_id}:{gene_symbol}" if panel_id else gene_symbol,
            "url": self._panel_gene_url(panel_id, gene_symbol),
            "gene": gene_symbol,
            "gene_name": _clean_cell(gene_data.get("gene_name")),
            "hgnc_id": _clean_cell(gene_data.get("hgnc_id")),
            "omim_gene": self._clean_list(gene_data.get("omim_gene")),
            "aliases": self._clean_list(gene_data.get("alias")),
            "biotype": _clean_cell(gene_data.get("biotype")),
            "panel_id": panel_id,
            "panel_name": panel_name,
            "panel_status": _clean_cell(panel.get("status")),
            "panel_version": _clean_cell(panel.get("version")),
            "panel_version_created": _clean_cell(panel.get("version_created")),
            "panel_disease_group": _clean_cell(panel.get("disease_group")),
            "panel_disease_sub_group": _clean_cell(panel.get("disease_sub_group")),
            "panel_types": self._names(panel.get("types")),
            "relevant_disorders": self._clean_list(panel.get("relevant_disorders")),
            "confidence_level": confidence_level,
            "confidence_label": confidence_label,
            "phenotypes": self._clean_list(item.get("phenotypes")),
            "mode_of_inheritance": _clean_cell(item.get("mode_of_inheritance")),
            "mode_of_pathogenicity": _clean_cell(item.get("mode_of_pathogenicity")),
            "penetrance": _clean_cell(item.get("penetrance")),
            "evidence": self._clean_list(item.get("evidence")),
            "publications": self._clean_list(item.get("publications")),
            "tags": self._clean_list(item.get("tags")),
            "transcripts": self._clean_list(item.get("transcript")),
            "genomic_locations": self._locations(gene_data.get("ensembl_genes")),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _summary(
        self,
        item: dict[str, Any],
        gene_data: dict[str, Any],
        panel: dict[str, Any],
        gene_symbol: str,
        panel_name: str,
        confidence_level: str,
        confidence_label: str,
    ) -> str:
        panel_status = _clean_cell(panel.get("status"))
        panel_version = _clean_cell(panel.get("version"))
        panel_text = panel_name
        if panel_status or panel_version:
            details = ", ".join(part for part in (panel_status, f"v{panel_version}" if panel_version else "") if part)
            panel_text = f"{panel_text} ({details})"

        parts = [f"{gene_symbol} is listed on PanelApp panel {panel_text}"]
        if confidence_level:
            confidence_text = f"confidence {confidence_level}"
            if confidence_label:
                confidence_text = f"{confidence_text} ({confidence_label})"
            parts.append(confidence_text)
        phenotypes = self._clean_list(item.get("phenotypes"))
        if phenotypes:
            parts.append(f"phenotypes: {', '.join(phenotypes[:3])}")
        inheritance = _clean_cell(item.get("mode_of_inheritance"))
        if inheritance:
            parts.append(f"inheritance: {inheritance}")
        penetrance = _clean_cell(item.get("penetrance"))
        if penetrance:
            parts.append(f"penetrance: {penetrance}")
        evidence = self._clean_list(item.get("evidence"))
        if evidence:
            parts.append(f"evidence: {', '.join(evidence[:4])}")
        publications = self._clean_list(item.get("publications"))
        if publications:
            parts.append(f"publications: {', '.join(publications[:4])}")
        panel_types = self._names(panel.get("types"))
        if panel_types:
            parts.append(f"panel type: {', '.join(panel_types[:3])}")
        locations = self._locations(gene_data.get("ensembl_genes"))
        if locations:
            first_location = locations[0]
            location_text = f"{first_location.get('assembly')} {first_location.get('location')}".strip()
            ensembl_id = first_location.get("ensembl_id")
            if ensembl_id:
                location_text = f"{location_text} ({ensembl_id})"
            if location_text:
                parts.append(f"gene location: {location_text}")
        return "; ".join(parts) + "."

    def _confidence_label(self, value: str) -> str:
        labels = {"3": "green/high evidence", "2": "amber/moderate evidence", "1": "red/low evidence"}
        return labels.get(value, "")

    def _panel_gene_url(self, panel_id: str, gene_symbol: str) -> str:
        if panel_id:
            return f"https://panelapp.genomicsengland.co.uk/panels/{quote(panel_id)}/gene/{quote(gene_symbol)}/"
        return f"https://panelapp.genomicsengland.co.uk/api/v1/genes/?entity_name={quote(gene_symbol)}"

    def _clean_list(self, values: Any) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str):
            return [values.strip()] if values.strip() else []
        if not isinstance(values, list):
            values = [values]
        cleaned: list[str] = []
        for value in values:
            text = _clean_cell(value)
            if text and text not in cleaned:
                cleaned.append(text)
        return cleaned

    def _names(self, values: Any) -> list[str]:
        if values is None:
            return []
        if not isinstance(values, list):
            values = [values]
        names: list[str] = []
        for value in values:
            text = _clean_cell(value.get("name") if isinstance(value, dict) else value)
            if text and text not in names:
                names.append(text)
        return names

    def _locations(self, ensembl_genes: Any) -> list[dict[str, str]]:
        if not isinstance(ensembl_genes, dict):
            return []
        rows: list[dict[str, str]] = []
        for assembly, releases in ensembl_genes.items():
            if not isinstance(releases, dict):
                continue
            for release, payload in releases.items():
                if not isinstance(payload, dict):
                    continue
                row = {
                    "assembly": _clean_cell(assembly),
                    "release": _clean_cell(release),
                    "location": _clean_cell(payload.get("location")),
                    "ensembl_id": _clean_cell(payload.get("ensembl_id")),
                }
                if row["location"] or row["ensembl_id"]:
                    rows.append({key: value for key, value in row.items() if value})
        return rows


class MaveDbConnector(BaseConnector):
    """MaveDB connector for published multiplexed functional assay score sets."""

    MAX_SCORE_SETS = 5
    GENE_API_BASE_URL = "https://api.mavedb.org/api/v1/genes"
    SCORE_SET_WEB_BASE_URL = "https://www.mavedb.org/score-sets"

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        gene = context.gene.strip()
        url = f"{self.GENE_API_BASE_URL}/{quote(gene, safe='')}"
        params = {"limit": self.MAX_SCORE_SETS, "offset": 0}
        try:
            payload = self.client.get_json(
                url,
                params=params,
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=[url], started=started)

        if not isinstance(payload, dict):
            message = "MaveDB gene response was not a JSON object."
            return SourceResult(
                self.spec.key,
                "failed",
                message,
                [],
                errors=[message],
                queried_urls=[url],
                elapsed_ms=_elapsed_ms(started),
            )

        raw_score_sets = payload.get("scoreSets")
        score_sets = [item for item in raw_score_sets if isinstance(item, dict)] if isinstance(raw_score_sets, list) else []
        records = [self._record_from_score_set(item, payload, context) for item in score_sets[: self.MAX_SCORE_SETS]]
        records = [record for record in records if record]
        gene_symbol = _clean_cell(payload.get("symbol")) or gene
        if not records:
            return SourceResult(
                self.spec.key,
                "ok",
                f"Queried MaveDB; {gene_symbol} has no published MAVE score sets.",
                [],
                queried_urls=[url],
                elapsed_ms=_elapsed_ms(started),
            )

        total = _clean_cell(payload.get("total")) or str(len(records))
        message = f"Queried MaveDB; {len(records)} published score set record(s) returned for {gene_symbol}."
        if total and total != str(len(records)):
            message = f"{message} {total} total score set(s) are available."
        return SourceResult(
            self.spec.key,
            "ok",
            message,
            records,
            queried_urls=[url],
            elapsed_ms=_elapsed_ms(started),
        )

    def _record_from_score_set(
        self,
        score_set: dict[str, Any],
        gene_payload: dict[str, Any],
        context: KnowledgeQuery,
    ) -> dict[str, Any]:
        urn = _clean_cell(score_set.get("urn"))
        title = _clean_cell(score_set.get("title")) or urn or "MaveDB score set"
        gene_symbol = _clean_cell(gene_payload.get("symbol")) or context.gene
        experiment = score_set.get("experiment") if isinstance(score_set.get("experiment"), dict) else {}
        target_genes = self._target_genes(score_set.get("targetGenes"))
        publications = self._publications(
            list(self._as_list(score_set.get("primaryPublicationIdentifiers")))
            + list(self._as_list(score_set.get("secondaryPublicationIdentifiers")))
        )
        record = {
            "category": "functional_assay_score_set",
            "source": self.spec.name,
            "label": title,
            "summary": self._summary(score_set, gene_payload, experiment, target_genes, publications, urn, title),
            "source_id": urn,
            "url": self._score_set_url(urn) if urn else self.spec.homepage,
            "gene": gene_symbol,
            "gene_name": _clean_cell(gene_payload.get("name")),
            "hgnc_id": _clean_cell(gene_payload.get("hgncId")),
            "omim_id": _clean_cell(gene_payload.get("omimId")),
            "gene_location": _clean_cell(gene_payload.get("location")),
            "locus_group": _clean_cell(gene_payload.get("locusGroup")),
            "score_set_urn": urn,
            "title": title,
            "short_description": _clean_cell(score_set.get("shortDescription")),
            "published_date": _clean_cell(score_set.get("publishedDate")),
            "num_variants": score_set.get("numVariants"),
            "total_gene_score_sets": gene_payload.get("total"),
            "total_gene_scored_variants": gene_payload.get("totalScoredVariants"),
            "experiment_urn": _clean_cell(experiment.get("urn")),
            "experiment_title": _clean_cell(experiment.get("title")),
            "experiment_short_description": _clean_cell(experiment.get("shortDescription")),
            "target_genes": target_genes,
            "publications": publications,
            "license": self._license_text(score_set.get("license")),
            "private": score_set.get("private"),
            "record_type": _clean_cell(score_set.get("recordType")),
        }
        return {key: value for key, value in record.items() if value not in ("", [], {}, None)}

    def _summary(
        self,
        score_set: dict[str, Any],
        gene_payload: dict[str, Any],
        experiment: dict[str, Any],
        target_genes: list[dict[str, Any]],
        publications: list[str],
        urn: str,
        title: str,
    ) -> str:
        gene_symbol = _clean_cell(gene_payload.get("symbol"))
        gene_name = _clean_cell(gene_payload.get("name"))
        gene_text = gene_symbol or "queried gene"
        if gene_name and gene_name.lower() != gene_text.lower():
            gene_text = f"{gene_text} ({gene_name})"
        lead = f"MaveDB score set {urn} for {gene_text}: {title}" if urn else f"MaveDB score set for {gene_text}: {title}"

        parts: list[str] = []
        num_variants = self._format_int(score_set.get("numVariants"))
        if num_variants:
            parts.append(f"{num_variants} scored variants")
        published_date = _clean_cell(score_set.get("publishedDate"))
        if published_date:
            parts.append(f"published {published_date}")
        experiment_title = _clean_cell(experiment.get("title"))
        if experiment_title and experiment_title.lower() != title.lower():
            parts.append(f"experiment: {experiment_title}")
        short_description = _clean_cell(score_set.get("shortDescription"))
        if short_description and short_description.lower() != title.lower():
            parts.append(f"assay summary: {self._clip(short_description, 220)}")
        target_names = [row["name"] for row in target_genes if row.get("name")]
        if target_names:
            parts.append(f"target genes: {', '.join(target_names[:4])}")
        if publications:
            parts.append(f"publications: {', '.join(publications[:4])}")
        license_text = self._license_text(score_set.get("license"))
        if license_text:
            parts.append(f"license {license_text}")
        total_scored = self._format_int(gene_payload.get("totalScoredVariants"))
        if total_scored and total_scored != num_variants:
            parts.append(f"{gene_symbol or 'gene'} total scored variants in MaveDB: {total_scored}")
        return f"{lead}; {'; '.join(parts)}." if parts else f"{lead}."

    def _target_genes(self, values: Any) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for value in self._as_list(values):
            if isinstance(value, dict):
                name = (
                    _clean_cell(value.get("name"))
                    or _clean_cell(value.get("symbol"))
                    or _clean_cell(value.get("mappedHgncName"))
                    or _clean_cell(value.get("hgncSymbol"))
                )
                row = {
                    "name": name,
                    "category": _clean_cell(value.get("category")),
                    "mapped_hgnc_name": _clean_cell(value.get("mappedHgncName")),
                    "uniprot_id": _clean_cell(value.get("uniprotIdFromMappedMetadata") or value.get("uniprotId")),
                    "external_identifiers": self._external_identifiers(value.get("externalIdentifiers")),
                }
            else:
                row = {"name": _clean_cell(value)}
            cleaned = {key: item for key, item in row.items() if item not in ("", [], {}, None)}
            if cleaned and cleaned not in rows:
                rows.append(cleaned)
        return rows

    def _external_identifiers(self, values: Any) -> list[str]:
        identifiers: list[str] = []
        for value in self._as_list(values):
            if isinstance(value, dict):
                identifier_payload = value.get("identifier") if isinstance(value.get("identifier"), dict) else value
                db_name = _clean_cell(identifier_payload.get("dbName") or identifier_payload.get("db_name"))
                identifier = _clean_cell(
                    identifier_payload.get("identifier")
                    or identifier_payload.get("id")
                    or identifier_payload.get("accession")
                )
                if db_name and identifier and not identifier.lower().startswith(f"{db_name.lower()}:"):
                    text = f"{db_name}:{identifier}"
                else:
                    text = identifier or db_name
            else:
                text = _clean_cell(value)
            if text and text not in identifiers:
                identifiers.append(text)
        return identifiers

    def _publications(self, values: list[Any]) -> list[str]:
        publications: list[str] = []
        for value in values:
            if isinstance(value, dict):
                db_name = _clean_cell(value.get("dbName") or value.get("db_name") or value.get("source"))
                identifier = _clean_cell(
                    value.get("identifier")
                    or value.get("id")
                    or value.get("pmid")
                    or value.get("doi")
                    or value.get("accession")
                )
                if db_name.lower() == "pubmed" and identifier:
                    text = f"PMID {identifier}"
                elif db_name and identifier:
                    text = f"{db_name} {identifier}"
                else:
                    text = identifier or _clean_cell(value.get("title"))
            else:
                text = _clean_cell(value)
            if text and text not in publications:
                publications.append(text)
        return publications

    def _license_text(self, value: Any) -> str:
        if isinstance(value, dict):
            short_name = _clean_cell(value.get("shortName") or value.get("short_name") or value.get("name"))
            version = _clean_cell(value.get("version"))
            if short_name and version and version not in short_name:
                return f"{short_name} {version}"
            return short_name or _clean_cell(value.get("longName") or value.get("long_name"))
        return _clean_cell(value)

    def _score_set_url(self, urn: str) -> str:
        return f"{self.SCORE_SET_WEB_BASE_URL}/{quote(urn, safe=':')}"

    def _as_list(self, values: Any) -> list[Any]:
        if values is None:
            return []
        return values if isinstance(values, list) else [values]

    def _format_int(self, value: Any) -> str:
        try:
            return f"{int(value):,}"
        except (TypeError, ValueError):
            return _clean_cell(value)

    def _clip(self, value: str, limit: int) -> str:
        if len(value) <= limit:
            return value
        return value[: limit - 1].rstrip() + "..."


class StringConnector(BaseConnector):
    """Version-pinned STRING functional-association connector for human genes."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        started = time.monotonic()
        url = STRING_V12_NETWORK_URL
        try:
            payload = self.client.get_json(
                url,
                params={
                    "identifiers": context.gene,
                    "species": 9606,
                    "required_score": STRING_REQUIRED_SCORE,
                    "add_nodes": STRING_DIRECT_PARTNER_LIMIT,
                    "caller_identity": "NophiGene",
                },
                headers={"Accept": "application/json"},
                rate_limit_per_second=self.spec.rate_limit_per_second,
            )
        except KnowledgeRequestError as exc:
            return _request_failure_result(self.spec, exc, queried_urls=[url], started=started)
        if not isinstance(payload, list):
            message = "STRING v12 network response was not a JSON list."
            return SourceResult(
                self.spec.key,
                "failed",
                message,
                errors=[message],
                queried_urls=[url],
                elapsed_ms=_elapsed_ms(started),
            )
        records = self._direct_records(payload, context)
        message = (
            f"Queried STRING v12.0; {len(records)} direct functional association(s) returned for {context.gene}."
            if records
            else f"Queried STRING v12.0; no direct functional associations met the score threshold for {context.gene}."
        )
        return SourceResult(
            self.spec.key,
            "ok",
            message,
            records,
            queried_urls=[url],
            elapsed_ms=_elapsed_ms(started),
        )

    def _direct_records(self, payload: list[Any], context: KnowledgeQuery) -> list[dict[str, Any]]:
        query_gene = context.gene.upper()
        records: list[dict[str, Any]] = []
        seen: set[str] = set()
        for item in payload:
            if not isinstance(item, dict):
                continue
            gene_a = _clean_cell(item.get("preferredName_A")).upper()
            gene_b = _clean_cell(item.get("preferredName_B")).upper()
            if gene_a == query_gene:
                partner = gene_b
                source_string_id = _clean_cell(item.get("stringId_A"))
                target_string_id = _clean_cell(item.get("stringId_B"))
            elif gene_b == query_gene:
                partner = gene_a
                source_string_id = _clean_cell(item.get("stringId_B"))
                target_string_id = _clean_cell(item.get("stringId_A"))
            else:
                continue
            if not partner or partner == query_gene or partner in seen:
                continue
            score = self._score(item.get("score"))
            channels = {
                "neighborhood": self._score(item.get("nscore")),
                "fusion": self._score(item.get("fscore")),
                "cooccurrence": self._score(item.get("pscore")),
                "coexpression": self._score(item.get("ascore")),
                "experiments": self._score(item.get("escore")),
                "databases": self._score(item.get("dscore")),
                "text_mining": self._score(item.get("tscore")),
            }
            channel_text = ", ".join(
                f"{name.replace('_', ' ')} {value:.3f}"
                for name, value in channels.items()
                if value is not None and value > 0
            )
            score_text = f"{score:.3f}" if score is not None else "not reported"
            summary = (
                f"STRING v12.0 functional association between {query_gene} and {partner}; combined score "
                f"{score_text}{'; evidence channels: ' + channel_text if channel_text else ''}. "
                "STRING associations do not necessarily indicate direct physical binding."
            )
            record_id = f"string-v12:{source_string_id}:{target_string_id}"
            records.append(
                {
                    "record_id": record_id,
                    "source_id": record_id,
                    "category": "protein_functional_association",
                    "source": self.spec.name,
                    "source_key": self.spec.key,
                    "source_release": "12.0",
                    "label": f"{query_gene}–{partner} STRING functional association",
                    "summary": summary,
                    "gene": query_gene,
                    "source_gene": query_gene,
                    "partner_gene": partner,
                    "target_gene": partner,
                    "edge_type": "functional_association",
                    "interaction_type": "functional_association",
                    "directed": False,
                    "native_score": score,
                    "score": score,
                    "score_type": "STRING combined score (0-1)",
                    "evidence_channels": channels,
                    "association_scope": "functional association; not necessarily direct physical binding",
                    "string_id_source": source_string_id,
                    "string_id_target": target_string_id,
                    "url": f"https://version-12-0.string-db.org/network/{source_string_id}",
                }
            )
            seen.add(partner)
            if len(records) >= STRING_DIRECT_PARTNER_LIMIT:
                break
        records.sort(key=lambda record: (-(record.get("native_score") or 0.0), record["partner_gene"]))
        return records

    @staticmethod
    def _score(value: Any) -> float | None:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None


class LinkoutOfficialConnector(BaseConnector):
    """Official-source linkout connector for open resources without a compact stable API path."""

    def query(self, context: KnowledgeQuery) -> SourceResult:
        record = self._metadata_record()
        term = _literature_query(context)
        record["label"] = f"{self.spec.name} query context for {term}"
        record["summary"] = (
            f"{self.spec.name} is included in the dynamic KB registry. This v1 connector records "
            "source metadata and uses direct linkout/provenance until a compact official API shape is added."
        )
        return SourceResult(
            source_key=self.spec.key,
            status="metadata_only",
            message="Official source registered; linkout metadata recorded.",
            records=[record],
        )


CONNECTOR_CLASSES = {
    "metadata": BaseConnector,
    "auth_metadata": AuthMetadataConnector,
    "licensed_metadata": LicensedMetadataConnector,
    "clinvar": ClinVarConnector,
    "medgen": MedGenConnector,
    "dbsnp": NcbiEutilsConnector,
    "ncbi_gene": NcbiEutilsConnector,
    "pubmed": NcbiEutilsConnector,
    "pmc": NcbiEutilsConnector,
    "geo": NcbiEutilsConnector,
    "litvar": NcbiEutilsConnector,
    "clingen": ClinGenConnector,
    "ensembl": EnsemblConnector,
    "ucsc": UcscConnector,
    "europe_pmc": LiteratureSearchConnector,
    "openalex": LiteratureSearchConnector,
    "crossref": LiteratureSearchConnector,
    "semantic_scholar": LiteratureSearchConnector,
    "biorxiv": LiteratureSearchConnector,
    "medrxiv": LiteratureSearchConnector,
    "gwas_catalog": GwasCatalogConnector,
    "pgs_catalog": PgsCatalogConnector,
    "igsr": IgsrConnector,
    "encode": EncodeConnector,
    "ewas_catalog": EwasCatalogConnector,
    "ewas_atlas": EwasAtlasConnector,
    "screen": ScreenConnector,
    "biostudies": BioStudiesConnector,
    "gnomad": GraphqlConnector,
    "civic": GraphqlConnector,
    "pharmgkb": LinkoutOfficialConnector,
    "mavedb": MaveDbConnector,
    "panelapp": PanelAppConnector,
    "pharmvar": LinkoutOfficialConnector,
    "cpic": LinkoutOfficialConnector,
    "fda_pgx": LinkoutOfficialConnector,
    "dgidb": LinkoutOfficialConnector,
    "string": StringConnector,
}


def connector_for(
    spec: SourceSpec,
    client: RequestClient,
    credential: ResolvedCredential,
) -> BaseConnector:
    """Return a connector instance for one source spec."""
    connector_class = CONNECTOR_CLASSES.get(spec.connector_kind, BaseConnector)
    return connector_class(spec, client, credential)
