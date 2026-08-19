"""Purpose-based workflow presets for dynamic variant knowledge bases."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .local_articles import LOCAL_ARTICLE_SOURCE_KEY, LOCAL_ARTICLE_WORKFLOW_KEY
from .registry import LANE_LABELS, get_source_spec

LOCAL_WORKFLOW_SOURCE_KEYS = {LOCAL_ARTICLE_SOURCE_KEY}


CORE_SAFETY_WORKFLOW_KEYS = (
    "clinical_variant_triage",
    "population_frequency_association",
    "regulatory_epigenomic_context",
    "gene_interaction_network",
    "pharmacogenomics_drug_response",
    "literature_dataset_discovery",
)


@dataclass(frozen=True)
class WorkflowSpec:
    """A deterministic source sequence for one interpretation purpose."""

    key: str
    label: str
    purpose: str
    default_selected: bool
    ordered_source_keys: tuple[str, ...]
    evidence_lanes: tuple[str, ...]
    report_section: str
    requires_vcf: bool = True
    requires_manifest: bool = False
    licensed_notes: tuple[str, ...] = ()

    def valid_source_keys(self) -> tuple[str, ...]:
        """Return source keys that exist in the current source registry."""
        return tuple(
            key for key in self.ordered_source_keys if get_source_spec(key) is not None or key in LOCAL_WORKFLOW_SOURCE_KEYS
        )

    def missing_source_keys(self) -> tuple[str, ...]:
        """Return configured source keys that are not currently registered."""
        return tuple(
            key for key in self.ordered_source_keys if get_source_spec(key) is None and key not in LOCAL_WORKFLOW_SOURCE_KEYS
        )

    def to_card(self, *, selected: bool | None = None) -> dict[str, Any]:
        """Return a compact, redacted UI/API card."""
        valid_source_keys = self.valid_source_keys()
        missing_source_keys = self.missing_source_keys()
        return {
            "key": self.key,
            "label": self.label,
            "purpose": self.purpose,
            "default_selected": self.default_selected,
            "selected": self.default_selected if selected is None else bool(selected),
            "ordered_source_keys": list(valid_source_keys),
            "configured_source_count": len(self.ordered_source_keys),
            "source_count": len(valid_source_keys),
            "missing_source_keys": list(missing_source_keys),
            "evidence_lanes": list(self.evidence_lanes),
            "evidence_lane_labels": [LANE_LABELS.get(lane, lane.replace("_", " ").title()) for lane in self.evidence_lanes],
            "report_section": self.report_section,
            "requires_vcf": self.requires_vcf,
            "requires_manifest": self.requires_manifest,
            "licensed_notes": list(self.licensed_notes),
        }


WORKFLOW_SPECS: tuple[WorkflowSpec, ...] = (
    WorkflowSpec(
        key="clinical_variant_triage",
        label="Clinical Variant Triage",
        purpose=(
            "Prioritize clinical assertions, gene-disease validity, transcript/variant context, "
            "cancer evidence, diagnostic panels, and functional assay evidence."
        ),
        default_selected=True,
        ordered_source_keys=(
            "clinvar",
            "clingen",
            "medgen",
            "ensembl",
            "opentargets",
            "dbsnp",
            "civic",
            "panelapp",
            "mavedb",
            "omim",
            "oncokb",
            "hgmd",
            "varsome",
            "franklin",
        ),
        evidence_lanes=("clinical", "population", "licensed"),
        report_section="Clinical Variant Triage",
        requires_vcf=True,
        requires_manifest=False,
        licensed_notes=(
            "OMIM/OncoKB require configured official credentials where applicable.",
            "HGMD, VarSome, and Franklin use permitted exports or linkout metadata when no licensed API is configured.",
        ),
    ),
    WorkflowSpec(
        key="population_frequency_association",
        label="Population Frequency and Association",
        purpose=(
            "Check population frequency, rsID context, GWAS/PGS associations, and 1000 Genomes/IGSR data-file context "
            "to separate rare-disease plausibility from common-trait context."
        ),
        default_selected=True,
        ordered_source_keys=("gnomad", "dbsnp", "gwas_catalog", "pgs_catalog", "igsr"),
        evidence_lanes=("population",),
        report_section="Population Frequency and Association",
        requires_vcf=True,
        requires_manifest=False,
    ),
    WorkflowSpec(
        key="regulatory_epigenomic_context",
        label="Regulatory and Epigenomic Context",
        purpose=(
            "Combine genome-browser annotations, ENCODE/SCREEN regulatory records, public functional-genomics "
            "studies, EWAS resources, and manifest-derived methylation loci."
        ),
        default_selected=True,
        ordered_source_keys=(
            "ucsc",
            "encode",
            "screen",
            "gtex",
            "eqtl_catalogue",
            "geo",
            "biostudies",
            "ewas_catalog",
            "ewas_atlas",
            "roadmap_epigenomics",
            "ihec",
            "methbank",
            "diseasemeth",
            "epifactors",
            "chip_atlas",
        ),
        evidence_lanes=("regulatory", "literature"),
        report_section="Regulatory and Epigenomic Context",
        requires_vcf=False,
        requires_manifest=True,
        licensed_notes=("Roadmap, IHEC, MethBank, and related resources may return metadata/linkouts when no live connector is configured.",),
    ),
    WorkflowSpec(
        key="gene_interaction_network",
        label="Gene and Protein Interaction Network",
        purpose=(
            "Retrieve version-pinned, source-backed functional association partners for the queried gene and "
            "preserve the source's native confidence and evidence-channel scores."
        ),
        default_selected=True,
        ordered_source_keys=("string",),
        evidence_lanes=("interactions",),
        report_section="Interactions",
        requires_vcf=False,
        requires_manifest=False,
        licensed_notes=(
            "STRING functional associations do not necessarily represent direct physical binding or causal regulation.",
        ),
    ),
    WorkflowSpec(
        key="pharmacogenomics_drug_response",
        label="Pharmacogenomics and Drug Response",
        purpose=(
            "Collect gene-drug, haplotype, FDA label, guideline, and drug-gene interaction context for medication "
            "response and safety review."
        ),
        default_selected=True,
        ordered_source_keys=(
            "pharmgkb",
            "cpic",
            "pharmvar",
            "fda_pharmacogenetic_associations",
            "fda_pharmacogenomic_biomarkers",
            "dgidb",
            "drugbank",
        ),
        evidence_lanes=("pharmacogenomics", "licensed"),
        report_section="Pharmacogenomics and Drug Response",
        requires_vcf=True,
        requires_manifest=False,
        licensed_notes=("DrugBank uses official credentials, permitted exports, or linkout metadata only.",),
    ),
    WorkflowSpec(
        key="literature_dataset_discovery",
        label="Literature and Dataset Discovery",
        purpose=(
            "Search open literature, citation metadata, preprints, variant-literature tools, and dataset catalogs "
            "for publications and study artifacts around the queried gene and variants."
        ),
        default_selected=True,
        ordered_source_keys=(
            "pubmed",
            "pmc",
            "europe_pmc",
            "openalex",
            "crossref",
            "semantic_scholar",
            "biorxiv",
            "medrxiv",
            "litvar",
            "google_scholar",
            "mastermind",
            "embase",
            "scopus",
            "web_of_science",
        ),
        evidence_lanes=("literature", "licensed"),
        report_section="Literature and Dataset Discovery",
        requires_vcf=False,
        requires_manifest=False,
        licensed_notes=(
            "Google Scholar, Mastermind, Embase, Scopus, and Web of Science are API/import/linkout only; no scraping is performed.",
        ),
    ),
    WorkflowSpec(
        key=LOCAL_ARTICLE_WORKFLOW_KEY,
        label="Local PDF Article Evidence",
        purpose=(
            "Extract short gene-relevant evidence snippets from a user-provided folder of legally obtained "
            "scientific article PDFs."
        ),
        default_selected=False,
        ordered_source_keys=(LOCAL_ARTICLE_SOURCE_KEY,),
        evidence_lanes=("literature",),
        report_section="Local Article Evidence",
        requires_vcf=False,
        requires_manifest=False,
        licensed_notes=(
            "Only local user-provided PDFs are parsed. The app stores snippets and provenance, not full article text.",
            "No publisher login automation, DRM bypass, or web scraping is performed.",
        ),
    ),
    WorkflowSpec(
        key="nutrition_exposome_context",
        label="Nutrition and Exposome Context",
        purpose=(
            "Gather nutrition, food-chemistry, phytochemical, environmental-exposure, and nutrigenomics metadata "
            "for exploratory gene-environment interpretation."
        ),
        default_selected=False,
        ordered_source_keys=(
            "foodb",
            "exposome_explorer",
            "phenol_explorer",
            "phytohub",
            "nutrigenomedb",
            "dbnp",
            "nutrichem",
            "fsta",
        ),
        evidence_lanes=("nutrition_exposome", "literature"),
        report_section="Nutrition and Exposome Context",
        requires_vcf=False,
        requires_manifest=False,
        licensed_notes=("FSTA and other subscription indexes use metadata, licensed APIs, or permitted exports only.",),
    ),
    WorkflowSpec(
        key="licensed_aggregator_review",
        label="Licensed Aggregator Review",
        purpose=(
            "Review configured licensed aggregators and subscription indexes as a compliance-safe second pass using "
            "official APIs, user-provided exports, or source linkouts."
        ),
        default_selected=False,
        ordered_source_keys=(
            "genecards",
            "hgmd",
            "varsome",
            "franklin",
            "mastermind",
            "embase",
            "scopus",
            "web_of_science",
        ),
        evidence_lanes=("licensed", "literature", "clinical"),
        report_section="Licensed Aggregator Review",
        requires_vcf=False,
        requires_manifest=False,
        licensed_notes=("Licensed aggregators never use scraping, login automation, CAPTCHA handling, or paywall bypass.",),
    ),
)


def _normalize_key(value: Any) -> str:
    return "_".join(str(value or "").strip().lower().replace("-", "_").split())


def list_workflow_specs() -> tuple[WorkflowSpec, ...]:
    """Return all known workflow presets in deterministic display order."""
    return WORKFLOW_SPECS


def get_workflow_spec(workflow_key: str) -> WorkflowSpec | None:
    """Return one workflow preset by key."""
    normalized = _normalize_key(workflow_key)
    for spec in WORKFLOW_SPECS:
        if spec.key == normalized:
            return spec
    return None


def default_workflow_keys() -> list[str]:
    """Return the Core safety workflow preset keys."""
    return [spec.key for spec in WORKFLOW_SPECS if spec.default_selected]


def select_workflow_specs(workflow_keys: list[str] | tuple[str, ...] | None) -> tuple[WorkflowSpec, ...]:
    """Return selected workflows, defaulting to Core safety presets when omitted."""
    if workflow_keys is None:
        requested = set(default_workflow_keys())
    else:
        requested = {_normalize_key(key) for key in workflow_keys if _normalize_key(key)}
    selected: list[WorkflowSpec] = []
    for spec in WORKFLOW_SPECS:
        if spec.key in requested:
            selected.append(spec)
    return tuple(selected)


def source_keys_for_workflows(workflows: list[WorkflowSpec] | tuple[WorkflowSpec, ...]) -> list[str]:
    """Return the ordered source union for selected workflow presets."""
    ordered: list[str] = []
    seen: set[str] = set()
    for workflow in workflows:
        for source_key in workflow.valid_source_keys():
            if source_key not in seen:
                seen.add(source_key)
                ordered.append(source_key)
    return ordered


def default_workflow_source_keys() -> list[str]:
    """Return the ordered source union for default Core safety workflows."""
    return source_keys_for_workflows(select_workflow_specs(None))


def list_workflow_cards(*, selected_keys: list[str] | tuple[str, ...] | None = None) -> list[dict[str, Any]]:
    """Return UI/API cards for workflow presets."""
    if selected_keys is None:
        selected = set(default_workflow_keys())
    else:
        selected = {_normalize_key(key) for key in selected_keys if _normalize_key(key)}
    return [spec.to_card(selected=spec.key in selected) for spec in WORKFLOW_SPECS]
