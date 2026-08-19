"""Source registry derived from ``resources.txt`` with connector metadata."""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from .imports import ACCEPTED_IMPORT_FORMATS, CANONICAL_IMPORT_FIELDS
from .models import SourceSpec

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESOURCES_PATH = PROJECT_ROOT / "resources.txt"

LANE_LABELS = {
    "clinical": "Clinical variant evidence",
    "population": "Population genetics",
    "regulatory": "Regulatory and epigenomic evidence",
    "literature": "Literature discovery",
    "pharmacogenomics": "Pharmacogenomics and drug response",
    "nutrition_exposome": "Nutrition and exposome",
    "licensed": "Licensed aggregators",
    "interactions": "Gene and protein interactions",
}

DEFAULT_LICENSE_NOTE = (
    "Use according to the source terms. This app records provenance and does not redistribute "
    "licensed full-text or proprietary records."
)

METADATA_ONLY_NOTE = (
    "No safe open programmatic connector is configured. The card is shown for completeness and "
    "returns metadata/linkout status until a legal endpoint is configured."
)

LICENSED_NOTE = (
    "License-gated source. Do not query or ingest proprietary records unless the user has a valid "
    "license and an official endpoint."
)

LIVE_CONNECTOR_KINDS = {
    "clinvar",
    "dbsnp",
    "ncbi_gene",
    "medgen",
    "pubmed",
    "pmc",
    "geo",
    "litvar",
    "ensembl",
    "ucsc",
    "europe_pmc",
    "openalex",
    "crossref",
    "semantic_scholar",
    "biorxiv",
    "medrxiv",
    "gwas_catalog",
    "pgs_catalog",
    "igsr",
    "encode",
    "screen",
    "biostudies",
    "gnomad",
    "civic",
    "pharmgkb",
    "mavedb",
    "panelapp",
    "pharmvar",
    "cpic",
    "fda_pgx",
    "dgidb",
    "string",
}


def _slugify(value: str) -> str:
    text = value.encode("ascii", "ignore").decode("ascii")
    text = text.lower()
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return text or "source"


def _clean_entry_line(line: str) -> str:
    return " ".join(line.strip().split())


def _split_name_description(line: str) -> tuple[str, str]:
    if "\t" in line:
        name, description = line.split("\t", 1)
        return _clean_entry_line(name), _clean_entry_line(description)
    for delimiter in (" - ", " — ", " – "):
        if delimiter in line:
            name, description = line.split(delimiter, 1)
            return _clean_entry_line(name), _clean_entry_line(description)
    return _clean_entry_line(line), ""


def _canonical_key(name: str) -> str:
    lower = name.lower()
    aliases = [
        ("1000 genomes", "igsr"),
        ("pharmgkb", "pharmgkb"),
        ("clinpgx", "pharmgkb"),
        ("encode", "encode"),
        ("civic", "civic"),
        ("database of genomic variants", "dgv"),
        ("dgv", "dgv"),
        ("ncbi gene", "ncbi_gene"),
        ("string", "string"),
        ("omim", "omim"),
        ("genereviews", "genereviews"),
        ("orphanet", "orphanet"),
        ("genecards", "genecards"),
        ("disgenet", "disgenet"),
        ("clinvar", "clinvar"),
        ("clingen", "clingen"),
        ("medgen", "medgen"),
        ("hgmd", "hgmd"),
        ("dbsnp", "dbsnp"),
        ("dbvar", "dbvar"),
        ("gnomad", "gnomad"),
        ("lovd", "lovd"),
        ("decipher", "decipher"),
        ("mavedb", "mavedb"),
        ("gwas catalog", "gwas_catalog"),
        ("pgs catalog", "pgs_catalog"),
        ("gtex", "gtex"),
        ("eqtl catalogue", "eqtl_catalogue"),
        ("panelapp", "panelapp"),
        ("open targets", "opentargets"),
        ("oncokb", "oncokb"),
        ("vicc", "vicc"),
        ("ewas catalog", "ewas_catalog"),
        ("ewas atlas", "ewas_atlas"),
        ("screen", "screen"),
        ("ihec", "ihec"),
        ("roadmap", "roadmap_epigenomics"),
        ("geo", "geo"),
        ("biostudies", "biostudies"),
        ("arrayexpress", "biostudies"),
        ("methbank", "methbank"),
        ("diseasemeth", "diseasemeth"),
        ("epifactors", "epifactors"),
        ("chip-atlas", "chip_atlas"),
        ("europe pmc", "europe_pmc"),
        ("pubmed central", "pmc"),
        ("pmc", "pmc"),
        ("pubmed", "pubmed"),
        ("medline", "pubmed"),
        ("embase", "embase"),
        ("scopus", "scopus"),
        ("web of science", "web_of_science"),
        ("dimensions", "dimensions"),
        ("openalex", "openalex"),
        ("semantic scholar", "semantic_scholar"),
        ("google scholar", "google_scholar"),
        ("crossref", "crossref"),
        ("doaj", "doaj"),
        ("core", "core"),
        ("openaire", "openaire"),
        ("biorxiv", "biorxiv"),
        ("medrxiv", "medrxiv"),
        ("cochrane", "cochrane"),
        ("agricola", "agricola"),
        ("cab abstracts", "cab_abstracts"),
        ("fsta", "fsta"),
        ("litvar", "litvar"),
        ("mastermind", "mastermind"),
        ("varsome", "varsome"),
        ("franklin", "franklin"),
        ("ensembl", "ensembl"),
        ("ucsc", "ucsc"),
        ("nutrigenomedb", "nutrigenomedb"),
        ("dbnp", "dbnp"),
        ("nutrichem", "nutrichem"),
        ("exposome-explorer", "exposome_explorer"),
        ("phenol-explorer", "phenol_explorer"),
        ("foodb", "foodb"),
        ("phytohub", "phytohub"),
        ("cpic", "cpic"),
        ("pharmvar", "pharmvar"),
        ("pharmcat", "pharmcat"),
        ("fda table of pharmacogenetic associations", "fda_pharmacogenetic_associations"),
        ("pharmacogenomic biomarkers", "fda_pharmacogenomic_biomarkers"),
        ("dgidb", "dgidb"),
        ("base", "base"),
    ]
    for needle, key in aliases:
        if needle in lower:
            return key
    return _slugify(name)


def _infer_lane(key: str, name: str) -> str:
    lower = f"{key} {name}".lower()
    if any(token in lower for token in ("drug", "pharm", "cpic", "fda", "dgidb")):
        return "pharmacogenomics"
    if any(token in lower for token in ("food", "nutri", "exposome", "phenol", "phyto")):
        return "nutrition_exposome"
    if any(
        token in lower
        for token in (
            "pubmed",
            "pmc",
            "literature",
            "scholar",
            "crossref",
            "openalex",
            "biorxiv",
            "medrxiv",
            "embase",
            "scopus",
            "cochrane",
            "agricola",
            "abstract",
            "doaj",
            "openaire",
            "dimensions",
            "mastermind",
        )
    ):
        return "literature"
    if any(
        token in lower
        for token in (
            "encode",
            "ewas",
            "eqtl",
            "gtex",
            "screen",
            "ihec",
            "roadmap",
            "geo",
            "biostudies",
            "arrayexpress",
            "meth",
            "epifactor",
            "chip",
            "regulatory",
        )
    ):
        return "regulatory"
    if any(token in lower for token in ("1000", "igsr", "gnomad", "population", "gwas", "pgs", "dbsnp")):
        return "population"
    if any(token in lower for token in ("genecards", "varsome", "franklin", "hgmd", "disgenet")):
        return "licensed"
    if any(token in lower for token in ("string", "reactome", "intact", "biogrid", "interaction")):
        return "interactions"
    return "clinical"


SOURCE_OVERRIDES: dict[str, dict[str, Any]] = {
    "igsr": {
        "connector_kind": "igsr",
        "access_type": "open_download",
        "homepage": "https://www.internationalgenome.org/",
        "license_note": "Open population data; this connector records IGSR/1000 Genomes file provenance and extraction context.",
        "lane": "population",
        "rate_limit_per_second": 1.0,
    },
    "pharmgkb": {
        "connector_kind": "pharmgkb",
        "access_type": "open_api",
        "homepage": "https://www.pharmgkb.org/",
        "license_note": "PharmGKB content is CC BY-SA 4.0; preserve attribution and ShareAlike terms.",
        "lane": "pharmacogenomics",
        "rate_limit_per_second": 2.0,
    },
    "encode": {
        "connector_kind": "encode",
        "access_type": "open_api",
        "homepage": "https://www.encodeproject.org/",
        "lane": "regulatory",
        "rate_limit_per_second": 5.0,
    },
    "civic": {
        "connector_kind": "civic",
        "access_type": "open_api",
        "homepage": "https://civicdb.org/",
        "lane": "clinical",
    },
    "string": {
        "connector_kind": "string",
        "access_type": "open_api",
        "homepage": "https://version-12-0.string-db.org/",
        "license_note": (
            "STRING v12.0 API results are used with source/version provenance. Functional associations do not "
            "necessarily represent direct physical binding."
        ),
        "lane": "interactions",
        "rate_limit_per_second": 1.0,
        "supports_variant": False,
        "supports_gene": True,
        "supports_region": False,
    },
    "ncbi_gene": {
        "connector_kind": "ncbi_gene",
        "access_type": "open_api",
        "homepage": "https://www.ncbi.nlm.nih.gov/gene/",
    },
    "omim": {
        "connector_kind": "auth_metadata",
        "access_type": "auth_api",
        "homepage": "https://www.omim.org/",
        "env_var": "NOPHIGENE_SOURCE_OMIM_TOKEN",
        "license_note": "OMIM API requires registration and commercial use may require a license.",
    },
    "disgenet": {
        "connector_kind": "auth_metadata",
        "access_type": "auth_api",
        "homepage": "https://www.disgenet.org/",
        "env_var": "NOPHIGENE_SOURCE_DISGENET_TOKEN",
        "license_note": "DISGENET access conditions include commercial plans.",
    },
    "genecards": {
        "connector_kind": "auth_metadata",
        "access_type": "auth_api",
        "homepage": "https://www.genecards.org/",
        "env_var": "NOPHIGENE_SOURCE_GENECARDS_TOKEN",
        "license_note": LICENSED_NOTE,
        "lane": "licensed",
        "ingestion_modes": ("official_api", "user_export", "linkout_only"),
    },
    "clinvar": {
        "connector_kind": "clinvar",
        "access_type": "open_api",
        "homepage": "https://www.ncbi.nlm.nih.gov/clinvar/",
    },
    "clingen": {
        "connector_kind": "clingen",
        "access_type": "open_api",
        "homepage": "https://clinicalgenome.org/",
    },
    "medgen": {
        "connector_kind": "medgen",
        "access_type": "open_api",
        "homepage": "https://www.ncbi.nlm.nih.gov/medgen/",
        "lane": "clinical",
    },
    "hgmd": {
        "connector_kind": "licensed_metadata",
        "access_type": "licensed",
        "homepage": "https://www.hgmd.cf.ac.uk/",
        "license_note": LICENSED_NOTE,
    },
    "dbsnp": {
        "connector_kind": "dbsnp",
        "access_type": "open_api",
        "homepage": "https://www.ncbi.nlm.nih.gov/snp/",
        "lane": "population",
    },
    "dbvar": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://www.ncbi.nlm.nih.gov/dbvar/",
    },
    "gnomad": {
        "connector_kind": "gnomad",
        "access_type": "open_api",
        "homepage": "https://gnomad.broadinstitute.org/",
        "lane": "population",
    },
    "lovd": {
        "connector_kind": "metadata",
        "access_type": "mixed",
        "homepage": "https://www.lovd.nl/",
    },
    "mavedb": {
        "connector_kind": "mavedb",
        "access_type": "open_api",
        "homepage": "https://www.mavedb.org/",
    },
    "gwas_catalog": {
        "connector_kind": "gwas_catalog",
        "access_type": "open_api",
        "homepage": "https://www.ebi.ac.uk/gwas/",
        "lane": "population",
    },
    "pgs_catalog": {
        "connector_kind": "pgs_catalog",
        "access_type": "open_api",
        "homepage": "https://www.pgscatalog.org/",
        "lane": "population",
    },
    "gtex": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://gtexportal.org/",
        "lane": "regulatory",
        "license_note": "GTEx is an open contextual source. Until a versioned official query connector is configured, this report records it as not assessed rather than negative evidence.",
    },
    "eqtl_catalogue": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://www.ebi.ac.uk/eqtl/",
        "lane": "regulatory",
        "license_note": "eQTL Catalogue is an open contextual source. Until a versioned official query connector is configured, this report records it as not assessed rather than negative evidence.",
    },
    "panelapp": {
        "connector_kind": "panelapp",
        "access_type": "open_api",
        "homepage": "https://panelapp.genomicsengland.co.uk/",
    },
    "opentargets": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://platform.opentargets.org/",
        "lane": "clinical",
        "license_note": "Open Targets is an open contextual source. Until a versioned official query connector is configured, this report records it as not assessed rather than negative evidence.",
    },
    "oncokb": {
        "connector_kind": "auth_metadata",
        "access_type": "auth_api",
        "homepage": "https://www.oncokb.org/",
        "env_var": "NOPHIGENE_SOURCE_ONCOKB_TOKEN",
        "license_note": "OncoKB API requires registration; commercial/clinical use may require a license.",
    },
    "vicc": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://search.cancervariants.org/",
    },
    "ewas_catalog": {
        "connector_kind": "ewas_catalog",
        "access_type": "open_api",
        "homepage": "https://www.ewascatalog.org/",
        "lane": "regulatory",
        "rate_limit_per_second": 1.0,
    },
    "ewas_atlas": {
        "connector_kind": "ewas_atlas",
        "access_type": "open_api",
        "homepage": "https://ngdc.cncb.ac.cn/ewas/",
        "lane": "regulatory",
        "rate_limit_per_second": 1.0,
    },
    "screen": {
        "connector_kind": "screen",
        "access_type": "open_api",
        "homepage": "https://screen.encodeproject.org/",
        "lane": "regulatory",
        "rate_limit_per_second": 5.0,
    },
    "geo": {
        "connector_kind": "geo",
        "access_type": "open_api",
        "homepage": "https://www.ncbi.nlm.nih.gov/geo/",
        "lane": "regulatory",
    },
    "biostudies": {
        "connector_kind": "biostudies",
        "access_type": "open_api",
        "homepage": "https://www.ebi.ac.uk/biostudies/",
        "lane": "regulatory",
    },
    "chip_atlas": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://chip-atlas.org/",
        "lane": "regulatory",
    },
    "pubmed": {
        "connector_kind": "pubmed",
        "access_type": "open_api",
        "homepage": "https://pubmed.ncbi.nlm.nih.gov/",
        "lane": "literature",
        "supports_literature": True,
    },
    "pmc": {
        "connector_kind": "pmc",
        "access_type": "open_api",
        "homepage": "https://pmc.ncbi.nlm.nih.gov/",
        "lane": "literature",
        "supports_literature": True,
    },
    "europe_pmc": {
        "connector_kind": "europe_pmc",
        "access_type": "open_api",
        "homepage": "https://europepmc.org/",
        "lane": "literature",
        "supports_literature": True,
    },
    "embase": {
        "connector_kind": "licensed_metadata",
        "access_type": "licensed",
        "homepage": "https://www.elsevier.com/products/embase",
        "license_note": LICENSED_NOTE,
        "lane": "literature",
    },
    "scopus": {
        "connector_kind": "licensed_metadata",
        "access_type": "licensed",
        "homepage": "https://www.scopus.com/",
        "license_note": LICENSED_NOTE,
        "lane": "literature",
    },
    "web_of_science": {
        "connector_kind": "licensed_metadata",
        "access_type": "licensed",
        "homepage": "https://www.webofscience.com/",
        "license_note": LICENSED_NOTE,
        "lane": "literature",
    },
    "dimensions": {
        "connector_kind": "metadata",
        "access_type": "mixed",
        "homepage": "https://www.dimensions.ai/",
        "lane": "literature",
    },
    "openalex": {
        "connector_kind": "openalex",
        "access_type": "open_api",
        "homepage": "https://openalex.org/",
        "lane": "literature",
        "supports_literature": True,
    },
    "semantic_scholar": {
        "connector_kind": "semantic_scholar",
        "access_type": "open_api",
        "homepage": "https://www.semanticscholar.org/product/api",
        "lane": "literature",
        "supports_literature": True,
    },
    "google_scholar": {
        "connector_kind": "licensed_metadata",
        "access_type": "no_open_api",
        "homepage": "https://scholar.google.com/",
        "license_note": "No official public API is configured; scraping is intentionally not implemented.",
        "lane": "literature",
    },
    "crossref": {
        "connector_kind": "crossref",
        "access_type": "open_api",
        "homepage": "https://www.crossref.org/",
        "lane": "literature",
        "supports_literature": True,
    },
    "doaj": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://doaj.org/",
        "lane": "literature",
    },
    "core": {
        "connector_kind": "auth_metadata",
        "access_type": "auth_api",
        "homepage": "https://core.ac.uk/",
        "env_var": "NOPHIGENE_SOURCE_CORE_TOKEN",
        "lane": "literature",
    },
    "openaire": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://explore.openaire.eu/",
        "lane": "literature",
    },
    "base": {
        "connector_kind": "metadata",
        "access_type": "open_api",
        "homepage": "https://www.base-search.net/",
        "lane": "literature",
    },
    "biorxiv": {
        "connector_kind": "biorxiv",
        "access_type": "open_api",
        "homepage": "https://www.biorxiv.org/",
        "lane": "literature",
        "supports_literature": True,
    },
    "medrxiv": {
        "connector_kind": "medrxiv",
        "access_type": "open_api",
        "homepage": "https://www.medrxiv.org/",
        "lane": "literature",
        "supports_literature": True,
    },
    "litvar": {
        "connector_kind": "litvar",
        "access_type": "open_api",
        "homepage": "https://www.ncbi.nlm.nih.gov/research/litvar2/",
        "lane": "literature",
        "supports_literature": True,
    },
    "mastermind": {
        "connector_kind": "licensed_metadata",
        "access_type": "licensed",
        "homepage": "https://mastermind.genomenon.com/",
        "license_note": LICENSED_NOTE,
        "lane": "licensed",
    },
    "varsome": {
        "connector_kind": "auth_metadata",
        "access_type": "auth_api",
        "homepage": "https://varsome.com/",
        "env_var": "NOPHIGENE_SOURCE_VARSOME_TOKEN",
        "license_note": LICENSED_NOTE,
        "lane": "licensed",
        "ingestion_modes": ("official_api", "user_export", "linkout_only"),
    },
    "franklin": {
        "connector_kind": "licensed_metadata",
        "access_type": "licensed",
        "homepage": "https://franklin.genoox.com/",
        "license_note": LICENSED_NOTE,
        "lane": "licensed",
    },
    "ensembl": {
        "connector_kind": "ensembl",
        "access_type": "open_api",
        "homepage": "https://rest.ensembl.org/",
        "lane": "clinical",
    },
    "ucsc": {
        "connector_kind": "ucsc",
        "access_type": "open_api",
        "homepage": "https://api.genome.ucsc.edu/",
        "lane": "regulatory",
        "rate_limit_per_second": 1.0,
    },
    "drugbank": {
        "connector_kind": "auth_metadata",
        "access_type": "auth_api",
        "homepage": "https://go.drugbank.com/",
        "env_var": "NOPHIGENE_SOURCE_DRUGBANK_TOKEN",
        "license_note": LICENSED_NOTE,
        "lane": "pharmacogenomics",
        "ingestion_modes": ("official_api", "user_export", "linkout_only"),
        "requires_export": True,
    },
    "cpic": {
        "connector_kind": "cpic",
        "access_type": "open_api",
        "homepage": "https://cpicpgx.org/",
        "lane": "pharmacogenomics",
    },
    "pharmvar": {
        "connector_kind": "pharmvar",
        "access_type": "open_api",
        "homepage": "https://www.pharmvar.org/",
        "lane": "pharmacogenomics",
    },
    "pharmcat": {
        "connector_kind": "metadata",
        "access_type": "open_tool",
        "homepage": "https://pharmcat.org/",
        "lane": "pharmacogenomics",
    },
    "fda_pharmacogenetic_associations": {
        "connector_kind": "fda_pgx",
        "access_type": "open_download",
        "homepage": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
        "lane": "pharmacogenomics",
    },
    "fda_pharmacogenomic_biomarkers": {
        "connector_kind": "fda_pgx",
        "access_type": "open_download",
        "homepage": "https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling",
        "lane": "pharmacogenomics",
    },
    "dgidb": {
        "connector_kind": "dgidb",
        "access_type": "open_api",
        "homepage": "https://dgidb.org/",
        "lane": "pharmacogenomics",
    },
}


def _default_access_type(key: str, name: str) -> str:
    lower = f"{key} {name}".lower()
    if any(token in lower for token in ("hgmd", "embase", "scopus", "web_of_science", "mastermind", "franklin")):
        return "licensed"
    if any(token in lower for token in ("genecards", "drugbank", "omim", "varsome", "oncokb", "disgenet")):
        return "auth_api"
    return "metadata_only"


def _default_ingestion_modes(access_type: str, connector_kind: str) -> tuple[str, ...]:
    if connector_kind == "licensed_metadata" or access_type in {"licensed", "no_open_api"}:
        return ("user_export", "linkout_only")
    if access_type == "auth_api":
        return ("official_api", "linkout_only")
    if connector_kind in LIVE_CONNECTOR_KINDS:
        return ("official_api", "linkout_only")
    return ("linkout_only",)


def _default_requires_export(access_type: str, connector_kind: str, ingestion_modes: tuple[str, ...]) -> bool:
    if "user_export" not in ingestion_modes:
        return False
    return connector_kind == "licensed_metadata" or access_type in {"licensed", "no_open_api"}


def _build_spec(name: str, description: str) -> SourceSpec:
    key = _canonical_key(name)
    override = dict(SOURCE_OVERRIDES.get(key, {}))
    lane = override.pop("lane", _infer_lane(key, name))
    access_type = override.pop("access_type", _default_access_type(key, name))
    connector_kind = override.pop(
        "connector_kind",
        "licensed_metadata" if access_type == "licensed" else "metadata",
    )
    license_note = override.pop(
        "license_note",
        LICENSED_NOTE if access_type == "licensed" else METADATA_ONLY_NOTE if connector_kind == "metadata" else DEFAULT_LICENSE_NOTE,
    )
    env_var = override.pop("env_var", "")
    if access_type in {"auth_api", "licensed"} and not env_var and connector_kind != "licensed_metadata":
        env_var = f"NOPHIGENE_SOURCE_{key.upper()}_TOKEN"
    ingestion_modes = tuple(override.pop("ingestion_modes", _default_ingestion_modes(access_type, connector_kind)))
    requires_export = bool(
        override.pop(
            "requires_export",
            _default_requires_export(access_type, connector_kind, ingestion_modes),
        )
    )
    return SourceSpec(
        key=key,
        name=name,
        description=description,
        lane=lane,
        access_type=access_type,
        connector_kind=connector_kind,
        license_note=license_note,
        env_var=env_var,
        ingestion_modes=ingestion_modes,
        requires_export=requires_export,
        accepted_import_formats=tuple(override.pop("accepted_import_formats", ACCEPTED_IMPORT_FORMATS)),
        import_schema=tuple(override.pop("import_schema", CANONICAL_IMPORT_FIELDS)),
        **override,
    )


@lru_cache(maxsize=1)
def list_source_specs() -> tuple[SourceSpec, ...]:
    """Return one source spec for every active line in ``resources.txt``."""
    if not RESOURCES_PATH.exists():
        return ()
    specs: list[SourceSpec] = []
    seen: set[str] = set()
    for raw_line in RESOURCES_PATH.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        name, description = _split_name_description(line)
        spec = _build_spec(name, description)
        if spec.key in seen:
            continue
        seen.add(spec.key)
        specs.append(spec)
    return tuple(specs)


def get_source_spec(source_key: str) -> SourceSpec | None:
    """Return one source spec by key."""
    normalized = _slugify(source_key)
    for spec in list_source_specs():
        if spec.key == normalized:
            return spec
    return None


def select_source_specs(source_keys: list[str] | tuple[str, ...] | None) -> tuple[SourceSpec, ...]:
    """Return selected sources, defaulting to all resources."""
    specs = list_source_specs()
    if not source_keys:
        return specs
    lookup = {spec.key: spec for spec in specs}
    selected: list[SourceSpec] = []
    seen: set[str] = set()
    for key in source_keys:
        normalized = _slugify(str(key))
        spec = lookup.get(normalized)
        if spec is None or spec.key in seen:
            continue
        seen.add(spec.key)
        selected.append(spec)
    return tuple(selected)


def list_source_cards(
    *,
    selected_keys: list[str] | tuple[str, ...] | None = None,
    credential_statuses: dict[str, str] | None = None,
    import_statuses: dict[str, str] | None = None,
    import_paths: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Return grouped-card payloads for the UI and API."""
    selected = {spec.key for spec in select_source_specs(selected_keys)} if selected_keys else {
        spec.key for spec in list_source_specs()
    }
    credential_statuses = credential_statuses or {}
    import_statuses = import_statuses or {}
    import_paths = import_paths or {}
    return [
        spec.to_card(
            credential_status=credential_statuses.get(spec.key, "not_required" if not spec.env_var else "missing"),
            import_status=import_statuses.get(spec.key, ""),
            import_path=import_paths.get(spec.key, ""),
            selected=spec.key in selected,
        )
        for spec in list_source_specs()
    ]
