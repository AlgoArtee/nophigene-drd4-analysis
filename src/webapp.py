"""Flask-based web UI for the gene-focused analysis workflow."""

from __future__ import annotations

import io
import hashlib
import hmac
import json
import os
import re
import uuid
import traceback
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
from flask import Flask, jsonify, redirect, render_template, request, send_from_directory, session, url_for
from werkzeug.utils import secure_filename

try:
    from .analysis import (
        ANALYSIS_SCOPE_OPTIONS,
        AnalysisError,
        DEFAULT_ANALYSIS_SCOPE,
        DEFAULT_GENE_NAME,
        DEFAULT_REGION,
        DEFAULT_REPORT_NAME,
        GENE_DATA_BUNDLE_PATH,
        GENERAL_ANALYSIS_DATABASE_COLUMNS,
        _prepare_methylation_table_for_output,
        _prepare_variant_table_for_output,
        get_analysis_scope_label,
        get_analysis_scope_slug,
        load_gene_interpretation_database,
        load_variants,
        normalize_analysis_scope,
        run_analysis,
    )
    from .bam_extraction import (
        HG38_FASTA,
        HG38_REFERENCE_DIR,
        ExtractionError,
        default_extracted_vcf_path,
        extract_region_vcf,
        get_extraction_tool_status,
        get_hg38_reference_status,
        prepare_hg38_reference,
    )
    from .gene_region_extraction import find_gene_region
    from .helper_functions.filter_manifest_region import (
        sanitize_gene_name_for_filename,
        save_filtered_manifest,
    )
    from .human_protein_catalog import FEATURED_HUMAN_PROTEIN_QUERIES, get_human_protein_catalog
    from .variant_knowledge.credentials import credential_status_for_specs
    from .variant_knowledge.merger import load_dynamic_knowledge_base
    from .variant_knowledge.orchestrator import build_dynamic_knowledge_base
    from .variant_knowledge.registry import LANE_LABELS, list_source_cards, list_source_specs
    from .variant_knowledge.workflows import (
        default_workflow_keys,
        default_workflow_source_keys,
        list_workflow_cards,
        list_workflow_specs,
        select_workflow_specs,
        source_keys_for_workflows,
    )
    from .workflow import (
        build_scope_regions as build_shared_scope_regions,
        format_region_with_padding as format_shared_region_with_padding,
        genome_build_from_knowledge_base,
        knowledge_base_matches_build,
    )
    from .workbench.reporting import build_canonical_report
    from .workbench.evidence import CORE_SOURCE_KEYS, list_core_source_metadata
except ImportError:
    from analysis import (
        ANALYSIS_SCOPE_OPTIONS,
        AnalysisError,
        DEFAULT_ANALYSIS_SCOPE,
        DEFAULT_GENE_NAME,
        DEFAULT_REGION,
        DEFAULT_REPORT_NAME,
        GENE_DATA_BUNDLE_PATH,
        GENERAL_ANALYSIS_DATABASE_COLUMNS,
        _prepare_methylation_table_for_output,
        _prepare_variant_table_for_output,
        get_analysis_scope_label,
        get_analysis_scope_slug,
        load_gene_interpretation_database,
        load_variants,
        normalize_analysis_scope,
        run_analysis,
    )
    from bam_extraction import (
        HG38_FASTA,
        HG38_REFERENCE_DIR,
        ExtractionError,
        default_extracted_vcf_path,
        extract_region_vcf,
        get_extraction_tool_status,
        get_hg38_reference_status,
        prepare_hg38_reference,
    )
    from gene_region_extraction import find_gene_region
    from helper_functions.filter_manifest_region import sanitize_gene_name_for_filename, save_filtered_manifest
    from human_protein_catalog import FEATURED_HUMAN_PROTEIN_QUERIES, get_human_protein_catalog
    from variant_knowledge.credentials import credential_status_for_specs
    from variant_knowledge.merger import load_dynamic_knowledge_base
    from variant_knowledge.orchestrator import build_dynamic_knowledge_base
    from variant_knowledge.registry import LANE_LABELS, list_source_cards, list_source_specs
    from variant_knowledge.workflows import (
        default_workflow_keys,
        default_workflow_source_keys,
        list_workflow_cards,
        list_workflow_specs,
        select_workflow_specs,
        source_keys_for_workflows,
    )
    from workflow import (
        build_scope_regions as build_shared_scope_regions,
        format_region_with_padding as format_shared_region_with_padding,
        genome_build_from_knowledge_base,
        knowledge_base_matches_build,
    )
    from workbench.reporting import build_canonical_report
    from workbench.evidence import CORE_SOURCE_KEYS, list_core_source_metadata

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
PREPROCESSED_MANIFEST_DIR = Path(__file__).resolve().parent / "gene_data"
SESSION_PREPROCESS_KEY = "preprocess_state"
SESSION_EXTRACTION_KEY = "extraction_state"
SESSION_KNOWLEDGE_SOURCES_KEY = "knowledge_sources_state"
SESSION_KNOWLEDGE_CREDENTIAL_ID_KEY = "knowledge_sources_credential_id"
KNOWLEDGE_SOURCE_CREDENTIAL_STORE: dict[str, dict[str, str]] = {}
VARIANT_RAW_PAGE_SIZE = 25
GENERAL_ANALYSIS_DATABASE_FILENAME = "general_gene_analysis_database.csv"
REPORT_HISTORY_SUFFIXES = {".html", ".json", ".csv"}
GENERIC_PROCESSED_GENE_TOKENS = {
    "analysis",
    "catalog",
    "database",
    "gene",
    "genes",
    "general",
    "lookup",
    "manifest",
    "manifests",
    "methylation",
    "output",
    "preprocessed",
    "processed",
    "protein",
    "proteins",
    "report",
    "reports",
    "result",
    "results",
    "summary",
}
FUNCTIONAL_FAMILY_DEFINITIONS = [
    {
        "key": "longevity",
        "label": "Longevity & Healthy Aging",
        "description": "Aging, senescence, telomeres, genome maintenance, and resilience pathways.",
        "keywords": {
            "longevity": 12,
            "healthy aging": 12,
            "centenarian": 12,
            "aging": 7,
            "senescence": 6,
            "telomere": 6,
            "genome maintenance": 6,
            "progeroid": 6,
            "sirtuin": 5,
        },
    },
    {
        "key": "senses",
        "label": "Senses & Sensory Signaling",
        "description": "Vision, hearing, smell, taste, touch, pain, itch, and sensory transduction.",
        "keywords": {
            "human sensory biology": 14,
            "human senses": 14,
            "sensory transduction": 12,
            "phototransduction": 10,
            "retina": 7,
            "hearing": 7,
            "auditory": 7,
            "olfactory": 7,
            "gustatory": 7,
            "taste receptor": 7,
            "nocicept": 6,
            "pain": 4,
            "itch": 4,
        },
    },
    {
        "key": "asthma_allergy",
        "label": "Asthma, Allergy & Airways",
        "description": "Airway biology, type 2 inflammation, mast cells, eosinophils, and epithelial barriers.",
        "keywords": {
            "asthma": 14,
            "allergic": 12,
            "allergy": 12,
            "atopic": 10,
            "airway": 8,
            "eosinophil": 8,
            "mast cell": 8,
            "ige": 7,
            "type 2 inflammation": 8,
        },
    },
    {
        "key": "cancer",
        "label": "Cancer & Genome Integrity",
        "description": "Tumor biology, DNA repair, cell-cycle control, and cancer predisposition.",
        "keywords": {
            "cancer": 10,
            "tumor": 9,
            "oncology": 9,
            "carcinoma": 8,
            "melanoma": 8,
            "leukemia": 8,
            "dna repair": 7,
            "tumor suppressor": 9,
            "oncogene": 8,
            "malign": 7,
        },
    },
    {
        "key": "metabolism",
        "label": "Metabolism & Endocrinology",
        "description": "Glucose, lipid, energy, hormone, obesity, and diabetes pathways.",
        "keywords": {
            "metabolic": 8,
            "metabolism": 7,
            "diabetes": 9,
            "glucose": 7,
            "insulin": 7,
            "lipid": 7,
            "obesity": 8,
            "endocrine": 7,
            "adipose": 6,
            "cholesterol": 6,
        },
    },
    {
        "key": "neurobiology",
        "label": "Brain, Behavior & Sleep",
        "description": "Neurodevelopment, neurotransmission, cognition, behavior, and sleep biology.",
        "keywords": {
            "neurobiology": 10,
            "neuronal": 8,
            "neuron": 7,
            "brain": 7,
            "schizophrenia": 9,
            "parkinson": 9,
            "alzheimer": 9,
            "neurodevelopment": 8,
            "neurotrans": 7,
            "synaptic": 7,
            "sleep": 8,
            "circadian": 8,
            "behavior": 6,
        },
    },
    {
        "key": "immunity",
        "label": "Immunity & Inflammation",
        "description": "Innate and adaptive immunity, cytokines, infection, and inflammatory signaling.",
        "keywords": {
            "immune": 7,
            "immunity": 7,
            "inflammation": 7,
            "inflammatory": 7,
            "cytokine": 7,
            "interferon": 7,
            "autoimmune": 8,
            "antiviral": 7,
            "pathogen": 6,
            "inflammasome": 8,
        },
    },
    {
        "key": "cardiovascular",
        "label": "Heart, Circulation & Blood",
        "description": "Cardiac rhythm, muscle, vessels, blood pressure, coagulation, and hematology.",
        "keywords": {
            "cardiovascular": 9,
            "cardiac": 8,
            "heart": 7,
            "vascular": 7,
            "blood pressure": 8,
            "coagulation": 8,
            "thromb": 7,
            "cardiomyopathy": 9,
            "arrhythm": 8,
            "hematolog": 7,
        },
    },
    {
        "key": "digestion",
        "label": "Digestion, Liver & Gut",
        "description": "Digestion, intestinal barriers, liver, pancreas, nutrients, and microbiome-facing biology.",
        "keywords": {
            "digestion": 9,
            "digestive": 8,
            "intestinal": 8,
            "gut": 8,
            "gastro": 7,
            "pancrea": 7,
            "liver": 7,
            "hepatic": 7,
            "bowel": 7,
            "microbiome": 7,
            "nutrient": 5,
        },
    },
    {
        "key": "pigmentation_skin",
        "label": "Pigmentation, Skin & Hair",
        "description": "Melanin, melanocytes, epidermal barriers, hair, and ectodermal traits.",
        "keywords": {
            "pigmentation": 11,
            "melanin": 10,
            "melanocyte": 10,
            "skin": 7,
            "epiderm": 8,
            "hair": 7,
            "keratin": 8,
            "ectoderm": 8,
        },
    },
    {
        "key": "mitochondria_stress",
        "label": "Mitochondria & Cellular Stress",
        "description": "Mitochondrial function, redox balance, proteostasis, and cellular stress responses.",
        "keywords": {
            "mitochond": 10,
            "oxidative stress": 8,
            "redox": 8,
            "proteostasis": 7,
            "ferropt": 8,
            "respiratory chain": 8,
            "heat shock": 7,
        },
    },
    {
        "key": "development_reproduction",
        "label": "Development & Reproduction",
        "description": "Embryonic development, cilia, fertility, meiosis, and tissue patterning.",
        "keywords": {
            "development": 6,
            "developmental": 7,
            "embry": 8,
            "fertility": 8,
            "reproduct": 8,
            "meiosis": 9,
            "spermi": 9,
            "ciliary": 8,
            "cilia": 7,
            "laterality": 8,
            "morphogenesis": 7,
        },
    },
    {
        "key": "pharmacogenomics",
        "label": "Pharmacogenomics & Detoxification",
        "description": "Drug response, metabolism, toxicity, and treatment-modifying variation.",
        "keywords": {
            "pharmacogen": 12,
            "drug response": 10,
            "treatment response": 9,
            "medication": 8,
            "drug metabolism": 10,
            "toxicity": 8,
            "detox": 8,
            "adverse": 5,
        },
    },
    {
        "key": "systems",
        "label": "Rare & Systems Biology",
        "description": "Cross-system, rare-disease, and emerging functional biology not captured above.",
        "keywords": {},
    },
]

app = Flask(__name__, template_folder=str(Path(__file__).resolve().parent / "templates"))


def _read_runtime_secret(env_name: str) -> str:
    path_text = os.environ.get(env_name, "").strip()
    if not path_text:
        return ""
    path = Path(path_text)
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        raise RuntimeError(f"Configured runtime secret could not be read: {env_name}") from exc


RUNTIME_SESSION_TOKEN = _read_runtime_secret("NOPHIGENE_SESSION_TOKEN_FILE")
app.secret_key = os.environ.get("NOPHIGENE_SECRET_KEY") or (
    hashlib.sha256(RUNTIME_SESSION_TOKEN.encode("utf-8")).hexdigest()
    if RUNTIME_SESSION_TOKEN
    else "nophigene-local-dev"
)
try:
    from .api import register_api
except ImportError:
    from api import register_api

register_api(app)
app.config.setdefault("NOPHIGENE_BROWSER_TOKEN_USED", False)


@app.before_request
def require_local_authenticated_session():
    """Exchange the launcher token for a session or accept it as an API bearer token."""
    if not RUNTIME_SESSION_TOKEN:
        return None
    if request.path == "/api/v2/health":
        return None
    supplied = str(request.args.get("access_token") or "")
    authorization = str(request.headers.get("Authorization") or "")
    bearer = authorization[7:].strip() if authorization.lower().startswith("bearer ") else ""
    if (
        supplied
        and not app.config["NOPHIGENE_BROWSER_TOKEN_USED"]
        and hmac.compare_digest(supplied, RUNTIME_SESSION_TOKEN)
    ):
        app.config["NOPHIGENE_BROWSER_TOKEN_USED"] = True
        session["nophigene_authenticated"] = True
        clean_args = request.args.to_dict(flat=False)
        clean_args.pop("access_token", None)
        return redirect(request.path, code=303)
    if bearer and hmac.compare_digest(bearer, RUNTIME_SESSION_TOKEN):
        return None
    if session.get("nophigene_authenticated"):
        return None
    return jsonify({"error": {"code": "authentication_required", "message": "Open NophiGene through the authenticated launcher."}}), 401


def _build_app_structure_qa_items() -> list[dict[str, object]]:
    """Return general app-structure Q&A content for the static architecture page."""
    return [
        {
            "question": "How are gene regions retrieved, and how are promoter-only, gene-only, and promoter+gene scopes computed?",
            "answer_lines": [
                (
                    "The preprocessing `Find Region from Gene Name` action starts with public annotation lookup, but the public lookup is used to identify the gene-body interval, not to invent a promoter. "
                    "The app asks NCBI Entrez Gene, Ensembl GRCh37, and UCSC hg19 for candidate intervals for the submitted HGNC-style gene symbol."
                ),
                (
                    "NCBI contributes the Gene `genomicinfo` interval, Ensembl contributes the GRCh37 `lookup/symbol/homo_sapiens/{gene}` interval, and UCSC now uses the assembly-wide `/search` endpoint instead of the older DRD4-only `chr11` prototype. "
                    "UCSC results are filtered to exact leading symbol matches so similarly named or merely related genes do not contaminate the region."
                ),
                (
                    "The public source candidates are still recorded in preprocessing as `candidate_regions`, and the app picks the widest usable public interval as the provisional `selected_region`. "
                    "This widest interval is a practical reconciliation step across RefSeq, Ensembl, UCSC, HGNC, GENCODE, and transcript-track differences."
                ),
                (
                    "Biologically, those public database intervals should be read as gene/transcript-body coordinates. "
                    "In common annotation formats, fields such as UCSC genePred `txStart` and `txEnd`, Ensembl gene `start` and `end`, and NCBI Gene genomic coordinates describe the transcribed locus; the promoter is not automatically part of that interval."
                ),
                (
                    "Because promoter definitions vary by assay and biological question, this app computes the operational promoter locally. "
                    "For curated genes, the local `{gene}_interpretation_db.json` stores `gene_region`, `promoter_review_region`, and `recommended_promoter_plus_gene_region` under `gene_context`."
                ),
                (
                    "The local standard promoter window is strand-aware. "
                    "For plus-strand genes, promoter-only is the 1 kb window immediately before the gene start. "
                    "For reverse-strand genes, promoter-only is the 1 kb window immediately after the gene end, because the transcription start is at the higher genomic coordinate."
                ),
                (
                    "Promoter+gene is therefore not made by simply writing `promoter_start-gene_end`. "
                    "It is the coordinate union of the promoter interval and gene-body interval: `min(all starts and ends)` through `max(all starts and ends)` on the same chromosome. "
                    "This is what fixes reverse-strand genes such as SIRT6, where the valid standard region is `19:4174106-4183560`, not the invalid `19:4182561-4182560`."
                ),
                (
                    "When a curated local interpretation database exists, its validated scope regions override the raw public `selected_region` for the UI. "
                    "The Run Analysis selector then maps `Promoter + gene` to the standard union region, `Promoter only` to `promoter_review_region`, and `Gene only` to `gene_region`."
                ),
                (
                    "When no curated local database exists for a gene, the app falls back conservatively: gene-only is the public selected interval, promoter+gene is the selected interval padded 1 kb upstream from the lower coordinate, and promoter-only is left unavailable because the app does not yet know the strand-specific TSS with enough confidence."
                ),
                (
                    "The filtered methylation manifest is built for the standard promoter+gene region during preprocessing. "
                    "Focused promoter-only or gene-only reports then narrow the loaded manifest rows again to the selected report-focus interval, while the central analysis database remains tied to the standard promoter+gene run."
                ),
                (
                    "The intended lookup policy is therefore: use public databases for robust gene-body coordinate discovery, use the local curated knowledge database for biologically explicit promoter/TSS scope, and store promoter+gene as a validated coordinate union so both plus- and reverse-strand genes behave correctly."
                ),
            ],
        },
        {
            "question": "How are the local databases built from literature, and how does a probe get mapped to a locus or variant?",
            "answer_lines": [
                (
                    "The bundled interpretation and population databases are curated manually from literature first. "
                    "For each supported gene, the app stores a local interpretation JSON, a population JSON, and a "
                    "filtered EPIC manifest CSV in `src/gene_data/`."
                ),
                (
                    "During curation, papers are reviewed and converted into structured entries such as gene-level "
                    "context, curated variant records, bundled evidence links, and a gene-level methylation review area."
                ),
                (
                    "Important clarification: the current whitelist probes are usually not stored as 'these exact EPIC probe IDs were individually named in a paper.' "
                    "Instead, the current generator builds a gene-level shortlist from the bundled EPIC manifest subset for that gene."
                ),
                (
                    "That shortlist is produced in `scripts/generate_curated_gene_knowledge_bases.py` by taking the gene's filtered manifest subset, preferring probes annotated as "
                    "`TSS`, `5'UTR`, or `1stExon`, ranking them by distance to the transcription start site, and then keeping up to the first 10 probe IDs."
                ),
                (
                    "Concrete example: in `src/gene_data/sirt6_interpretation_db.json`, the SIRT6 whitelist contains "
                    "`cg15635336` and `cg09936839`. Those probes come from the bundled SIRT6 EPIC subset and are promoter-proximal/TSS-facing probes chosen by that ranking rule, "
                    "not because the JSON currently records that each exact probe ID was explicitly cited in a publication."
                ),
                (
                    "For SIRT6 on the reverse strand, the transcription start is near genomic position 4,182,560. "
                    "`cg15635336` is at 4,182,521 with a distance of 39 bp and a `1stExon;TSS1500` annotation, while `cg09936839` is at 4,181,854 with a distance of 706 bp and a `Body;TSS1500` annotation. "
                    "That is why those two probes rise to the top of the bundled whitelist."
                ),
                (
                    "After the shortlist exists, the app opens the bundled manifest subset CSV for the same gene. For SIRT6, it reads "
                    "`src/gene_data/SIRT6_epigenetics_hg19.csv`, normalizes `IlmnID`, `CHR`, and `MAPINFO` to "
                    "`probe_id`, `chrom`, and `pos`, and uses those fields to format the probe locus."
                ),
                (
                    "Concrete locus example: the `cg15635336` manifest row carries `CHR=19` and `MAPINFO=4182521`, "
                    "so the UI shows that probe at `chr19:4,182,521`."
                ),
                (
                    "The same manifest row may also carry nearby SNP annotations. For `cg15635336`, the bundled SIRT6 "
                    "manifest subset lists `SNP_ID=rs201182672` and `SNP_DISTANCE=50`, so the UI shows that as a nearby "
                    "manifest locus. That nearby-locus column is purely manifest-derived proximity annotation, not a literature claim."
                ),
                (
                    "There is a second important clarification for the current code: when the local bundles are generated, that same gene-level whitelist is copied into each curated variant record for the gene. "
                    "So the current 'probe to variant' table should be read as 'this variant is being discussed against the same bundled gene-level methylation hotspot' rather than 'this paper proved this exact probe belongs to this exact variant.'"
                ),
                (
                    "The paper column is gathered only from the linked curated variant records' `literature_findings` and `evidence` entries. "
                    "In the SIRT6 example, the table can show Li et al., 2016 for `rs350846` and Simon et al., 2022 for `centSIRT6`, but those papers support the curated variant entries, not necessarily a direct publication-specific assay of `cg15635336` or `cg09936839`."
                ),
                (
                    "In short: literature creates the gene record and the variant records; the whitelist probes are currently a bundled promoter/TSS-focused EPIC subset chosen from the manifest for that gene; "
                    "nearby SNPs come from manifest proximity fields; and the current probe-to-variant links are gene-level bundled associations unless a future bundle adds probe-specific evidence explicitly."
                ),
            ],
        },
        {
            "question": "How do I use the Knowledge Sources tab and dynamic variant knowledge-base preprocessing?",
            "answer_lines": [
                (
                    "The Knowledge Sources tab controls which external databases are used when the app builds a run-specific dynamic variant knowledge base. "
                    "Workflow cards sit above the database cards: the Core safety workflows are checked by default, and those presets select clinical, population, regulatory, pharmacogenomic, and literature/dataset evidence lanes for a normal run."
                ),
                (
                    "To customize a run, open Knowledge Sources, check or uncheck one or more purpose workflows, then fine-tune individual database cards. "
                    "Workflow presets define the source union, while source checkboxes can add or remove individual databases before you click Save Source Settings."
                ),
                (
                    "You can also enter session-only API credentials for sources that advertise an environment variable such as `NOPHIGENE_SOURCE_OMIM_TOKEN`. "
                    "The selected workflow keys and source keys are stored in the Flask session, but submitted tokens are kept only in the server's in-memory credential store for that browser session."
                ),
                (
                    "The cards also show the source access type, connector kind, credential status, import status, and ingestion modes. "
                    "`official_api` means the app can use an official programmatic connector when available, `user_export` means a legally obtained CSV or JSON export can be uploaded, and `linkout_only` means the dynamic KB records source metadata and search context without ingesting proprietary records."
                ),
                (
                    "After source settings are saved, return to Preprocessing. The dynamic KB step requires a resolved gene region and a filtered methylation manifest. "
                    "It also needs a VCF source: the app reuses the Extraction tab's last regional VCF when available, otherwise you can choose an existing VCF in the VCF Source for Dynamic KB field."
                ),
                (
                    "Optional local article evidence is configured in the same preprocessing step. "
                    "Check `Use local PDF article folder`, enter a folder containing legally obtained scientific article PDFs, and choose whether subfolders should be scanned. "
                    "The app extracts only short snippets that mention the queried gene or supplied aliases, then stores normalized summaries and file/path hashes rather than full article text."
                ),
                (
                    "Click Build Variant Knowledge Base after region resolution, methylation subset creation, and variant source selection. "
                    "The builder reads the observed VCF variants in the selected interval, combines them with the filtered EPIC methylation loci, runs selected workflows sequentially, de-duplicates shared database calls with an in-run source-result cache, optionally extracts local PDF article snippets, and writes `variant_kb.json` under `results/dynamic_knowledge_bases/`."
                ),
                (
                    "When analysis runs after a dynamic KB has been built, the app merges the dynamic variant records into the local curated gene bundle for that analysis only. "
                    "This means bundled curated knowledge still works as before, while current run evidence can add extra variant records, workflow runs, provider statuses, literature records, population records, epigenetic locus records, provenance, and license notes."
                ),
                (
                    "Local PDF evidence is also a run artifact. "
                    "It appears in `local_article_evidence`, `literature_records`, workflow summaries, report tables, and sidecar files under the dynamic KB output folder so it is visible in the queried gene's results without becoming a permanent curated bundle."
                ),
                (
                    "The dynamic KB is a run artifact, not a source-controlled curated bundle. "
                    "If a provider fails, needs credentials, needs an export, or only supports linkout metadata, the provider status is recorded and the rest of the selected sources continue; partial failures should not block the whole preprocessing step."
                ),
            ],
        },
        {
            "question": "How do licensed source exports, API credentials, and the API/CLI dynamic KB workflow work?",
            "answer_lines": [
                (
                    "The app does not scrape Google Scholar, HGMD, GeneCards, VarSome, Franklin, Mastermind, Embase, Scopus, Web of Science, or other subscription indexes without an official licensed endpoint. "
                    "For those sources, the supported paths are official APIs when legally configured, user-provided exports, or linkout/status metadata."
                ),
                (
                    "For licensed export ingestion, upload a permitted CSV or JSON file on the relevant Knowledge Sources card. "
                    "The accepted normalized fields are `source_key`, `record_id`, `gene`, `variant`, `rsid`, `title`, `summary`, `assertion`, `phenotype`, `drug`, `score`, `url`, `citation`, `evidence_level`, and `license_note`."
                ),
                (
                    "The import parser intentionally discards unknown columns from uploaded exports. "
                    "Dynamic KB artifacts store evidence summaries and provenance, not raw licensed database dumps. Provenance records include the source key/name, filename hash, file hash, row count, normalized record count, license note, timestamp, warnings, and parse errors."
                ),
                (
                    "Credential and export precedence is conservative: an available official connector with usable records wins first; otherwise a matching user export is used; otherwise the source returns `needs_credentials`, `needs_export`, or `metadata_only` with linkout provenance. "
                    "Secrets are never written to cookies, job manifests, cache files, reports, or `variant_kb.json`."
                ),
                (
                    "The local API exposes `GET /api/v1/knowledge-sources` to list cards, ingestion modes, import schema, and readiness metadata, and `POST /api/v1/knowledge-sources/test` to check whether selected sources are queryable, need credentials, need an export, or are import-ready. "
                    "`GET /api/v1/knowledge-workflows` lists purpose workflow cards, default status, source keys, evidence lanes, and caveats. Workflow jobs can enable dynamic KB building with operation `build_knowledge_bases` or with `options.use_dynamic_knowledge_base` on `analyze` and `full_workflow`."
                ),
                (
                    "API jobs accept `options.knowledge_workflows` and `options.knowledge_sources` as lists or comma-separated strings, and `options.knowledge_source_imports` as a source-key to path map. "
                    "For local article folders, API jobs also accept `options.use_local_article_evidence`, `options.article_pdf_folder`, `options.article_pdf_recursive`, and `options.max_article_pdfs`. "
                    "The command-line builder mirrors that with `scripts/build_dynamic_variant_knowledge_base.py --selected-workflows ... --selected-sources ... --source-import SOURCE_KEY=PATH --use-local-article-evidence --article-pdf-folder PATH`, repeatable for multiple import files."
                ),
                (
                    "The expected practical order is: resolve the gene region, prepare the filtered manifest, extract or choose the VCF, select Knowledge Sources and uploads, build the dynamic KB, then run analysis. "
                    "Analysis reports show dynamic KB availability, Dynamic Workflow Summary sections, provider details per workflow, and the JSON report includes the dynamic KB path/status and workflow metadata so the evidence trail remains auditable."
                ),
            ],
        }
    ]


def _as_relative_display(path: Path) -> str:
    """Render a project-relative path for form fields and UI display."""
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def _resolve_user_path(raw_path: str) -> Path:
    """Resolve a user-entered path against the project root when needed."""
    candidate = Path(raw_path.strip())
    if candidate.is_absolute():
        return candidate
    return PROJECT_ROOT / candidate


def _safe_nonnegative_int(value: Any, default: int, *, maximum: int | None = None) -> int:
    """Parse optional numeric UI fields without throwing away the surrounding form state."""
    try:
        parsed = int(str(value).strip())
    except (TypeError, ValueError):
        parsed = default
    parsed = max(0, parsed)
    if maximum is not None:
        parsed = min(parsed, maximum)
    return parsed


def _format_interval_from_record(record: dict[str, Any] | None, *, default_chrom: str = "") -> str:
    """Format a local knowledge-base interval record as a region string."""
    if not record:
        return ""
    chrom = str(record.get("chrom") or record.get("chromosome") or default_chrom).strip().removeprefix("chr")
    start = record.get("start")
    end = record.get("end")
    if not chrom or start is None or end is None:
        return ""
    try:
        return f"{chrom}:{int(start)}-{int(end)}"
    except (TypeError, ValueError):
        return ""


def _parse_interval_text(region: str) -> tuple[str, int, int] | None:
    """Parse ``chrom:start-end`` text without raising for optional UI state."""
    match = re.fullmatch(
        r"(?:chr)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)",
        str(region).replace(",", "").strip(),
    )
    if match is None:
        return None
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        return None
    return match.group("chrom"), start, end


def _region_covers(candidate_region: str, required_regions: list[str]) -> bool:
    """Return whether one interval covers all required same-chromosome intervals."""
    candidate = _parse_interval_text(candidate_region)
    if candidate is None:
        return False
    candidate_chrom, candidate_start, candidate_end = candidate
    for required_region in required_regions:
        required = _parse_interval_text(required_region)
        if required is None:
            return False
        required_chrom, required_start, required_end = required
        if str(required_chrom).removeprefix("chr") != str(candidate_chrom).removeprefix("chr"):
            return False
        if required_start < candidate_start or required_end > candidate_end:
            return False
    return True


def _format_region_union(regions: list[str]) -> str:
    """Return the coordinate union of valid same-chromosome intervals."""
    parsed_regions = [_parse_interval_text(region) for region in regions if region]
    parsed_regions = [region for region in parsed_regions if region is not None]
    if not parsed_regions:
        return ""
    chrom = str(parsed_regions[0][0]).removeprefix("chr")
    if any(str(region[0]).removeprefix("chr") != chrom for region in parsed_regions):
        return ""
    start = min(region[1] for region in parsed_regions)
    end = max(region[2] for region in parsed_regions)
    return f"{chrom}:{start}-{end}"


def _format_region_with_padding(region: str, upstream_bp: int = 1000) -> str:
    """Return a conservative promoter+gene fallback when no curated scope exists."""
    try:
        return format_shared_region_with_padding(region, upstream_bp=upstream_bp)
    except ValueError:
        return region


def _build_analysis_scope_regions(gene_name: str, selected_gene_region: str) -> dict[str, str]:
    """Build promoter+gene, promoter-only, and gene-only regions for the UI."""
    normalized_gene_name = gene_name.strip().upper() or DEFAULT_GENE_NAME
    knowledge_base = load_gene_interpretation_database(normalized_gene_name)
    genome_build = genome_build_from_knowledge_base(knowledge_base) or "hg19"
    return build_shared_scope_regions(
        normalized_gene_name,
        selected_gene_region,
        genome_build=genome_build,
        knowledge_base=knowledge_base,
    )


def _knowledge_base_uses_genome_build(gene_name: str, genome_build: str) -> bool:
    """Return whether a bundled gene knowledge base explicitly matches a build."""
    knowledge_base = load_gene_interpretation_database(gene_name.strip().upper() or DEFAULT_GENE_NAME)
    try:
        return knowledge_base_matches_build(knowledge_base, genome_build)
    except ValueError:
        return False


def _build_extraction_scope_regions(
    gene_name: str,
    selected_gene_region: str,
    genome_build: str,
) -> dict[str, str]:
    """Build extraction scope regions without reusing hg19-only local intervals."""
    knowledge_base = load_gene_interpretation_database(gene_name.strip().upper() or DEFAULT_GENE_NAME)
    return build_shared_scope_regions(
        gene_name,
        selected_gene_region,
        genome_build=genome_build,
        knowledge_base=knowledge_base,
    )


def _allows_empty_methylation_subset(gene_name: str) -> bool:
    """Return whether preprocessing may unlock analysis with a zero-probe subset."""
    knowledge_base = load_gene_interpretation_database(gene_name.strip().upper() or DEFAULT_GENE_NAME)
    if knowledge_base is None:
        return False

    gene_context = knowledge_base.get("gene_context", {})
    chromosome = str(gene_context.get("chromosome", "")).strip().upper()
    relevant_probe_ids = gene_context.get("relevant_methylation_probe_ids")
    return chromosome in {"M", "MT"} or relevant_probe_ids == []


def _build_analysis_scope_options(preprocess_state: dict[str, Any]) -> list[dict[str, str]]:
    """Return report-focus options with the current gene's regions attached."""
    scope_regions = dict(preprocess_state.get("scope_regions") or {})
    selected_scope = normalize_analysis_scope(str(preprocess_state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE)))
    options: list[dict[str, str]] = []
    for scope_key, scope_config in ANALYSIS_SCOPE_OPTIONS.items():
        region = str(scope_regions.get(scope_key, "")).strip()
        options.append(
            {
                "key": scope_key,
                "label": scope_config["label"],
                "description": scope_config["description"],
                "region": region,
                "output": _default_report_path_for_scope(
                    str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)),
                    scope_key,
                ),
                "selected": "true" if scope_key == selected_scope else "false",
                "disabled": "true" if not region else "false",
            }
        )
    return options


def _default_report_path_for_scope(gene_name: str, analysis_scope: str = DEFAULT_ANALYSIS_SCOPE) -> str:
    """Return the default report path for a gene and report focus."""
    sanitized_gene_name = sanitize_gene_name_for_filename(gene_name).lower()
    scope_slug = get_analysis_scope_slug(analysis_scope)
    return f"results/{sanitized_gene_name}_{scope_slug}_report.html"


def discover_vcf_files() -> list[str]:
    """Find VCF candidates under the mounted data directory."""
    if not DATA_DIR.exists():
        return []
    matches = sorted(DATA_DIR.rglob("*.vcf.gz")) + sorted(DATA_DIR.rglob("*.vcf"))
    return [_as_relative_display(path) for path in matches]


def discover_bam_files() -> list[str]:
    """Find BAM candidates under the mounted data directory."""
    if not DATA_DIR.exists():
        return []
    matches = sorted(DATA_DIR.rglob("*.bam"))
    return [_as_relative_display(path) for path in matches]


def _merge_path_options(*groups: list[str]) -> list[str]:
    """Merge path option lists while preserving order."""
    merged: list[str] = []
    seen: set[str] = set()
    for group in groups:
        for raw_path in group:
            path = str(raw_path).strip()
            if not path or path in seen:
                continue
            merged.append(path)
            seen.add(path)
    return merged


def search_bam_files(search_root: str, *, max_results: int = 200) -> list[str]:
    """Search a user-selected folder tree for BAM files."""
    raw_root = Path(str(search_root).strip()).expanduser()
    root = raw_root if raw_root.is_absolute() else PROJECT_ROOT / raw_root
    if not root.exists():
        raise AnalysisError(f"BAM search folder does not exist: {root}")
    if not root.is_dir():
        raise AnalysisError(f"BAM search path must be a folder: {root}")

    matches: list[str] = []
    for current_root, dir_names, file_names in os.walk(root, onerror=lambda _error: None):
        dir_names[:] = [
            name
            for name in dir_names
            if name not in {".git", ".venv", "__pycache__", ".pytest_cache", ".docker-local"}
        ]
        for file_name in file_names:
            if not file_name.lower().endswith(".bam"):
                continue
            matches.append(_as_relative_display(Path(current_root) / file_name))
            if len(matches) >= max_results:
                return sorted(matches)
    return sorted(matches)


def browse_bam_file(initial_path: str = "") -> str:
    """Open a native file picker and return the selected BAM path."""
    raw_initial_path = Path(str(initial_path or "").strip()).expanduser()
    if not raw_initial_path:
        initial_dir = DATA_DIR if DATA_DIR.exists() else PROJECT_ROOT
    elif raw_initial_path.is_file():
        initial_dir = raw_initial_path.parent
    elif raw_initial_path.is_dir():
        initial_dir = raw_initial_path
    else:
        initial_dir = raw_initial_path.parent if raw_initial_path.parent.exists() else PROJECT_ROOT

    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as exc:
        raise AnalysisError(
            "The native BAM file picker is unavailable in this Python runtime. "
            "Use folder search or enter the BAM path manually."
        ) from exc

    root = tk.Tk()
    root.withdraw()
    try:
        root.attributes("-topmost", True)
    except Exception:
        pass
    try:
        selected_path = filedialog.askopenfilename(
            title="Select GRCh38 BAM file",
            initialdir=str(initial_dir),
            filetypes=(("BAM files", "*.bam"), ("All files", "*.*")),
        )
    finally:
        root.destroy()

    if not selected_path:
        raise AnalysisError("No BAM file was selected.")

    selected = Path(selected_path)
    if selected.suffix.lower() != ".bam":
        raise AnalysisError("Select a file ending in .bam.")
    return _as_relative_display(selected)


def discover_population_stats_files() -> list[str]:
    """Find optional CSV and JSON population statistics files in ``data/``."""
    if not DATA_DIR.exists():
        return []

    include_tokens = (
        "population",
        "popstats",
        "allele",
        "frequency",
        "gnomad",
        "topmed",
        "ancestry",
        "cohort",
        "reference",
    )
    exclude_tokens = (
        "manifest",
        "sample_sheet",
        "samplesheet",
        "processed",
        "methylation",
        "beta_values",
        "control_probes",
        "noob_",
    )

    matches: list[Path] = []
    for path in sorted(DATA_DIR.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".csv", ".json"}:
            continue

        relative_display = _as_relative_display(path).lower()
        file_name = path.name.lower()
        if any(token in relative_display or token in file_name for token in exclude_tokens):
            continue
        if not any(token in relative_display or token in file_name for token in include_tokens):
            continue
        matches.append(path)

    return [_as_relative_display(path) for path in matches]


def discover_manifest_files() -> list[str]:
    """Find EPIC or manifest-like CSV inputs that can be prefiltered for one gene."""
    if not DATA_DIR.exists():
        return []

    matches: list[Path] = []
    for path in DATA_DIR.rglob("*"):
        lower_name = path.name.lower()
        if not path.is_file():
            continue
        if not (
            lower_name.endswith(".csv")
            or lower_name.endswith(".csv.gz")
            or lower_name.endswith(".txt")
        ):
            continue
        if "manifest" not in lower_name and "epic" not in lower_name:
            continue
        matches.append(path)

    return [_as_relative_display(path) for path in sorted(matches)]


def discover_idat_prefixes() -> list[str]:
    """Discover IDAT sample prefixes that have both green and red channels."""
    if not DATA_DIR.exists():
        return []

    prefixes: list[str] = []
    seen: set[str] = set()
    for green_file in sorted(DATA_DIR.rglob("*_Grn.idat")):
        red_file = green_file.with_name(green_file.name.replace("_Grn.idat", "_Red.idat"))
        if not red_file.exists():
            continue

        prefix = green_file.with_name(green_file.name[: -len("_Grn.idat")])
        display = _as_relative_display(prefix)
        if display not in seen:
            seen.add(display)
            prefixes.append(display)
    return prefixes


def _format_artifact_size(num_bytes: int) -> str:
    """Render a compact file-size string for history cards."""
    thresholds = (
        (1024**3, "GB"),
        (1024**2, "MB"),
        (1024, "KB"),
    )
    for threshold, suffix in thresholds:
        if num_bytes >= threshold:
            return f"{num_bytes / threshold:.1f} {suffix}"
    return f"{num_bytes} B"


def _infer_report_label(report_path: Path) -> str:
    """Derive a short, user-facing label from a generated report filename."""
    stem = report_path.stem
    if stem.endswith("_report"):
        gene_name = stem[: -len("_report")]
        if gene_name:
            return gene_name.upper()
    return stem.replace("_", " ").strip() or report_path.name


def _is_report_history_artifact(path: Path) -> bool:
    """Return whether an artifact should appear in report history."""
    if not path.is_file():
        return False
    if path.suffix.lower() not in REPORT_HISTORY_SUFFIXES:
        return False
    if path.name == GENERAL_ANALYSIS_DATABASE_FILENAME:
        return False
    if path.name.lower().endswith("_methylation.csv"):
        return False
    return True


def _normalize_processed_gene_symbol(value: Any) -> str:
    """Normalize a user- or filename-derived gene symbol for matching."""
    normalized = re.sub(r"[^A-Za-z0-9.-]+", "", str(value or "").strip()).upper()
    if not normalized:
        return ""
    if normalized.lower() in GENERIC_PROCESSED_GENE_TOKENS:
        return ""
    return normalized


def _infer_processed_gene_symbol(report_path: Path) -> str:
    """Infer the processed gene symbol represented by a saved artifact path."""
    stem = report_path.stem.strip()
    if not stem:
        return ""

    stem_lower = stem.lower()
    for scope_config in ANALYSIS_SCOPE_OPTIONS.values():
        scope_slug = str(scope_config.get("slug", "")).lower()
        scope_suffix = f"_{scope_slug}_report"
        if scope_slug and stem_lower.endswith(scope_suffix):
            return _normalize_processed_gene_symbol(stem[: -len(scope_suffix)])

    if stem_lower.endswith("_report"):
        return _normalize_processed_gene_symbol(stem[: -len("_report")])

    return _normalize_processed_gene_symbol(stem.split("_", 1)[0])


def _general_analysis_database_path() -> Path:
    """Return the central one-row-per-observed-variant database path used by the UI."""
    return RESULTS_DIR / GENERAL_ANALYSIS_DATABASE_FILENAME


def load_general_analysis_database() -> dict[str, Any]:
    """Load the central analysis database for the dedicated UI tab."""
    database_path = _general_analysis_database_path()
    payload: dict[str, Any] = {
        "path": _as_relative_display(database_path),
        "url": "",
        "exists": database_path.exists(),
        "columns": list(GENERAL_ANALYSIS_DATABASE_COLUMNS),
        "rows": [],
        "row_count": 0,
        "modified_display": "",
        "size_display": "",
        "error": "",
    }

    if not database_path.exists():
        return payload

    try:
        stats = database_path.stat()
        database = pd.read_csv(database_path, dtype=object, keep_default_na=False)
    except Exception as exc:
        payload["error"] = f"Could not read the central database: {exc}"
        return payload

    for column in GENERAL_ANALYSIS_DATABASE_COLUMNS:
        if column not in database.columns:
            database[column] = ""

    primary_columns = list(GENERAL_ANALYSIS_DATABASE_COLUMNS)
    extra_columns = [column for column in database.columns if column not in primary_columns]
    database = database[primary_columns + extra_columns].fillna("")

    payload.update(
        {
            "url": url_for(
                "result_artifact",
                artifact_path=database_path.relative_to(RESULTS_DIR).as_posix(),
            ),
            "columns": primary_columns + extra_columns,
            "rows": database.to_dict(orient="records"),
            "row_count": int(len(database)),
            "modified_display": datetime.fromtimestamp(stats.st_mtime).strftime("%Y-%m-%d %H:%M"),
            "size_display": _format_artifact_size(stats.st_size),
        }
    )
    return payload


def discover_report_history() -> list[dict[str, str]]:
    """List previously generated reports under ``results/`` for the History tab."""
    if not RESULTS_DIR.exists():
        return []

    sortable_entries: list[tuple[float, dict[str, str]]] = []
    for path in RESULTS_DIR.rglob("*"):
        if not _is_report_history_artifact(path):
            continue

        stats = path.stat()
        methylation_output = path.with_name(f"{path.stem}_methylation.csv")
        sortable_entries.append(
            (
                stats.st_mtime,
                {
                    "label": _infer_report_label(path),
                    "report_name": path.name,
                    "report_path": _as_relative_display(path),
                    "report_url": url_for("result_artifact", artifact_path=path.relative_to(RESULTS_DIR).as_posix()),
                    "report_type": path.suffix.removeprefix(".").upper() or "FILE",
                    "processed_gene_symbol": _infer_processed_gene_symbol(path),
                    "modified_display": datetime.fromtimestamp(stats.st_mtime).strftime("%Y-%m-%d %H:%M"),
                    "size_display": _format_artifact_size(stats.st_size),
                    "methylation_path": (
                        _as_relative_display(methylation_output) if methylation_output.exists() else ""
                    ),
                    "methylation_url": (
                        url_for(
                            "result_artifact",
                            artifact_path=methylation_output.relative_to(RESULTS_DIR).as_posix(),
                        )
                        if methylation_output.exists()
                        else ""
                    ),
                },
            )
        )

    sortable_entries.sort(key=lambda item: item[0], reverse=True)
    return [entry for _, entry in sortable_entries]


def _classify_functional_family(search_text: str) -> str:
    """Assign one primary functional family from curated biological context."""
    normalized_text = " ".join(str(search_text or "").lower().split())
    best_key = "systems"
    best_score = 0

    for family in FUNCTIONAL_FAMILY_DEFINITIONS:
        score = sum(
            weight * normalized_text.count(keyword)
            for keyword, weight in family["keywords"].items()
        )
        if score > best_score:
            best_key = str(family["key"])
            best_score = score

    return best_key


@lru_cache(maxsize=1)
def _load_functional_gene_catalog() -> tuple[dict[str, str], ...]:
    """Load concise metadata for every bundled interpretation knowledge base."""
    payloads_by_gene: dict[str, dict[str, Any]] = {}

    if GENE_DATA_BUNDLE_PATH.exists():
        try:
            with zipfile.ZipFile(GENE_DATA_BUNDLE_PATH) as bundle:
                for member_name in bundle.namelist():
                    if not member_name.endswith("_interpretation_db.json"):
                        continue
                    try:
                        payload = json.loads(bundle.read(member_name))
                    except (KeyError, UnicodeDecodeError, json.JSONDecodeError):
                        continue
                    gene_context = payload.get("gene_context", {})
                    gene_name = str(gene_context.get("gene_name", "")).strip().upper()
                    if gene_name:
                        payloads_by_gene[gene_name] = payload
        except (OSError, zipfile.BadZipFile):
            pass

    if PREPROCESSED_MANIFEST_DIR.exists():
        for database_path in PREPROCESSED_MANIFEST_DIR.glob("*_interpretation_db.json"):
            try:
                payload = json.loads(database_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            gene_context = payload.get("gene_context", {})
            gene_name = str(gene_context.get("gene_name", "")).strip().upper()
            if gene_name:
                payloads_by_gene[gene_name] = payload

    catalog: list[dict[str, str]] = []
    for gene_name, payload in payloads_by_gene.items():
        gene_context = payload.get("gene_context", {})
        condition_overview = gene_context.get("condition_research_overview", [])
        if not isinstance(condition_overview, list):
            condition_overview = []
        summary = " ".join(str(gene_context.get("gene_summary", "")).split())
        classification_text = " ".join(
            [
                str(payload.get("database_name", "")),
                summary,
                str(gene_context.get("clinical_context", "")),
                " ".join(str(item) for item in condition_overview),
                " ".join(str(item) for item in gene_context.get("variant_effect_overview", []) or []),
            ]
        )
        gene_region = _format_interval_from_record(
            gene_context.get("gene_region"),
            default_chrom=str(gene_context.get("chromosome", "")),
        )
        catalog.append(
            {
                "gene_name": gene_name,
                "summary": summary or "Curated local interpretation knowledge base.",
                "family_key": _classify_functional_family(classification_text),
                "assembly": str(gene_context.get("assembly", "")).strip() or "Genome build not stated",
                "gene_region": gene_region,
                "database_name": str(payload.get("database_name", "")).strip(),
            }
        )

    catalog.sort(key=lambda item: item["gene_name"])
    return tuple(catalog)


def _build_functional_map(report_history: list[dict[str, str]]) -> dict[str, Any]:
    """Join knowledge-base genes to functional families and generated reports."""
    newest_html_reports: dict[str, dict[str, str]] = {}
    for report in report_history:
        if str(report.get("report_type", "")).upper() != "HTML":
            continue
        gene_name = _normalize_processed_gene_symbol(report.get("processed_gene_symbol", ""))
        if gene_name and gene_name not in newest_html_reports:
            newest_html_reports[gene_name] = report

    genes_by_family: dict[str, list[dict[str, Any]]] = {
        str(family["key"]): [] for family in FUNCTIONAL_FAMILY_DEFINITIONS
    }
    processed_gene_count = 0
    for catalog_gene in _load_functional_gene_catalog():
        gene = dict(catalog_gene)
        report = newest_html_reports.get(gene["gene_name"])
        gene["report_url"] = str(report.get("report_url", "")) if report else ""
        gene["report_name"] = str(report.get("report_name", "")) if report else ""
        gene["processed"] = bool(gene["report_url"])
        if gene["processed"]:
            processed_gene_count += 1
        genes_by_family.setdefault(str(gene["family_key"]), []).append(gene)

    families: list[dict[str, Any]] = []
    for family_definition in FUNCTIONAL_FAMILY_DEFINITIONS:
        family_genes = genes_by_family.get(str(family_definition["key"]), [])
        if not family_genes:
            continue
        families.append(
            {
                "key": family_definition["key"],
                "label": family_definition["label"],
                "description": family_definition["description"],
                "genes": family_genes,
                "gene_count": len(family_genes),
                "processed_count": sum(1 for gene in family_genes if gene["processed"]),
            }
        )

    return {
        "families": families,
        "gene_count": sum(family["gene_count"] for family in families),
        "family_count": len(families),
        "processed_gene_count": processed_gene_count,
    }


def discover_processed_gene_symbols(
    report_history: list[dict[str, str]] | None = None,
    general_database: dict[str, Any] | None = None,
) -> list[str]:
    """Return gene symbols already represented by generated outputs."""
    gene_symbols: set[str] = set()

    if report_history is None:
        if RESULTS_DIR.exists():
            for path in RESULTS_DIR.rglob("*"):
                if _is_report_history_artifact(path):
                    symbol = _infer_processed_gene_symbol(path)
                    if symbol:
                        gene_symbols.add(symbol)
    else:
        for item in report_history:
            symbol = _normalize_processed_gene_symbol(item.get("processed_gene_symbol", ""))
            if not symbol and item.get("report_name"):
                symbol = _infer_processed_gene_symbol(Path(str(item["report_name"])))
            if symbol:
                gene_symbols.add(symbol)

    if general_database is None:
        database_path = _general_analysis_database_path()
        if database_path.exists():
            try:
                database = pd.read_csv(database_path, dtype=object, keep_default_na=False, usecols=["gene"])
            except Exception:
                database = None
            if database is not None:
                for value in database["gene"].dropna().unique():
                    symbol = _normalize_processed_gene_symbol(value)
                    if symbol:
                        gene_symbols.add(symbol)
    else:
        for row in general_database.get("rows", []):
            symbol = _normalize_processed_gene_symbol(row.get("gene") if isinstance(row, dict) else "")
            if symbol:
                gene_symbols.add(symbol)

    return sorted(gene_symbols)


def _protein_matches_processed_gene_symbols(protein: dict[str, Any], processed_gene_symbols: set[str]) -> bool:
    """Return whether a UniProt card belongs to an already processed gene."""
    candidates = {_normalize_processed_gene_symbol(protein.get("gene_name", ""))}
    gene_synonyms = protein.get("gene_synonyms", [])
    if isinstance(gene_synonyms, (list, tuple, set)):
        candidates.update(_normalize_processed_gene_symbol(item) for item in gene_synonyms)
    candidates.discard("")
    return bool(candidates.intersection(processed_gene_symbols))


def _biorender_visuals_for_protein(protein: dict[str, Any]) -> dict[str, Any] | None:
    """Return bundled BioRender figure-starter metadata for a UniProt protein card."""
    candidates = [_normalize_processed_gene_symbol(protein.get("gene_name", ""))]
    gene_synonyms = protein.get("gene_synonyms", [])
    if isinstance(gene_synonyms, (list, tuple, set)):
        candidates.extend(_normalize_processed_gene_symbol(item) for item in gene_synonyms)

    for candidate in candidates:
        if not candidate:
            continue
        knowledge_base = load_gene_interpretation_database(candidate)
        if knowledge_base is None:
            continue
        visuals = knowledge_base.get("gene_context", {}).get("biorender_visuals")
        if isinstance(visuals, dict) and visuals:
            return visuals
    return None


def _render_table(df: pd.DataFrame, rows: int = 12) -> str:
    """Render a compact preview table for the result cards."""
    return df.head(rows).to_html(index=False, classes="preview-table", border=0)


def _prepare_variant_preview_table(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize the visible variant preview so missing source IDs stay readable."""
    return _prepare_variant_table_for_output(df)


def _prepare_methylation_preview_table(df: pd.DataFrame) -> pd.DataFrame:
    """Apply stable methylation column ordering before rendering the preview table."""
    return _prepare_methylation_table_for_output(df)


def _serialize_table_rows(df: pd.DataFrame) -> list[dict[str, Any]]:
    """Convert dataframe rows into JSON-safe dictionaries for client-side pagination."""
    serialized_rows: list[dict[str, Any]] = []

    def _serialize_value(value: Any) -> Any:
        if isinstance(value, (list, tuple, set)):
            return [_serialize_value(item) for item in value]
        if isinstance(value, dict):
            return {str(key): _serialize_value(item) for key, item in value.items()}
        if isinstance(value, (bytes, bytearray)):
            return value.decode("utf-8", errors="replace")
        if hasattr(value, "item"):
            try:
                return value.item()
            except Exception:
                return value
        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
            return value
        return value

    for record in df.to_dict(orient="records"):
        serialized_record: dict[str, Any] = {}
        for key, value in record.items():
            serialized_record[str(key)] = _serialize_value(value)
        serialized_rows.append(serialized_record)
    return serialized_rows


DATA_SOURCE_GROUP_LABELS = {
    "curated": "Bundled curated evidence",
    "variant": "Genetic variant evidence",
    "methylation": "Methylation evidence",
    "population": "Population sources",
    "dynamic": "Dynamic knowledge sources",
    "literature": "Literature and local PDFs",
    "licensed": "Licensed/import-only notices",
}


def _clean_source_text(value: Any, *, limit: int = 700) -> str:
    """Return compact plain text for source cards."""
    text = " ".join(str(value or "").split())
    return text[:limit]


def _as_dict_list(value: Any) -> list[dict[str, Any]]:
    """Return only dictionary items from a possibly mixed list."""
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, dict)]


def _append_source_link(
    links: list[dict[str, str]],
    *,
    label: Any = "",
    url: Any = "",
) -> None:
    """Append one source link when it has either label or URL."""
    clean_label = _clean_source_text(label, limit=180)
    clean_url = _clean_source_text(url, limit=500)
    if not clean_label and not clean_url:
        return
    links.append(
        {
            "label": clean_label or clean_url,
            "url": clean_url,
        }
    )


def _dedupe_source_links(links: list[dict[str, str]]) -> list[dict[str, str]]:
    """Deduplicate links by URL when possible and then by label."""
    deduped: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for link in links:
        label = _clean_source_text(link.get("label", ""), limit=180)
        url = _clean_source_text(link.get("url", ""), limit=500)
        key = (url.casefold(), label.casefold() if not url else "")
        if key in seen:
            continue
        seen.add(key)
        deduped.append({"label": label or url, "url": url})
    return deduped


def _source_links_from_records(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    """Extract commonly shaped evidence links from normalized records."""
    links: list[dict[str, str]] = []
    for record in records:
        _append_source_link(
            links,
            label=record.get("label") or record.get("title") or record.get("paper") or record.get("source"),
            url=record.get("url"),
        )
        for evidence_item in _as_dict_list(record.get("evidence")):
            _append_source_link(
                links,
                label=evidence_item.get("label") or evidence_item.get("title") or evidence_item.get("source"),
                url=evidence_item.get("url"),
            )
        for finding in _as_dict_list(record.get("literature_findings")):
            _append_source_link(
                links,
                label=finding.get("paper") or finding.get("title") or finding.get("label"),
                url=finding.get("url"),
            )
        for link in _as_dict_list(record.get("research_links")):
            _append_source_link(
                links,
                label=link.get("label") or link.get("title"),
                url=link.get("url"),
            )
        for paper in _as_dict_list(record.get("papers")):
            _append_source_link(
                links,
                label=paper.get("label") or paper.get("title") or paper.get("source_variant"),
                url=paper.get("url"),
            )
    return _dedupe_source_links(links)


def _source_card(
    *,
    group: str,
    source_key: str,
    source_name: str,
    status: str,
    summary: str,
    record_count: int = 0,
    findings: list[Any] | None = None,
    links: list[dict[str, str]] | None = None,
    warnings: list[Any] | None = None,
    errors: list[Any] | None = None,
    license_note: str = "",
) -> dict[str, Any]:
    """Return one normalized card for the Data Sources result tab."""
    return {
        "group": group,
        "source_key": source_key,
        "source_name": source_name,
        "status": status,
        "record_count": int(record_count or 0),
        "summary": _clean_source_text(summary),
        "findings": [
            _clean_source_text(finding)
            for finding in (findings or [])
            if _clean_source_text(finding)
        ][:12],
        "links": _dedupe_source_links(links or [])[:20],
        "warnings": [
            _clean_source_text(warning)
            for warning in (warnings or [])
            if _clean_source_text(warning)
        ][:8],
        "errors": [
            _clean_source_text(error)
            for error in (errors or [])
            if _clean_source_text(error)
        ][:8],
        "license_note": _clean_source_text(license_note, limit=500),
    }


def _findings_from_variant_record(record: dict[str, Any]) -> list[str]:
    """Summarize one curated or dynamic variant evidence record."""
    label = _clean_source_text(record.get("display_name") or record.get("variant") or record.get("label"))
    text = _clean_source_text(
        record.get("clinical_significance")
        or record.get("clinical_interpretation")
        or record.get("summary")
        or record.get("assertion")
    )
    if label and text:
        return [f"{label}: {text}"]
    return [label or text] if (label or text) else []


def _source_text_contains(haystack: str, needle: str) -> bool:
    """Return whether compact source text already contains a value."""
    return bool(needle and needle.casefold() in haystack.casefold())


def _clinvar_identity_has_source_id(identity: str, source_id: str) -> bool:
    """Return whether a ClinVar identity already includes the variation ID."""
    if not source_id:
        return False
    normalized = identity.casefold()
    normalized_id = source_id.casefold()
    return (
        normalized == normalized_id
        or f"variation id {normalized_id}" in normalized
        or f"variation {normalized_id}" in normalized
    )


def _clinvar_source_record_identity(record: dict[str, Any]) -> str:
    """Build a compact ClinVar variant identity for source-card findings."""
    label = _clean_source_text(record.get("label") or record.get("title") or record.get("name"), limit=220)
    variant = _clean_source_text(record.get("variant"), limit=120)
    rsid = _clean_source_text(record.get("rsid"), limit=120)
    source_id = _clean_source_text(record.get("source_id"), limit=120)
    identity = label or variant or rsid or (f"Variation ID {source_id}" if source_id else "")
    if not identity:
        return ""

    details: list[str] = []
    if source_id and not _clinvar_identity_has_source_id(identity, source_id):
        details.append(f"Variation ID {source_id}")

    rsid_detail = rsid or (variant if variant.lower().startswith("rs") else "")
    if rsid_detail and not _source_text_contains(identity, rsid_detail):
        details.append(rsid_detail)

    deduped_details: list[str] = []
    for detail in details:
        if any(
            _source_text_contains(existing, detail) or _source_text_contains(detail, existing)
            for existing in deduped_details
        ):
            continue
        deduped_details.append(detail)
    if deduped_details:
        return f"{identity} ({'; '.join(deduped_details)})"
    return identity


def _finding_from_clinvar_source_record(record: dict[str, Any]) -> str:
    """Summarize one ClinVar source record with an explicit variant identity."""
    identity = _clinvar_source_record_identity(record)
    summary = _clean_source_text(
        record.get("summary")
        or record.get("clinical_significance")
        or record.get("clinical_interpretation")
        or record.get("assertion")
        or record.get("title")
        or record.get("label")
    )
    if identity and summary:
        if summary.casefold().startswith(identity.casefold()):
            return summary
        return f"{identity}: {summary}"
    return identity or summary


def _records_by_source(records: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Group normalized dynamic records by source key."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        source_key = _clean_source_text(record.get("source_key") or record.get("source") or "dynamic", limit=120)
        grouped.setdefault(source_key, []).append(record)
    return grouped


def _selected_source_specs(selected_source_keys: list[str]) -> list[Any]:
    """Return source specs in selected order, skipping unknown keys."""
    spec_lookup = {spec.key: spec for spec in list_source_specs()}
    return [spec_lookup[key] for key in selected_source_keys if key in spec_lookup]


def _build_data_sources_payload(
    *,
    knowledge_base: dict[str, Any],
    population_database: dict[str, Any],
    population_insights: dict[str, Any],
    methylation_insights: dict[str, Any],
    dynamic_payload: dict[str, Any] | None,
    selected_source_keys: list[str],
) -> dict[str, Any]:
    """Build compact source provenance cards for the Result Viewer."""
    cards: list[dict[str, Any]] = []
    knowledge_base = knowledge_base or {}
    population_database = population_database or {}
    population_insights = population_insights or {}
    methylation_insights = methylation_insights or {}
    gene_context = knowledge_base.get("gene_context", {}) if isinstance(knowledge_base.get("gene_context"), dict) else {}

    gene_evidence = _as_dict_list(gene_context.get("evidence"))
    gene_findings = [
        gene_context.get("gene_summary"),
        gene_context.get("clinical_context"),
        *list(gene_context.get("variant_effect_overview") or [])[:3],
        *list(gene_context.get("condition_research_overview") or [])[:3],
    ]
    cards.append(
        _source_card(
            group="curated",
            source_key="curated_gene_bundle",
            source_name=str(knowledge_base.get("database_name") or "Bundled gene interpretation database"),
            status="ok",
            record_count=len(knowledge_base.get("variant_records", []) or []),
            summary="Bundled curated gene context used as the stable interpretation base for this analysis.",
            findings=gene_findings,
            links=_source_links_from_records(gene_evidence),
        )
    )

    variant_records = _as_dict_list(knowledge_base.get("variant_records"))
    variant_links = _source_links_from_records(variant_records)
    variant_findings: list[str] = []
    for record in variant_records[:12]:
        variant_findings.extend(_findings_from_variant_record(record))
    cards.append(
        _source_card(
            group="variant",
            source_key="curated_variant_records",
            source_name="Curated variant records",
            status="ok" if variant_records else "metadata_only",
            record_count=len(variant_records),
            summary="Variant-level evidence bundled with the gene interpretation database.",
            findings=variant_findings,
            links=variant_links,
        )
    )

    methylation_records = _as_dict_list(methylation_insights.get("whitelist_probe_reference_rows"))
    methylation_links = _source_links_from_records(methylation_records)
    methylation_links.extend(_source_links_from_records(_as_dict_list(methylation_insights.get("evidence"))))
    methylation_findings = [
        methylation_insights.get("summary"),
        methylation_insights.get("clinical_context"),
        methylation_insights.get("whitelist_literature_context"),
        *list(methylation_insights.get("methylation_effects") or [])[:4],
        *list(methylation_insights.get("methylation_condition_research") or [])[:4],
    ]
    cards.append(
        _source_card(
            group="methylation",
            source_key="methylation_evidence",
            source_name="Methylation evidence",
            status="ok" if methylation_records or methylation_findings else "metadata_only",
            record_count=len(methylation_records),
            summary="Probe and methylation literature context used to interpret the selected gene subset.",
            findings=methylation_findings,
            links=methylation_links,
        )
    )

    population_sources = _as_dict_list(population_database.get("sources"))
    population_sources.extend(_as_dict_list(population_insights.get("sources")))
    population_findings = [
        population_database.get("coverage_note"),
        population_insights.get("summary"),
        population_insights.get("gene_population_patterns_intro"),
    ]
    population_record_count = len(population_insights.get("variant_population_records") or []) + len(
        population_insights.get("gene_population_patterns") or []
    )
    cards.append(
        _source_card(
            group="population",
            source_key="population_database",
            source_name=str(population_database.get("database_name") or "Population database"),
            status="ok" if population_sources or population_record_count else "metadata_only",
            record_count=population_record_count,
            summary="Population reference sources and cohort context available for the analyzed gene.",
            findings=population_findings,
            links=_source_links_from_records(population_sources),
        )
    )

    dynamic_source_records = _as_dict_list((dynamic_payload or {}).get("source_records"))
    dynamic_records_by_source = _records_by_source(dynamic_source_records)
    dynamic_statuses = _as_dict_list((dynamic_payload or {}).get("provider_statuses"))
    status_source_keys: set[str] = set()
    if dynamic_payload:
        for status in dynamic_statuses:
            source_key = _clean_source_text(status.get("source_key"), limit=120)
            if not source_key:
                continue
            status_source_keys.add(source_key)
            source_records = dynamic_records_by_source.get(source_key, [])
            if source_key == "clinvar":
                source_findings = [_finding_from_clinvar_source_record(record) for record in source_records[:8]]
            else:
                source_findings = [
                    record.get("summary") or record.get("title") or record.get("label")
                    for record in source_records[:8]
                ]
            links = _source_links_from_records(source_records)
            _append_source_link(links, label=status.get("name"), url=status.get("homepage"))
            for url in status.get("queried_urls", []) or []:
                _append_source_link(links, label=url, url=url)
            group = "licensed" if status.get("lane") == "licensed" else "dynamic"
            cards.append(
                _source_card(
                    group=group,
                    source_key=source_key,
                    source_name=str(status.get("name") or source_key),
                    status=str(status.get("status") or "metadata_only"),
                    record_count=int(status.get("record_count") or len(source_records)),
                    summary=str(status.get("message") or "Dynamic provider status recorded during preprocessing."),
                    findings=source_findings,
                    links=links,
                    warnings=status.get("warnings") or [],
                    errors=status.get("errors") or [],
                    license_note=str(status.get("license_note") or ""),
                )
            )

        literature_records = _as_dict_list(dynamic_payload.get("literature_records"))
        for source_key, records in _records_by_source(literature_records).items():
            source_name = source_key
            for status in dynamic_statuses:
                if status.get("source_key") == source_key:
                    source_name = str(status.get("name") or source_key)
                    break
            cards.append(
                _source_card(
                    group="literature",
                    source_key=f"{source_key}_literature",
                    source_name=f"{source_name} literature records",
                    status="ok" if records else "metadata_only",
                    record_count=len(records),
                    summary="Literature records returned by dynamic Knowledge Sources preprocessing.",
                    findings=[
                        record.get("summary") or record.get("title") or record.get("label")
                        for record in records[:10]
                    ],
                    links=_source_links_from_records(records),
                )
            )

        local_article_evidence = dynamic_payload.get("local_article_evidence")
        if isinstance(local_article_evidence, dict):
            local_records = _as_dict_list(local_article_evidence.get("records"))
            cards.append(
                _source_card(
                    group="literature",
                    source_key="local_pdf_articles",
                    source_name="Local PDF articles",
                    status=str(local_article_evidence.get("status") or "metadata_only"),
                    record_count=len(local_records),
                    summary=str(local_article_evidence.get("message") or "Local article evidence status."),
                    findings=[
                        record.get("snippet") or record.get("summary") or record.get("title")
                        for record in local_records[:10]
                    ],
                    links=_source_links_from_records(local_records),
                    warnings=(local_article_evidence.get("provenance") or {}).get("warnings", []),
                    errors=(local_article_evidence.get("provenance") or {}).get("errors", []),
                )
            )

    if not dynamic_payload:
        selected_specs = _selected_source_specs(selected_source_keys)
        for spec in selected_specs:
            group = "licensed" if spec.lane == "licensed" or spec.requires_export else "dynamic"
            cards.append(
                _source_card(
                    group=group,
                    source_key=spec.key,
                    source_name=spec.name,
                    status="not_run",
                    record_count=0,
                    summary=(
                        "Selected in Knowledge Sources but not queried for this analysis. "
                        "Run Build Variant Knowledge Base during preprocessing to populate this provider."
                    ),
                    links=[{"label": spec.name, "url": spec.homepage}] if spec.homepage else [],
                    license_note=spec.license_note,
                )
            )
    else:
        for spec in _selected_source_specs(selected_source_keys):
            if spec.key in status_source_keys:
                continue
            group = "licensed" if spec.lane == "licensed" or spec.requires_export else "dynamic"
            cards.append(
                _source_card(
                    group=group,
                    source_key=spec.key,
                    source_name=spec.name,
                    status="not_run",
                    record_count=0,
                    summary="Selected source did not appear in the dynamic KB provider statuses for this run.",
                    links=[{"label": spec.name, "url": spec.homepage}] if spec.homepage else [],
                    license_note=spec.license_note,
                )
            )

    groups: list[dict[str, Any]] = []
    for group_key, label in DATA_SOURCE_GROUP_LABELS.items():
        group_cards = [card for card in cards if card.get("group") == group_key]
        if group_cards:
            groups.append(
                {
                    "key": group_key,
                    "label": label,
                    "cards": group_cards,
                }
            )
    return {
        "summary": (
            f"{len(cards)} source card(s) assembled from bundled evidence"
            + (" and dynamic Knowledge Sources." if dynamic_payload else ". Dynamic Knowledge Sources were not run.")
        ),
        "groups": groups,
        "total_cards": len(cards),
        "dynamic_status": "available" if dynamic_payload else "not_run",
    }


def _empty_form_state() -> dict[str, str]:
    """Return the initial analysis form state used on first page load."""
    vcf_files = discover_vcf_files()
    idat_prefixes = discover_idat_prefixes()
    popstats_files = discover_population_stats_files()

    return {
        "vcf": vcf_files[0] if vcf_files else "",
        "idat": idat_prefixes[0] if idat_prefixes else "",
        "out": _default_report_path_for_scope(DEFAULT_GENE_NAME, DEFAULT_ANALYSIS_SCOPE),
        "region": DEFAULT_REGION,
        "analysis_scope": DEFAULT_ANALYSIS_SCOPE,
        "popstats": "",
        "manifest_file": "",
        "overwrite_general_database": "",
        "suggested_popstats": popstats_files[0] if popstats_files else "",
    }


def _empty_preprocess_state(manifest_files: list[str]) -> dict[str, Any]:
    """Return the default preprocessing state."""
    return {
        "gene_name": DEFAULT_GENE_NAME,
        "region": DEFAULT_REGION,
        "analysis_scope": DEFAULT_ANALYSIS_SCOPE,
        "scope_regions": {
            "promoter_plus_gene": DEFAULT_REGION,
            "promoter_only": "",
            "gene_only": DEFAULT_REGION,
        },
        "scope_region_source": "Default DRD4 promoter+gene scope",
        "manifest_source": manifest_files[0] if manifest_files else "",
        "filtered_manifest": "",
        "region_candidates": [],
        "selected_sources": [],
        "region_ready": False,
        "manifest_ready": False,
        "analysis_ready": False,
        "probe_count": 0,
        "build": "hg19",
        "logs": [],
        "region_recently_updated": False,
        "overwrite_filtered_manifest": False,
        "hg38_extraction_suggested": False,
        "hg38_extraction_message": "",
        "hg38_extraction_region": "",
        "knowledge_vcf_source": "",
        "use_local_article_evidence": False,
        "article_pdf_folder": "",
        "article_pdf_recursive": True,
        "max_article_pdfs": 100,
        "local_article_evidence_status": "",
        "local_article_evidence_record_count": 0,
        "local_article_pdf_count": 0,
        "dynamic_kb_ready": False,
        "dynamic_kb_path": "",
        "dynamic_kb_status": "",
        "dynamic_kb_provider_count": 0,
        "dynamic_kb_source_count": 0,
        "dynamic_kb_workflow_count": 0,
        "dynamic_kb_workflow_summary": [],
    }


def _load_preprocess_state(manifest_files: list[str]) -> dict[str, Any]:
    """Load the persisted preprocessing state from the Flask session."""
    state = _empty_preprocess_state(manifest_files)
    saved_state = session.get(SESSION_PREPROCESS_KEY)
    if isinstance(saved_state, dict):
        state.update(saved_state)

    if not state.get("scope_regions"):
        state["scope_regions"] = {
            "promoter_plus_gene": str(state.get("region", DEFAULT_REGION)),
            "promoter_only": "",
            "gene_only": str(state.get("region", DEFAULT_REGION)),
        }

    if not state.get("manifest_source") and manifest_files:
        state["manifest_source"] = manifest_files[0]
    state["max_article_pdfs"] = _safe_nonnegative_int(
        state.get("max_article_pdfs", 100),
        100,
        maximum=1000,
    )
    return state


def _store_preprocess_state(state: dict[str, Any]) -> None:
    """Persist preprocessing state back to the Flask session."""
    session[SESSION_PREPROCESS_KEY] = {
        "gene_name": str(state.get("gene_name", DEFAULT_GENE_NAME)),
        "region": str(state.get("region", DEFAULT_REGION)),
        "analysis_scope": normalize_analysis_scope(str(state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE))),
        "scope_regions": dict(state.get("scope_regions") or {}),
        "scope_region_source": str(state.get("scope_region_source", "")),
        "manifest_source": str(state.get("manifest_source", "")),
        "filtered_manifest": str(state.get("filtered_manifest", "")),
        "region_candidates": list(state.get("region_candidates", [])),
        "selected_sources": list(state.get("selected_sources", [])),
        "region_ready": bool(state.get("region_ready", False)),
        "manifest_ready": bool(state.get("manifest_ready", False)),
        "analysis_ready": bool(state.get("analysis_ready", False)),
        "probe_count": int(state.get("probe_count", 0)),
        "build": str(state.get("build", "hg19")),
        "logs": list(state.get("logs", []))[-160:],
        "region_recently_updated": bool(state.get("region_recently_updated", False)),
        "overwrite_filtered_manifest": bool(state.get("overwrite_filtered_manifest", False)),
        "hg38_extraction_suggested": bool(state.get("hg38_extraction_suggested", False)),
        "hg38_extraction_message": str(state.get("hg38_extraction_message", "")),
        "hg38_extraction_region": str(state.get("hg38_extraction_region", "")),
        "knowledge_vcf_source": str(state.get("knowledge_vcf_source", "")),
        "use_local_article_evidence": bool(state.get("use_local_article_evidence", False)),
        "article_pdf_folder": str(state.get("article_pdf_folder", "")),
        "article_pdf_recursive": bool(state.get("article_pdf_recursive", True)),
        "max_article_pdfs": _safe_nonnegative_int(state.get("max_article_pdfs", 100), 100, maximum=1000),
        "local_article_evidence_status": str(state.get("local_article_evidence_status", "")),
        "local_article_evidence_record_count": int(state.get("local_article_evidence_record_count", 0) or 0),
        "local_article_pdf_count": int(state.get("local_article_pdf_count", 0) or 0),
        "dynamic_kb_ready": bool(state.get("dynamic_kb_ready", False)),
        "dynamic_kb_path": str(state.get("dynamic_kb_path", "")),
        "dynamic_kb_status": str(state.get("dynamic_kb_status", "")),
        "dynamic_kb_provider_count": int(state.get("dynamic_kb_provider_count", 0) or 0),
        "dynamic_kb_source_count": int(state.get("dynamic_kb_source_count", 0) or 0),
        "dynamic_kb_workflow_count": int(state.get("dynamic_kb_workflow_count", 0) or 0),
        "dynamic_kb_workflow_summary": list(state.get("dynamic_kb_workflow_summary", []))[:20],
    }
    session.modified = True


def _default_extraction_output_display(
    gene_name: str,
    genome_build: str,
    analysis_scope: str,
) -> str:
    """Return the project-relative default extraction output path."""
    return _as_relative_display(default_extracted_vcf_path(gene_name, genome_build, analysis_scope))


def _empty_extraction_state(bam_files: list[str]) -> dict[str, Any]:
    """Return the default Extraction tab state."""
    scope_regions = {
        "promoter_plus_gene": DEFAULT_REGION,
        "promoter_only": "",
        "gene_only": DEFAULT_REGION,
    }
    return {
        "gene_name": DEFAULT_GENE_NAME,
        "genome_build": "hg38",
        "analysis_scope": DEFAULT_ANALYSIS_SCOPE,
        "scope_regions": scope_regions,
        "scope_region_source": "GRCh38 / hg38 extraction default",
        "region": scope_regions[DEFAULT_ANALYSIS_SCOPE],
        "bam_path": bam_files[0] if bam_files else "",
        "bam_search_root": _as_relative_display(DATA_DIR),
        "bam_search_results": [],
        "output_vcf": _default_extraction_output_display(DEFAULT_GENE_NAME, "hg38", DEFAULT_ANALYSIS_SCOPE),
        "reference_dir": _as_relative_display(HG38_REFERENCE_DIR),
        "region_candidates": [],
        "selected_sources": [],
        "region_ready": False,
        "reference_ready": False,
        "extraction_ready": False,
        "logs": [],
        "notice": "",
        "last_output_vcf": "",
        "last_resolved_region": "",
        "last_commands": {},
    }


def _load_extraction_state(bam_files: list[str]) -> dict[str, Any]:
    """Load persisted Extraction tab state from the Flask session."""
    state = _empty_extraction_state(bam_files)
    saved_state = session.get(SESSION_EXTRACTION_KEY)
    if isinstance(saved_state, dict):
        state.update(saved_state)
    if not state.get("bam_path") and bam_files:
        state["bam_path"] = bam_files[0]
    state["bam_search_results"] = _merge_path_options(list(state.get("bam_search_results", [])))
    if not state.get("scope_regions"):
        state["scope_regions"] = {
            "promoter_plus_gene": str(state.get("region", DEFAULT_REGION)),
            "promoter_only": "",
            "gene_only": str(state.get("region", DEFAULT_REGION)),
        }
    if not state.get("output_vcf"):
        state["output_vcf"] = _default_extraction_output_display(
            str(state.get("gene_name", DEFAULT_GENE_NAME)),
            str(state.get("genome_build", "hg38")),
            str(state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE)),
        )
    return state


def _store_extraction_state(state: dict[str, Any]) -> None:
    """Persist Extraction tab state back to the Flask session."""
    session[SESSION_EXTRACTION_KEY] = {
        "gene_name": str(state.get("gene_name", DEFAULT_GENE_NAME)),
        "genome_build": str(state.get("genome_build", "hg38")),
        "analysis_scope": normalize_analysis_scope(str(state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE))),
        "scope_regions": dict(state.get("scope_regions") or {}),
        "scope_region_source": str(state.get("scope_region_source", "")),
        "region": str(state.get("region", "")),
        "bam_path": str(state.get("bam_path", "")),
        "bam_search_root": str(state.get("bam_search_root", _as_relative_display(DATA_DIR))),
        "bam_search_results": list(state.get("bam_search_results", []))[:200],
        "output_vcf": str(state.get("output_vcf", "")),
        "reference_dir": str(state.get("reference_dir", _as_relative_display(HG38_REFERENCE_DIR))),
        "region_candidates": list(state.get("region_candidates", [])),
        "selected_sources": list(state.get("selected_sources", [])),
        "region_ready": bool(state.get("region_ready", False)),
        "reference_ready": bool(state.get("reference_ready", False)),
        "extraction_ready": bool(state.get("extraction_ready", False)),
        "logs": list(state.get("logs", []))[-180:],
        "notice": str(state.get("notice", "")),
        "last_output_vcf": str(state.get("last_output_vcf", "")),
        "last_resolved_region": str(state.get("last_resolved_region", "")),
        "last_commands": dict(state.get("last_commands") or {}),
    }
    session.modified = True


def _knowledge_credential_id() -> str:
    """Return a non-secret key for in-memory source credentials."""
    credential_id = session.get(SESSION_KNOWLEDGE_CREDENTIAL_ID_KEY)
    if not isinstance(credential_id, str) or not credential_id:
        credential_id = uuid.uuid4().hex
        session[SESSION_KNOWLEDGE_CREDENTIAL_ID_KEY] = credential_id
        session.modified = True
    return credential_id


def _session_knowledge_credentials() -> dict[str, str]:
    """Return credentials stored in server memory for this browser session."""
    return KNOWLEDGE_SOURCE_CREDENTIAL_STORE.setdefault(_knowledge_credential_id(), {})


def _empty_knowledge_sources_state() -> dict[str, Any]:
    return {
        "selected_workflows": [],
        "selected_sources": [],
        "external_consents": [],
        "source_imports": {},
        "notice": "",
    }


def _load_knowledge_sources_state() -> dict[str, Any]:
    state = _empty_knowledge_sources_state()
    saved_state = session.get(SESSION_KNOWLEDGE_SOURCES_KEY)
    if isinstance(saved_state, dict):
        state.update(saved_state)
    known_keys = {spec.key for spec in list_source_specs()}
    known_workflows = {workflow.key for workflow in list_workflow_specs()}
    selected_workflows = [
        key for key in state.get("selected_workflows", []) if key in known_workflows
    ]
    state["selected_workflows"] = selected_workflows
    workflow_source_keys = source_keys_for_workflows(select_workflow_specs(selected_workflows))
    selected = [key for key in state.get("selected_sources", []) if key in known_keys]
    state["selected_sources"] = selected or workflow_source_keys
    state["external_consents"] = [
        key for key in state.get("external_consents", []) if key in known_keys
    ]
    imports = state.get("source_imports") if isinstance(state.get("source_imports"), dict) else {}
    state["source_imports"] = {
        str(key): str(path)
        for key, path in imports.items()
        if key in known_keys and str(path).strip()
    }
    return state


def _store_knowledge_sources_state(state: dict[str, Any]) -> None:
    session[SESSION_KNOWLEDGE_SOURCES_KEY] = {
        "selected_workflows": list(state.get("selected_workflows", [])),
        "selected_sources": list(state.get("selected_sources", [])),
        "external_consents": list(state.get("external_consents", [])),
        "source_imports": dict(state.get("source_imports") or {}),
        "notice": str(state.get("notice", "")),
    }
    session.modified = True


def _knowledge_source_import_statuses(source_imports: dict[str, str]) -> dict[str, str]:
    statuses: dict[str, str] = {}
    for spec in list_source_specs():
        import_path = str(source_imports.get(spec.key, "")).strip()
        if import_path:
            statuses[spec.key] = "ready" if _resolve_user_path(import_path).exists() else "missing"
        elif "user_export" in spec.ingestion_modes:
            statuses[spec.key] = "needs_export"
    return statuses


def _build_knowledge_source_cards(state: dict[str, Any]) -> list[dict[str, Any]]:
    specs = [spec for spec in list_source_specs() if spec.key in CORE_SOURCE_KEYS]
    statuses = credential_status_for_specs(specs, _session_knowledge_credentials())
    source_imports = dict(state.get("source_imports") or {})
    cards = [
        card for card in list_source_cards(
        selected_keys=list(state.get("selected_sources", [])),
        credential_statuses=statuses,
        import_statuses=_knowledge_source_import_statuses(source_imports),
        import_paths=source_imports,
        ) if card.get("key") in CORE_SOURCE_KEYS
    ]
    selected = set(state.get("selected_sources", []))
    approved = set(state.get("external_consents", []))
    for card in cards:
        card["external_consent"] = card["key"] in approved
        card["selected"] = card["key"] in selected
        card["counts_as_assessed"] = card.get("connector_kind") not in {
            "metadata", "auth_metadata", "licensed_metadata"
        }
        card["available"] = card["counts_as_assessed"] or "user_export" in card.get("ingestion_modes", [])
    present = {card["key"] for card in cards}
    for metadata in list_core_source_metadata():
        if metadata["key"] in present:
            continue
        cards.append(
            {
                **metadata,
                "connector_kind": "not_installed",
                "selected": False,
                "external_consent": False,
                "counts_as_assessed": False,
                "available": False,
            }
        )
    return cards


def _build_knowledge_workflow_cards(state: dict[str, Any]) -> list[dict[str, Any]]:
    return list_workflow_cards(selected_keys=list(state.get("selected_workflows", [])))


def _save_uploaded_knowledge_source_imports(state: dict[str, Any]) -> int:
    """Persist uploaded licensed-source exports as local runtime artifacts."""
    source_imports = dict(state.get("source_imports") or {})
    upload_count = 0
    target_dir = RESULTS_DIR / "knowledge_source_imports" / _knowledge_credential_id()
    for spec in list_source_specs():
        if "user_export" not in spec.ingestion_modes:
            continue
        upload = request.files.get(f"source_import_{spec.key}")
        if upload is None or not upload.filename:
            continue
        safe_name = secure_filename(upload.filename) or f"{spec.key}.csv"
        target_dir.mkdir(parents=True, exist_ok=True)
        target_path = target_dir / f"{spec.key}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{safe_name}"
        upload.save(target_path)
        source_imports[spec.key] = _as_relative_display(target_path)
        upload_count += 1
    state["source_imports"] = source_imports
    return upload_count


def _group_knowledge_source_cards(cards: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: list[dict[str, Any]] = []
    for lane_key, label in LANE_LABELS.items():
        group_cards = [card for card in cards if card.get("lane") == lane_key]
        if group_cards:
            groups.append({"key": lane_key, "label": label, "cards": group_cards})
    remaining = [
        card
        for card in cards
        if card.get("lane") not in {group["key"] for group in groups}
    ]
    if remaining:
        groups.append({"key": "other", "label": "Other sources", "cards": remaining})
    return groups


def _default_knowledge_vcf_source(
    *,
    preprocess_state: dict[str, Any],
    extraction_state: dict[str, Any],
    vcf_files: list[str],
) -> str:
    """Return the preferred VCF for dynamic KB preprocessing."""
    if extraction_state.get("last_output_vcf"):
        return str(extraction_state["last_output_vcf"])
    if preprocess_state.get("knowledge_vcf_source"):
        return str(preprocess_state["knowledge_vcf_source"])
    return vcf_files[0] if vcf_files else ""


def _append_extraction_log(state: dict[str, Any], message: str, *, stream: str = "stdout") -> None:
    """Append a developer-facing Extraction log line."""
    normalized = str(message).rstrip()
    if not normalized:
        return
    logs = list(state.get("logs", []))
    for line in normalized.splitlines():
        logs.append(f"[{stream}] {line}")
    state["logs"] = logs[-180:]


def _build_extraction_scope_options(extraction_state: dict[str, Any]) -> list[dict[str, str]]:
    """Return extraction focus options with regions attached."""
    scope_regions = dict(extraction_state.get("scope_regions") or {})
    selected_scope = normalize_analysis_scope(str(extraction_state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE)))
    options: list[dict[str, str]] = []
    for scope_key, scope_config in ANALYSIS_SCOPE_OPTIONS.items():
        region = str(scope_regions.get(scope_key, "")).strip()
        options.append(
            {
                "key": scope_key,
                "label": scope_config["label"],
                "description": scope_config["description"],
                "region": region,
                "output_vcf": _default_extraction_output_display(
                    str(extraction_state.get("gene_name", DEFAULT_GENE_NAME)),
                    str(extraction_state.get("genome_build", "hg38")),
                    scope_key,
                ),
                "selected": "true" if scope_key == selected_scope else "false",
                "disabled": "true" if not region else "false",
            }
        )
    return options


def _gene_needs_hg38_extraction(gene_name: str) -> bool:
    """Return whether the local gene record is explicitly hg38/GRCh38-only."""
    knowledge_base = load_gene_interpretation_database(gene_name.strip().upper() or DEFAULT_GENE_NAME)
    if knowledge_base is None:
        return False
    gene_context = knowledge_base.get("gene_context", {})
    assembly = str(gene_context.get("assembly", "")).lower()
    if not ("hg38" in assembly or "grch38" in assembly):
        return False
    return "hg19" not in assembly and "grch37" not in assembly


def _prefill_extraction_state_from_hg38_preprocessing(
    extraction_state: dict[str, Any],
    *,
    gene_name: str,
    scope_regions: dict[str, str],
    selected_sources: list[str],
    region_candidates: list[dict[str, Any]],
) -> None:
    """Pre-populate Extraction after preprocessing detects an hg38-only gene."""
    region = str(scope_regions.get(DEFAULT_ANALYSIS_SCOPE) or scope_regions.get("gene_only") or "").strip()
    extraction_state.update(
        {
            "gene_name": gene_name,
            "genome_build": "hg38",
            "analysis_scope": DEFAULT_ANALYSIS_SCOPE,
            "scope_regions": dict(scope_regions),
            "scope_region_source": "Local GRCh38 curated promoter/gene intervals",
            "region": region,
            "output_vcf": _default_extraction_output_display(gene_name, "hg38", DEFAULT_ANALYSIS_SCOPE),
            "selected_sources": list(selected_sources),
            "region_candidates": list(region_candidates),
            "region_ready": bool(region),
            "extraction_ready": False,
            "last_output_vcf": "",
            "last_resolved_region": "",
            "last_commands": {},
        }
    )


def _append_preprocess_log(state: dict[str, Any], message: str, *, stream: str = "stdout") -> None:
    """Append a timestamp-free developer log line to the preprocessing state."""
    normalized = str(message).rstrip()
    if not normalized:
        return

    logs = list(state.get("logs", []))
    for line in normalized.splitlines():
        logs.append(f"[{stream}] {line}")
    state["logs"] = logs[-160:]


def _capture_preprocess_call(operation, *args, **kwargs):
    """Run a preprocessing helper while capturing any stdout/stderr it emits."""
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
        result = operation(*args, **kwargs)
    return result, stdout_buffer.getvalue(), stderr_buffer.getvalue()


def _expected_filtered_manifest_path(gene_name: str, genome_build: str) -> Path:
    """Return the canonical saved subset path for one gene/build combination."""
    return PREPROCESSED_MANIFEST_DIR / (
        f"{sanitize_gene_name_for_filename(gene_name)}_epigenetics_{genome_build}.csv"
    )


def _summarize_existing_filtered_manifest(output_path: Path) -> dict[str, object]:
    """Read an existing filtered manifest so it can be reused without rewriting it."""
    existing_df = pd.read_csv(output_path, low_memory=False)
    return {
        "output_path": output_path,
        "probe_count": int(len(existing_df)),
    }


def _filtered_manifest_metadata_path(output_path: Path) -> Path:
    """Return the small metadata sidecar for a filtered manifest subset."""
    return output_path.with_suffix(output_path.suffix + ".meta.json")


def _write_filtered_manifest_metadata(
    output_path: Path,
    *,
    gene_name: str,
    region: str,
    analysis_scope: str,
) -> None:
    """Persist the exact region used to create a filtered manifest subset."""
    metadata = {
        "gene_name": gene_name,
        "region": region,
        "analysis_scope": normalize_analysis_scope(analysis_scope),
    }
    _filtered_manifest_metadata_path(output_path).write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def _filtered_manifest_covers_region(output_path: Path, region: str) -> bool:
    """Return whether an existing subset appears to span the requested interval."""
    try:
        metadata = json.loads(_filtered_manifest_metadata_path(output_path).read_text(encoding="utf-8"))
    except Exception:
        return False

    return str(metadata.get("region", "")).replace(",", "") == str(region).replace(",", "")


def _apply_preprocessing_defaults(form: dict[str, str], preprocess_state: dict[str, Any]) -> None:
    """Propagate preprocessing results into the analysis form defaults."""
    analysis_scope = normalize_analysis_scope(str(preprocess_state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE)))
    scope_regions = dict(preprocess_state.get("scope_regions") or {})
    scoped_region = str(scope_regions.get(analysis_scope) or preprocess_state.get("region") or DEFAULT_REGION)
    form["analysis_scope"] = analysis_scope
    if scoped_region:
        form["region"] = scoped_region

    default_outputs = {
        f"results/{DEFAULT_REPORT_NAME}",
        _default_report_path_for_scope(DEFAULT_GENE_NAME, DEFAULT_ANALYSIS_SCOPE),
    }
    if form["out"] in default_outputs and preprocess_state.get("gene_name"):
        form["out"] = _default_report_path_for_scope(str(preprocess_state["gene_name"]), analysis_scope)
    elif not form["out"] and preprocess_state.get("gene_name"):
        form["out"] = _default_report_path_for_scope(str(preprocess_state["gene_name"]), analysis_scope)


def _build_field_info(
    form: dict[str, str],
    *,
    preprocess_state: dict[str, Any],
    vcf_files: list[str],
    idat_prefixes: list[str],
    popstats_files: list[str],
) -> dict[str, dict[str, str]]:
    """Build UI metadata for field examples and foldable biological explanations."""
    return {
        "vcf": {
            "example": form["vcf"] or (vcf_files[0] if vcf_files else "data/gene.vcf.gz"),
            "details": (
                "A VCF, or Variant Call Format file, lists sequence differences observed "
                "between your sample and the reference genome. In this app it provides the "
                "gene-region SNPs and other variant calls that will be summarized in the report."
            ),
        },
        "idat": {
            "example": form["idat"] or (idat_prefixes[0] if idat_prefixes else "data/202277800037_R01C01"),
            "details": (
                "An IDAT base path points to the paired red and green Illumina methylation "
                "array intensity files. Those raw signal files are what methylprep uses to "
                "calculate probe-level methylation beta values for the preprocessed gene subset."
            ),
        },
        "out": {
            "example": form["out"] or f"results/{DEFAULT_REPORT_NAME}",
            "details": (
                "This controls where the generated report artifact is written so you can review "
                "or share the analysis output after the run."
            ),
        },
        "region": {
            "example": form["region"] or str(preprocess_state.get("region", DEFAULT_REGION)),
            "details": (
                "The genomic region limits the analysis to the selected report focus. "
                "The default focus is promoter plus gene body; switch to promoter-only or gene-only "
                "only when you want a narrower report artifact."
            ),
        },
        "popstats": {
            "example": form["popstats"]
            or form["suggested_popstats"]
            or (popstats_files[0] if popstats_files else "data/reference_population.json"),
            "details": (
                "Population statistics files add cohort-level context such as allele frequency "
                "or prevalence summaries. They help interpret whether a gene-region variant looks rare, "
                "common, or enriched in reference populations."
            ),
        },
        "manifest_file": {
            "example": form["manifest_file"]
            or str(preprocess_state.get("manifest_source", "data/infinium-methylationepic-manifest-file.csv")),
            "details": (
                "Optional override. Leave this empty unless you specifically want to force methylprep to use "
                "a custom vendor manifest file. The preprocessing workflow already saves a gene-specific subset "
                "like src/gene_data/GENE_epigenetics_hg19.csv, and that saved subset is reused for the probe "
                "annotation join after methylprep finishes."
            ),
        },
    }


def _build_preprocess_field_info(
    preprocess_state: dict[str, Any],
    *,
    manifest_files: list[str],
) -> dict[str, dict[str, str]]:
    """Build field hints for the preprocessing tab."""
    return {
        "gene_name": {
            "example": str(preprocess_state.get("gene_name") or DEFAULT_GENE_NAME),
            "details": (
                "Enter an HGNC-style gene symbol. The preprocessing workflow will look up the genomic "
                "coordinates for that gene and use the resulting interval in the downstream analysis."
            ),
        },
        "region": {
            "example": str(preprocess_state.get("region") or DEFAULT_REGION),
            "details": (
                "This field is filled by the gene-region lookup step. The standard preprocessing span is "
                "promoter plus gene body, because promoter regulatory variants and methylation probes often sit "
                "just upstream of the transcribed interval."
            ),
        },
        "manifest_source": {
            "example": str(preprocess_state.get("manifest_source") or (manifest_files[0] if manifest_files else "data/infinium-methylationepic-manifest.csv")),
            "details": (
                "Use the full EPIC manifest here. The app will filter it down to the standard promoter+gene interval "
                "and save a much smaller CSV subset into src/gene_data for downstream methylation processing."
            ),
        },
    }


def _build_preprocess_result(preprocess_state: dict[str, Any]) -> dict[str, Any] | None:
    """Build the preprocessing result payload used by the lower status panel."""
    if not any(
        [
            preprocess_state.get("region_ready"),
            preprocess_state.get("manifest_ready"),
            preprocess_state.get("filtered_manifest"),
            preprocess_state.get("dynamic_kb_ready"),
            preprocess_state.get("dynamic_kb_path"),
            preprocess_state.get("logs"),
        ]
    ):
        return None

    filtered_manifest = str(preprocess_state.get("filtered_manifest", "")).strip()
    preview_html = None
    if filtered_manifest:
        preview_path = _resolve_user_path(filtered_manifest)
        if preview_path.exists():
            try:
                preview_df = pd.read_csv(preview_path).head(12)
                preview_html = _render_table(preview_df, rows=12)
            except Exception:
                preview_html = None

    return {
        "gene_name": str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)),
        "region": str(preprocess_state.get("region", "")),
        "analysis_scope": normalize_analysis_scope(str(preprocess_state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE))),
        "analysis_scope_label": get_analysis_scope_label(str(preprocess_state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE))),
        "scope_regions": dict(preprocess_state.get("scope_regions") or {}),
        "scope_region_source": str(preprocess_state.get("scope_region_source", "")),
        "manifest_source": str(preprocess_state.get("manifest_source", "")),
        "filtered_manifest": filtered_manifest,
        "region_ready": bool(preprocess_state.get("region_ready", False)),
        "manifest_ready": bool(preprocess_state.get("manifest_ready", False)),
        "dynamic_kb_ready": bool(preprocess_state.get("dynamic_kb_ready", False)),
        "dynamic_kb_path": str(preprocess_state.get("dynamic_kb_path", "")),
        "dynamic_kb_status": str(preprocess_state.get("dynamic_kb_status", "")),
        "dynamic_kb_provider_count": int(preprocess_state.get("dynamic_kb_provider_count", 0) or 0),
        "dynamic_kb_source_count": int(preprocess_state.get("dynamic_kb_source_count", 0) or 0),
        "dynamic_kb_workflow_count": int(preprocess_state.get("dynamic_kb_workflow_count", 0) or 0),
        "dynamic_kb_workflow_summary": list(preprocess_state.get("dynamic_kb_workflow_summary", [])),
        "use_local_article_evidence": bool(preprocess_state.get("use_local_article_evidence", False)),
        "article_pdf_folder": str(preprocess_state.get("article_pdf_folder", "")),
        "article_pdf_recursive": bool(preprocess_state.get("article_pdf_recursive", True)),
        "max_article_pdfs": int(preprocess_state.get("max_article_pdfs", 100) or 0),
        "local_article_evidence_status": str(preprocess_state.get("local_article_evidence_status", "")),
        "local_article_evidence_record_count": int(
            preprocess_state.get("local_article_evidence_record_count", 0) or 0
        ),
        "local_article_pdf_count": int(preprocess_state.get("local_article_pdf_count", 0) or 0),
        "analysis_ready": bool(preprocess_state.get("analysis_ready", False)),
        "probe_count": int(preprocess_state.get("probe_count", 0)),
        "selected_sources": list(preprocess_state.get("selected_sources", [])),
        "region_candidates": list(preprocess_state.get("region_candidates", [])),
        "build": str(preprocess_state.get("build", "hg19")),
        "preview_html": preview_html,
        "logs": list(preprocess_state.get("logs", [])),
        "region_recently_updated": bool(preprocess_state.get("region_recently_updated", False)),
        "hg38_extraction_suggested": bool(preprocess_state.get("hg38_extraction_suggested", False)),
        "hg38_extraction_message": str(preprocess_state.get("hg38_extraction_message", "")),
        "hg38_extraction_region": str(preprocess_state.get("hg38_extraction_region", "")),
        "progress_percent": (
            100
            if preprocess_state.get("dynamic_kb_ready")
            else 66
            if preprocess_state.get("manifest_ready")
            else 33
            if preprocess_state.get("region_ready")
            else 0
        ),
        "steps": [
            {
                "title": "Find Region from Gene Name",
                "status": "complete" if preprocess_state.get("region_ready") else "pending",
                "summary": (
                    f"Resolved standard promoter+gene scope to {preprocess_state.get('region', DEFAULT_REGION)}"
                    if preprocess_state.get("region_ready")
                    else "Waiting for a gene-symbol lookup."
                ),
            },
            {
                "title": "Select Methylation Data",
                "status": "complete" if preprocess_state.get("manifest_ready") else "pending",
                "summary": (
                    f"{preprocess_state.get('probe_count', 0)} probes saved to {filtered_manifest}"
                    if preprocess_state.get("manifest_ready")
                    else "Waiting for the filtered EPIC manifest subset."
                ),
            },
            {
                "title": "Build Variant Knowledge Base",
                "status": "complete" if preprocess_state.get("dynamic_kb_ready") else "optional",
                "summary": (
                    f"{preprocess_state.get('dynamic_kb_workflow_count', 0)} workflow(s), {preprocess_state.get('dynamic_kb_provider_count', 0)} provider statuses"
                    f" and {preprocess_state.get('local_article_evidence_record_count', 0)} local article snippet(s) saved to {preprocess_state.get('dynamic_kb_path', '')}"
                    if preprocess_state.get("dynamic_kb_ready")
                    else "Waiting for an observed-variant VCF and selected knowledge sources."
                ),
            },
        ],
    }


def _build_population_context_status(
    *,
    popstats: Any | None,
    population_database: dict[str, Any],
    population_insights: dict[str, Any],
) -> str:
    """Summarize which population context sources were available for the current run."""
    has_sidecar_file = popstats is not None
    has_curated_database = str(population_database.get("version", "")).lower() != "generic"
    if not has_curated_database:
        has_curated_database = bool(
            population_insights.get("variant_population_records")
            or population_insights.get("gene_population_patterns")
        )

    if has_sidecar_file and has_curated_database:
        return "File + DB"
    if has_curated_database:
        return "Database"
    if has_sidecar_file:
        return "File"
    return "None"


@app.get("/api/human-proteins")
def human_proteins_api() -> Any:
    """Return a page of the live human protein catalog for the UI tab."""
    query = request.args.get("q", "").strip()
    cursor = request.args.get("cursor", "").strip() or None
    longevity_page_raw = request.args.get("longevity_page", "1").strip()
    reviewed_only_raw = request.args.get("reviewed_only", "1").strip().lower()
    longevity_only_raw = request.args.get("longevity_only", "0").strip().lower()
    exclude_processed_raw = request.args.get("exclude_processed", "0").strip().lower()
    reviewed_only = reviewed_only_raw not in {"0", "false", "no"}
    longevity_only = longevity_only_raw in {"1", "true", "yes"}
    exclude_processed = exclude_processed_raw in {"1", "true", "yes"}
    try:
        longevity_page = max(1, int(longevity_page_raw))
    except ValueError:
        longevity_page = 1

    payload = get_human_protein_catalog(
        query=query,
        reviewed_only=reviewed_only,
        cursor=cursor,
        longevity_only=longevity_only,
        longevity_page=longevity_page,
    )

    processed_gene_symbols = discover_processed_gene_symbols() if exclude_processed else []
    excluded_processed_count = 0
    if exclude_processed:
        processed_gene_symbol_set = set(processed_gene_symbols)
        proteins = payload.get("proteins", [])
        if isinstance(proteins, list):
            filtered_proteins = [
                protein
                for protein in proteins
                if not isinstance(protein, dict)
                or not _protein_matches_processed_gene_symbols(protein, processed_gene_symbol_set)
            ]
            excluded_processed_count = len(proteins) - len(filtered_proteins)
            payload["proteins"] = filtered_proteins
            payload["records_returned"] = len(filtered_proteins)

    proteins = payload.get("proteins", [])
    if isinstance(proteins, list):
        for protein in proteins:
            if not isinstance(protein, dict):
                continue
            visuals = _biorender_visuals_for_protein(protein)
            if visuals:
                protein["biorender_visuals"] = visuals

    payload["exclude_processed"] = exclude_processed
    payload["processed_gene_symbols"] = processed_gene_symbols
    payload["processed_gene_count"] = len(processed_gene_symbols)
    payload["excluded_processed_count"] = excluded_processed_count
    return jsonify(payload)


@app.get("/results/<path:artifact_path>")
def result_artifact(artifact_path: str) -> Any:
    """Serve a generated artifact from the results directory."""
    return send_from_directory(str(RESULTS_DIR), artifact_path, as_attachment=False)


@app.get("/functional-map")
def functional_map_page() -> str:
    """Render the standalone functional-family map for curated genes."""
    report_history = discover_report_history()
    return render_template(
        "functional_map.html",
        functional_map=_build_functional_map(report_history),
    )


@app.route("/", methods=["GET", "POST"])
def index() -> str:
    """Render the landing page and handle preprocessing and analysis submissions."""
    result = None
    analysis_error = None
    preprocess_error = None
    preprocess_notice = None
    knowledge_sources_error = None
    knowledge_sources_notice = None
    extraction_error = None
    extraction_notice = None
    initial_tab = "overview"

    vcf_files = discover_vcf_files()
    bam_files = discover_bam_files()
    idat_prefixes = discover_idat_prefixes()
    popstats_files = discover_population_stats_files()
    manifest_files = discover_manifest_files()
    extraction_tool_status = get_extraction_tool_status()
    extraction_reference_status = get_hg38_reference_status()

    if request.method == "GET":
        session.pop(SESSION_PREPROCESS_KEY, None)
        session.pop(SESSION_EXTRACTION_KEY, None)
        preprocess_state = _empty_preprocess_state(manifest_files)
        extraction_state = _empty_extraction_state(bam_files)
        knowledge_sources_state = _load_knowledge_sources_state()
    else:
        preprocess_state = _load_preprocess_state(manifest_files)
        extraction_state = _load_extraction_state(bam_files)
        knowledge_sources_state = _load_knowledge_sources_state()
    extraction_state["reference_ready"] = bool(extraction_reference_status.get("ready", False))
    bam_files = _merge_path_options(bam_files, list(extraction_state.get("bam_search_results", [])))
    if not extraction_state.get("bam_path") and bam_files:
        extraction_state["bam_path"] = bam_files[0]
    analysis_unlocked = bool(preprocess_state.get("analysis_ready", False))

    form = _empty_form_state()
    _apply_preprocessing_defaults(form, preprocess_state)

    if request.method == "POST":
        workflow = request.form.get("workflow", "analysis").strip()

        if workflow == "knowledge_sources":
            initial_tab = "knowledge_sources"
            known_keys = {spec.key for spec in list_source_specs()}
            known_workflows = {workflow_spec.key for workflow_spec in list_workflow_specs()}
            selected_workflows = [
                key for key in request.form.getlist("knowledge_workflow") if key in known_workflows
            ]
            workflow_source_keys = source_keys_for_workflows(select_workflow_specs(selected_workflows))
            selected_sources = [
                key for key in request.form.getlist("knowledge_source") if key in known_keys
            ]
            if not selected_sources:
                selected_sources = workflow_source_keys
            knowledge_sources_state["selected_workflows"] = selected_workflows
            knowledge_sources_state["selected_sources"] = selected_sources
            knowledge_sources_state["external_consents"] = [
                key for key in request.form.getlist("external_consent") if key in selected_sources
            ]
            credential_store = _session_knowledge_credentials()
            for spec in list_source_specs():
                field_name = f"credential_{spec.key}"
                submitted_secret = request.form.get(field_name, "").strip()
                if submitted_secret:
                    credential_store[spec.key] = submitted_secret
            uploaded_imports = _save_uploaded_knowledge_source_imports(knowledge_sources_state)
            knowledge_sources_state["notice"] = (
                f"Saved {len(selected_workflows)} selected workflow(s), {len(selected_sources)} selected knowledge source(s)"
                f" and {uploaded_imports} licensed-source import upload(s) for preprocessing."
            )
            knowledge_sources_notice = knowledge_sources_state["notice"]
            _store_knowledge_sources_state(knowledge_sources_state)

        elif workflow == "functional_map":
            requested_gene_name = request.form.get("functional_gene_name", "").strip().upper()
            try:
                knowledge_base = load_gene_interpretation_database(requested_gene_name)
                if knowledge_base is None:
                    raise AnalysisError(
                        f"No curated interpretation knowledge base was found for {requested_gene_name or 'that gene'}."
                    )

                gene_context = knowledge_base.get("gene_context", {})
                selected_gene_region = _format_interval_from_record(
                    gene_context.get("gene_region"),
                    default_chrom=str(gene_context.get("chromosome", "")),
                )
                if not selected_gene_region:
                    raise AnalysisError(
                        f"The {requested_gene_name} knowledge base does not contain a usable gene interval."
                    )

                preprocess_state = _empty_preprocess_state(manifest_files)
                extraction_state = _empty_extraction_state(bam_files)
                scope_regions = _build_analysis_scope_regions(
                    requested_gene_name,
                    selected_gene_region,
                )
                selected_region = str(
                    scope_regions.get(DEFAULT_ANALYSIS_SCOPE) or selected_gene_region
                )
                preprocess_state.update(
                    {
                        "gene_name": requested_gene_name,
                        "region": selected_region,
                        "analysis_scope": DEFAULT_ANALYSIS_SCOPE,
                        "scope_regions": scope_regions,
                        "scope_region_source": "Functional Map curated knowledge base",
                        "region_candidates": [
                            {
                                "source": "Functional Map curated knowledge base",
                                "region": selected_gene_region,
                            }
                        ],
                        "selected_sources": ["Functional Map curated knowledge base"],
                        "region_ready": True,
                        "region_recently_updated": True,
                        "build": (
                            "hg38"
                            if _knowledge_base_uses_genome_build(requested_gene_name, "hg38")
                            else "hg19"
                        ),
                    }
                )
                _append_preprocess_log(
                    preprocess_state,
                    (
                        f"Selected {requested_gene_name} from Functional Map and loaded "
                        f"curated promoter+gene region {selected_region}."
                    ),
                )
                preprocess_notice = (
                    f"Selected {requested_gene_name} from Functional Map. "
                    f"Curated promoter+gene region {selected_region} is ready."
                )
                if _gene_needs_hg38_extraction(requested_gene_name):
                    hg38_message = (
                        f"{requested_gene_name} uses a GRCh38/hg38 knowledge-base interval. "
                        "Use the Extraction tab with a GRCh38-aligned BAM to create the regional VCF."
                    )
                    preprocess_state["hg38_extraction_suggested"] = True
                    preprocess_state["hg38_extraction_message"] = hg38_message
                    preprocess_state["hg38_extraction_region"] = selected_region
                    _prefill_extraction_state_from_hg38_preprocessing(
                        extraction_state,
                        gene_name=requested_gene_name,
                        scope_regions=scope_regions,
                        selected_sources=list(preprocess_state["selected_sources"]),
                        region_candidates=list(preprocess_state["region_candidates"]),
                    )
                    _append_preprocess_log(preprocess_state, hg38_message)
                    preprocess_notice = f"{preprocess_notice} {hg38_message}"

                _store_preprocess_state(preprocess_state)
                _store_extraction_state(extraction_state)
                initial_tab = "preprocessing"
                analysis_unlocked = False
                form = _empty_form_state()
                _apply_preprocessing_defaults(form, preprocess_state)
            except (AnalysisError, ValueError) as exc:
                preprocess_error = str(exc)
                initial_tab = "preprocessing"

        elif workflow == "preprocess":
            initial_tab = "preprocessing"
            preprocess_action = request.form.get("preprocess_action", "").strip()
            if preprocess_action == "reset_preprocessing":
                preprocess_state = _empty_preprocess_state(manifest_files)
                preprocess_notice = "Preprocessing was refreshed. Start again with Genetics."
                _append_preprocess_log(preprocess_state, preprocess_notice)
            else:
                previous_gene_name = str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)).strip().upper()
                requested_gene_name = (
                    request.form.get("gene_name", "").strip() or DEFAULT_GENE_NAME
                ).upper()
                preprocess_state.update(
                    {
                        "gene_name": requested_gene_name,
                        "region": request.form.get("preprocess_region", "").strip()
                        or str(preprocess_state.get("region", DEFAULT_REGION)),
                        "manifest_source": request.form.get("manifest_source", "").strip()
                        or str(preprocess_state.get("manifest_source", "")),
                        "overwrite_filtered_manifest": request.form.get("overwrite_filtered_manifest") == "1",
                        "knowledge_vcf_source": request.form.get("knowledge_vcf_source", "").strip()
                        or str(preprocess_state.get("knowledge_vcf_source", "")),
                        "use_local_article_evidence": request.form.get("use_local_article_evidence") == "1",
                        "article_pdf_folder": request.form.get("article_pdf_folder", "").strip()
                        or str(preprocess_state.get("article_pdf_folder", "")),
                        "article_pdf_recursive": request.form.get("article_pdf_recursive") == "1",
                        "max_article_pdfs": _safe_nonnegative_int(
                            request.form.get("max_article_pdfs", preprocess_state.get("max_article_pdfs", 100)),
                            100,
                            maximum=1000,
                        ),
                    }
                )
                if requested_gene_name != previous_gene_name:
                    preprocess_state["region_ready"] = False
                    preprocess_state["manifest_ready"] = False
                    preprocess_state["analysis_ready"] = False
                    preprocess_state["filtered_manifest"] = ""
                    preprocess_state["probe_count"] = 0
                    preprocess_state["selected_sources"] = []
                    preprocess_state["region_candidates"] = []
                    preprocess_state["analysis_scope"] = DEFAULT_ANALYSIS_SCOPE
                    preprocess_state["scope_regions"] = {}
                    preprocess_state["scope_region_source"] = ""
                    preprocess_state["logs"] = []
                    preprocess_state["region_recently_updated"] = False
                    preprocess_state["hg38_extraction_suggested"] = False
                    preprocess_state["hg38_extraction_message"] = ""
                    preprocess_state["hg38_extraction_region"] = ""
                    preprocess_state["knowledge_vcf_source"] = ""
                    preprocess_state["local_article_evidence_status"] = ""
                    preprocess_state["local_article_evidence_record_count"] = 0
                    preprocess_state["local_article_pdf_count"] = 0
                    preprocess_state["dynamic_kb_ready"] = False
                    preprocess_state["dynamic_kb_path"] = ""
                    preprocess_state["dynamic_kb_status"] = ""
                    preprocess_state["dynamic_kb_provider_count"] = 0
                    preprocess_state["dynamic_kb_source_count"] = 0
                    preprocess_state["dynamic_kb_workflow_count"] = 0
                    preprocess_state["dynamic_kb_workflow_summary"] = []

            try:
                if preprocess_action == "reset_preprocessing":
                    pass
                elif preprocess_action == "find_region":
                    _append_preprocess_log(
                        preprocess_state,
                        f"Starting gene-region lookup for {preprocess_state['gene_name']}.",
                    )
                    lookup, captured_stdout, captured_stderr = _capture_preprocess_call(
                        find_gene_region,
                        preprocess_state["gene_name"],
                    )
                    _append_preprocess_log(preprocess_state, captured_stdout, stream="stdout")
                    _append_preprocess_log(preprocess_state, captured_stderr, stream="stderr")
                    preprocess_state["gene_name"] = str(lookup["gene_name"])
                    scope_regions = _build_analysis_scope_regions(
                        str(lookup["gene_name"]),
                        str(lookup["selected_region"]),
                    )
                    preprocess_state["analysis_scope"] = DEFAULT_ANALYSIS_SCOPE
                    preprocess_state["scope_regions"] = scope_regions
                    preprocess_state["scope_region_source"] = (
                        "Local curated promoter/gene intervals"
                        if scope_regions.get("promoter_only")
                        else "Generic upstream promoter heuristic"
                    )
                    preprocess_state["region"] = str(
                        scope_regions.get(DEFAULT_ANALYSIS_SCOPE) or lookup["selected_region"]
                    )
                    preprocess_state["selected_sources"] = list(lookup["selected_sources"]) + [
                        preprocess_state["scope_region_source"]
                    ]
                    preprocess_state["region_candidates"] = list(lookup["candidate_regions"])
                    preprocess_state["region_ready"] = True
                    preprocess_state["manifest_ready"] = False
                    preprocess_state["analysis_ready"] = False
                    preprocess_state["filtered_manifest"] = ""
                    preprocess_state["probe_count"] = 0
                    preprocess_state["dynamic_kb_ready"] = False
                    preprocess_state["dynamic_kb_path"] = ""
                    preprocess_state["dynamic_kb_status"] = ""
                    preprocess_state["dynamic_kb_provider_count"] = 0
                    preprocess_state["dynamic_kb_source_count"] = 0
                    preprocess_state["dynamic_kb_workflow_count"] = 0
                    preprocess_state["dynamic_kb_workflow_summary"] = []
                    preprocess_state["local_article_evidence_status"] = ""
                    preprocess_state["local_article_evidence_record_count"] = 0
                    preprocess_state["local_article_pdf_count"] = 0
                    preprocess_state["region_recently_updated"] = True
                    _append_preprocess_log(
                        preprocess_state,
                        (
                            f"Resolved {preprocess_state['gene_name']} to standard promoter+gene region "
                            f"{preprocess_state['region']} using {', '.join(preprocess_state['selected_sources']) or 'the available sources'}."
                        ),
                    )
                    preprocess_notice = (
                        f"Resolved {preprocess_state['gene_name']} to standard promoter+gene region {preprocess_state['region']}."
                    )
                    if _gene_needs_hg38_extraction(str(preprocess_state["gene_name"])):
                        hg38_message = (
                            f"{preprocess_state['gene_name']} is annotated as a GRCh38/hg38 locus in the local "
                            "knowledge base. Use the Extraction tab with a GRCh38-aligned BAM to create the regional VCF."
                        )
                        preprocess_state["hg38_extraction_suggested"] = True
                        preprocess_state["hg38_extraction_message"] = hg38_message
                        preprocess_state["hg38_extraction_region"] = str(preprocess_state["region"])
                        _prefill_extraction_state_from_hg38_preprocessing(
                            extraction_state,
                            gene_name=str(preprocess_state["gene_name"]),
                            scope_regions=scope_regions,
                            selected_sources=list(preprocess_state["selected_sources"]),
                            region_candidates=list(preprocess_state["region_candidates"]),
                        )
                        _append_preprocess_log(preprocess_state, hg38_message)
                        _append_extraction_log(
                            extraction_state,
                            (
                                f"Preprocessing detected an hg38-only gene and prefilled Extraction with "
                                f"{extraction_state.get('region', '')}."
                            ),
                        )
                        preprocess_notice = f"{preprocess_notice} {hg38_message}"
                    else:
                        preprocess_state["hg38_extraction_suggested"] = False
                        preprocess_state["hg38_extraction_message"] = ""
                        preprocess_state["hg38_extraction_region"] = ""
                elif preprocess_action == "select_methylation":
                    if not preprocess_state.get("region"):
                        raise AnalysisError(
                            "Resolve the gene coordinates first, or enter a region before selecting methylation data."
                        )
                    if not preprocess_state.get("manifest_source"):
                        raise AnalysisError(
                            "Enter an EPIC manifest path before selecting methylation data."
                        )

                    manifest_source = str(_resolve_user_path(str(preprocess_state["manifest_source"])))
                    target_output_path = _expected_filtered_manifest_path(
                        str(preprocess_state["gene_name"]),
                        str(preprocess_state.get("build", "hg19")),
                    )
                    target_output_previously_existed = target_output_path.exists()
                    overwrite_filtered_manifest = bool(
                        preprocess_state.get("overwrite_filtered_manifest", False)
                    )
                    _append_preprocess_log(
                        preprocess_state,
                        (
                            f"Preparing filtered EPIC manifest for {preprocess_state['region']} from {manifest_source}. "
                            f"Overwrite existing subset: {'yes' if overwrite_filtered_manifest else 'no'}."
                        ),
                    )
                    existing_subset_matches_region = (
                        target_output_path.exists()
                        and _filtered_manifest_covers_region(target_output_path, str(preprocess_state["region"]))
                    )
                    if target_output_path.exists() and not overwrite_filtered_manifest and existing_subset_matches_region:
                        selection = _summarize_existing_filtered_manifest(target_output_path)
                        _append_preprocess_log(
                            preprocess_state,
                            (
                                f"Reusing existing filtered manifest at "
                                f"{_as_relative_display(target_output_path)} because overwrite is disabled."
                            ),
                        )
                    else:
                        selection, captured_stdout, captured_stderr = _capture_preprocess_call(
                            save_filtered_manifest,
                            gene_name=str(preprocess_state["gene_name"]),
                            manifest_path=manifest_source,
                            region=str(preprocess_state["region"]),
                            genome_build=str(preprocess_state.get("build", "hg19")),
                            output_dir=PREPROCESSED_MANIFEST_DIR,
                            allow_empty=_allows_empty_methylation_subset(
                                str(preprocess_state["gene_name"])
                            ),
                        )
                        _append_preprocess_log(preprocess_state, captured_stdout, stream="stdout")
                        _append_preprocess_log(preprocess_state, captured_stderr, stream="stderr")
                        _append_preprocess_log(
                            preprocess_state,
                            (
                                "Created a fresh filtered manifest subset."
                                if not target_output_previously_existed
                                else "Refreshed the filtered manifest subset for the current promoter+gene region."
                            ),
                        )
                        _write_filtered_manifest_metadata(
                            Path(selection["output_path"]),
                            gene_name=str(preprocess_state["gene_name"]),
                            region=str(preprocess_state["region"]),
                            analysis_scope=str(preprocess_state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE)),
                        )
                    preprocess_state["filtered_manifest"] = _as_relative_display(selection["output_path"])
                    preprocess_state["manifest_ready"] = True
                    preprocess_state["analysis_ready"] = True
                    preprocess_state["probe_count"] = int(selection["probe_count"])
                    preprocess_state["dynamic_kb_ready"] = False
                    preprocess_state["dynamic_kb_path"] = ""
                    preprocess_state["dynamic_kb_status"] = ""
                    preprocess_state["dynamic_kb_provider_count"] = 0
                    preprocess_state["dynamic_kb_source_count"] = 0
                    preprocess_state["dynamic_kb_workflow_count"] = 0
                    preprocess_state["dynamic_kb_workflow_summary"] = []
                    preprocess_state["local_article_evidence_status"] = ""
                    preprocess_state["local_article_evidence_record_count"] = 0
                    preprocess_state["local_article_pdf_count"] = 0
                    preprocess_state["region_recently_updated"] = False
                    _append_preprocess_log(
                        preprocess_state,
                        (
                            f"Prepared {selection['probe_count']} probes at "
                            f"{preprocess_state['filtered_manifest']}."
                        ),
                    )
                    preprocess_notice = (
                        f"Prepared {selection['probe_count']} probe(s) at {preprocess_state['filtered_manifest']}."
                    )
                elif preprocess_action == "build_variant_knowledge_base":
                    if not preprocess_state.get("region_ready"):
                        raise AnalysisError("Resolve the gene coordinates before building a dynamic knowledge base.")
                    if not preprocess_state.get("manifest_ready") or not preprocess_state.get("filtered_manifest"):
                        raise AnalysisError("Select methylation data before building a dynamic knowledge base.")
                    vcf_source = (
                        str(preprocess_state.get("knowledge_vcf_source", "")).strip()
                        or _default_knowledge_vcf_source(
                            preprocess_state=preprocess_state,
                            extraction_state=extraction_state,
                            vcf_files=vcf_files,
                        )
                    )
                    if not vcf_source:
                        raise AnalysisError(
                            "Choose a VCF source for dynamic knowledge-base generation, or extract one first."
                        )
                    vcf_path = _resolve_user_path(vcf_source)
                    if not vcf_path.exists():
                        raise AnalysisError(f"Dynamic KB VCF source was not found: {vcf_path}")
                    manifest_subset_path = _resolve_user_path(str(preprocess_state["filtered_manifest"]))
                    if not manifest_subset_path.exists():
                        raise AnalysisError(f"Filtered manifest was not found: {manifest_subset_path}")

                    article_pdf_folder = str(preprocess_state.get("article_pdf_folder", "")).strip()
                    article_pdf_folder_path = _resolve_user_path(article_pdf_folder) if article_pdf_folder else None
                    use_local_article_evidence = bool(
                        preprocess_state.get("use_local_article_evidence", False)
                        or article_pdf_folder
                    )
                    selected_sources = list(knowledge_sources_state.get("selected_sources", []))
                    if not selected_sources:
                        raise AnalysisError("Select at least one vetted evidence source before refreshing evidence.")
                    imported_sources = set((knowledge_sources_state.get("source_imports") or {}).keys())
                    approved_sources = set(knowledge_sources_state.get("external_consents", []))
                    unapproved_sources = [
                        key for key in selected_sources if key not in approved_sources and key not in imported_sources
                    ]
                    if unapproved_sources:
                        raise AnalysisError(
                            "Explicit per-source transfer consent is required for: "
                            + ", ".join(sorted(unapproved_sources))
                        )
                    selected_workflows = list(knowledge_sources_state.get("selected_workflows", []))
                    if not selected_workflows:
                        selected_workflows = default_workflow_keys()
                    output_dir = (
                        RESULTS_DIR
                        / "dynamic_knowledge_bases"
                        / f"{sanitize_gene_name_for_filename(str(preprocess_state['gene_name']).lower())}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                    )
                    _append_preprocess_log(
                        preprocess_state,
                        (
                            f"Building dynamic variant knowledge base for {preprocess_state['gene_name']} "
                            f"from {vcf_path} with {len(selected_workflows)} workflow(s) and "
                            f"{len(selected_sources)} selected source(s). "
                            f"Local PDF article evidence: {'enabled' if use_local_article_evidence else 'disabled'}."
                        ),
                    )
                    variants = load_variants(str(vcf_path), str(preprocess_state["region"]))
                    manifest_subset = pd.read_csv(manifest_subset_path)
                    dynamic_payload = build_dynamic_knowledge_base(
                        gene=str(preprocess_state["gene_name"]),
                        region=str(preprocess_state["region"]),
                        genome_build=str(preprocess_state.get("build", "hg19")),
                        variants=variants,
                        manifest_subset=manifest_subset,
                        selected_workflows=selected_workflows,
                        selected_sources=selected_sources,
                        credentials=_session_knowledge_credentials(),
                        source_imports=dict(knowledge_sources_state.get("source_imports") or {}),
                        use_local_article_evidence=use_local_article_evidence,
                        article_pdf_folder=article_pdf_folder_path,
                        article_pdf_recursive=bool(preprocess_state.get("article_pdf_recursive", True)),
                        max_article_pdfs=int(preprocess_state.get("max_article_pdfs", 100) or 0),
                        output_dir=output_dir,
                        cache_dir=PROJECT_ROOT / ".research-cache" / "variant_knowledge",
                    )
                    artifact_path = Path(str(dynamic_payload.get("artifact_path", output_dir / "variant_kb.json")))
                    selected_source_keys = list(
                        (dynamic_payload.get("provenance") or {}).get("selected_source_keys") or selected_sources
                    )
                    local_article_evidence = (
                        dynamic_payload.get("local_article_evidence")
                        if isinstance(dynamic_payload.get("local_article_evidence"), dict)
                        else {}
                    )
                    local_article_provenance = (
                        local_article_evidence.get("provenance")
                        if isinstance(local_article_evidence.get("provenance"), dict)
                        else {}
                    )
                    preprocess_state["knowledge_vcf_source"] = _as_relative_display(vcf_path)
                    preprocess_state["dynamic_kb_ready"] = True
                    preprocess_state["dynamic_kb_path"] = _as_relative_display(artifact_path)
                    preprocess_state["dynamic_kb_provider_count"] = len(dynamic_payload.get("provider_statuses", []))
                    preprocess_state["dynamic_kb_source_count"] = len(selected_source_keys)
                    preprocess_state["dynamic_kb_workflow_count"] = len(dynamic_payload.get("workflow_runs", []))
                    preprocess_state["local_article_evidence_status"] = str(local_article_evidence.get("message", ""))
                    preprocess_state["local_article_evidence_record_count"] = int(
                        local_article_provenance.get(
                            "record_count",
                            len(local_article_evidence.get("records", []) or []),
                        )
                        or 0
                    )
                    preprocess_state["local_article_pdf_count"] = int(
                        local_article_provenance.get("pdf_count", 0) or 0
                    )
                    preprocess_state["dynamic_kb_workflow_summary"] = [
                        {
                            "label": str(workflow.get("label", "")),
                            "status": str(workflow.get("status", "")),
                            "summary": str(workflow.get("summary", "")),
                        }
                        for workflow in dynamic_payload.get("workflow_runs", [])
                        if isinstance(workflow, dict)
                    ]
                    preprocess_state["dynamic_kb_status"] = (
                        f"Built dynamic KB with {preprocess_state['dynamic_kb_workflow_count']} workflow(s) "
                        f"and {preprocess_state['dynamic_kb_provider_count']} provider status record(s). "
                        f"Local article snippets: {preprocess_state['local_article_evidence_record_count']}."
                    )
                    preprocess_state["analysis_ready"] = True
                    _append_preprocess_log(
                        preprocess_state,
                        f"Saved dynamic knowledge base at {preprocess_state['dynamic_kb_path']}.",
                    )
                    if local_article_evidence:
                        _append_preprocess_log(
                            preprocess_state,
                            (
                                f"Local PDF article evidence: "
                                f"{preprocess_state['local_article_evidence_record_count']} snippet(s) "
                                f"from {preprocess_state['local_article_pdf_count']} PDF file(s)."
                            ),
                        )
                    preprocess_notice = preprocess_state["dynamic_kb_status"]
                else:
                    raise AnalysisError("Choose a preprocessing action before submitting the form.")
            except (AnalysisError, ValueError) as exc:
                preprocess_error = str(exc)
                _append_preprocess_log(preprocess_state, str(exc), stream="stderr")
            except Exception as exc:
                preprocess_error = str(exc)
                _append_preprocess_log(preprocess_state, traceback.format_exc(), stream="stderr")

            _store_preprocess_state(preprocess_state)
            _store_extraction_state(extraction_state)
            analysis_unlocked = bool(preprocess_state.get("analysis_ready", False))
            form = _empty_form_state()
            _apply_preprocessing_defaults(form, preprocess_state)

        elif workflow == "extraction":
            initial_tab = "extraction"
            previous_gene_name = str(extraction_state.get("gene_name", DEFAULT_GENE_NAME)).strip().upper()
            requested_gene_name = (
                request.form.get("extraction_gene_name", "").strip() or DEFAULT_GENE_NAME
            ).upper()
            requested_scope = normalize_analysis_scope(
                request.form.get(
                    "extraction_scope",
                    str(extraction_state.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE)),
                )
            )
            requested_build = request.form.get("extraction_genome_build", "hg38").strip() or "hg38"
            extraction_state.update(
                {
                    "gene_name": requested_gene_name,
                    "genome_build": requested_build,
                    "analysis_scope": requested_scope,
                    "region": request.form.get("extraction_region", "").strip()
                    or str(extraction_state.get("region", "")),
                    "bam_path": request.form.get("bam_path", "").strip()
                    or str(extraction_state.get("bam_path", "")),
                    "bam_search_root": request.form.get("bam_search_root", "").strip()
                    or str(extraction_state.get("bam_search_root", _as_relative_display(DATA_DIR))),
                    "output_vcf": request.form.get("output_vcf", "").strip()
                    or _default_extraction_output_display(requested_gene_name, requested_build, requested_scope),
                    "notice": "",
                }
            )
            if requested_gene_name != previous_gene_name:
                extraction_state["region_ready"] = False
                extraction_state["extraction_ready"] = False
                extraction_state["region_candidates"] = []
                extraction_state["selected_sources"] = []
                extraction_state["scope_regions"] = {}
                extraction_state["scope_region_source"] = ""
                extraction_state["last_output_vcf"] = ""
                extraction_state["last_resolved_region"] = ""
                extraction_state["last_commands"] = {}
                extraction_state["logs"] = []

            extraction_action = request.form.get("extraction_action", "").strip()

            try:
                if extraction_action == "search_bam_files":
                    search_root = str(extraction_state.get("bam_search_root") or DATA_DIR)
                    _append_extraction_log(
                        extraction_state,
                        f"Searching for BAM files under {search_root}.",
                    )
                    search_results = search_bam_files(search_root)
                    extraction_state["bam_search_results"] = search_results
                    bam_files = _merge_path_options(discover_bam_files(), search_results)
                    if search_results and (
                        not extraction_state.get("bam_path")
                        or str(extraction_state.get("bam_path")) not in bam_files
                    ):
                        extraction_state["bam_path"] = search_results[0]
                    message = (
                        f"Found {len(search_results)} BAM file(s) under {search_root}."
                        if search_results
                        else f"No BAM files were found under {search_root}."
                    )
                    _append_extraction_log(extraction_state, message)
                    extraction_notice = message
                    extraction_state["notice"] = extraction_notice
                elif extraction_action == "browse_bam_file":
                    if extraction_tool_status.get("docker_runtime"):
                        raise AnalysisError(
                            "A native Windows BAM picker cannot open from inside Docker. "
                            "Use folder search inside the mounted data folder, or run the local UI with local extraction enabled."
                        )
                    selected_bam = browse_bam_file(
                        str(extraction_state.get("bam_path") or extraction_state.get("bam_search_root") or DATA_DIR)
                    )
                    extraction_state["bam_path"] = selected_bam
                    selected_parent = Path(selected_bam).parent
                    extraction_state["bam_search_root"] = selected_parent.as_posix()
                    extraction_state["bam_search_results"] = _merge_path_options(
                        [selected_bam],
                        list(extraction_state.get("bam_search_results", [])),
                    )
                    bam_files = _merge_path_options(discover_bam_files(), list(extraction_state["bam_search_results"]))
                    message = f"Selected BAM file {selected_bam}."
                    _append_extraction_log(extraction_state, message)
                    extraction_notice = message
                    extraction_state["notice"] = extraction_notice
                elif extraction_action == "find_region":
                    _append_extraction_log(
                        extraction_state,
                        f"Starting GRCh38/hg38 gene-region lookup for {extraction_state['gene_name']}.",
                    )
                    lookup, captured_stdout, captured_stderr = _capture_preprocess_call(
                        find_gene_region,
                        extraction_state["gene_name"],
                        genome_build="hg38",
                    )
                    _append_extraction_log(extraction_state, captured_stdout, stream="stdout")
                    _append_extraction_log(extraction_state, captured_stderr, stream="stderr")
                    extraction_state["gene_name"] = str(lookup["gene_name"])
                    scope_regions = _build_extraction_scope_regions(
                        str(lookup["gene_name"]),
                        str(lookup["selected_region"]),
                        "hg38",
                    )
                    extraction_state["scope_regions"] = scope_regions
                    extraction_state["scope_region_source"] = (
                        "Local GRCh38 curated promoter/gene intervals"
                        if _knowledge_base_uses_genome_build(str(lookup["gene_name"]), "hg38")
                        and scope_regions.get("promoter_only")
                        else "GRCh38 public lookup with generic upstream promoter heuristic"
                    )
                    extraction_state["analysis_scope"] = DEFAULT_ANALYSIS_SCOPE
                    extraction_state["region"] = str(
                        scope_regions.get(DEFAULT_ANALYSIS_SCOPE) or lookup["selected_region"]
                    )
                    extraction_state["output_vcf"] = _default_extraction_output_display(
                        str(extraction_state["gene_name"]),
                        "hg38",
                        DEFAULT_ANALYSIS_SCOPE,
                    )
                    extraction_state["selected_sources"] = list(lookup["selected_sources"]) + [
                        extraction_state["scope_region_source"]
                    ]
                    extraction_state["region_candidates"] = list(lookup["candidate_regions"])
                    extraction_state["region_ready"] = True
                    extraction_state["extraction_ready"] = False
                    _append_extraction_log(
                        extraction_state,
                        (
                            f"Resolved {extraction_state['gene_name']} to hg38 extraction region "
                            f"{extraction_state['region']} using {', '.join(extraction_state['selected_sources']) or 'the available sources'}."
                        ),
                    )
                    extraction_notice = (
                        f"Resolved {extraction_state['gene_name']} to hg38 extraction region {extraction_state['region']}."
                    )
                    extraction_state["notice"] = extraction_notice
                elif extraction_action == "prepare_reference":
                    if not extraction_tool_status.get("available"):
                        raise AnalysisError(str(extraction_tool_status.get("message", "")))
                    _append_extraction_log(
                        extraction_state,
                        (
                            "Preparing UCSC hg38 analysis-set reference under "
                            f"{_as_relative_display(HG38_REFERENCE_DIR)}."
                        ),
                    )
                    prepared_status, captured_stdout, captured_stderr = _capture_preprocess_call(
                        prepare_hg38_reference,
                    )
                    _append_extraction_log(extraction_state, captured_stdout, stream="stdout")
                    _append_extraction_log(extraction_state, captured_stderr, stream="stderr")
                    extraction_reference_status = dict(prepared_status)
                    extraction_state["reference_ready"] = bool(prepared_status.get("ready", False))
                    _append_extraction_log(
                        extraction_state,
                        str(prepared_status.get("message", "Reference preparation finished.")),
                    )
                    extraction_notice = str(prepared_status.get("message", "Reference preparation finished."))
                    extraction_state["notice"] = extraction_notice
                elif extraction_action == "extract_vcf":
                    if not extraction_tool_status.get("available"):
                        raise AnalysisError(str(extraction_tool_status.get("message", "")))
                    extraction_reference_status = get_hg38_reference_status()
                    if not extraction_reference_status.get("ready"):
                        raise AnalysisError(
                            "Prepare the hg38 reference before extracting a regional VCF."
                        )
                    if not extraction_state.get("bam_path"):
                        raise AnalysisError("Choose a BAM file from data/ before extracting variants.")
                    if not extraction_state.get("region"):
                        raise AnalysisError("Resolve or enter an hg38 region before extracting variants.")

                    bam_path = _resolve_user_path(str(extraction_state["bam_path"]))
                    output_vcf = _resolve_user_path(str(extraction_state["output_vcf"]))
                    _append_extraction_log(
                        extraction_state,
                        (
                            f"Extracting {extraction_state['region']} from {bam_path} "
                            f"to {output_vcf}."
                        ),
                    )
                    extraction_result, captured_stdout, captured_stderr = _capture_preprocess_call(
                        extract_region_vcf,
                        bam_path=bam_path,
                        region=str(extraction_state["region"]),
                        output_vcf=output_vcf,
                        reference_fasta=HG38_FASTA,
                    )
                    _append_extraction_log(extraction_state, captured_stdout, stream="stdout")
                    _append_extraction_log(extraction_state, captured_stderr, stream="stderr")
                    extracted_vcf = Path(extraction_result["output_vcf"])
                    resolved_region = str(extraction_result["resolved_region"])
                    extracted_vcf_display = _as_relative_display(extracted_vcf)
                    extraction_state["output_vcf"] = extracted_vcf_display
                    extraction_state["last_output_vcf"] = extracted_vcf_display
                    extraction_state["last_resolved_region"] = resolved_region
                    extraction_state["last_commands"] = dict(extraction_result.get("commands") or {})
                    extraction_state["extraction_ready"] = True
                    extraction_state["region_ready"] = True
                    _append_extraction_log(
                        extraction_state,
                        f"Extracted and indexed regional VCF at {extracted_vcf_display}.",
                    )

                    previous_preprocess_gene = str(
                        preprocess_state.get("gene_name", DEFAULT_GENE_NAME)
                    ).strip().upper()
                    if previous_preprocess_gene != requested_gene_name:
                        preprocess_state["manifest_ready"] = False
                        preprocess_state["filtered_manifest"] = ""
                        preprocess_state["probe_count"] = 0
                    scope_regions = dict(extraction_state.get("scope_regions") or {})
                    scope_regions[requested_scope] = resolved_region
                    preprocess_state.update(
                        {
                            "gene_name": requested_gene_name,
                            "region": resolved_region,
                            "analysis_scope": requested_scope,
                            "scope_regions": scope_regions,
                            "scope_region_source": "GRCh38 BAM Extraction regional VCF",
                            "selected_sources": list(extraction_state.get("selected_sources", [])),
                            "region_candidates": list(extraction_state.get("region_candidates", [])),
                            "region_ready": True,
                            "analysis_ready": True,
                            "build": "hg38",
                            "region_recently_updated": True,
                        }
                    )
                    _append_preprocess_log(
                        preprocess_state,
                        (
                            f"Extraction populated the Run Analysis VCF with {extracted_vcf_display} "
                            f"and active hg38 region {resolved_region}."
                        ),
                    )
                    _store_preprocess_state(preprocess_state)
                    analysis_unlocked = True
                    vcf_files = discover_vcf_files()
                    form = _empty_form_state()
                    _apply_preprocessing_defaults(form, preprocess_state)
                    form["vcf"] = extracted_vcf_display
                    form["region"] = resolved_region
                    form["analysis_scope"] = requested_scope
                    extraction_notice = (
                        f"Extracted {requested_gene_name} regional VCF to {extracted_vcf_display}."
                    )
                    extraction_state["notice"] = extraction_notice
                else:
                    raise AnalysisError("Choose an extraction action before submitting the form.")
            except (AnalysisError, ExtractionError, ValueError) as exc:
                extraction_error = str(exc)
                extraction_state["notice"] = ""
                _append_extraction_log(extraction_state, str(exc), stream="stderr")
            except Exception:
                extraction_error = "Extraction failed unexpectedly. Review the Extraction console."
                extraction_state["notice"] = ""
                _append_extraction_log(extraction_state, traceback.format_exc(), stream="stderr")

            extraction_reference_status = get_hg38_reference_status()
            extraction_state["reference_ready"] = bool(extraction_reference_status.get("ready", False))
            bam_files = _merge_path_options(discover_bam_files(), list(extraction_state.get("bam_search_results", [])))
            _store_extraction_state(extraction_state)
            analysis_unlocked = bool(preprocess_state.get("analysis_ready", False))
            if not form.get("vcf"):
                form = _empty_form_state()
                _apply_preprocessing_defaults(form, preprocess_state)

        else:
            if not analysis_unlocked:
                preprocess_error = "Complete preprocessing before running the analysis workflow."
                initial_tab = "preprocessing"
            else:
                initial_tab = "analysis"
                form = _empty_form_state()
                _apply_preprocessing_defaults(form, preprocess_state)
                form.update(
                    {
                        "vcf": request.form.get("vcf", "").strip(),
                        "idat": request.form.get("idat", "").strip(),
                        "out": request.form.get("out", "").strip() or form["out"],
                        "analysis_scope": normalize_analysis_scope(
                            request.form.get("analysis_scope", form.get("analysis_scope", DEFAULT_ANALYSIS_SCOPE))
                        ),
                        "region": request.form.get("region", "").strip() or form["region"],
                        "popstats": request.form.get("popstats", "").strip(),
                        "manifest_file": request.form.get("manifest_file", "").strip() or form["manifest_file"],
                        "overwrite_general_database": (
                            "1" if request.form.get("overwrite_general_database") else ""
                        ),
                    }
                )

                try:
                    analysis_result = run_analysis(
                        vcf_path=str(_resolve_user_path(form["vcf"])),
                        idat_base=str(_resolve_user_path(form["idat"])),
                        output_path=str(_resolve_user_path(form["out"])),
                        gene_name=str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)),
                        region=form["region"],
                        analysis_scope=form["analysis_scope"],
                        popstats_source=str(_resolve_user_path(form["popstats"])) if form["popstats"] else None,
                        manifest_filepath=(
                            str(_resolve_user_path(form["manifest_file"])) if form["manifest_file"] else None
                        ),
                        overwrite_general_database=bool(form["overwrite_general_database"]),
                        general_database_path=_general_analysis_database_path(),
                        dynamic_knowledge_base_path=(
                            _resolve_user_path(str(preprocess_state.get("dynamic_kb_path", "")))
                            if preprocess_state.get("dynamic_kb_ready") and preprocess_state.get("dynamic_kb_path")
                            else None
                        ),
                        genome_build=str(preprocess_state.get("build", "")) or None,
                        interpretation_mode="dual",
                    )

                    methylation_probe_preview = analysis_result.methylation_insights.get("probe_preview")
                    variant_preview = _prepare_variant_preview_table(analysis_result.variants)
                    methylation_preview = _prepare_methylation_preview_table(analysis_result.methylation)
                    variant_rows = _serialize_table_rows(variant_preview)
                    dynamic_knowledge_base_path_raw = getattr(
                        analysis_result,
                        "dynamic_knowledge_base_path",
                        None,
                    )
                    dynamic_knowledge_base_path = (
                        Path(dynamic_knowledge_base_path_raw)
                        if dynamic_knowledge_base_path_raw
                        else None
                    )
                    dynamic_payload = load_dynamic_knowledge_base(dynamic_knowledge_base_path)
                    data_sources = _build_data_sources_payload(
                        knowledge_base=analysis_result.knowledge_base,
                        population_database=analysis_result.population_database,
                        population_insights=analysis_result.population_insights,
                        methylation_insights=analysis_result.methylation_insights,
                        dynamic_payload=dynamic_payload,
                        selected_source_keys=list(knowledge_sources_state.get("selected_sources", [])),
                    )
                    result = build_canonical_report(
                        {
                            "gene": str(preprocess_state.get("gene_name", DEFAULT_GENE_NAME)),
                            "genome_build": str(preprocess_state.get("build", "")),
                            "region": form["region"],
                            "analysis_scope": form["analysis_scope"],
                            "variants": analysis_result.variants,
                            "methylation": analysis_result.methylation,
                            "population_statistics": analysis_result.popstats,
                            "knowledge_base": analysis_result.knowledge_base,
                            "dynamic_knowledge_base": dynamic_payload or {},
                            "interpretation": getattr(analysis_result, "interpretation", {}),
                            "source_provenance": {
                                "vcf": form["vcf"],
                                "idat": form["idat"],
                                "manifest": form["manifest_file"],
                            },
                            "artifacts": {
                                "report": _as_relative_display(analysis_result.report_path),
                                "methylation": _as_relative_display(analysis_result.methylation_output_path),
                            },
                        }
                    )
                except AnalysisError as exc:
                    analysis_error = str(exc)

    preprocess_result = _build_preprocess_result(preprocess_state)
    report_history = discover_report_history()
    general_database = {"exists": False, "rows": [], "columns": [], "row_count": 0, "path": "retired"}
    processed_gene_symbols = discover_processed_gene_symbols(report_history, None)
    analysis_scope_options = _build_analysis_scope_options(preprocess_state)
    extraction_scope_options = _build_extraction_scope_options(extraction_state)
    if not preprocess_state.get("knowledge_vcf_source"):
        preprocess_state["knowledge_vcf_source"] = _default_knowledge_vcf_source(
            preprocess_state=preprocess_state,
            extraction_state=extraction_state,
            vcf_files=vcf_files,
        )
    knowledge_source_cards = _build_knowledge_source_cards(knowledge_sources_state)
    knowledge_source_groups = _group_knowledge_source_cards(knowledge_source_cards)
    knowledge_workflow_cards = _build_knowledge_workflow_cards(knowledge_sources_state)
    task_initial = {
        "analysis": "results" if result else "run",
        "history": "history",
        "knowledge_sources": "settings",
        "proteins": "data_explorer",
        "structure": "data_explorer",
        "overview": "run",
        "preprocessing": "run",
        "extraction": "run",
    }.get(initial_tab, "results" if result else "run")

    return render_template(
        "v2/index.html",
        form=form,
        error=analysis_error,
        preprocess_error=preprocess_error,
        preprocess_notice=preprocess_notice,
        knowledge_sources_error=knowledge_sources_error,
        knowledge_sources_notice=knowledge_sources_notice,
        knowledge_sources_state=knowledge_sources_state,
        knowledge_workflow_cards=knowledge_workflow_cards,
        knowledge_source_groups=knowledge_source_groups,
        knowledge_source_cards=knowledge_source_cards,
        knowledge_source_lane_labels=LANE_LABELS,
        extraction_error=extraction_error,
        extraction_notice=extraction_notice,
        preprocess_state=preprocess_state,
        preprocess_result=preprocess_result,
        extraction_state=extraction_state,
        extraction_scope_options=extraction_scope_options,
        extraction_tool_status=extraction_tool_status,
        extraction_reference_status=extraction_reference_status,
        analysis_unlocked=analysis_unlocked,
        result=result,
        initial_tab=task_initial,
        field_info=_build_field_info(
            form,
            preprocess_state=preprocess_state,
            vcf_files=vcf_files,
            idat_prefixes=idat_prefixes,
            popstats_files=popstats_files,
        ),
        preprocess_field_info=_build_preprocess_field_info(
            preprocess_state,
            manifest_files=manifest_files,
        ),
        data_dir=_as_relative_display(DATA_DIR),
        results_dir=_as_relative_display(RESULTS_DIR),
        vcf_files=vcf_files,
        bam_files=bam_files,
        idat_prefixes=idat_prefixes,
        popstats_files=popstats_files,
        manifest_files=manifest_files,
        report_history=report_history,
        general_database=general_database,
        processed_gene_symbols=processed_gene_symbols,
        analysis_scope_options=analysis_scope_options,
        featured_protein_queries=FEATURED_HUMAN_PROTEIN_QUERIES,
    )


def run_server(host: str = "127.0.0.1", port: int = 8766, debug: bool = False) -> None:
    """Start the web server."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PREPROCESSED_MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    app.run(host=host, port=port, debug=debug)
