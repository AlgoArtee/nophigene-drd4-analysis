# NophiGene Version 2

NophiGene is a local, single-user, evidence-first bioinformatics workbench. It keeps measured observations, exploratory statistics, scientific literature, established medical evidence, interaction networks, and computational predictions in separate data types and result sections.

It is research-use-only and does not provide diagnoses or treatment recommendations.

## What changed in Version 2

- Canonical report schema `3.0` with Summary plus Objective Data, Statistics, Scientific Literature, Medical Information, Interactions, and Predictions.
- Normalized SQLAlchemy 2 persistence backed by Community SQLCipher in the supported Linux container.
- Versioned Alembic migrations, encrypted pre-migration backups, seven-daily/four-weekly rotation, and a hash-chained audit log.
- `/api/v2` for new runs; `/api/v1` remains read-only for legacy jobs and artifacts.
- Explicit external-provider consent with an exact payload preview and hash confirmation.
- Typed, progressively expandable one-to-three-hop interaction graphs capped at 150 genes.
- Independent model manifests and input gates for AlphaGenome, AlphaMissense, EVE/popEVE, AlphaFold DB, Borzoi, Enformer, Sei, DeepSEA, ESM, MethylBERT, Methyl-GP, MethylProphet, Melody, and disease-specific adapters.
- SeSAMe is the default single-sample EPIC preprocessing contract; minfi/noob remains available for compatible cohort workflows.
- DANDELION 0.1.0 cohort analysis with immutable dataset manifests, both gene→gene and SNP→gene exposure modes, per-exposure BH-FDR, signed offline R jobs, cancellation, and source-backed statistical interaction hypotheses.

The unsupported Predictive Theses system, hard-coded HERC2 eye-colour rules, synthetic cases, active central CSV store, and monolithic legacy result template have been removed. Historical files can be inventoried, imported without synthetic fields, and placed in an encrypted read-only archive.

## Supported runtime

The supported production-like runtime is Windows 11 with WSL2 and Docker Desktop. Flask and SQLCipher run in the pinned Linux container and bind only to `127.0.0.1`.

The Windows launcher retrieves or creates the database key in Windows Credential Manager, creates access-restricted Compose secret files for the container lifetime, generates a one-time browser-session token, and opens the local UI. Secrets are not passed in command-line arguments or environment variables and are removed by the stop launcher.

Start:

```powershell
.\scripts\start-v2.ps1
```

Stop and remove runtime secret files:

```powershell
.\scripts\stop-v2.ps1
```

The UI is served at [http://127.0.0.1:8766](http://127.0.0.1:8766). The launcher can accept `-Port` and `-NoBrowser`.

The start and stop launchers print numbered stages, resolved paths, Compose state, health information, elapsed time, and failure diagnostics without printing secret values. Useful options include `-SkipBuild`, `-StartupTimeoutSeconds`, `-ShutdownTimeoutSeconds`, and `-DryRun`. `Start NophiGene UI.cmd` and `Stop NophiGene UI.cmd` are the only root CMD launchers; retired launcher variants are preserved under `version1/launchers` and are unsupported.

## Development setup

Python 3.12 is used by the app image and supported local test environment. Plain SQLite is allowed only for explicit development/tests; the container sets `NOPHIGENE_REQUIRE_ENCRYPTION=1` and fails closed without SQLCipher.

Only `.venv-v2` is used. The broken Python 3.10 Version 1 environment was removed; its interpreter metadata and recreation notes are archived under `version1/environment`.

```powershell
py -3.12 -m venv .venv-v2
.\.venv-v2\Scripts\python.exe -m pip install --upgrade pip
.\.venv-v2\Scripts\python.exe -m pip install -r requirements-app.txt
.\.venv-v2\Scripts\python.exe -m pytest -q tests\test_workbench_v2.py
```

Run a plaintext local development server only when working with non-sensitive fixtures:

```powershell
.\.venv-v2\Scripts\python.exe src\app.py web --host 127.0.0.1 --port 8766
```

## API

- API index: [http://127.0.0.1:8766/api/v2](http://127.0.0.1:8766/api/v2)
- OpenAPI: [http://127.0.0.1:8766/api/v2/openapi.json](http://127.0.0.1:8766/api/v2/openapi.json)
- Health: [http://127.0.0.1:8766/api/v2/health](http://127.0.0.1:8766/api/v2/health)

`POST /api/v2/runs` accepts one gene in the UI and at most 100 unique genes through the API. A reusable sample profile is required for processing operations.

Example request:

```json
{
  "genes": ["DRD4"],
  "profile_id": "sample-profile-id",
  "genome_build": "hg38",
  "sample_context": {
    "tissue": "blood",
    "platform": "Illumina EPIC",
    "normalization": "SeSAMe",
    "ancestry": "not declared",
    "phenotype_terms": []
  },
  "source_set": [],
  "models": []
}
```

External evidence refresh is a two-step operation:

1. Read `/api/v2/runs/{id}/evidence/payload-preview?sources=pubmed,clinvar`.
2. Submit the returned `payload_sha256`, selected sources, and per-source consent to `/api/v2/runs/{id}/evidence/refresh`.

Imported licensed evidence does not require external-transfer consent. A linkout or unavailable adapter never counts as an assessed source.

## Scientific safeguards

- Primary variant rows require PASS plus GQ ≥20 and DP ≥10 when those fields exist. Missing QC values are retained and reported rather than invented.
- Multi-allelic alleles must be split and normalized; native GRCh37 and GRCh38 identities remain distinct.
- Ambiguous liftover, invalid REF mapping, or incomplete chain/tool provenance blocks the affected adapter.
- Methylation beta values are displayed for interpretation. Compatible tests can use M-values.
- Single-sample comparisons require at least 30 raw values from an exactly compatible tissue/platform/normalization/build reference cohort.
- Public and user cohorts are never silently pooled. Raw p-values, effects, uncertainty, and within-family BH q-values remain separate.
- Medical records require an authoritative clinical source and release/effective date. GWAS, trials, preprints, case reports, and adverse-event signals stay in Literature.
- PGx diplotypes resolve only when every defining locus is QC-covered and exactly one solution remains; otherwise phenotype is not assessed.
- Model outputs remain independent. No NophiGene consensus score is calculated.

## Models and workers

Model manifests are allowlisted and checksummed. A definition implements the `inspect_inputs`, `estimate_resources`, `prepare`, `execute`, `normalize`, and `validate` contract. The browser-facing container cannot install or launch model containers.

Installation confirmation returns a signed-runner requirement; a privilege-separated runner must validate assets, licenses, hashes, disk, VRAM, and container digests. Local model containers must run without network access after asset installation. GPU models remain unavailable until WSL2/NVIDIA preflight succeeds.

Wave status:

- Wave 1: API/precomputed/structure manifests and normalized job/result contracts.
- Wave 2: offline GPU-container contracts for regulatory and ESM models; assets and adapters must be installed explicitly.
- Wave 3: input-gated experimental methylation adapters. MethylBERT requires Bismark BAM/SAM with XM tags; MethylProphet requires matched expression; Methyl-GP is blocked for ordinary human EPIC 5mC; Melody is blocked pending a verified official release.

The R methylation worker is under `docker/methylation`. It has no network at execution time and writes `measurements.csv` plus `qc.json`; failed and flagged probes remain available outside primary results.

### DANDELION cohort statistics

DANDELION is integrated as a statistical method, not as an AI model. Put cohort files below `data/dandelion`, register a JSON manifest in Run, and then queue an offline analysis. Registration validates CSV/TSV headers (RDS is validated by the worker), records structured phenotype/build/namespace/context, and hashes every input without copying it. The optional “managed copy” action creates a separate immutable dataset backed by AES-256 encrypted archives using a key distinct from the SQLCipher and runner-signing keys.

The `dandelion-worker` uses pinned R 4.6.1 and DANDELION 0.1.0. It has no network, no Docker socket, no elevated capabilities, and accepts only HMAC-signed manifests. It executes one CPU job at a time, streams delimited exposure columns in configurable chunks, validates p-values, preserves raw RDS results, and checks that package significance calls equal Benjamini–Hochberg q-values for the same exposure family. See [docs/DANDELION.md](docs/DANDELION.md) for the file contract and API examples.

DANDELION results populate Statistics and typed directed Interactions. Objective Data is not applicable for this cohort workflow; Literature and Medical remain not assessed; Predictions remains not requested. No result is presented as causal or medically established.

## Database and migrations

Initialize a development database:

```powershell
.\.venv-v2\Scripts\python.exe src\app.py db init
```

The container entrypoint creates an encrypted backup before `alembic upgrade head`. It creates daily backups and a weekly backup on Sundays, then rotates to seven daily and four weekly archives.

Useful database commands:

```powershell
.\.venv-v2\Scripts\python.exe src\app.py db inventory
.\.venv-v2\Scripts\python.exe src\app.py db migrate-legacy
.\.venv-v2\Scripts\python.exe src\app.py db migrate-legacy --apply
.\.venv-v2\Scripts\python.exe src\app.py db verify-migration
.\.venv-v2\Scripts\python.exe src\app.py db verify-audit
.\.venv-v2\Scripts\python.exe src\app.py db archive-legacy --apply --password-file .\path\to\password.txt
```

Legacy import is dry-run-first. Synthetic synthesis fields are excluded. The encrypted archive contains a manifest and originals; the originals are not deleted automatically.

## Repository layout

- `src/workbench/`: normalized schema, persistence, reports, statistics, evidence gates, graphs, models, PGx, audit, backups, exports, and migration.
- `src/api/v2_routes.py`: Version 2 API and data-disclosure gates.
- `src/templates/v2/` and `src/static/`: task-based UI and vendored runtime assets.
- `migrations/`: Alembic environment and schema revision.
- `docker/methylation/`: offline SeSAMe/minfi preprocessing worker.
- `docker/dandelion/`: offline signed DANDELION R worker and JSON contract.
- `tests/test_workbench_v2.py`: Version 2 scientific/security acceptance tests.
- `version1/`: retired Version 1 launchers, template, predictive generator, notes, exploratory material, environment metadata, and local loose-store archive.

Large inputs are never silently copied. Input paths and checksums are stored; managed encrypted copies are always an explicit user action.
