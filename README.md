# nophigene-drd4-analysis

Containerized DRD4 analysis with a browser-based workbench and a CLI fallback.

The project now supports two ways to run the same Python pipeline:

- `web` mode: a lightweight browser UI for selecting mounted input files and running the analysis visually
- `cli` mode: the original command-line workflow for scripted or batch-style execution

The Docker image starts in `web` mode by default.

## Implementation approach

The UI was added as a thin layer on top of the existing project structure instead of replacing the pipeline:

- [src/analysis.py](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/src/analysis.py:1) now contains reusable analysis functions plus `run_analysis()` so both the browser and CLI call the same code path
- [src/webapp.py](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/src/webapp.py:1) provides the Flask UI and scans the mounted `data/` directory for likely inputs
- [src/templates/index.html](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/src/templates/index.html:1) holds the browser UI template
- [src/app.py](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/src/app.py:1) is the Docker entrypoint and dispatches to either `web` or `cli`
- [Dockerfile](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Dockerfile:1) now exposes port `8000` and starts the browser app by default

This keeps the analysis logic centralized and makes Docker usage straightforward:

- the container still reads inputs from mounted folders
- the UI only orchestrates file selection and execution
- outputs are still written back to the mounted `results/` folder on the host

## What the web UI does

The browser workbench is designed for the current repo shape and mounted-volume workflow.

It:

- discovers `*.vcf` and `*.vcf.gz` files under `data/`
- discovers IDAT pairs by looking for matching `_Grn.idat` and `_Red.idat` files
- optionally discovers CSV and JSON files that can be used as population statistics sidecars
- runs the same Python analysis pipeline as the CLI
- writes a report artifact plus a companion methylation CSV to `results/`
- shows a preview of the variant and methylation tables after each run

## Prerequisites

You need:

- Docker Engine 20.10 or newer
- a host-side `data/` directory containing your input files
- a host-side `results/` directory for generated outputs

Create those folders from the project root if they do not exist yet:

```powershell
New-Item -ItemType Directory -Force data, results
```

Or on macOS/Linux:

```bash
mkdir -p data results
```

## Expected input files

The current UI and CLI expect:

- a VCF or bgzipped VCF for the DRD4 region, for example `data/drd4.vcf.gz`
- an IDAT sample prefix, which means the actual files should exist as:
  - `data/202277800037_R01C01_Grn.idat`
  - `data/202277800037_R01C01_Red.idat`
- optionally, a population statistics sidecar file such as:
  - `data/gnomad.csv`
  - `data/gnomad.json`

Important:

- in the UI, the IDAT input is the shared prefix only, for example `data/202277800037_R01C01`
- the workflow currently uses the curated DRD4 methylation manifest already checked into the repo

## Build the Docker image

From the project root:

```bash
docker build -t nophigene-drd4-analysis:latest .
```

## Run the web UI with Docker

This is the main way to use the app now.

### Windows PowerShell

```powershell
docker run --rm -it `
  -p 8000:8000 `
  -v "${PWD}\data:/home/appuser/app/data" `
  -v "${PWD}\results:/home/appuser/app/results" `
  nophigene-drd4-analysis:latest
```

### macOS/Linux

```bash
docker run --rm -it \
  -p 8000:8000 \
  -v "${PWD}/data":/home/appuser/app/data \
  -v "${PWD}/results":/home/appuser/app/results \
  nophigene-drd4-analysis:latest
```

Then open:

- [http://localhost:8000](http://localhost:8000)

### What this command does

- `--rm` removes the container after you stop it
- `-it` keeps logs visible in the terminal
- `-p 8000:8000` publishes the UI to your host machine
- the `data/` mount makes your input files visible inside the container
- the `results/` mount makes generated files persist on your host

## Using the web UI

After opening the browser interface:

1. Pick or enter the VCF path.
2. Pick or enter the IDAT base path.
3. Choose the report output path, usually something like `results/drd4_report.html`.
4. Optionally set a different genomic region.
5. Optionally add a population statistics CSV or JSON file.
6. Click `Run DRD4 analysis`.

When the run finishes, the page will show:

- PASS variant count
- methylation probe count
- the generated report path
- the generated methylation CSV path
- a table preview of both outputs

## Output files

The pipeline writes:

- the main report file at the path you choose in the UI or CLI
- a companion methylation CSV next to that report

Example:

- requested report: `results/drd4_report.html`
- generated methylation file: `results/drd4_report_methylation.csv`

If you use an HTML report, you can open it directly in your browser from the host `results/` folder.

## CLI mode inside Docker

The browser UI is the default container behavior, but the original CLI flow is still available.

Use `cli` after the image name to override the default startup mode.

### Windows PowerShell

```powershell
docker run --rm -it `
  -v "${PWD}\data:/home/appuser/app/data" `
  -v "${PWD}\results:/home/appuser/app/results" `
  nophigene-drd4-analysis:latest `
  cli `
  --vcf /home/appuser/app/data/drd4.vcf.gz `
  --idat /home/appuser/app/data/202277800037_R01C01 `
  --out /home/appuser/app/results/drd4_report.html `
  --region 11:63671737-63677367
```

### macOS/Linux

```bash
docker run --rm -it \
  -v "${PWD}/data":/home/appuser/app/data \
  -v "${PWD}/results":/home/appuser/app/results \
  nophigene-drd4-analysis:latest \
  cli \
  --vcf /home/appuser/app/data/drd4.vcf.gz \
  --idat /home/appuser/app/data/202277800037_R01C01 \
  --out /home/appuser/app/results/drd4_report.html \
  --region 11:63671737-63677367
```

You can also pass:

- `--popstats /home/appuser/app/data/gnomad.json`
- `--manifest-file /home/appuser/app/data/custom_manifest.csv.gz`

## Run the web UI on a different port

If port `8000` is already taken, change both the Docker port mapping and the launcher port:

```bash
docker run --rm -it \
  -p 8050:8050 \
  -v "${PWD}/data":/home/appuser/app/data \
  -v "${PWD}/results":/home/appuser/app/results \
  nophigene-drd4-analysis:latest \
  web --host 0.0.0.0 --port 8050
```

Then open:

- [http://localhost:8050](http://localhost:8050)

## Run locally without Docker

Docker is the intended workflow, but for local development you can also run the UI directly once your Python environment is healthy.

```bash
python src/app.py web --host 0.0.0.0 --port 8000
```

Or the CLI:

```bash
python src/app.py cli --vcf data/drd4.vcf.gz --idat data/202277800037_R01C01 --out results/drd4_report.html
```

## Development notes

The UI is intentionally lightweight.

That means:

- the browser does not upload large genomics files into the container
- instead, you mount `data/` and `results/` and point the form at files that already exist there
- this keeps the Docker workflow simple and avoids building a separate file-upload backend before the pipeline itself is more mature

## Troubleshooting

### The browser opens but the suggestion lists are empty

Check that:

- you mounted the host `data/` directory into `/home/appuser/app/data`
- your files are actually inside that folder on the host
- the IDAT sample has both `_Grn.idat` and `_Red.idat`

### The page says a VCF or IDAT file is missing

The UI validates the exact paths you submit. Double-check the path text in the form and make sure the mounted file exists inside the container path space.

### The report file was created but I expected more analysis content

The report generator is intentionally lightweight right now. It produces a working HTML, JSON, or CSV summary with table previews rather than a fully designed biological interpretation report.

### `pytest` does not run locally

The repository currently contains a broken local virtual environment reference on this machine. Recreate the environment before relying on local test commands outside Docker.

## Current limitations

- the pipeline is still DRD4-specific
- the methylation workflow still relies on the checked-in DRD4 region manifest
- population statistics are loaded but not yet deeply merged into the variant table
- the generated report is currently a structured summary, not a full scientific narrative

## Suggested next steps

- add asynchronous job execution so long-running runs do not block a single web request
- add richer report sections with figures and domain-specific interpretation
- add stronger validation around VCF indexing and IDAT naming
- add real integration tests that exercise both `cli` and `web` flows in Docker
