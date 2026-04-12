# nophigene-drd4-analysis

Local-first DRD4 analysis workbench with an optional Docker path.

The project now has two supported run modes:

- local mode: the default and recommended workflow, launched from the repaired `.venv`
- Docker mode: a slimmer secondary option for reproducible container runs

## Why the workflow changed

This app is currently a single-process Python web app plus a CLI pipeline:

- [src/webapp.py](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/src/webapp.py:1) provides the Flask UI
- [src/analysis.py](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/src/analysis.py:1) contains the reusable analysis workflow
- [src/app.py](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/src/app.py:1) dispatches to either web or CLI mode

Because the current app does not need multiple services, a database, or orchestration, Docker was adding more local overhead than value. The biggest pain points were:

- large build context
- long image build times
- heavy scientific dependencies getting pulled into every app build
- Docker Desktop startup latency for simple local runs

The repo is now structured so:

- local launch is the default
- Docker is still available, but slimmer
- app runtime dependencies are separated from optional research extras

## Dependency layout

The dependency files are now split by purpose:

- [requirements-app.txt](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/requirements-app.txt:1)
  - minimal runtime set for the UI, CLI, and tests
- [requirements-research.txt](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/requirements-research.txt:1)
  - optional heavier packages for exploratory or future workflows
- [requirements.txt](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/requirements.txt:1)
  - convenience alias to the app requirements

At the moment, the app runtime uses:

- `Flask`
- `pandas`
- `numpy`
- `scikit-allel`
- `methylprep`
- `requests`

Moved out of the default app runtime:

- `deepchem`
- `biomart`
- `matplotlib`
- `pysam`

Important note:

- `pysam` remains listed as an optional research dependency, but it still does not install cleanly on this Windows setup

## Launchers

### Default local launchers

These are the main “starter icon” files for day-to-day use on Windows:

- [Start NophiGene UI.cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Start%20NophiGene%20UI.cmd)
- [Stop NophiGene UI.cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Stop%20NophiGene%20UI.cmd)

They call:

- [scripts/start_nophigene_ui_local.ps1](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/scripts/start_nophigene_ui_local.ps1:1)
- [scripts/stop_nophigene_ui_local.ps1](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/scripts/stop_nophigene_ui_local.ps1:1)

What the local start launcher does:

- checks that `.venv\Scripts\python.exe` exists
- checks that key app dependencies can be imported
- creates `data/` and `results/` if needed
- starts the UI from the local environment
- waits for the server to respond on `http://127.0.0.1:8000`
- opens the browser automatically
- tracks the running process in a local PID file

What the local stop launcher does:

- stops the tracked local UI process
- removes the PID file
- exits cleanly if nothing is running

### Secondary Docker launchers

These are still available if you want a containerized run:

- [Start NophiGene UI (Docker).cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Start%20NophiGene%20UI%20%28Docker%29.cmd)
- [Stop NophiGene UI (Docker).cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Stop%20NophiGene%20UI%20%28Docker%29.cmd)

They call:

- [scripts/start_nophigene_ui.ps1](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/scripts/start_nophigene_ui.ps1:1)
- [scripts/stop_nophigene_ui.ps1](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/scripts/stop_nophigene_ui.ps1:1)

## Local setup

### Recommended Python version

Use Python `3.10`.

The current project `.venv` has already been repaired to point to Python `3.10.11`.

### Install app dependencies

If you need to recreate the local environment from scratch:

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements-app.txt
```

### Install optional research extras

Only do this if you need the non-runtime stack:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-research.txt
```

## Local-first usage

### Fastest path

Double-click:

- [Start NophiGene UI.cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Start%20NophiGene%20UI.cmd)

Then open:

- [http://127.0.0.1:8000](http://127.0.0.1:8000)

When finished, double-click:

- [Stop NophiGene UI.cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Stop%20NophiGene%20UI.cmd)

### Manual local launch

You can also run the UI directly:

```powershell
.\.venv\Scripts\python.exe src\app.py web --host 127.0.0.1 --port 8000
```

Or run the CLI path:

```powershell
.\.venv\Scripts\python.exe src\app.py cli --vcf data/drd4.vcf.gz --idat data/202277800037_R01C01 --out results/drd4_report.html
```

## Expected input layout

The UI and CLI both assume a local project structure like this:

```text
data/
  drd4.vcf.gz
  202277800037_R01C01_Grn.idat
  202277800037_R01C01_Red.idat
results/
```

Important:

- the IDAT argument or form field uses the shared prefix only
- example: `data/202277800037_R01C01`

## What the UI writes

Each run creates:

- a report file at the path you choose
- a companion methylation CSV beside the report

Example:

- requested report: `results/drd4_report.html`
- generated methylation file: `results/drd4_report_methylation.csv`

## Docker is now secondary

Docker still works, but it is no longer the recommended local default.

### What changed to make Docker lighter

- the image now installs from [requirements-app.txt](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/requirements-app.txt:1) instead of the full research stack
- [Dockerfile](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Dockerfile:1) now copies only `src/` and the app requirements into the image
- [.dockerignore](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/.dockerignore:1) now excludes:
  - `data/`
  - `results/`
  - `.venv/`
  - `.docker-local/`
  - `.pytest_cache/`
  - notebooks, tests, and launcher scripts

That should materially reduce Docker build context size and image churn.

### Build the Docker image manually

```bash
docker build -t nophigene-drd4-analysis:latest .
```

### Run Docker manually

```bash
docker run --rm -it \
  -p 8000:8000 \
  -v "${PWD}/data":/home/appuser/app/data \
  -v "${PWD}/results":/home/appuser/app/results \
  nophigene-drd4-analysis:latest
```

Then open:

- [http://127.0.0.1:8000](http://127.0.0.1:8000)

### Use the Docker launcher

If you still want the automated Docker flow, double-click:

- [Start NophiGene UI (Docker).cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Start%20NophiGene%20UI%20%28Docker%29.cmd)

When done, use:

- [Stop NophiGene UI (Docker).cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Stop%20NophiGene%20UI%20%28Docker%29.cmd)

## VS Code

VS Code is now aligned with the local-first setup:

- [settings.json](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/.vscode/settings.json:1) points to `.venv\Scripts\python.exe`
- [launch.json](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/.vscode/launch.json:1) includes:
  - a local UI launch config
  - a CLI launch config

## Troubleshooting

### Double-clicking the local launcher says dependencies are missing

Reinstall the app runtime:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-app.txt
```

### The browser does not open automatically

Open:

- [http://127.0.0.1:8000](http://127.0.0.1:8000)

### The local server failed to start

Check the local launcher logs in the repo root:

- `.nophigene-ui.log`
- `.nophigene-ui.err.log`

### Docker is still slow

That is now expected to be less severe than before, but local `.venv` launch is still the recommended path for iterative work.

### `pysam` still fails on Windows

That package remains optional and is not required for the current UI or CLI path.

## Recommended workflow now

For daily use:

1. Put your input files in `data/`.
2. Double-click [Start NophiGene UI.cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Start%20NophiGene%20UI.cmd).
3. Run the analysis from the browser.
4. Open outputs from `results/`.
5. Double-click [Stop NophiGene UI.cmd](/C:/Users/Mewxy/Desktop/YouTopy/NophiGene/nophigene-drd4-analysis/Stop%20NophiGene%20UI.cmd) when finished.
