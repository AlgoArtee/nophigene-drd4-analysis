# nophigene-drd4-analysis

Pipeline for analyzing DRD4 genetic & epigenetic data.  


## Prerequisites

- **Docker Engine ≥ 20.10** installed and running  
- A `data/` directory at the project root containing your input files:
  ```bash
  mkdir -p data
  cp /path/to/drd4.vcf.gz    data/drd4.vcf.gz
  cp /path/to/drd4_meth.bed  data/drd4_meth.bed
  # (Optional) population stats JSON
  cp /path/to/gnomad.json     data/gnomad.json
  ```
- An empty `results/` directory for outputs:
  ```bash
  mkdir -p results
  ```

> **Note:** `data/` and `results/` are listed in `.gitignore` to avoid committing large files.

---

## Build the Docker Image

From the project root (where `Dockerfile` lives), run:

```bash
docker build -t nophigene-drd4-analysis:latest .
```

- `-t nophigene-drd4-analysis:latest` tags the image for easy reference.

---

## Run the Analysis

Mount your `data/` and `results/` directories into the container and invoke the analysis:

```bash
docker run --rm -it \
  -v "${PWD}/data":/home/appuser/app/data \
  -v "${PWD}/results":/home/appuser/app/results \
  nophigene-drd4-analysis:latest \
    --vcf data/drd4.vcf.gz \
    --bed data/drd4_meth.bed \
    --out results/drd4_report.html
```

- `--rm` removes the container after it finishes  
- `-it` runs in interactive mode (displays logs)  
- `-v` flags mount host folders into the container  
- The final flags (`--vcf`, `--bed`, `--out`) correspond to the CLI in `src/analysis.py`

---

## Verify Success

After the container exits, list your results:

```bash
ls results/
# expect: drd4_report.html
```

Open the report in your browser.

```bash
open results/drd4_report.html      # macOS
xdg-open results/drd4_report.html   # Linux
```

---

## Next Steps

1. **Development:**  
   Reopen this project in VS Code with the Dev Container (`.devcontainer/`) for interactive debugging, notebook support, and consistent environment.  
2. **Testing:**  
   Inside the container, run:
   ```bash
   pytest
   ```
   to validate all imports and any unit tests you add under `tests/`.  
3. **CI/CD:**  
   Automate your workflow by adding a GitHub Actions pipeline that:
   - Builds the production Docker image
   - Runs `pytest`
   - Scans for vulnerabilities (e.g., via Trivy or Snyk)

Happy analyzing!
