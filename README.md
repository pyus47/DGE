# Differential Gene Expression Explorer

This project is a Streamlit app for bulk RNA-seq differential gene expression analysis. It provides one interface for three widely used analysis packages:

- `DESeq2`
- `edgeR`
- `limma-voom`

The app accepts a raw counts matrix plus sample metadata, lets you define a contrast, then returns ranked genes, significance calls, plots, and downloadable normalized expression tables.

## Features

- Upload counts and metadata as CSV, TSV, or TXT
- Run one-click analysis with `DESeq2`, `edgeR`, or `limma-voom`
- Optional batch or covariate term in the design formula
- Volcano plot, MA plot, and p-value histogram
- Export full results, significant hits, and normalized counts
- Example dataset included for smoke testing

## Expected Input Format

### Counts matrix

The counts table should contain one row per gene and one column per sample.

| gene_id | ctrl_1 | ctrl_2 | treated_1 |
| --- | ---: | ---: | ---: |
| ENSG000001 | 125 | 130 | 218 |
| ENSG000002 | 40 | 51 | 39 |

Notes:

- The first identifier column defaults to `gene_id`
- Sample columns should contain raw counts, not TPM or FPKM
- Gene IDs must be unique

### Metadata table

The metadata table should contain one row per sample.

| sample | condition | batch |
| --- | --- | --- |
| ctrl_1 | control | A |
| ctrl_2 | control | A |
| treated_1 | treated | B |

Notes:

- The sample column defaults to `sample`
- The primary contrast column defaults to `condition`
- The optional adjustment column defaults to `batch`
- Sample names must match the count matrix columns exactly

## Setup

You need both Python and R available on your machine.

### 1. Create a Python environment

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2. Install required R packages

```powershell
Rscript .\install_r_packages.R
```

This installs:

- `DESeq2`
- `edgeR`
- `limma`
- `tibble`

If a method still shows as unavailable in the app, the usual causes are:

- `rpy2` is missing from the Python environment running Streamlit
- R is installed, but Streamlit is using a different R than the one where Bioconductor packages were installed
- the required Bioconductor packages were not installed into the R library visible to `rpy2`

### 3. Launch the app

```powershell
streamlit run .\app.py
```

## Deploy With Docker

If you want a reproducible environment with Python, R, and the required Bioconductor packages together, you can run the app in Docker.

### Build the image

```powershell
docker build -t dge-explorer .
```

### Start the container

```powershell
docker run --rm -p 8501:8501 dge-explorer
```

Then open `http://localhost:8501`.

### Start with Docker Compose (recommended)

```powershell
docker compose up --build -d app
```

Stop it with:

```powershell
docker compose down
```

### Live reload during development

If you want the app to reload when you edit local files, run the container with a bind mount and polling file watcher:

```powershell
docker run --rm -p 8501:8501 -v "${PWD}:/app" -e STREAMLIT_SERVER_FILE_WATCHER_TYPE=poll dge-explorer
```

Or using Docker Compose dev profile:

```powershell
docker compose --profile dev up --build -d app-dev
```

Stop dev mode:

```powershell
docker compose --profile dev down
```

Notes:

- Keep the first `-p 8501:8501` so the app is reachable at `http://localhost:8501`
- The `-v "${PWD}:/app"` mount makes code changes visible inside the container
- Polling watcher is more reliable on Docker Desktop file mounts

## Analysis Notes

- `DESeq2` is often a strong default for standard two-group bulk RNA-seq workflows.
- `edgeR` is commonly used for negative binomial GLM-based count analysis.
- `limma-voom` is popular when you want limma's modeling workflow on RNA-seq data.
- The app currently assumes a single primary contrast and one optional covariate.

## Quality Control Included

- Library size per sample
- Detected genes per sample
- PCA on log-CPM values
- Sample-to-sample distance heatmap

## Project Structure

- `app.py`: Streamlit user interface
- `dge_app/analysis.py`: R-backed differential expression execution
- `dge_app/qc.py`: pre-analysis QC metrics and transforms
- `dge_app/validation.py`: input checks and design validation
- `dge_app/plots.py`: result visualizations
- `install_r_packages.R`: helper script for Bioconductor dependencies
- `Dockerfile`: containerized deployment recipe

## Limitations

- This is intended for bulk RNA-seq count matrices, not single-cell workflows
- The local machine must have a working R installation discoverable by `rpy2`
- Complex multi-factor or interaction designs are not yet exposed in the UI
