# CellFlow

<div align="center">
  <img src="docs/assets/logo.png" alt="CellFlow Logo" width="200" />
</div>

## Overview

A reproducible Nextflow pipeline for end-to-end single-cell RNA-seq analysis.
Takes a complete Seurat `.rds` or Scanpy `.h5ad` object and walks it through
doublet removal (DoubletFinder), QC filtering, normalization, highly-variable-gene
selection, PCA, optional Harmony batch correction, UMAP, and reference-based
cell-type annotation (SingleR / `HumanPrimaryCellAtlasData`). Outputs an
annotated AnnData or Seurat object plus diagnostic UMAP plots.

## Pipeline Workflow

```mermaid
flowchart TD
    Input([Input: .rds or .h5ad]) --> Ingest
    
    subgraph Preprocessing
        Ingest --> DoubletFinder[DoubletFinder<br/><i>(Optional)</i>]
        DoubletFinder --> QC[QC Filter]
        QC --> Normalize[Normalize, HVG & PCA]
    end
    
    subgraph Embedding & Integration
        Normalize --> UMAP_Pre[Pre-integration UMAP]
        UMAP_Pre --> HarmonyCheck{Skip Harmony?}
        HarmonyCheck -- No --> HarmonyRun[Harmony Batch Correction]
        HarmonyRun --> UMAP_Post[Post-integration UMAP]
        HarmonyCheck -- Yes --> Annotate
        UMAP_Post --> Annotate
    end
    
    subgraph Annotation
        Annotate[SingleR Cell Type Annotation] --> Plot[Plot Final UMAP]
        Plot --> Export[Export .rds / .h5ad]
    end
    
    Export --> Outputs([Outputs:<br/>final.rds/h5ad, UMAPs, manifest.json])
```

Every run is replayable. A single `--seed` is threaded through every stochastic
step (PCA, neighbors, UMAP, Harmony, DoubletFinder, SingleR); BLAS/OpenMP threads
are pinned for bit-exact numerics; the conda env is captured in `conda.lock.txt`;
DoubletFinder is installed at a pinned commit SHA; and each run emits
`manifest.json` recording SHA-256 hashes of all inputs, outputs, and references.
Same seed + same lockfile + same architecture → byte-identical outputs. Runs on
Linux x86_64 natively and on Apple Silicon Macs under Rosetta 2.

## Quick start

```bash
# 1. Clone
git clone https://github.com/Rengoku02/scrna-seq-pipeline.git
cd scrna-seq-pipeline

# 2. One-time conda env setup (installs Nextflow + Java + R + Python deps).
#    Apple Silicon Mac? Install Rosetta 2 first:
#       softwareupdate --install-rosetta --agree-to-license
bash setup.sh
conda activate scrna-demo

# 3. Smoke-test on the bundled pbmc_small dataset (~80 cells, ~1-2 min)
nextflow run main.nf -profile test

# 4. Run on your own data
nextflow run main.nf --input /path/to/your.rds --outdir results
```

After step 4 you should have `results_test/{final.rds, umap_pre.png,
umap_final.png, manifest.json, ...}`.

## Run on your own data

```bash
nextflow run main.nf --input /path/to/data.rds  --outdir results
nextflow run main.nf --input /path/to/data.h5ad --outdir results
```

The pipeline branches on the input file extension. Use `nextflow run main.nf
-resume` to skip already-cached stages on subsequent runs.

All Nextflow params are passed with `--<name> <value>` (double dash, snake_case).

**Flow control**

| Param                       | What it does                                                      |
|-----------------------------|-------------------------------------------------------------------|
| `--subsample N`             | cap object to N cells during ingest                               |
| `--skip_harmony`            | skip the Harmony batch-correction stage (default: Harmony runs)   |
| `--skip_doublet`            | bypass DoubletFinder                                              |
| `--output_format FMT`       | `rds` / `h5ad` / `both` (default: matches input type)             |
| `--keep_intermediates`      | keep per-stage .h5ad files (default: cleaned up at end)           |

**Doublet detection**

| Param               | What it does                                                |
|---------------------|-------------------------------------------------------------|
| `--doublet_rate F`  | expected doublet rate (default: auto from cell count)       |

**QC thresholds** (all optional — scripts have sensible defaults)

| Param               | What it does                                                |
|---------------------|-------------------------------------------------------------|
| `--min_genes N`     | min genes per cell                                          |
| `--max_genes N`     | max genes per cell                                          |
| `--min_umis N`      | min UMIs per cell                                           |
| `--max_mito_pct F`  | max mitochondrial %                                         |
| `--max_ribo_pct F`  | max ribosomal %                                             |

**Embedding / integration**

| Param               | What it does                                                |
|---------------------|-------------------------------------------------------------|
| `--hvg_flavor`      | `seurat` (default) / `seurat_v3` / `cell_ranger`            |
| `--n_top_genes N`   | number of HVGs to keep                                      |
| `--n_pcs N`         | number of principal components                              |
| `--harmony_theta F` | Harmony diversity penalty                                   |

**Annotation**

| Param                 | What it does                                              |
|-----------------------|-----------------------------------------------------------|
| `--singler_ref NAME`  | SingleR reference (e.g. `BlueprintEncodeData`)            |
| `--singler_labels C`  | reference label column (`label.main` / `label.fine`)      |

**Reproducibility**

| Param           | What it does                                                   |
|-----------------|----------------------------------------------------------------|
| `--seed N`      | random seed threaded through every stochastic step (default 42)|
| `--threads N`   | BLAS/OpenMP thread cap (default 1 for bit-exact determinism)   |

## The Harmony stage

By default Harmony always runs and SingleR annotates the corrected embedding.
To skip Harmony and annotate the uncorrected UMAP, pass `--skip_harmony`.
The pre-integration UMAP is written to `umap_pre.png` regardless, so you can
inspect for batch effects before deciding whether your next run should keep
or skip Harmony.

## Outputs

Default behavior: intermediate `.h5ad` files are **cleaned up** at the end.
Pass `--keep_intermediates` to keep them (useful for debugging).

After a successful run, `--outdir` contains:

```
final.rds  and/or  final.h5ad      # depends on --output_format
umap_pre.png                       # pre-integration UMAP
umap_post.png                      # only if Harmony ran
umap_final.png                     # colored by singler_labels
doublet/complete_stats.txt         # doublet rate, singlets retained, pK
qc/complete_stats.txt              # cells before/after, thresholds applied
manifest.json                      # run provenance + hashes
pipeline_trace.txt                 # Nextflow per-process trace
pipeline_report.html               # Nextflow execution report
```

With `--keep_intermediates`, you additionally get:

```
loaded/complete.h5ad
doublet/complete.h5ad
qc/complete.h5ad
normalized/complete.h5ad
merged_pre.h5ad
merged_harmony.h5ad                # only if Harmony ran
merged_umap.h5ad                   # only if Harmony ran
```

Quick sanity check on the labels:

```bash
python -c "
import anndata as ad
a = ad.read_h5ad('results/final.h5ad')
print(a.obs.filter(like='label').head())
print(a.obs.iloc[:, -1].value_counts().head(10))
"
```

## Files

| Path                         | What                                                  |
|------------------------------|-------------------------------------------------------|
| `main.nf`                    | Nextflow workflow entry point                         |
| `nextflow.config`            | Pipeline params, env vars, profiles                   |
| `modules/local/*.nf`         | One file per process (12 stages)                      |
| `conf/base.config`           | Default per-process settings                          |
| `conf/test.config`           | Smoke-test profile (uses `test/data/test.rds`)        |
| `test/data/test.rds`         | Bundled `pbmc_small` (~80 cells) for smoke test       |
| `environment.yml`            | Conda env spec (R + Python together)                  |
| `conda.lock.txt`             | Byte-exact env replay spec                            |
| `external_pins.txt`          | Pinned DoubletFinder commit SHA                       |
| `setup.sh`                   | Build env + install DoubletFinder from GitHub         |
| `bin/ingest_seurat_h5ad.R`   | complete Seurat .rds → one h5ad                       |
| `bin/ingest_h5ad_full.py`    | complete h5ad → one pipeline-ready h5ad               |
| `bin/run_doubletfinder.R`    | DoubletFinder doublet removal                         |
| `bin/cell_qc_filter.py`      | QC thresholds on genes/UMIs/mito%                     |
| `bin/normalize_hvg.py`       | Full-gene `normalize_total` + `log1p`; HVG work copy for PCA |
| `bin/umap_checkpoint.py`     | PCA + pre-integration UMAP                            |
| `bin/run_harmony.py`         | Harmony on PCA → `obsm['X_pca_harmony']`              |
| `bin/post_harmony_umap.py`   | Neighbors + UMAP on Harmony embedding + PNG          |
| `bin/run_singler.R`          | SingleR annotation (HumanPrimaryCellAtlasData)        |
| `bin/plot_final_umap.py`     | Final UMAP colored by SingleR labels                  |
| `bin/export_h5ad_to_rds.R`   | Rebuild a Seurat .rds from final.h5ad (for `--output-format rds`) |

All `bin/` scripts take `--input <h5ad> --output <h5ad>` (or equivalent) —
each is runnable standalone for debugging.

## Reproducing a run

Every run writes `manifest.json` to `--outdir` capturing:

- input file path + SHA-256
- each output file's SHA-256
- resolved parameters (every CLI flag)
- seed, thread counts, environment variables
- conda env name, Python / R versions, key package versions
- SingleR reference and celldex version
- DoubletFinder commit SHA (from `external_pins.txt`)
- pipeline git SHA (if under version control)

To reproduce a prior run byte-for-byte:

1. **Rebuild the env from the lockfile** (not from `environment.yml`).
   `setup.sh` writes `demo/conda.lock.txt` automatically; use:
   ```bash
   conda create --name scrna-demo --file demo/conda.lock.txt
   ```
   This pins exact versions + build hashes.

2. **Replay with the same seed and thread count** from the old manifest:
   ```bash
   nextflow run main.nf \
     --input   "$(jq -r .input.path results/manifest.json)" \
     --outdir  rerun \
     --seed    "$(jq -r .environment.seed    results/manifest.json)" \
     --threads "$(jq -r .environment.threads results/manifest.json)"
   ```

3. **Verify output parity:**
   ```bash
   diff <(jq -r '.outputs[] | "\(.sha256)  \(.path)"' results/manifest.json | sort) \
        <(jq -r '.outputs[] | "\(.sha256)  \(.path)"' rerun/manifest.json   | sort)
   ```
   Empty diff → byte-identical outputs.

`--threads 1` makes runs reproducible across machines; using more threads
speeds things up but BLAS reordering makes embeddings drift slightly.

## Troubleshooting

- **`SingleR` is slow the first time** — `HumanPrimaryCellAtlasData()` downloads
  ~200 MB from ExperimentHub on first call. Pre-warm it:
  `Rscript -e 'celldex::HumanPrimaryCellAtlasData()'`.
- **`scikit-misc` / `seurat_v3` HVG fails on ARM64** — stick to the default
  `--hvg_flavor seurat`.
- **DoubletFinder install failed** — re-run `setup.sh`, or use `--skip_doublet`.
- **Apple Silicon: `setup.sh` fails inside conda** — install Rosetta 2 first:
  `softwareupdate --install-rosetta --agree-to-license`. The conda env runs
  `osx-64` builds under Rosetta because Bioconductor packages have no native
  `osx-arm64` binaries.
- **Want to start over?** `nextflow run` keeps a `work/` cache. Delete it
  with `rm -rf work/ .nextflow*` to force a clean re-run.
- **`anndata` (R) can't find Python** — make sure the `scrna-demo` conda env is
  active; `r-anndata` uses `reticulate` and needs the co-installed Python.
