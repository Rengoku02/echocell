# scRNA-seq Analysis Pipeline

End-to-end single-cell RNA-seq analysis from a complete Seurat `.rds` or
Scanpy `.h5ad` object through to cell-type-annotated output. Runs locally
inside a single conda environment.

## Stages

```
ingest (complete object)  ->  doublet (DoubletFinder)  ->  QC filter
   ->  normalize + HVG + PCA  ->  pre-integration UMAP
   ->  [checkpoint: view UMAP, confirm Harmony?]
   ->  Harmony + post-integration UMAP   (optional)
   ->  SingleR annotation  ->  final.h5ad + umap_final.png
```

## Setup (one-time)

```bash
cd demo
bash setup.sh                 # builds the conda env; installs DoubletFinder
conda activate scrna-demo
```

`setup.sh` uses `mamba` if present, else `conda`. Takes ~5-10 min.

## Run

Two entry points — both consume one complete object:

```bash
# Seurat .rds
bash run.sh --from-seurat /path/to/merged.rds --outdir results --subsample 1000

# Scanpy .h5ad
bash run.sh --from-h5ad   /path/to/data.h5ad  --outdir results --subsample 1000
```

Useful flags:

**Flow control**

| Flag                        | What it does                                                      |
|-----------------------------|-------------------------------------------------------------------|
| `--subsample N`             | cap object to N cells during ingest                               |
| `--no-prompt`               | skip the Harmony y/n checkpoint (always run Harmony)              |
| `--skip-doublet`            | bypass DoubletFinder (faster demo, less realistic)                |
| `--output-format FMT`       | `rds` / `h5ad` / `both` (default: matches input type)             |
| `--keep-intermediates`      | keep per-stage .h5ad files (default: cleaned up at end)           |

**Doublet detection**

| Flag                | What it does                                                |
|---------------------|-------------------------------------------------------------|
| `--doublet-rate F`  | expected doublet rate (default: auto from cell count)       |

**QC thresholds** (all optional — scripts have sensible defaults)

| Flag                | What it does                                                |
|---------------------|-------------------------------------------------------------|
| `--min-genes N`     | min genes per cell                                          |
| `--max-genes N`     | max genes per cell                                          |
| `--min-umis N`      | min UMIs per cell                                           |
| `--max-mito-pct F`  | max mitochondrial %                                         |
| `--max-ribo-pct F`  | max ribosomal %                                             |

**Embedding / integration**

| Flag                | What it does                                                |
|---------------------|-------------------------------------------------------------|
| `--hvg-flavor`      | `seurat` (default) / `seurat_v3` / `cell_ranger`            |
| `--n-top-genes N`   | number of HVGs to keep                                      |
| `--n-pcs N`         | number of principal components                              |
| `--harmony-theta F` | Harmony diversity penalty                                   |

**Annotation**

| Flag                  | What it does                                              |
|-----------------------|-----------------------------------------------------------|
| `--singler-ref NAME`  | SingleR reference (e.g. `BlueprintEncodeData`)            |
| `--singler-labels C`  | reference label column (`label.main` / `label.fine`)      |

**Reproducibility**

| Flag            | What it does                                                   |
|-----------------|----------------------------------------------------------------|
| `--seed N`      | random seed threaded through every stochastic step (default 42)|
| `--threads N`   | BLAS/OpenMP thread cap (default 1 for bit-exact determinism)   |

## The Harmony checkpoint

After the pre-integration UMAP is written, the script opens `umap_pre.png`
and asks:

```
Pre-integration UMAP: results/umap_pre.png
Run Harmony? [Y/n]
```

- Answer **Y** (or just Enter) → Harmony runs, post-integration UMAP is produced,
  SingleR is run on the corrected embedding.
- Answer **n** → Harmony is skipped, SingleR runs on the uncorrected UMAP.

Batch correction is a decision, not a default — inspect the UMAP first.

## Outputs

Default behavior: intermediate `.h5ad` files are **cleaned up** at the end.
Pass `--keep-intermediates` to keep them (useful for debugging).

After a successful run, `--outdir` contains:

```
final.rds  and/or  final.h5ad      # depends on --output-format
umap_pre.png                       # pre-integration UMAP
umap_post.png                      # only if Harmony ran
umap_final.png                     # colored by singler_labels
doublet/complete_stats.txt         # doublet rate, singlets retained, pK
qc/complete_stats.txt              # cells before/after, thresholds applied
log.txt                            # full pipeline stdout/stderr
manifest.json                      # run provenance + hashes
```

With `--keep-intermediates`, you additionally get:

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
| `environment.yml`            | Conda env spec (one env, R + Python together)         |
| `setup.sh`                   | Build env + install DoubletFinder from GitHub         |
| `run.sh`                     | Main pipeline driver                                  |
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
   bash run.sh \
     --from-seurat "$(jq -r .input.path results/manifest.json)" \
     --outdir rerun \
     --seed    "$(jq -r .environment.seed results/manifest.json)" \
     --threads "$(jq -r .environment.threads results/manifest.json)" \
     --no-prompt
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
  `--hvg-flavor seurat`.
- **DoubletFinder install failed** — re-run `setup.sh`, or use `--skip-doublet`.
- **`anndata` (R) can't find Python** — make sure the `scrna-demo` conda env is
  active; `r-anndata` uses `reticulate` and needs the co-installed Python.
