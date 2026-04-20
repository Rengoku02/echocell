#!/usr/bin/env python3
"""
UMAP checkpoint: merge per-sample AnnDatas and compute UMAP on uncorrected PCA.

This is a critical checkpoint that:
1. Merges all per-sample normalized AnnDatas (concatenate with batch_key='sample_id')
2. Recomputes a merged PCA from HVGs, then computes neighbors + UMAP
3. Generates a UMAP plot colored by sample_id
4. Saves the merged AnnData and the plot
"""

import argparse
import os
import random

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="UMAP checkpoint — merge and visualize")
    parser.add_argument("--inputs", required=True, help="Comma-separated list of input h5ad files")
    parser.add_argument("--output", required=True, help="Output merged h5ad file")
    parser.add_argument("--plot_output", required=True, help="Output UMAP plot (PNG)")
    parser.add_argument("--checkpoint_label", default="pre_integration", help="Label for this checkpoint")
    parser.add_argument("--n_top_genes", type=int, default=2000, help="Number of merged HVGs for PCA")
    parser.add_argument("--hvg_flavor", default="seurat", help="Merged HVG selection method")
    parser.add_argument("--n_pcs", type=int, default=30, help="Number of PCs for neighbors/UMAP")
    parser.add_argument("--scale_max_value", type=float, default=10, help="Max value after scaling")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for PCA/neighbors/UMAP")
    return parser.parse_args()


def main():
    args = parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    # Load all per-sample AnnDatas
    input_files = [f.strip() for f in args.inputs.split(",") if f.strip()]
    adatas = []
    for fpath in input_files:
        a = sc.read_h5ad(fpath)
        # Ensure sample_id is in .obs for batch coloring
        if "sample_id" not in a.obs.columns:
            # Try to infer from filename
            sample_id = fpath.replace("_normalized.h5ad", "").split("/")[-1]
            a.obs["sample_id"] = sample_id
        adatas.append(a)

    # Merge / concatenate
    if len(adatas) == 1:
        adata = adatas[0].copy()
    else:
        adata = ad.concat(adatas, join="inner", merge="same")

    print(f"[umap_checkpoint] Merged {len(adatas)} samples: {adata.n_obs} cells, {adata.n_vars} genes")

    # Recompute PCA on the merged expression matrix. Per-sample PCA embeddings
    # are not comparable after concatenation.
    for key in ("X_pca", "X_pca_harmony"):
        if key in adata.obsm:
            del adata.obsm[key]
    if "pca" in adata.uns:
        del adata.uns["pca"]
    if "PCs" in adata.varm:
        del adata.varm["PCs"]

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(args.n_top_genes, adata.n_vars),
        flavor=args.hvg_flavor,
    )
    hvg_mask = adata.var["highly_variable"].to_numpy()
    work = adata[:, hvg_mask].copy()
    sc.pp.scale(work, max_value=args.scale_max_value)
    sc.tl.pca(work, n_comps=min(args.n_pcs, work.n_obs - 1, work.n_vars - 1),
              random_state=args.seed)
    adata.obsm["X_pca"] = work.obsm["X_pca"]
    adata.uns["pca"] = work.uns["pca"]

    # Compute neighbors and UMAP on uncorrected PCA
    sc.pp.neighbors(adata, n_pcs=min(args.n_pcs, adata.obsm["X_pca"].shape[1]),
                    use_rep="X_pca", random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)

    # Generate UMAP plot colored by sample_id
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="sample_id", ax=ax, show=False, title=f"UMAP — {args.checkpoint_label}")
    plt.tight_layout()
    fig.savefig(args.plot_output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Save merged AnnData
    adata.write_h5ad(args.output)

    print(f"[umap_checkpoint] Saved {args.output} and {args.plot_output}")


if __name__ == "__main__":
    main()
