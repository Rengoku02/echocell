#!/usr/bin/env python3
"""
Normalization, HVG selection, regression, scaling, and PCA for scRNA-seq data.

Takes a QC-passed AnnData object and produces a full-gene normalized object
with HVGs marked and PCA computed from a scaled HVG working copy.
"""

import argparse
import os
import random

# Pin BLAS/OpenMP threads for deterministic numerics. Must happen before
# numpy/scipy are imported by anything.
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="Normalize, select HVGs, scale, PCA")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--output", required=True, help="Output h5ad file")
    parser.add_argument("--sample_id", required=True, help="Sample identifier")
    parser.add_argument("--target_sum", type=float, default=10000, help="Library-size normalization target")
    parser.add_argument("--n_top_genes", type=int, default=2000, help="Number of HVGs to select")
    parser.add_argument("--hvg_flavor", default="seurat", help="HVG selection method")
    parser.add_argument("--regress_out", default="pct_counts_mt", help="Comma-separated variables to regress out")
    parser.add_argument("--scale_max_value", type=float, default=10, help="Max value after scaling")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for PCA / HVG")
    return parser.parse_args()


def main():
    args = parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    adata = sc.read_h5ad(args.input)

    # Library-size normalization
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)

    # Preserve full log-normalized expression for annotation. Do not subset or
    # scale adata.X; downstream SingleR consumes this matrix as logcounts.
    adata.raw = adata

    # Highly variable gene selection
    if args.hvg_flavor == "seurat_v3":
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=args.n_top_genes,
            flavor=args.hvg_flavor,
            subset=False,
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=args.n_top_genes,
            flavor=args.hvg_flavor,
        )
    # seurat_v3 has no random_state arg; np.random.seed() above covers it.

    # Compute PCA on a temporary HVG-only object. This avoids densifying the
    # full expression matrix during regression/scaling.
    hvg_mask = adata.var["highly_variable"].to_numpy()
    work = adata[:, hvg_mask].copy()

    # Regress out confounders
    regress_vars = [v.strip() for v in args.regress_out.split(",") if v.strip()]
    if regress_vars:
        # Only regress variables that exist in .obs
        valid_vars = [v for v in regress_vars if v in work.obs.columns]
        if valid_vars:
            sc.pp.regress_out(work, valid_vars)

    # Scale
    sc.pp.scale(work, max_value=args.scale_max_value)

    # PCA
    sc.tl.pca(work, random_state=args.seed)
    adata.obsm["X_pca"] = work.obsm["X_pca"]
    adata.uns["pca"] = work.uns["pca"]
    pcs = np.zeros((adata.n_vars, work.varm["PCs"].shape[1]), dtype=work.varm["PCs"].dtype)
    pcs[hvg_mask, :] = work.varm["PCs"]
    adata.varm["PCs"] = pcs

    # Write output
    adata.write_h5ad(args.output)

    print(
        f"[normalize_hvg] {args.sample_id}: "
        f"{adata.n_obs} cells, {int(hvg_mask.sum())} HVGs, full matrix retained"
    )


if __name__ == "__main__":
    main()
