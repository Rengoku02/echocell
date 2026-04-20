#!/usr/bin/env python3
"""Compute neighbors + UMAP on the Harmony-corrected embedding and plot by sample."""

import argparse
import os
import random

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="h5ad with obsm['X_pca_harmony']")
    p.add_argument("--output", required=True)
    p.add_argument("--plot_output", required=True)
    p.add_argument("--use_rep", default="X_pca_harmony")
    p.add_argument("--n_neighbors", type=int, default=15)
    p.add_argument("--seed", type=int, default=42, help="Random seed for neighbors/UMAP")
    args = p.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    adata = sc.read_h5ad(args.input)
    if args.use_rep not in adata.obsm:
        raise SystemExit(f"ERROR: {args.use_rep} not in .obsm. Present: {list(adata.obsm.keys())}")

    sc.pp.neighbors(adata, use_rep=args.use_rep, n_neighbors=args.n_neighbors,
                    random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)

    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="sample_id", ax=ax, show=False,
               title="Post-Harmony UMAP (colored by sample)", frameon=False)
    fig.tight_layout()
    fig.savefig(args.plot_output, dpi=150)
    plt.close(fig)

    adata.write_h5ad(args.output)
    print(f"[post_harmony_umap] wrote {args.output} and {args.plot_output}")


if __name__ == "__main__":
    main()
