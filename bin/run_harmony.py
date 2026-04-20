#!/usr/bin/env python3
"""
Batch integration using Harmony.

Corrects the PCA embedding for batch effects and stores the corrected PCs
in adata.obsm['X_pca_harmony'].
"""

import argparse
import os
import random
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import scanpy as sc
import harmonypy as hm


def parse_args():
    parser = argparse.ArgumentParser(description="Harmony batch integration")
    parser.add_argument("--input", required=True, help="Input h5ad (merged, with PCA)")
    parser.add_argument("--output", required=True, help="Output h5ad with corrected embedding")
    parser.add_argument("--batch_key", default="sample_id", help="Batch key in .obs")
    parser.add_argument("--theta", type=float, default=2.0, help="Harmony theta (diversity penalty)")
    parser.add_argument("--max_iter", type=int, default=20, help="Max Harmony iterations")
    parser.add_argument("--n_pcs", type=int, default=30, help="Number of PCs to use")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for Harmony")
    return parser.parse_args()


def main():
    args = parse_args()

    # harmonypy does not expose a seed arg; seeding NumPy globally is the
    # only way to make its k-means init deterministic.
    random.seed(args.seed)
    np.random.seed(args.seed)

    adata = sc.read_h5ad(args.input)

    # Ensure PCA is computed
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, n_comps=args.n_pcs, random_state=args.seed)

    # Call harmonypy directly. The scanpy.external wrapper mishandles
    # Z_corr's shape across harmonypy versions and raises a shape error.
    pcs = adata.obsm["X_pca"][:, : args.n_pcs]
    ho = hm.run_harmony(
        pcs,
        adata.obs,
        vars_use=[args.batch_key],
        theta=args.theta,
        max_iter_harmony=args.max_iter,
    )
    z = np.asarray(ho.Z_corr)
    if z.shape[0] == pcs.shape[1] and z.shape[1] == pcs.shape[0]:
        z = z.T
    adata.obsm["X_pca_harmony"] = z

    print(
        f"[run_harmony] Integrated {adata.n_obs} cells across "
        f"{adata.obs[args.batch_key].nunique()} batches. "
        f"Corrected embedding stored in obsm['X_pca_harmony']."
    )

    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
