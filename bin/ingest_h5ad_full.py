#!/usr/bin/env python3
"""Prepare a complete AnnData .h5ad for the pipeline without splitting."""

from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Path to input h5ad")
    parser.add_argument("--output", required=True, help="Output h5ad")
    parser.add_argument("--sample_id", default="complete", help="Fallback sample_id")
    parser.add_argument("--subsample", type=int, default=None, help="Subsample N cells from the complete object")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    print(f"[ingest_h5ad_full] Loading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    if args.subsample and adata.n_obs > args.subsample:
        idx = rng.choice(adata.n_obs, args.subsample, replace=False)
        adata = adata[idx, :].copy()
        print(f"[ingest_h5ad_full] Subsampled to {adata.n_obs} cells")

    # Downstream doublet/QC expects counts in X. Prefer an explicit counts layer,
    # then raw counts, then the existing X matrix.
    if "counts" in adata.layers:
        x = adata.layers["counts"]
        var = adata.var.copy()
    elif adata.raw is not None:
        x = adata.raw.X
        var = adata.raw.var.copy()
    else:
        x = adata.X
        var = adata.var.copy()

    adata = ad.AnnData(
        X=x.copy() if sparse.issparse(x) else np.asarray(x).copy(),
        obs=adata.obs.copy(),
        var=var,
    )
    if "sample_id" not in adata.obs.columns:
        adata.obs["sample_id"] = args.sample_id

    adata.write_h5ad(args.output)
    print(f"[ingest_h5ad_full] Wrote {args.output} ({adata.n_obs} cells x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()
