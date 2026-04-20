#!/usr/bin/env python3
"""Plot the final UMAP colored by SingleR labels."""

import argparse
import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="final.h5ad with singler_labels in .obs and UMAP in .obsm")
    p.add_argument("--plot_output", required=True)
    p.add_argument("--color", default="singler_labels",
                   help=".obs column to color by (default: singler_labels)")
    args = p.parse_args()

    adata = sc.read_h5ad(args.input)

    if "X_umap" not in adata.obsm:
        raise SystemExit("ERROR: X_umap missing in .obsm — run UMAP before annotation plot.")
    if args.color not in adata.obs.columns:
        raise SystemExit(f"ERROR: {args.color} not in .obs. Present: {list(adata.obs.columns)}")

    fig, ax = plt.subplots(figsize=(10, 7))
    sc.pl.umap(adata, color=args.color, ax=ax, show=False,
               title=f"Final UMAP (colored by {args.color})",
               frameon=False, legend_loc="right margin", legend_fontsize=7)
    fig.tight_layout()
    fig.savefig(args.plot_output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"[plot_final_umap] wrote {args.plot_output}")


if __name__ == "__main__":
    main()
