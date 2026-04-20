#!/usr/bin/env python3
"""
Cell-level QC filtering for scRNA-seq data.

Computes QC metrics (mitochondrial %, ribosomal %, gene counts, UMI counts)
and filters cells based on parameterized thresholds. All thresholds are
passed via CLI arguments — nothing is hardcoded.
"""

import argparse
import sys

import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="Cell-level QC filtering")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--output", required=True, help="Output h5ad file (QC-passed cells)")
    parser.add_argument("--stats_output", required=True, help="Output text file with QC statistics")
    parser.add_argument("--sample_id", required=True, help="Sample identifier")
    parser.add_argument("--min_genes", type=int, default=200, help="Min genes per cell")
    parser.add_argument("--max_genes", type=int, default=5000, help="Max genes per cell")
    parser.add_argument("--min_umis", type=int, default=500, help="Min UMIs per cell")
    parser.add_argument("--max_umis", type=int, default=None, help="Max UMIs per cell (null=no limit)")
    parser.add_argument("--max_mito_pct", type=float, default=20.0, help="Max mitochondrial %")
    parser.add_argument("--max_ribo_pct", type=float, default=None, help="Max ribosomal % (null=no filter)")
    parser.add_argument("--min_cells_per_gene", type=int, default=3, help="Min cells expressing a gene")
    return parser.parse_args()


def main():
    args = parse_args()

    # Validate input
    adata = sc.read_h5ad(args.input)
    cells_before = adata.n_obs
    genes_before = adata.n_vars

    # Annotate mitochondrial and ribosomal genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )

    # Build cell-level filter mask
    keep = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes)
        & (adata.obs["n_genes_by_counts"] <= args.max_genes)
        & (adata.obs["total_counts"] >= args.min_umis)
        & (adata.obs["pct_counts_mt"] <= args.max_mito_pct)
    )

    if args.max_umis is not None:
        keep &= adata.obs["total_counts"] <= args.max_umis

    if args.max_ribo_pct is not None:
        keep &= adata.obs["pct_counts_ribo"] <= args.max_ribo_pct

    # Apply cell filter
    adata = adata[keep].copy()

    # Gene-level filter: remove genes expressed in too few cells
    sc.pp.filter_genes(adata, min_cells=args.min_cells_per_gene)

    cells_after = adata.n_obs
    genes_after = adata.n_vars

    # Write filtered AnnData
    adata.write_h5ad(args.output)

    with open(args.stats_output, "w") as f:
        f.write("=== Cell QC filter stats ===\n")
        f.write(f"sample_id                  : {args.sample_id}\n")
        f.write(f"cells_before               : {cells_before}\n")
        f.write(f"cells_after                : {cells_after}\n")
        f.write(f"cells_removed              : {cells_before - cells_after}\n")
        f.write(f"genes_before               : {genes_before}\n")
        f.write(f"genes_after                : {genes_after}\n")
        f.write("--- filters applied ---\n")
        f.write(f"min_genes                  : {args.min_genes}\n")
        f.write(f"max_genes                  : {args.max_genes}\n")
        f.write(f"min_umis                   : {args.min_umis}\n")
        f.write(f"max_umis                   : {args.max_umis}\n")
        f.write(f"max_mito_pct               : {args.max_mito_pct}\n")
        f.write(f"max_ribo_pct               : {args.max_ribo_pct}\n")
        f.write(f"min_cells_per_gene         : {args.min_cells_per_gene}\n")
        f.write("--- post-filter medians ---\n")
        f.write(f"median_genes_post_filter   : {float(adata.obs['n_genes_by_counts'].median()):.1f}\n")
        f.write(f"median_umis_post_filter    : {float(adata.obs['total_counts'].median()):.1f}\n")
        f.write(f"median_mito_pct_post_filter: {float(adata.obs['pct_counts_mt'].median()):.3f}\n")

    print(
        f"[cell_qc_filter] {args.sample_id}: {cells_before} → {cells_after} cells "
        f"({cells_before - cells_after} removed), "
        f"{genes_before} → {genes_after} genes"
    )


if __name__ == "__main__":
    main()
