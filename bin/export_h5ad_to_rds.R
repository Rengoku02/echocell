#!/usr/bin/env Rscript

# =============================================================================
#  Export final.h5ad to a Seurat .rds
# =============================================================================
#  The pipeline's lingua franca is AnnData. When the user wants Seurat back
#  out, we rebuild a Seurat object from the h5ad's X matrix, .obs, and .obsm
#  embeddings. Note: X at this stage is log-normalized, not raw counts —
#  it goes into the 'data' slot, with the same matrix duplicated into 'counts'
#  so the object is self-consistent for downstream Seurat ops.

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(anndata)
    library(Matrix)
})

option_list <- list(
    make_option("--input", type = "character", help = "Input h5ad file (e.g. final.h5ad)"),
    make_option("--output", type = "character", help = "Output .rds file"),
    make_option("--seed", type = "integer", default = 42, help = "Random seed")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

cat(sprintf("[export_h5ad_to_rds] Reading %s\n", opt$input))
adata <- read_h5ad(opt$input)

# genes x cells; Seurat expects this orientation
expr <- as(t(adata$X), "CsparseMatrix")
colnames(expr) <- rownames(adata$obs)
rownames(expr) <- rownames(adata$var)

sobj <- CreateSeuratObject(counts = expr, meta.data = adata$obs, project = "scrna_demo")
sobj <- SetAssayData(sobj, assay = "RNA", layer = "data", new.data = expr)

add_reduction <- function(sobj, key_in_obsm, reduction_name, key_prefix) {
    if (!(key_in_obsm %in% names(adata$obsm))) return(sobj)
    emb <- as.matrix(adata$obsm[[key_in_obsm]])
    rownames(emb) <- rownames(adata$obs)
    colnames(emb) <- paste0(key_prefix, seq_len(ncol(emb)))
    sobj[[reduction_name]] <- CreateDimReducObject(
        embeddings = emb, key = paste0(key_prefix, "_"), assay = "RNA"
    )
    sobj
}

sobj <- add_reduction(sobj, "X_umap", "umap", "UMAP")
sobj <- add_reduction(sobj, "X_pca", "pca", "PC")
sobj <- add_reduction(sobj, "X_pca_harmony", "harmony", "harmony")

saveRDS(sobj, opt$output)
cat(sprintf("[export_h5ad_to_rds] Wrote %s (%d cells, %d genes, reductions: %s)\n",
            opt$output, ncol(sobj), nrow(sobj),
            paste(Reductions(sobj), collapse = ", ")))
