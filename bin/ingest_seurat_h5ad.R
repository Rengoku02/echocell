#!/usr/bin/env Rscript

# Convert a complete Seurat object to one h5ad without splitting by sample.

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(anndata)
})

option_list <- list(
    make_option("--rds", type = "character", help = "Path to a Seurat .rds file"),
    make_option("--output", type = "character", help = "Output h5ad file"),
    make_option("--sample_id", type = "character", default = "complete",
                help = "Fallback sample_id if metadata has no sample_id column"),
    make_option("--subsample", type = "integer", default = NULL,
                help = "Subsample N cells from the complete object"),
    make_option("--seed", type = "integer", default = 42,
                help = "Random seed for subsampling")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rds)) stop("--rds is required")
if (is.null(opt$output)) stop("--output is required")
set.seed(opt$seed)

cat(sprintf("[ingest_seurat_h5ad] Loading %s\n", opt$rds))
sobj <- readRDS(opt$rds)
sobj <- tryCatch(UpdateSeuratObject(sobj), error = function(e) sobj)
sobj <- tryCatch(JoinLayers(sobj), error = function(e) sobj)

if (!is.null(opt$subsample) && ncol(sobj) > opt$subsample) {
    cells <- sample(colnames(sobj), opt$subsample)
    sobj <- subset(sobj, cells = cells)
    cat(sprintf("[ingest_seurat_h5ad] Subsampled to %d cells\n", ncol(sobj)))
}

counts <- tryCatch(
    GetAssayData(sobj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(sobj, assay = "RNA", slot = "counts")
)
if (all(dim(counts) == 0)) {
    counts <- GetAssayData(sobj, assay = "RNA", slot = "counts")
}

obs <- sobj@meta.data
if (!("sample_id" %in% colnames(obs))) {
    obs$sample_id <- opt$sample_id
}
obs[] <- lapply(obs, function(x) {
    if (is.factor(x)) as.character(x) else x
})

adata <- AnnData(
    X = t(counts),
    obs = obs,
    var = data.frame(row.names = rownames(counts))
)

write_h5ad(adata, opt$output)
cat(sprintf("[ingest_seurat_h5ad] Wrote %s (%d cells x %d genes)\n",
            opt$output, nrow(obs), nrow(counts)))
