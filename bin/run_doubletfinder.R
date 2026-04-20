#!/usr/bin/env Rscript

# =============================================================================
#  DoubletFinder — doublet detection & removal
# =============================================================================
#  Reads an h5ad, converts to Seurat, runs DoubletFinder, removes predicted
#  doublets, and writes the singlets back as h5ad.

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(DoubletFinder)
    library(anndata)
})

# ---- CLI arguments ----
option_list <- list(
    make_option("--input", type = "character", help = "Input h5ad file"),
    make_option("--output", type = "character", help = "Output h5ad (singlets only)"),
    make_option("--stats_output", type = "character", help = "Output text file with doublet stats"),
    make_option("--sample_id", type = "character", help = "Sample identifier"),
    make_option("--pcs", type = "integer", default = 30, help = "Number of PCs"),
    make_option("--doublet_rate", type = "double", default = NULL,
                help = "Expected doublet rate (NULL = auto from 10x guidelines)"),
    make_option("--sct", action = "store_true", default = FALSE,
                help = "Use SCTransform normalization"),
    make_option("--seed", type = "integer", default = 42,
                help = "Random seed for paramSweep / PCA / clustering / doubletFinder")
)

opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)
cat(sprintf("[run_doubletfinder] Processing sample: %s (seed=%d)\n",
            opt$sample_id, opt$seed))

# ---- Load h5ad and convert to Seurat ----
adata <- read_h5ad(opt$input)
counts <- t(adata$X)  # genes x cells
colnames(counts) <- rownames(adata$obs)
rownames(counts) <- rownames(adata$var)

sobj <- CreateSeuratObject(counts = counts, project = opt$sample_id)

# Copy over obs metadata
for (col in colnames(adata$obs)) {
    sobj@meta.data[[col]] <- adata$obs[[col]]
}

total_cells <- ncol(sobj)

# ---- Preprocessing for DoubletFinder ----
if (opt$sct) {
    sobj <- SCTransform(sobj, verbose = FALSE)
} else {
    sobj <- NormalizeData(sobj, verbose = FALSE)
    sobj <- FindVariableFeatures(sobj, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)
    sobj <- ScaleData(sobj, verbose = FALSE)
}

n_pcs <- min(opt$pcs, ncol(sobj) - 1)
set.seed(opt$seed)
sobj <- RunPCA(sobj, npcs = n_pcs, verbose = FALSE, seed.use = opt$seed)
set.seed(opt$seed)
sobj <- FindNeighbors(sobj, dims = 1:n_pcs, verbose = FALSE)
sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE, random.seed = opt$seed)

# ---- Estimate doublet rate ----
if (is.null(opt$doublet_rate)) {
    # 10x Genomics guideline: ~0.8% per 1000 cells loaded
    doublet_rate <- 0.008 * (total_cells / 1000)
    # Cap at reasonable bounds
    doublet_rate <- min(max(doublet_rate, 0.01), 0.25)
    cat(sprintf("[run_doubletfinder] Auto-estimated doublet rate: %.3f (%d cells)\n",
                doublet_rate, total_cells))
} else {
    doublet_rate <- opt$doublet_rate
}

nExp_poi <- round(doublet_rate * total_cells)

# ---- Run DoubletFinder ----
# Parameter sweep for pK
set.seed(opt$seed)
sweep.res <- paramSweep(sobj, PCs = 1:n_pcs, sct = opt$sct)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Select optimal pK
optimal_pk <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]
))

cat(sprintf("[run_doubletfinder] Optimal pK: %.3f\n", optimal_pk))

# Run DoubletFinder
set.seed(opt$seed)
sobj <- doubletFinder(sobj,
                       PCs = 1:n_pcs,
                       pN = 0.25,
                       pK = optimal_pk,
                       nExp = nExp_poi,
                       sct = opt$sct)

# Get the classification column (name varies by run)
df_col <- grep("^DF.classifications", colnames(sobj@meta.data), value = TRUE)[1]
doublet_class <- sobj@meta.data[[df_col]]

doublets_detected <- sum(doublet_class == "Doublet")
singlets_retained <- sum(doublet_class == "Singlet")

cat(sprintf("[run_doubletfinder] %d doublets detected, %d singlets retained\n",
            doublets_detected, singlets_retained))

# ---- Filter to singlets ----
sobj_singlets <- subset(sobj, cells = colnames(sobj)[doublet_class == "Singlet"])

# ---- Convert back to AnnData and write ----
counts_mat <- tryCatch(
    GetAssayData(sobj_singlets, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(sobj_singlets, assay = "RNA", slot = "counts")
)
counts_out <- t(counts_mat)
obs_out <- sobj_singlets@meta.data
# Clean up DoubletFinder columns from obs
obs_out <- obs_out[, !grepl("^(DF\\.|pANN)", colnames(obs_out))]

adata_out <- AnnData(
    X = counts_out,
    obs = obs_out,
    var = data.frame(row.names = rownames(counts_mat))
)

write_h5ad(adata_out, opt$output)

# ---- Write stats (plain text) ----
stats_lines <- c(
    "=== DoubletFinder stats ===",
    sprintf("sample_id            : %s", opt$sample_id),
    sprintf("total_cells          : %d", total_cells),
    sprintf("doublets_detected    : %d", doublets_detected),
    sprintf("doublet_rate         : %.4f", doublets_detected / total_cells),
    sprintf("singlets_retained    : %d", singlets_retained),
    sprintf("optimal_pk           : %.3f", optimal_pk),
    sprintf("expected_doublet_rate: %.4f", doublet_rate)
)
writeLines(stats_lines, opt$stats_output)
cat(sprintf("[run_doubletfinder] Written: %s, %s\n", opt$output, opt$stats_output))
