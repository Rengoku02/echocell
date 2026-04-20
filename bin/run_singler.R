#!/usr/bin/env Rscript

# =============================================================================
#  SingleR cell type annotation
# =============================================================================
#  Reads an h5ad, runs SingleR with a specified reference dataset,
#  adds annotation labels to .obs, and writes back to h5ad.

suppressPackageStartupMessages({
    library(optparse)
    library(SingleR)
    library(celldex)
    library(anndata)
    library(SingleCellExperiment)
    library(scuttle)
})

# ---- CLI arguments ----
option_list <- list(
    make_option("--input", type = "character", help = "Input h5ad file"),
    make_option("--output", type = "character", help = "Output h5ad file"),
    make_option("--ref", type = "character", default = "HumanPrimaryCellAtlasData",
                help = "SingleR reference dataset name"),
    make_option("--labels", type = "character", default = "label.main",
                help = "Label column from reference (label.main or label.fine)"),
    make_option("--fine_tune", action = "store_true", default = FALSE,
                help = "Enable fine-tuning step"),
    make_option("--seed", type = "integer", default = 42,
                help = "Random seed for SingleR")
)

opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)
cat(sprintf("[run_singler] Loading input: %s (seed=%d)\n", opt$input, opt$seed))

# ---- Load h5ad ----
adata <- read_h5ad(opt$input)

# Create a SingleCellExperiment from the full log-normalized expression matrix.
# The Python normalization stage intentionally keeps adata$X full-gene and
# unscaled so SingleR receives valid logcounts.
logcounts <- t(adata$X)  # genes x cells
colnames(logcounts) <- rownames(adata$obs)
rownames(logcounts) <- rownames(adata$var)

sce <- SingleCellExperiment(list(logcounts = logcounts))

# ---- Load reference ----
cat(sprintf("[run_singler] Loading reference: %s\n", opt$ref))

ref <- switch(opt$ref,
    "HumanPrimaryCellAtlasData" = HumanPrimaryCellAtlasData(),
    "BlueprintEncodeData"       = BlueprintEncodeData(),
    "MonacoImmuneData"          = MonacoImmuneData(),
    "DatabaseImmuneCellExpressionData" = DatabaseImmuneCellExpressionData(),
    "NovershternHematopoieticData" = NovershternHematopoieticData(),
    "MouseRNAseqData"           = MouseRNAseqData(),
    "ImmGenData"                = ImmGenData(),
    {
        cat(sprintf("[run_singler] Unknown reference '%s', defaulting to HumanPrimaryCellAtlasData\n", opt$ref))
        HumanPrimaryCellAtlasData()
    }
)

celldex_version <- as.character(packageVersion("celldex"))
singler_version <- as.character(packageVersion("SingleR"))
cat(sprintf("[run_singler] celldex=%s SingleR=%s ref=%s dim=%dx%d\n",
            celldex_version, singler_version, opt$ref, nrow(ref), ncol(ref)))

# ---- Run SingleR ----
cat("[run_singler] Running SingleR annotation...\n")

labels_col <- opt$labels
if (!(labels_col %in% colnames(colData(ref)))) {
    stop(sprintf(
        "Reference '%s' does not contain label column '%s'. Available columns: %s",
        opt$ref,
        labels_col,
        paste(colnames(colData(ref)), collapse = ", ")
    ))
}
ref_labels <- ref[[labels_col]]

set.seed(opt$seed)
predictions <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels,
    fine.tune = opt$fine_tune
)

cat(sprintf("[run_singler] Annotated %d cells\n", nrow(predictions)))

# ---- Add labels to AnnData .obs ----
adata$obs[["singler_labels"]] <- predictions$labels
adata$obs[["singler_pruned_labels"]] <- predictions$pruned.labels
adata$obs[["singler_scores_max"]] <- apply(predictions$scores, 1, max)

# ---- Record provenance in .uns so it survives into final.h5ad ----
adata$uns[["singler_ref"]]         <- opt$ref
adata$uns[["singler_ref_version"]] <- celldex_version
adata$uns[["singler_version"]]     <- singler_version
adata$uns[["celldex_version"]]     <- celldex_version
adata$uns[["seed"]]                <- opt$seed

# ---- Write output ----
write_h5ad(adata, opt$output)
cat(sprintf("[run_singler] Written: %s\n", opt$output))

# Print label distribution
label_table <- sort(table(predictions$labels), decreasing = TRUE)
cat("[run_singler] Top cell types:\n")
print(head(label_table, 10))
