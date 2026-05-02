#!/usr/bin/env Rscript

# Regenerate test/data/test.rds — Seurat's built-in pbmc_small (~80 cells x
# 230 genes), ships with the Seurat package and is BSD-licensed for
# redistribution. Tiny enough to commit, sufficient to exercise every
# stage of the pipeline end-to-end.

suppressPackageStartupMessages(library(Seurat))

data(pbmc_small)
pbmc_small$sample_id <- "pbmc_small"
pbmc_small <- UpdateSeuratObject(pbmc_small)

dir.create(file.path(dirname(sys.frame(1)$ofile %||% "."), "data"),
           showWarnings = FALSE)

out <- file.path("test", "data", "test.rds")
saveRDS(pbmc_small, out)
cat(sprintf("Wrote %s (%d cells x %d genes)\n",
            out, ncol(pbmc_small), nrow(pbmc_small)))
