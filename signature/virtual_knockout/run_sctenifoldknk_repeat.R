#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[[1]] else "CCDC80"
seed <- if (length(args) >= 2) as.integer(args[[2]]) else 20260906
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
.libPaths(c(file.path(root, "virtual_knockout/Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(scTenifoldKnk)
})

set.seed(seed)
obj <- readRDS(file.path(root, "external_CAF_data/integrated/GSE135337_GSE190888_GSE130001_GSE129845_CAF_Harmony_integrated.rds"))
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stopifnot(target %in% rownames(counts))
cells <- sample(colnames(counts), min(500, ncol(counts)))
vf <- intersect(VariableFeatures(obj), rownames(counts))
genes <- unique(c(target, vf[seq_len(min(999, length(vf)))]))
X <- as.matrix(counts[genes, cells, drop = FALSE])
storage.mode(X) <- "numeric"

outdir <- file.path(root, "analysis_revisions_20260906", "virtual_knockout_repeat", target, paste0("seed_", seed))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
res <- scTenifoldKnk(countMatrix = X, gKO = target,
                     qc_minLSize = 0, nc_nNet = 5,
                     nc_nCells = min(100, ncol(X)), nc_nComp = 3,
                     td_K = 3, nCores = 1, verbose = TRUE)
saveRDS(res, file.path(outdir, paste0(target, "_result.rds")))
write.csv(res$diffRegulation, file.path(outdir, paste0(target, "_diffRegulation.csv")), row.names = FALSE)
writeLines(c(paste0("target=", target), paste0("seed=", seed),
             paste0("n_genes_input=", nrow(X)), paste0("n_cells_input=", ncol(X)),
             "method=scTenifoldKnk_2.1.0", "repeat_analysis=TRUE",
             "observational scRNA-seq; virtual perturbation is not causal proof"),
           file.path(outdir, "RUN_STATUS.txt"))
