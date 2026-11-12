#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[[1]] else "CRISPLD2"
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
.libPaths(c(file.path(root, "virtual_knockout/Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(scTenifoldKnk)
})

set.seed(20260906)
obj <- readRDS(file.path(root, "external_CAF_data/integrated/GSE135337_GSE190888_GSE130001_GSE129845_CAF_Harmony_integrated.rds"))
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stopifnot(target %in% rownames(counts))

# Pilot design: 500 CAF cells, 1000 variable genes plus the target.
cells <- sample(colnames(counts), min(500, ncol(counts)))
vf <- VariableFeatures(obj)
vf <- intersect(vf, rownames(counts))
vf <- vf[seq_len(min(999, length(vf)))]
genes <- unique(c(target, vf))
X <- counts[genes, cells, drop = FALSE]
X <- as.matrix(X)
storage.mode(X) <- "numeric"

outdir <- file.path(root, "analysis_revisions_20260906", "virtual_knockout_pilot")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(target = target, n_genes = nrow(X), n_cells = ncol(X),
                     mean_count = mean(X[target, ]), detected_fraction = mean(X[target, ] > 0)),
          file.path(outdir, paste0(target, "_input_audit.csv")), row.names = FALSE)

res <- scTenifoldKnk(countMatrix = X, gKO = target,
                     qc_minLSize = 0, nc_nNet = 5, nc_nCells = min(100, ncol(X)),
                     nc_nComp = 3, td_K = 3, nCores = 1, verbose = TRUE)
saveRDS(res, file.path(outdir, paste0(target, "_scTenifoldKnk_result.rds")))
write.csv(res$diffRegulation, file.path(outdir, paste0(target, "_diffRegulation.csv")), row.names = FALSE)
writeLines(c(paste0("target=", target), paste0("n_genes=", nrow(X)),
             paste0("n_cells=", ncol(X)), "method=scTenifoldKnk_2.1.0",
             "pilot_only=TRUE; observational scRNA-seq in-silico perturbation; not causal proof"),
           file.path(outdir, paste0(target, "_RUN_STATUS.txt")))
