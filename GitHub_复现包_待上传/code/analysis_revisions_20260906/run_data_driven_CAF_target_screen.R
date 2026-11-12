#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

set.seed(20260906)
root <- normalizePath(".")
infile <- file.path(root, "external_CAF_data/integrated/GSE135337_GSE190888_GSE130001_GSE129845_CAF_Harmony_integrated.rds")
outdir <- file.path(root, "analysis_revisions_20260906/NNMT_data_driven_target_screen")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(infile)
md <- as.data.frame(x[[]])
stopifnot("sample_id" %in% colnames(md), "dataset" %in% colnames(md))
md$sample_id <- as.character(md$sample_id)
md$dataset <- as.character(md$dataset)

# These are predefined CAF-state axes in the existing object; no GDF15 variable is used.
module_cols <- c("ECM_CAF_module1", "Activated_myCAF_module1", "iCAF_module1",
                 "Hypoxia_CAF_module1", "TGFb_CAF_module1", "Vascular_CAF_module1")
module_cols <- intersect(module_cols, colnames(md))
if (length(module_cols) < 3) stop("Too few CAF module columns available")

# Sample-level CAF state table. Cells are aggregated before association testing.
state <- as.data.table(md[, c("sample_id", "dataset", module_cols), drop = FALSE])
state <- state[, lapply(.SD, function(z) mean(as.numeric(z), na.rm = TRUE)),
               by = .(sample_id, dataset), .SDcols = module_cols]
state <- state[complete.cases(state[, ..module_cols])]
state_mat <- scale(as.matrix(state[, ..module_cols]))
pc <- prcomp(state_mat, center = FALSE, scale. = FALSE)
state$CAF_state_PC1 <- as.numeric(pc$x[, 1])

# Orient PC1 so positive values broadly reflect activated/ECM/TGF-beta features.
activation_cols <- intersect(c("ECM_CAF_module1", "Activated_myCAF_module1", "TGFb_CAF_module1"), module_cols)
if (length(activation_cols) > 0 && cor(state$CAF_state_PC1, rowMeans(state[, ..activation_cols]), use = "complete.obs") < 0)
  state$CAF_state_PC1 <- -state$CAF_state_PC1

write.csv(state, file.path(outdir, "sample_level_CAF_state_axes.csv"), row.names = FALSE)

# Use normalized RNA values only for screening/visualization; no integrated assay values.
expr <- GetAssayData(x, assay = "RNA", layer = "data")
common <- intersect(colnames(expr), rownames(md))
expr <- expr[, common, drop = FALSE]
md <- md[common, , drop = FALSE]

# Aggregate normalized expression to the sample level.
sample_ids <- md$sample_id
sample_levels <- unique(sample_ids)
sample_expr <- matrix(NA_real_, nrow = nrow(expr), ncol = length(sample_levels),
                      dimnames = list(rownames(expr), sample_levels))
for (s in sample_levels) {
  ix <- which(sample_ids == s)
  sample_expr[, s] <- Matrix::rowMeans(expr[, ix, drop = FALSE])
}

sample_anno <- unique(md[, c("sample_id", "dataset"), drop = FALSE])
sample_anno <- sample_anno[match(colnames(sample_expr), sample_anno$sample_id), , drop = FALSE]
state <- state[match(colnames(sample_expr), state$sample_id), , drop = FALSE]
keep_samples <- complete.cases(state[, c("CAF_state_PC1", module_cols), drop = FALSE])
sample_expr <- sample_expr[, keep_samples, drop = FALSE]
state <- state[keep_samples, , drop = FALSE]
sample_anno <- sample_anno[keep_samples, , drop = FALSE]

# Exclude genes that are poor target candidates for this screen.
gene <- rownames(sample_expr)
is_housekeeping <- gene %in% c("GAPDH", "ACTB", "B2M", "MALAT1", "RPLP0")
is_ribo_mito <- grepl("^(RPL|RPS|MT-)", gene)
candidate_genes <- gene[!(is_housekeeping | is_ribo_mito)]

# Candidate-level sample associations with independent CAF-state axes.
targets <- c("CAF_state_PC1", module_cols)
res <- rbindlist(lapply(candidate_genes, function(g) {
  y <- as.numeric(sample_expr[g, ])
  rows <- lapply(targets, function(t) {
    z <- as.numeric(state[[t]])
    ok <- is.finite(y) & is.finite(z)
    if (sum(ok) < 6 || sd(y[ok]) == 0 || sd(z[ok]) == 0)
      return(data.table(gene = g, axis = t, n_samples = sum(ok), rho = NA_real_, p_value = NA_real_))
    ct <- suppressWarnings(cor.test(y[ok], z[ok], method = "spearman", exact = FALSE))
    data.table(gene = g, axis = t, n_samples = sum(ok), rho = unname(ct$estimate), p_value = ct$p.value)
  })
  rbindlist(rows)
}), fill = TRUE)
res[, FDR := p.adjust(p_value, method = "BH"), by = axis]

# Summarize stability across axes. This is ranking, not a claim of independent replication.
ranked <- res[, .(
  n_axes_tested = sum(is.finite(rho)),
  mean_abs_rho = mean(abs(rho), na.rm = TRUE),
  max_abs_rho = max(abs(rho), na.rm = TRUE),
  n_positive = sum(rho > 0, na.rm = TRUE),
  n_negative = sum(rho < 0, na.rm = TRUE),
  min_FDR = min(FDR, na.rm = TRUE),
  axes_FDR_lt_0.10 = sum(FDR < 0.10, na.rm = TRUE)
), by = gene]
ranked[!is.finite(min_FDR), min_FDR := NA_real_]
ranked[, direction_consistency := pmax(n_positive, n_negative) / n_axes_tested]
ranked[, data_driven_score := 0.45 * mean_abs_rho +
                              0.25 * max_abs_rho +
                              0.20 * direction_consistency +
                              0.10 * pmin(axes_FDR_lt_0.10 / 2, 1)]
setorder(ranked, -data_driven_score, -mean_abs_rho)

# Track biologically actionable candidates without using them to choose the ranking.
watchlist <- c("NNMT", "CXCL12", "FAP", "MFAP4", "PDGFRB", "PDGFRA", "TGFB1", "TGFBR1",
               "TGFBR2", "VCAN", "FGF7", "POSTN", "TGFBI", "COL11A1", "CCL2", "CXCL14",
               "CXCL1", "CXCL8", "PI16", "CFD", "C7", "SMOC2", "GAS1", "CRISPLD2")
ranked[, watchlist := gene %in% watchlist]
fwrite(ranked, file.path(outdir, "all_gene_CAF_state_target_ranking.csv"))
fwrite(res[gene %in% watchlist], file.path(outdir, "watchlist_axis_associations.csv"))
fwrite(ranked[watchlist == TRUE], file.path(outdir, "watchlist_ranked_candidates.csv"))

# Basic audit for reproducibility and interpretation.
audit <- data.table(
  metric = c("cells", "genes", "samples", "datasets", "state_axes", "NNMT_present", "NNMT_rank"),
  value = c(nrow(md), nrow(expr), nrow(state), length(unique(sample_anno$dataset)),
            length(targets), "NNMT" %in% rownames(sample_expr),
            match("NNMT", ranked$gene))
)
fwrite(audit, file.path(outdir, "screen_audit.csv"))
writeLines(c(
  "Data-driven CAF target screen completed",
  paste("Input cells:", nrow(md)),
  paste("Input genes:", nrow(expr)),
  paste("Sample-level replicates:", nrow(state)),
  paste("Datasets:", length(unique(sample_anno$dataset))),
  paste("State axes:", paste(targets, collapse = ", ")),
  "Primary association unit: sample-level mean normalized RNA expression",
  "GDF15 was not used as a phenotype or ranking variable",
  "Cell-level p-values were not used as the primary evidence"
), file.path(outdir, "RUN_STATUS.txt"))

message("Completed: ", outdir)
