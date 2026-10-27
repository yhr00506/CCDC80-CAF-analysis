suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/external_CAF_data"
raw <- file.path(base, "raw")
out <- file.path(base, "new_cohorts_20260906")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "GSE222315_CAF"), showWarnings = FALSE)
dir.create(file.path(out, "GSE293189_CAF"), showWarnings = FALSE)

marker_sets <- list(
  CAF = c("COL1A1","COL1A2","DCN","LUM","COL3A1","PDGFRA","PDGFRB","C7","CFD","COL6A1","COL6A2","COL6A3","FAP","ACTA2"),
  epithelial = c("EPCAM","KRT8","KRT18","KRT19","KRT7","KRT17","KRT14"),
  immune = c("PTPRC","CD3D","CD3E","MS4A1","CD79A","LYZ","FCER1G"),
  endothelial = c("PECAM1","VWF","KDR","EMCN")
)

score_cells <- function(obj) {
  counts <- GetAssayData(obj, layer = "counts")
  for (nm in names(marker_sets)) {
    g <- intersect(marker_sets[[nm]], rownames(counts))
    obj[[paste0(nm, "_score")]] <- if (length(g) >= 2) Matrix::colMeans(counts[g, , drop = FALSE]) else 0
  }
  md <- obj[[]]
  obj$CAF_candidate <- md$CAF_score > md$epithelial_score &
    md$CAF_score > md$immune_score & md$CAF_score > md$endothelial_score
  obj
}

read_expression <- function(f, dataset, sample, has_name = FALSE) {
  x <- fread(cmd = paste("gzip -cd", shQuote(f)), data.table = FALSE, showProgress = FALSE)
  if (has_name) {
    gene_id <- x[[1]]
    gene_name <- x[[2]]
    gene <- ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
    x <- x[, -(1:2), drop = FALSE]
  } else {
    gene <- x[[1]]
    x <- x[, -1, drop = FALSE]
  }
  rownames(x) <- make.unique(as.character(gene))
  colnames(x) <- paste(dataset, sample, colnames(x), sep = "_")
  mat <- Matrix(as.matrix(x), sparse = TRUE)
  obj <- CreateSeuratObject(mat, project = dataset, min.cells = 0, min.features = 0)
  obj$dataset <- dataset
  obj$sample_id <- sample
  obj$tissue <- ifelse(grepl("NAT|_N_|_N[0-9]|-N", sample, ignore.case = TRUE), "adjacent_normal", "BLCA_tumor")
  obj
}

process_one <- function(f, dataset, sample, has_name, dest) {
  message("Processing ", dataset, " / ", sample)
  z <- score_cells(read_expression(f, dataset, sample, has_name))
  md <- z[[]]
  n_caf <- sum(md$CAF_candidate)
  n_tumor <- sum(md$tissue == "BLCA_tumor")
  if (n_caf > 0 && n_tumor > 0) {
    caf <- subset(z, subset = CAF_candidate & tissue == "BLCA_tumor")
    g <- intersect("GDF15", rownames(caf))
    z_counts <- GetAssayData(z, layer = "counts")
    gdf15_values <- if (length(g)) as.numeric(z_counts[g, colnames(caf), drop = FALSE]) else rep(0, ncol(caf))
    caf$GDF15_score <- gdf15_values
    caf$GDF15_group <- ifelse(caf$GDF15_score > 0, "GDF15_detected", "GDF15_not_detected")
    saveRDS(caf, file.path(dest, paste0(sample, "_CAF.rds")), compress = "xz")
    fwrite(as.data.table(caf[[]], keep.rownames = "cell_id"), file.path(dest, paste0(sample, "_CAF_metadata.csv")))
  }
  data.table(dataset = dataset, sample = sample, n_total = ncol(z), n_tumor = n_tumor,
             n_CAF_candidate = n_caf, n_CAF_tumor = sum(md$CAF_candidate & md$tissue == "BLCA_tumor"),
             GDF15_detected_CAF = if (n_caf > 0 && n_tumor > 0) sum(caf$GDF15_group == "GDF15_detected") else 0)
}

run_cohort <- function(dataset, pattern, has_name, dest, tumor_only = FALSE) {
  fs <- list.files(file.path(raw, dataset), pattern = pattern, full.names = TRUE)
  if (tumor_only) fs <- fs[grepl("T", basename(fs), ignore.case = TRUE)]
  if (!length(fs)) stop("No files found for ", dataset)
  audits <- lapply(fs, function(f) {
    sample <- sub("_.*$", "", basename(f))
    if (dataset == "GSE293189") sample <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(f)))
    process_one(f, dataset, sample, has_name, dest)
  })
  audit <- rbindlist(audits, fill = TRUE)
  fwrite(audit, file.path(dest, "sample_CAF_audit.csv"))
  writeLines(c(paste("dataset", dataset), paste("files_processed", nrow(audit)),
               paste("tumor_samples_with_CAF", sum(audit$n_CAF_tumor > 0)),
               paste("total_CAF", sum(audit$n_CAF_tumor)),
               paste("total_GDF15_detected_CAF", sum(audit$GDF15_detected_CAF))),
             file.path(dest, "RUN_STATUS.txt"))
  audit
}

audit_222315 <- run_cohort("GSE222315", "_expression\\.txt\\.gz$", TRUE,
                           file.path(out, "GSE222315_CAF"), tumor_only = FALSE)
audit_293189 <- run_cohort("GSE293189", "_dge\\.txt\\.gz$", FALSE,
                           file.path(out, "GSE293189_CAF"), tumor_only = TRUE)
fwrite(rbindlist(list(audit_222315, audit_293189), fill = TRUE), file.path(out, "new_cohort_CAF_audit.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
