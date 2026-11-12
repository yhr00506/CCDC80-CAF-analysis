suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(harmony)
  library(data.table)
})

out <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/external_CAF_data"
dir.create(file.path(out, "integrated"), recursive = TRUE, showWarnings = FALSE)

read_external <- function(f, dataset, sample) {
  x <- data.table::fread(cmd = paste("zcat", shQuote(f)), data.table = FALSE, showProgress = FALSE)
  genes <- x[[1]]
  x <- x[, -1, drop = FALSE]
  if (length(genes) != nrow(x)) stop("gene-name/row dimension mismatch in ", f)
  rownames(x) <- make.unique(as.character(genes))
  colnames(x) <- paste(dataset, sample, colnames(x), sep = "_")
  x <- Matrix(as.matrix(x), sparse = TRUE)
  obj <- CreateSeuratObject(x, project = dataset, min.cells = 0, min.features = 0)
  obj$dataset <- dataset
  obj$sample_id <- sample
  obj$tissue <- "BLCA_tumor"
  obj
}

marker_sets <- list(
  CAF = c("COL1A1","COL1A2","DCN","LUM","COL3A1","PDGFRA","PDGFRB","C7","CFD","COL6A1","COL6A2","COL6A3","FAP","ACTA2"),
  epithelial = c("EPCAM","KRT8","KRT18","KRT19","KRT7","KRT17","KRT14"),
  immune = c("PTPRC","CD3D","CD3E","MS4A1","CD79A","LYZ","FCER1G"),
  endothelial = c("PECAM1","VWF","KDR","EMCN")
)

score_cells <- function(obj) {
  for (nm in names(marker_sets)) {
    g <- intersect(marker_sets[[nm]], rownames(obj))
    obj[[paste0(nm, "_score")]] <- if (length(g) >= 2) Matrix::colMeans(GetAssayData(obj, layer = "counts")[g, , drop = FALSE]) else 0
  }
  md <- obj[[]]
  keep <- md$CAF_score > md$epithelial_score & md$CAF_score > md$immune_score & md$CAF_score > md$endothelial_score
  obj$CAF_candidate <- keep
  obj
}

files <- list.files(file.path(out, "processed", "GSE190888"), pattern = "expression_matrix\\.csv\\.gz$", full.names = TRUE)
ext_caf <- lapply(files, function(f) {
  bn <- basename(f); gsm <- sub("_.*", "", bn)
  z <- read_external(f, "GSE190888", gsm)
  z <- score_cells(z)
  z <- subset(z, subset = CAF_candidate)
  gc()
  z
})
names(ext_caf) <- basename(files)

old_path <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/01_scRNA_CAF_discovery/10_GDF15_high_low_CAF_FIXED/CAF_GDF15_high_low.rds"
old <- readRDS(old_path)
old <- JoinLayers(old)
old$dataset <- "GSE135337"
old$sample_id <- as.character(old$sample)
old$tissue <- as.character(old$tissue)
if ("tissue" %in% colnames(old[[]])) old <- subset(old, subset = tissue == "Tumor")
old <- score_cells(old)

audit <- rbindlist(c(list(data.table(dataset="GSE135337", sample="GSE135337_CAF_object", n_total=ncol(old), n_CAF=sum(old$CAF_candidate))),
                     lapply(ext_caf, function(z) data.table(dataset=unique(z$dataset), sample=unique(z$sample_id), n_total=NA_integer_, n_CAF=ncol(z)))), fill=TRUE)
fwrite(audit, file.path(out, "integrated", "CAF_candidate_audit.csv"))

old_caf <- subset(old, subset = CAF_candidate)
all_caf <- Reduce(function(a,b) merge(a,b), c(list(old_caf), ext_caf))
all_caf <- JoinLayers(all_caf)
all_caf <- NormalizeData(all_caf, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
all_caf <- FindVariableFeatures(all_caf, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
all_caf <- ScaleData(all_caf, features = VariableFeatures(all_caf), verbose = FALSE)
all_caf <- RunPCA(all_caf, features = VariableFeatures(all_caf), npcs = 30, verbose = FALSE)
all_caf <- RunHarmony(all_caf, group.by.vars = "dataset", reduction.use = "pca", dims.use = 1:30, assay.use = "RNA", plot_convergence = FALSE)
all_caf <- FindNeighbors(all_caf, reduction = "harmony", dims = 1:30, verbose = FALSE)
all_caf <- FindClusters(all_caf, resolution = 0.4, verbose = FALSE)
all_caf <- RunUMAP(all_caf, reduction = "harmony", dims = 1:30, verbose = FALSE)

gdf15 <- intersect("GDF15", rownames(all_caf))
if (length(gdf15)) all_caf$GDF15_score <- as.numeric(GetAssayData(all_caf, layer="data")[gdf15, ]) else all_caf$GDF15_score <- 0
all_caf$GDF15_group <- ifelse(all_caf$GDF15_score > 0, "GDF15_detected", "GDF15_not_detected")

fwrite(as.data.table(all_caf[[]], keep.rownames = "cell_id"), file.path(out, "integrated", "integrated_CAF_metadata.csv"))
saveRDS(all_caf, file.path(out, "integrated", "GSE135337_GSE190888_CAF_Harmony_integrated.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(out, "integrated", "sessionInfo.txt"))
writeLines(c("Integration completed", paste("cells", ncol(all_caf)), paste("datasets", paste(unique(all_caf$dataset), collapse=",")), paste("GDF15_detected", sum(all_caf$GDF15_group=="GDF15_detected"))), file.path(out, "integrated", "RUN_STATUS.txt"))
