library(Seurat)
base_dir <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
input_file <- file.path(base_dir, "external_CAF_data/integrated/GSE135337_GSE190888_GSE130001_GSE129845_CAF_Harmony_integrated.rds")
out_dir <- file.path(base_dir, "final_analysis_20260907/Figure1_CAF_atlas")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(out_dir, "GSE135337_GSE190888_GSE130001_GSE129845_CAF_global_annotated.rds")
obj <- readRDS(input_file)
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
genes <- rownames(expr)
marker_sets <- list(ECM_CAF = c("COL1A1", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "FBN1", "MFAP5", "DCN", "LUM", "SPARC"), TGFb_CAF = c("ACTA2", "TAGLN", "TGFB1", "PDGFRB", "POSTN", "CTHRC1", "THBS2", "MYL9"), iCAF = c("CXCL12", "IL6", "CFD", "C7", "CCL2", "FBLN1", "PI16", "DPT"), Vascular_CAF = c("RGS5", "MCAM", "CSPG4", "PECAM1", "VWF", "EMCN", "KDR"), Hypoxia_CAF = c("VEGFA", "CA9", "HILPDA", "ADM", "LDHA", "SLC2A1", "ENO1"), apCAF = c("HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "HLA-DQA1", "HLA-DQB1"))
score_mat <- matrix(NA_real_, nrow = ncol(expr), ncol = length(marker_sets), dimnames = list(colnames(expr), names(marker_sets)))
for (nm in names(marker_sets)) { g <- intersect(marker_sets[[nm]], genes); if (!length(g)) next; x <- as.matrix(expr[g, , drop = FALSE]); s <- apply(x, 1, sd); x <- sweep(x, 1, rowMeans(x), "-"); x <- sweep(x, 1, ifelse(s > 0, s, 1), "/"); score_mat[, nm] <- colMeans(x, na.rm = TRUE) }
best <- apply(score_mat, 1, function(x) names(which.max(x)))
sorted_scores <- t(apply(score_mat, 1, sort, decreasing = TRUE))
margin <- sorted_scores[, 1] - sorted_scores[, 2]
confidence <- ifelse(margin >= 0.25, "high", ifelse(margin >= 0.10, "moderate", "low"))
obj$CAF_subtype_global <- factor(best, levels = names(marker_sets))
obj$CAF_annotation_confidence <- confidence
obj$CAF_annotation_margin <- as.numeric(margin)
obj$CAF_annotation_source <- "global_marker_score_reannotation_20260907"
for (nm in names(marker_sets)) obj[[paste0("CAF_score_global_", nm)]] <- score_mat[, nm]
old <- obj[[]]
audit <- data.frame(total_cells = nrow(old), old_subtype_annotated = sum(!is.na(old$CAF_subtype_clean)), new_global_annotated = sum(!is.na(obj$CAF_subtype_global)), high_confidence = sum(confidence == "high"), moderate_confidence = sum(confidence == "moderate"), low_confidence = sum(confidence == "low"))
write.csv(audit, file.path(out_dir, "global_annotation_audit.csv"), row.names = FALSE)
write.csv(as.data.frame(table(dataset = old$dataset, CAF_subtype_global = obj$CAF_subtype_global)), file.path(out_dir, "global_annotation_by_cohort.csv"), row.names = FALSE)
write.csv(as.data.frame(table(old_subtype = old$CAF_subtype_clean, new_global_subtype = obj$CAF_subtype_global)), file.path(out_dir, "global_annotation_old_new_confusion.csv"), row.names = FALSE)
write.csv(data.frame(cell_id = colnames(obj), score_mat, global_subtype = obj$CAF_subtype_global, confidence = confidence, margin = margin), file.path(out_dir, "global_annotation_cell_scores.csv"), row.names = FALSE)
saveRDS(obj, output_file)
writeLines(c("Global CAF subtype re-annotation completed.", "Original RDS was preserved.", paste0("Cells=", nrow(obj), "; old annotated=", audit$old_subtype_annotated, "; global annotated=", audit$new_global_annotated), "Subtypes were assigned by standardized marker-module scores across the integrated object."), file.path(out_dir, "GLOBAL_ANNOTATION_STATUS.txt"))
print(audit)
print(table(obj$dataset, obj$CAF_subtype_global))
