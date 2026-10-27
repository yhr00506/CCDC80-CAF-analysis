library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

base_dir <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
object_file <- file.path(base_dir, "external_CAF_data/integrated/GSE135337_GSE190888_GSE130001_GSE129845_CAF_Harmony_integrated.rds")
out_dir <- file.path(base_dir, "final_analysis_20260907/Figure1_CAF_atlas")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
width_mm <- 180; height_mm <- 230; text_size_pt <- 8; point_size <- 0.18; line_width <- 0.35
wong <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#000000")
continuous_map <- "viridis"
theme_fig <- theme_classic(base_family = "Times", base_size = text_size_pt) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), axis.title = element_text(size = text_size_pt, family = "Times"), axis.text = element_text(size = text_size_pt, family = "Times", colour = "black"), legend.title = element_text(size = text_size_pt, family = "Times"), legend.text = element_text(size = text_size_pt, family = "Times"), strip.text = element_text(size = text_size_pt, family = "Times"), panel.grid.major = element_line(colour = "grey90", linewidth = 0.2), panel.grid.minor = element_blank(), plot.margin = margin(3, 4, 3, 4, unit = "mm"))

obj <- readRDS(object_file)
md <- obj[[]]
umap <- as.data.frame(obj@reductions$umap@cell.embeddings[, 1:2, drop = FALSE])
umap$cell_id <- rownames(umap)
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
dat <- md %>% tibble::rownames_to_column("cell_id") %>% left_join(umap, by = "cell_id")
dat$dataset <- factor(dat$dataset, levels = c("GSE135337", "GSE190888", "GSE130001", "GSE129845"))
dat$CAF_subtype_clean <- factor(dat$CAF_subtype_clean, levels = c("ECM_CAF", "TGFb_CAF", "Vascular_CAF", "Hypoxia_CAF"))

pA <- ggplot(dat, aes(UMAP_1, UMAP_2, colour = dataset)) + geom_point(size = point_size, alpha = 0.65, stroke = 0) + scale_colour_manual(values = wong[seq_len(nlevels(dat$dataset))], drop = FALSE) + labs(x = "UMAP 1", y = "UMAP 2", colour = "Cohort") + theme_fig
pB <- ggplot(dat, aes(UMAP_1, UMAP_2, colour = CAF_subtype_clean)) + geom_point(size = point_size, alpha = 0.65, stroke = 0) + scale_colour_manual(values = wong[seq_len(nlevels(dat$CAF_subtype_clean))], drop = FALSE, na.value = "grey80") + labs(x = "UMAP 1", y = "UMAP 2", colour = "CAF subtype") + theme_fig

marker_sets <- list(ECM = c("COL1A1", "COL3A1", "COL5A1"), TGFb = c("TGFB1", "ACTA2", "TAGLN"), Vascular = c("RGS5", "MCAM", "CSPG4"), Hypoxia = c("VEGFA", "CA9", "HILPDA"))
marker_df <- bind_rows(lapply(names(marker_sets), function(g) tibble::tibble(program = g, gene = marker_sets[[g]])))
present <- intersect(unique(marker_df$gene), rownames(obj))
expr <- GetAssayData(obj, assay = "RNA", slot = "data")[present, , drop = FALSE]
expr_df <- as.data.frame(as.matrix(expr)) %>% tibble::rownames_to_column("gene") %>% pivot_longer(-gene, names_to = "cell_id", values_to = "expression")
dot_df <- expr_df %>% left_join(md %>% tibble::rownames_to_column("cell_id") %>% select(cell_id, CAF_subtype_clean), by = "cell_id") %>% group_by(CAF_subtype_clean, gene) %>% summarise(mean_expression = mean(expression), pct_detected = mean(expression > 0) * 100, .groups = "drop") %>% left_join(marker_df, by = "gene") %>% mutate(gene = factor(gene, levels = rev(present)))
write.csv(dot_df, file.path(out_dir, "Figure1C_marker_dotplot_source.csv"), row.names = FALSE)
pC <- ggplot(dot_df, aes(gene, CAF_subtype_clean, size = pct_detected, colour = mean_expression)) + geom_point(alpha = 0.9) + scale_size_area(max_size = 5, name = "% detected") + scale_colour_viridis_c(option = continuous_map, name = "Mean expression") + labs(x = NULL, y = NULL) + theme_fig + theme(axis.text.x = element_text(angle = 45, hjust = 1))

score_cols <- c("ECM_CAF_module1", "Activated_myCAF_module1", "iCAF_module1", "Hypoxia_CAF_module1", "TGFb_CAF_module1", "Vascular_CAF_module1")
score_long <- dat %>% select(CAF_subtype_clean, all_of(score_cols)) %>% pivot_longer(-CAF_subtype_clean, names_to = "program", values_to = "score")
score_long$program <- factor(score_long$program, levels = score_cols, labels = c("ECM", "Activated myCAF", "iCAF", "Hypoxia", "TGFβ", "Vascular"))
write.csv(score_long, file.path(out_dir, "Figure1D_program_scores_source.csv"), row.names = FALSE)
pD <- ggplot(score_long, aes(CAF_subtype_clean, score, fill = CAF_subtype_clean)) + geom_violin(scale = "width", trim = TRUE, colour = "grey25", linewidth = line_width) + geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = line_width) + facet_wrap(~program, ncol = 3, scales = "free_y") + scale_fill_manual(values = wong[seq_len(nlevels(dat$CAF_subtype_clean))], drop = FALSE) + labs(x = NULL, y = "Program score") + theme_fig + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

count_df <- dat %>% count(dataset, CAF_subtype_clean) %>% group_by(dataset) %>% mutate(fraction = n / sum(n)) %>% ungroup()
write.csv(count_df, file.path(out_dir, "Figure1E_cohort_composition_source.csv"), row.names = FALSE)
pE <- ggplot(count_df, aes(dataset, fraction, fill = CAF_subtype_clean)) + geom_col(width = 0.72, colour = "white", linewidth = 0.2) + scale_y_continuous(labels = percent_format(accuracy = 1)) + scale_fill_manual(values = wong[seq_len(nlevels(dat$CAF_subtype_clean))], drop = FALSE) + labs(x = NULL, y = "CAF fraction", fill = "CAF subtype") + theme_fig + theme(axis.text.x = element_text(angle = 35, hjust = 1))

program_cols <- c("ECM_CAF_module1", "Activated_myCAF_module1", "iCAF_module1", "apCAF_module1", "Hypoxia_CAF_module1", "TGFb_CAF_module1", "Vascular_CAF_module1", "CAF_score")
sample_df <- dat %>% group_by(dataset, sample_id) %>% summarise(across(all_of(program_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")
cor_mat <- cor(sample_df[, program_cols], use = "pairwise.complete.obs", method = "spearman")
cor_df <- as.data.frame(as.table(cor_mat)) %>% rename(program_x = Var1, program_y = Var2, rho = Freq)
write.csv(sample_df, file.path(out_dir, "Figure1F_sample_program_source.csv"), row.names = FALSE)
write.csv(cor_df, file.path(out_dir, "Figure1F_program_correlation_source.csv"), row.names = FALSE)
pF <- ggplot(cor_df, aes(program_x, program_y, fill = rho)) + geom_tile(colour = "white", linewidth = 0.25) + geom_text(aes(label = sprintf("%.2f", rho)), size = 2.3, family = "Times") + scale_fill_viridis_c(option = continuous_map, limits = c(-1, 1), name = "Spearman ρ") + labs(x = NULL, y = NULL) + theme_fig + theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig1 <- (pA | pB) / (pC | pD) / (pE | pF) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(family = "Times", size = text_size_pt, face = "bold"))
ggsave(file.path(out_dir, "Figure1_CAF_atlas.pdf"), fig1, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "Figure1_CAF_atlas.png"), fig1, width = width_mm, height = height_mm, units = "mm", dpi = 300)
writeLines(c("Figure 1 CAF atlas generated from the existing integrated Seurat object.", paste0("Cells=", nrow(dat), "; genes=", ncol(obj), "; cohorts=", nlevels(dat$dataset)), "Continuous colors use the viridis family."), file.path(out_dir, "Figure1_STATUS.txt"))
