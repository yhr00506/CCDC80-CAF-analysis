#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

analysis_dir <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/final_analysis_20260907"
upstream_dir <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/analysis_revisions_20260906"
out_dir <- file.path(analysis_dir, "Figure4_CAF_myeloid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

patient_file <- file.path(upstream_dir, "CCDC80_meta_GSE293189_GSE135337", "combined_patient_level_CCDC80_TME.csv")
corr_file <- file.path(upstream_dir, "CCDC80_meta_GSE293189_GSE135337", "within_cohort_CCDC80_TME_correlations.csv")
meta_file <- file.path(upstream_dir, "CCDC80_meta_GSE293189_GSE135337", "fisher_effect_size_meta_analysis.csv")
lr_file <- file.path(upstream_dir, "GSE135337_CAF_myeloid_LR", "GSE135337_LR_ranked_pairs.csv")

patient <- read.csv(patient_file, check.names = FALSE, stringsAsFactors = FALSE)
corr <- read.csv(corr_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
lr <- read.csv(lr_file, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(nrow(patient) == 23, all(c("CCDC80_CAF_z", "CAF_Myeloid_score", "cohort") %in% names(patient)))

patient$cohort_label <- factor(patient$cohort,
  levels = c("GSE293189", "GSE135337"),
  labels = c("A  GSE293189\n(n = 16)", "B  GSE135337\n(n = 7)"))

rho_rows <- corr[corr$variable == "CAF_Myeloid_score", c("cohort", "n", "rho", "p", "FDR")]
rho_rows$label <- ifelse(rho_rows$cohort == "GSE293189", "GSE293189", "GSE135337")
rho_rows$lo <- NA_real_
rho_rows$hi <- NA_real_
for (i in seq_len(nrow(rho_rows))) {
  z <- atanh(rho_rows$rho[i])
  se <- 1 / sqrt(rho_rows$n[i] - 3)
  rho_rows$lo[i] <- tanh(z - 1.96 * se)
  rho_rows$hi[i] <- tanh(z + 1.96 * se)
}
meta_row <- meta[meta$variable == "CAF_Myeloid_score", ]

fmt_p <- function(x) {
  ifelse(x < 0.001, "P < 0.001", paste0("P = ", formatC(x, format = "f", digits = 3)))
}

palette <- c("A  GSE293189\n(n = 16)" = "#356A8A", "B  GSE135337\n(n = 7)" = "#C98245")
theme_fig <- theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9, hjust = 0, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 7, colour = "#444444", hjust = 0, margin = margin(b = 5)),
    axis.title = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 7, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    panel.grid.major = element_line(linewidth = 0.18, colour = "#E9E9E9"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 6, 5.5, 6)
  )

# (a) Patient-level association, shown separately by cohort.
annot <- merge(
  data.frame(cohort = c("GSE293189", "GSE135337"), stringsAsFactors = FALSE),
  rho_rows[, c("cohort", "rho", "p")], by = "cohort", all.x = TRUE
)
annot$cohort_label <- factor(annot$cohort,
  levels = c("GSE293189", "GSE135337"),
  labels = c("A  GSE293189\n(n = 16)", "B  GSE135337\n(n = 7)"))
annot$x <- -Inf
annot$y <- Inf
annot$txt <- paste0("Spearman ρ = ", sprintf("%.2f", annot$rho), "\n", fmt_p(annot$p))

p_a <- ggplot(patient, aes(x = CCDC80_CAF_z, y = CAF_Myeloid_score, colour = cohort_label)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_point(size = 2.2, alpha = 0.92, stroke = 0.25) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.65, linetype = "solid") +
  facet_wrap(~ cohort_label, nrow = 1) +
  scale_colour_manual(values = palette) +
  labs(title = "Patient-level association is cohort-dependent",
       subtitle = "Each point represents one patient; lines are visual linear fits",
       x = "CCDC80-associated CAF score (z)", y = "CAF–myeloid score") +
  theme_fig +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
        panel.spacing = unit(1.1, "lines"))

# Add annotations separately to avoid inheriting facet data.

# (b) Effect-size summary across the two cohorts and prespecified TME readouts.
forest <- corr[corr$variable %in% c("Myeloid_fraction", "Myeloid_score", "CAF_Myeloid_score", "Exhaustion_score", "Cytotoxicity_score"), ]
forest$metric <- factor(forest$variable,
  levels = c("CAF_Myeloid_score", "Myeloid_score", "Myeloid_fraction", "Exhaustion_score", "Cytotoxicity_score"),
  labels = c("CAF–myeloid score", "Myeloid score", "Myeloid fraction", "T-cell exhaustion", "Cytotoxicity"))
forest$lo <- NA_real_; forest$hi <- NA_real_
for (i in seq_len(nrow(forest))) {
  z <- atanh(forest$rho[i]); se <- 1 / sqrt(forest$n[i] - 3)
  forest$lo[i] <- tanh(z - 1.96 * se); forest$hi[i] <- tanh(z + 1.96 * se)
}
forest$cohort_label <- factor(forest$cohort, levels = c("GSE293189", "GSE135337"))
forest$y <- as.numeric(forest$metric) + ifelse(forest$cohort == "GSE293189", 0.13, -0.13)

meta_plot <- meta[meta$variable %in% c("Myeloid_fraction", "Myeloid_score", "CAF_Myeloid_score", "Exhaustion_score", "Cytotoxicity_score"), ]
meta_labels <- c(Myeloid_fraction = "Myeloid fraction", Myeloid_score = "Myeloid score", CAF_Myeloid_score = "CAF–myeloid score", Exhaustion_score = "T-cell exhaustion", Cytotoxicity_score = "Cytotoxicity")
meta_plot$metric <- factor(unname(meta_labels[meta_plot$variable]), levels = levels(forest$metric))
meta_plot$y <- as.numeric(meta_plot$metric)

p_b <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#777777") +
  geom_segment(data = forest, aes(x = lo, xend = hi, y = y, yend = y, colour = cohort_label), linewidth = 0.55) +
  geom_point(data = forest, aes(x = rho, y = y, colour = cohort_label), size = 2.0) +
  geom_point(data = meta_plot, aes(x = rho_fisher, y = y), shape = 23, size = 2.7,
             fill = "#252525", colour = "#252525") +
  scale_colour_manual(values = c(GSE293189 = palette[[1]], GSE135337 = palette[[2]])) +
  scale_y_continuous(breaks = seq_along(levels(forest$metric)), labels = levels(forest$metric),
                     limits = c(0.4, 5.6), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-1.05, 1.05), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(title = "(C) Effect-size summary across cohorts",
       subtitle = paste0("Diamonds: Fisher-combined estimate; n = 16 and 7 patients per cohort; pooled CAF–myeloid FDR = ",
                          formatC(meta_row$FDR, format = "f", digits = 4)),
       x = "Spearman ρ (95% CI)", y = NULL) +
  theme_fig +
  theme(axis.text.y = element_text(size = 7), plot.subtitle = element_text(size = 6.5))

# (c) LR ranking is explicitly framed as hypothesis generation.
lr_show <- head(lr[is.finite(lr$mean_LR) & lr$mean_LR > 0, ][order(lr[is.finite(lr$mean_LR) & lr$mean_LR > 0, ]$mean_LR, decreasing = TRUE), ], 8)
lr_show$pair <- factor(lr_show$pair, levels = rev(lr_show$pair))
p_c <- ggplot(lr_show, aes(x = mean_LR, y = pair)) +
  geom_segment(aes(x = min(mean_LR) * 0.8, xend = mean_LR, y = pair, yend = pair), linewidth = 0.65, colour = "#6B8797") +
  geom_point(size = 2.1, colour = "#356A8A") +
  scale_x_continuous(trans = "log10", breaks = c(1e-4, 1e-3, 1e-2), labels = c("10⁻⁴", "10⁻³", "10⁻²")) +
  labs(title = "(F) Candidate CAF–myeloid interactions",
       subtitle = "Prioritization score in GSE135337 (n = 7); not causal communication",
       x = "Mean ligand × receptor score (log scale)", y = NULL) +
  theme_fig + theme(axis.text.y = element_text(size = 6.6), plot.subtitle = element_text(size = 6.5))

p_d <- ggplot(patient, aes(x = CCDC80_CAF_z, y = Myeloid_score, colour = cohort_label)) + geom_point(size = 2.0, alpha = 0.92) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) + scale_colour_manual(values = palette) + labs(title = "(D) Myeloid score", x = "CCDC80-associated CAF score (z)", y = "Myeloid score") + theme_fig
p_e <- ggplot(patient, aes(x = CCDC80_CAF_z, y = Myeloid_fraction, colour = cohort_label)) + geom_point(size = 2.0, alpha = 0.92) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) + scale_colour_manual(values = palette) + scale_y_continuous(labels = scales::percent) + labs(title = "(E) Myeloid fraction", x = "CCDC80-associated CAF score (z)", y = "Myeloid fraction") + theme_fig
fig <- p_a / (p_b | p_d) / (p_e | p_c) + plot_layout(nrow = 3, ncol = 1, heights = c(1.15, 1, 1), guides = "keep")

pdf_file <- file.path(out_dir, "Figure4_CAF_myeloid_revised.pdf")
svg_file <- file.path(out_dir, "Figure4_CAF_myeloid_revised.svg")
tif_file <- file.path(out_dir, "Figure4_CAF_myeloid_revised.tif")
ggsave(pdf_file, fig, width = 180, height = 260, units = "mm", device = cairo_pdf, bg = "white")
ggsave(svg_file, fig, width = 180, height = 260, units = "mm", device = svglite, bg = "white")
ggsave(tif_file, fig, width = 180, height = 260, units = "mm", dpi = 600, compression = "lzw", device = "tiff", bg = "white")

# Traceability bundle.
write.csv(patient, file.path(out_dir, "Figure4_source_patient_level.csv"), row.names = FALSE, quote = TRUE)
write.csv(forest, file.path(out_dir, "Figure4_source_cohort_correlations.csv"), row.names = FALSE, quote = TRUE)
write.csv(meta_plot, file.path(out_dir, "Figure4_source_meta_effects.csv"), row.names = FALSE, quote = TRUE)
write.csv(lr_show, file.path(out_dir, "Figure4_source_LR_top8.csv"), row.names = FALSE, quote = TRUE)
writeLines(c(
  "Figure 4 revised: patient-level CCDC80-associated CAF / myeloid association.",
  "Panel a: CCDC80_CAF_z versus CAF_Myeloid_score, patient-level, stratified by cohort.",
  "Panel b: within-cohort Spearman correlations with Fisher-combined estimates; confidence intervals are Fisher-z approximations.",
  "Panel c: ligand-receptor prioritization in GSE135337; mean ligand x receptor score, hypothesis-generating only.",
  paste0("Source patient table: ", patient_file),
  paste0("Source correlation table: ", corr_file),
  paste0("Source meta-analysis table: ", meta_file),
  paste0("Source LR table: ", lr_file)
), file.path(out_dir, "Figure4_REVISED_README.txt"))

message("Wrote: ", pdf_file)
message("Wrote: ", svg_file)
message("Wrote: ", tif_file)
