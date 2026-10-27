#!/usr/bin/env Rscript

# Figure 8: quantitative robustness and incremental prediction value
# Continuous colour is not used; discrete series use the Wong palette.
# All tunable values are kept here. No schematic, workflow, or conceptual diagram.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(survival)
  library(riskRegression)
  library(svglite)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggrepel)
})

# ----------------------------- parameter block -----------------------------
text_size_pt <- 8
line_width <- 0.45
data_line_width <- 0.75
point_size <- 1.8
width_mm <- 240
height_mm <- 280
time_points <- c(1095, 1825, 2920) # 3, 5, and 8 years
bootstrap_n <- 1000
optimism_bootstrap_n <- 200 # computationally stable optimism correction; Panel A uses the existing 1,000-run audit
wong <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#000000")
continuous_map <- "viridis" # unused here; if continuous colour is added, use equal-count quantile breaks
series_cols <- c(
  "Stage" = wong[1],
  "Stage + age" = wong[2],
  "Stage + CAF" = wong[3],
  "Stage + CAF + ECM" = wong[4],
  "Full CCDC80" = wong[5]
)

root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
base <- file.path(root, "analysis_revisions_20260906")
out <- file.path(root, "final_analysis_20260907", "Figure8_robustness_incremental")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_classic(base_family = "Times", base_size = text_size_pt) +
  theme(
    text = element_text(family = "Times", size = text_size_pt),
    plot.title = element_blank(), plot.subtitle = element_blank(),
    axis.title = element_text(family = "Times", size = text_size_pt),
    axis.text = element_text(family = "Times", size = text_size_pt),
    legend.title = element_text(family = "Times", size = text_size_pt),
    legend.text = element_text(family = "Times", size = text_size_pt),
    strip.text = element_text(family = "Times", size = text_size_pt),
    panel.grid.major = element_line(linewidth = line_width, colour = "#E5E5E5"),
    panel.grid.minor = element_blank(), plot.margin = margin(4, 5, 4, 5),
    legend.position = "bottom"
  )

dat <- read_csv(file.path(base, "independent_ECM_comparator/TCGA_incremental_value_input.csv"), show_col_types = FALSE) %>%
  mutate(
    stage_simpleIII_IV = factor(stage_simple, levels = c("I_II", "III_IV")),
    time = as.numeric(os_time), event = as.numeric(os_event),
    stage_num = as.numeric(stage_simpleIII_IV == "III_IV")
  ) %>%
  filter(complete.cases(time, event, age, CAF_score, ECM_signature, CCDC80_ECM_program, stage_simpleIII_IV))

dat$stage_simpleIII_IV <- relevel(dat$stage_simpleIII_IV, ref = "I_II")

fits <- list(
  "Stage" = coxph(Surv(time, event) ~ stage_simpleIII_IV, data = dat, x = TRUE, y = TRUE),
  "Stage + age" = coxph(Surv(time, event) ~ stage_simpleIII_IV + age, data = dat, x = TRUE, y = TRUE),
  "Stage + CAF" = coxph(Surv(time, event) ~ stage_simpleIII_IV + age + CAF_score, data = dat, x = TRUE, y = TRUE),
  "Stage + CAF + ECM" = coxph(Surv(time, event) ~ stage_simpleIII_IV + age + CAF_score + ECM_signature, data = dat, x = TRUE, y = TRUE),
  "Full CCDC80" = coxph(Surv(time, event) ~ stage_simpleIII_IV + age + CAF_score + ECM_signature + CCDC80_ECM_program, data = dat, x = TRUE, y = TRUE)
)

# A: use the precomputed 1,000 bootstrap Cox fits from the project audit.
boot <- read_csv(file.path(base, "independent_ECM_comparator/TCGA_CCDC80_program_bootstrap_1000.csv"), show_col_types = FALSE)
boot_summary <- read_csv(file.path(base, "independent_ECM_comparator/TCGA_CCDC80_program_bootstrap_summary.csv"), show_col_types = FALSE)
boot_ci <- quantile(boot$HR, probs = c(.025, .5, .975), na.rm = TRUE)
boot_df <- boot %>% mutate(logHR = log(HR), density_group = "Bootstrap")
pA <- ggplot(boot_df, aes(x = HR)) +
  geom_density(fill = wong[2], colour = wong[7], alpha = .35, linewidth = data_line_width) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = line_width, colour = "#555555") +
  geom_vline(xintercept = boot_ci[2], colour = wong[7], linewidth = data_line_width) +
  annotate("text", x = boot_ci[2], y = Inf, label = sprintf("median %.2f; 95%% CI %.2f–%.2f", boot_ci[2], boot_ci[1], boot_ci[3]), vjust = 1.7, hjust = -.05, family = "Times", size = text_size_pt / 2.85) +
  labs(x = "Bootstrap hazard ratio", y = "Density") + theme_pub

# B: precomputed cohort omission estimates plus pooled estimate from full model.
loo <- read_csv(file.path(base, "independent_ECM_comparator/leave_one_cohort_out_CCDC80_incremental_Cox.csv"), show_col_types = FALSE) %>%
  transmute(model = paste0("Omit ", omitted), HR = HR, lo = CI_low, hi = CI_high, p = p)
full_coef <- broom::tidy(fits[["Full CCDC80"]], exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "CCDC80_ECM_program") %>% transmute(model = "All cohorts", HR = estimate, lo = conf.low, hi = conf.high, p = p.value)
loo_df <- bind_rows(full_coef, loo) %>% mutate(model = factor(model, levels = rev(model)))
pB <- ggplot(loo_df, aes(x = HR, y = model)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = line_width, colour = "#555555") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .18, linewidth = data_line_width, colour = wong[7]) +
  geom_point(size = point_size, colour = wong[7]) +
  geom_text_repel(aes(label = sprintf("%.2f [%.2f–%.2f]; P=%.3g", HR, lo, hi, p)), direction = "y", hjust = 0, nudge_x = .12, size = text_size_pt / 3.2, family = "Times", min.segment.length = 0, box.padding = .25, seed = 20260908) +
  scale_x_log10(expand = expansion(mult = c(.05, .42))) + labs(x = "CCDC80 program HR (95% CI)", y = NULL) + theme_pub

# C: effect path of the program as covariates are added.
program_fits <- list(
  "Program only" = coxph(Surv(time, event) ~ CCDC80_ECM_program, data = dat, x = TRUE, y = TRUE),
  "Stage + program" = coxph(Surv(time, event) ~ stage_simpleIII_IV + CCDC80_ECM_program, data = dat, x = TRUE, y = TRUE),
  "Stage + age + CAF + program" = coxph(Surv(time, event) ~ stage_simpleIII_IV + age + CAF_score + CCDC80_ECM_program, data = dat, x = TRUE, y = TRUE),
  "Full CCDC80" = fits[["Full CCDC80"]]
)
path_df <- bind_rows(lapply(names(program_fits), function(nm) {
  z <- broom::tidy(program_fits[[nm]], exponentiate = TRUE, conf.int = TRUE)
  z <- z[z$term == "CCDC80_ECM_program", , drop = FALSE]
  data.frame(model = nm, HR = z$estimate, lo = z$conf.low, hi = z$conf.high, p = z$p.value)
})) %>% mutate(model = factor(model, levels = names(program_fits)))
pC <- ggplot(path_df, aes(x = model, y = HR, group = 1)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = line_width, colour = "#555555") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = data_line_width, colour = wong[7]) +
  geom_point(size = point_size, colour = wong[7]) +
  geom_text_repel(aes(label = sprintf("HR %.2f; P=%.3g", HR, p)), direction = "both", box.padding = .65, point.padding = .45, force = 2, min.segment.length = 0, seed = 20260908, size = text_size_pt / 3.2, family = "Times") +
  scale_y_log10() + labs(x = NULL, y = "CCDC80 program HR (95% CI)") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# D: apparent and optimism-corrected C-index via bootstrap optimism correction.
concordance_index <- function(fit, newdata) {
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp"))
  z <- data.frame(time = newdata$time, event = newdata$event, lp = -lp)
  if (sum(z$event) == 0 || length(unique(z$time)) < 2) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(time, event) ~ lp, data = z)$concordance), error = function(e) NA_real_)
}
set.seed(20260908)
ci_rows <- lapply(names(fits), function(nm) {
  fit <- fits[[nm]]
  apparent <- concordance_index(fit, dat)
  opt <- numeric(bootstrap_n)
  for (i in seq_len(optimism_bootstrap_n)) {
    idx <- sample.int(nrow(dat), replace = TRUE)
    bdat <- dat[idx, , drop = FALSE]
    bfit <- tryCatch(update(fit, data = bdat), error = function(e) NULL)
    if (is.null(bfit)) { opt[i] <- NA_real_; next }
    opt[i] <- concordance_index(bfit, bdat) - concordance_index(bfit, dat)
  }
  corrected <- apparent - mean(opt, na.rm = TRUE)
  data.frame(model = nm, apparent = apparent, corrected = corrected,
             lo = quantile(opt, .025, na.rm = TRUE), hi = quantile(opt, .975, na.rm = TRUE))
}) %>% bind_rows()
ci_long <- ci_rows %>% select(model, apparent, corrected) %>% pivot_longer(-model, names_to = "metric", values_to = "cindex") %>%
  mutate(model = factor(model, levels = names(fits)), metric = factor(metric, levels = c("apparent", "corrected")))
pD <- ggplot(ci_long, aes(x = cindex, y = model, colour = metric)) +
  geom_point(position = position_dodge(width = .45), size = point_size) +
  geom_text(aes(label = sprintf("%.3f", cindex)), position = position_dodge(width = .45), hjust = -.15, size = text_size_pt / 3.2, family = "Times") +
  geom_segment(data = ci_rows, aes(x = apparent, xend = corrected, y = model, yend = model), inherit.aes = FALSE, linewidth = line_width, colour = "#777777") +
  scale_colour_manual(values = c("apparent" = wong[1], "corrected" = wong[7]), labels = c("Apparent", "Bootstrap optimism-corrected")) +
  labs(x = "Harrell's C-index", y = NULL, colour = NULL) + theme_pub

# E/F: riskRegression implementation for time-dependent AUC and Brier scores.
score <- Score(fits, formula = Surv(time, event) ~ 1, data = dat, times = time_points,
               metrics = c("AUC", "Brier"), summary = "IPA", conf.int = TRUE, plots = NULL)
score_df <- as.data.frame(score$AUC$score) %>% filter(model != "Null model") %>%
  transmute(model = factor(model, levels = names(fits)), time = times / 365.25, value = AUC, lo = lower, hi = upper)
brier_df <- as.data.frame(score$Brier$score) %>% filter(model != "Null model") %>%
  transmute(model = factor(model, levels = names(fits)), time = times / 365.25, value = Brier, lo = lower, hi = upper)
auc_heat <- score_df %>% mutate(time_label = paste0(round(time), " y"))
pE <- ggplot(auc_heat, aes(x = time_label, y = model, fill = value)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = sprintf("%.2f", value)), family = "Times", size = text_size_pt / 3.2) +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1), breaks = c(.4, .6, .8), name = "AUC") +
  labs(x = "Prediction horizon", y = NULL) + theme_pub + theme(legend.position = "right")
calibration_time <- 1825
calibration_df <- bind_rows(lapply(names(fits), function(nm) {
  dd <- dat
  bh <- basehaz(fits[[nm]], centered = FALSE)
  H0 <- approx(bh$time, bh$hazard, xout = calibration_time, rule = 2)$y
  lp <- as.numeric(predict(fits[[nm]], newdata = dd, type = "lp"))
  dd$predicted_risk <- 1 - exp(-H0 * exp(lp))
  dd$bin <- dplyr::ntile(dd$predicted_risk, 10)
  bind_rows(lapply(sort(unique(dd$bin)), function(bb) {
    x <- dd[dd$bin == bb, , drop = FALSE]
    sf <- survfit(Surv(time, event) ~ 1, data = x)
    obs <- 1 - summary(sf, times = calibration_time, extend = TRUE)$surv
    data.frame(model = nm, bin = bb, predicted = mean(x$predicted_risk), observed = obs, n = nrow(x))
  }))
}))
calibration_plot_df <- calibration_df %>% filter(model %in% c("Stage", "Full CCDC80"))
pF <- ggplot(calibration_plot_df, aes(x = predicted, y = observed, colour = model, group = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = line_width, colour = "#555555") +
  geom_line(linewidth = data_line_width) + geom_point(size = point_size) +
  scale_colour_manual(values = series_cols) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted 5-year risk", y = "Observed 5-year risk", colour = NULL) + theme_pub

fig <- (pA | pB) / (pC | pD) / (pE | pF) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(family = "Times", size = text_size_pt, face = "bold"))
ggsave(file.path(out, "Figure8_robustness_incremental.pdf"), fig, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "transparent")
ggsave(file.path(out, "Figure8_robustness_incremental.svg"), fig, width = width_mm, height = height_mm, units = "mm", device = svglite, bg = "transparent")
ggsave(file.path(out, "Figure8_robustness_incremental.tif"), fig, width = width_mm, height = height_mm, units = "mm", dpi = 600, compression = "lzw", device = "tiff", bg = "white")

write_csv(boot_df, file.path(out, "Figure8_source_bootstrap.csv"))
write_csv(loo_df, file.path(out, "Figure8_source_leave_one_cohort_out.csv"))
write_csv(path_df, file.path(out, "Figure8_source_nested_effect_path.csv"))
write_csv(ci_rows, file.path(out, "Figure8_source_cindex_bootstrap_optimism.csv"))
write_csv(score_df, file.path(out, "Figure8_source_time_dependent_AUC.csv"))
write_csv(brier_df, file.path(out, "Figure8_source_Brier.csv"))
write_csv(calibration_df, file.path(out, "Figure8_source_calibration_5y.csv"))
writeLines(c(
  "Figure 8 quantitative robustness and incremental prediction value completed.",
  paste0("TCGA complete-case N=", nrow(dat), "; events=", sum(dat$event), "; bootstrap=", bootstrap_n),
  "No schematic, workflow, or conceptual diagram was generated.",
  paste0("AUC/Brier were calculated with riskRegression::Score; C-index optimism was estimated by ", optimism_bootstrap_n, " bootstrap resamples."),
  sprintf("CCDC80 bootstrap median HR=%.3f; 95%% interval=%.3f-%.3f", boot_ci[2], boot_ci[1], boot_ci[3])
), file.path(out, "Figure8_STATUS.txt"))
