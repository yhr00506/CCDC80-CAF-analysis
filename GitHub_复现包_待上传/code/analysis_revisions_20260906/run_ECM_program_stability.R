base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
rev <- file.path(base, "analysis_revisions_20260906")
out <- file.path(rev, "independent_ECM_comparator")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
library(survival)
set.seed(20260907)

tcga <- read.csv(file.path(out, "TCGA_incremental_value_input.csv"), check.names = FALSE)
tcga$stage_simple <- factor(tcga$stage_simple)

fit_boot <- function(dat, idx) {
  b <- dat[idx, , drop = FALSE]
  fit <- try(coxph(Surv(os_time, os_event) ~ stage_simple + age + CAF_score + ECM_signature + CCDC80_ECM_program, data = b), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(HR = NA, p = NA, coef = NA))
  sm <- summary(fit)$coefficients
  if (!"CCDC80_ECM_program" %in% rownames(sm)) return(c(HR = NA, p = NA, coef = NA))
  c(HR = exp(sm["CCDC80_ECM_program", "coef"]), p = sm["CCDC80_ECM_program", "Pr(>|z|)"], coef = sm["CCDC80_ECM_program", "coef"])
}

nboot <- 1000
boot <- do.call(rbind, lapply(seq_len(nboot), function(i) fit_boot(tcga, sample.int(nrow(tcga), nrow(tcga), replace = TRUE))))
boot <- as.data.frame(boot)
boot$iteration <- seq_len(nrow(boot))
write.csv(boot, file.path(out, "TCGA_CCDC80_program_bootstrap_1000.csv"), row.names = FALSE)
valid_hr <- boot$HR[is.finite(boot$HR)]
valid_p <- boot$p[is.finite(boot$p)]
boot_summary <- data.frame(
  n_boot = nboot,
  n_valid = length(valid_hr),
  HR_median = median(valid_hr),
  HR_q025 = quantile(valid_hr, 0.025),
  HR_q975 = quantile(valid_hr, 0.975),
  fraction_HR_gt_1 = mean(valid_hr > 1),
  fraction_p_lt_005 = mean(valid_p < 0.05),
  fraction_coef_positive = mean(boot$coef[is.finite(boot$coef)] > 0)
)
write.csv(boot_summary, file.path(out, "TCGA_CCDC80_program_bootstrap_summary.csv"), row.names = FALSE)

load_cohort <- function(name, program_file, signature_file) {
  p <- read.csv(file.path(rev, "CCDC80_ECM_program_validation", program_file))
  names(p)[names(p) == "program_score"] <- "CCDC80_ECM_program"
  s <- read.csv(file.path(rev, "ECM_CAF_signature_validation", signature_file))
  s <- s[, intersect(c("sample_id", "os_time", "os_event", "stage", "CAF_score", "ECM_signature"), names(s))]
  d <- merge(p, s, by = "sample_id")
  d$cohort <- name
  d
}
g13507 <- load_cohort("GSE13507", "GSE13507_ECM_program_scores.csv", "GSE13507_ECM_signature_input.csv")
g32894 <- load_cohort("GSE32894", "GSE32894_ECM_program_scores.csv", "GSE32894_ECM_signature_input.csv")
tcga_simple <- tcga[, c("sample_id", "os_time", "os_event", "ECM_signature", "CCDC80_ECM_program")]
tcga_simple$cohort <- "TCGA"
pooled <- rbind(tcga_simple, g13507[, names(tcga_simple)], g32894[, names(tcga_simple)])
pooled <- pooled[complete.cases(pooled), ]
pooled$cohort <- factor(pooled$cohort)

fit_loo <- function(omit) {
  d <- pooled[pooled$cohort != omit, ]
  f <- try(coxph(Surv(os_time, os_event) ~ strata(cohort) + ECM_signature + CCDC80_ECM_program, data = d), silent = TRUE)
  if (inherits(f, "try-error")) return(data.frame(omitted = omit, n = nrow(d), events = sum(d$os_event), HR = NA, CI_low = NA, CI_high = NA, p = NA))
  sm <- summary(f)$coefficients
  ci <- summary(f)$conf.int
  data.frame(omitted = omit, n = nrow(d), events = sum(d$os_event), HR = ci["CCDC80_ECM_program", "exp(coef)"], CI_low = ci["CCDC80_ECM_program", "lower .95"], CI_high = ci["CCDC80_ECM_program", "upper .95"], p = sm["CCDC80_ECM_program", "Pr(>|z|)"])
}
loo <- do.call(rbind, lapply(levels(pooled$cohort), fit_loo))
write.csv(loo, file.path(out, "leave_one_cohort_out_CCDC80_incremental_Cox.csv"), row.names = FALSE)

writeLines(c("Completed CCDC80 program stability analysis.", paste0("TCGA bootstrap iterations=", nboot), "Leave-one-cohort-out model: stratified Cox with independent ECM signature + CCDC80 program.", "Interpretation is robustness of association, not causality."), file.path(out, "STABILITY_STATUS.txt"))
print(boot_summary)
print(loo)
