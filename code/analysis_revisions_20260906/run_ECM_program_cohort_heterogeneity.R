base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
rev <- file.path(base, "analysis_revisions_20260906")
out <- file.path(rev, "independent_ECM_comparator")
library(survival)

load_cohort <- function(name, program_file, signature_file) {
  p <- read.csv(file.path(rev, "CCDC80_ECM_program_validation", program_file))
  names(p)[names(p) == "program_score"] <- "CCDC80_ECM_program"
  s <- read.csv(file.path(rev, "ECM_CAF_signature_validation", signature_file))
  keep <- intersect(c("sample_id", "os_time", "os_event", "stage_simple", "age", "CAF_score", "ECM_signature"), names(s))
  d <- merge(p, s[, keep], by = "sample_id")
  d$cohort <- name
  d
}

tcga <- read.csv(file.path(out, "TCGA_incremental_value_input.csv"), check.names = FALSE)
tcga$cohort <- "TCGA"
g13507 <- load_cohort("GSE13507", "GSE13507_ECM_program_scores.csv", "GSE13507_ECM_signature_input.csv")
g32894 <- load_cohort("GSE32894", "GSE32894_ECM_program_scores.csv", "GSE32894_ECM_signature_input.csv")

standardize <- function(d) {
  d$CCDC80_ECM_program <- as.numeric(scale(d$CCDC80_ECM_program))
  d$ECM_signature <- as.numeric(scale(d$ECM_signature))
  d
}
tcga <- standardize(tcga)
g13507 <- standardize(g13507)
g32894 <- standardize(g32894)

fit_one <- function(d, model_name) {
  d <- d[complete.cases(d[, c("os_time", "os_event", "ECM_signature", "CCDC80_ECM_program")]), ]
  f <- coxph(Surv(os_time, os_event) ~ ECM_signature + CCDC80_ECM_program, data = d)
  sm <- summary(f)$coefficients
  ci <- summary(f)$conf.int
  data.frame(model = model_name, n = nrow(d), events = sum(d$os_event), term = rownames(sm), HR = ci[, "exp(coef)"], CI_low = ci[, "lower .95"], CI_high = ci[, "upper .95"], p = sm[, "Pr(>|z|)"], row.names = NULL)
}
by_cohort <- do.call(rbind, list(fit_one(tcga, "TCGA"), fit_one(g13507, "GSE13507"), fit_one(g32894, "GSE32894")))
write.csv(by_cohort, file.path(out, "cohort_specific_ECM_CCDC80_Cox.csv"), row.names = FALSE)

pooled <- rbind(tcga[, c("os_time", "os_event", "ECM_signature", "CCDC80_ECM_program", "cohort")], g13507[, c("os_time", "os_event", "ECM_signature", "CCDC80_ECM_program", "cohort")], g32894[, c("os_time", "os_event", "ECM_signature", "CCDC80_ECM_program", "cohort")])
pooled <- pooled[complete.cases(pooled), ]
pooled$cohort <- factor(pooled$cohort)
f_add <- coxph(Surv(os_time, os_event) ~ strata(cohort) + ECM_signature + CCDC80_ECM_program, data = pooled)
f_int <- coxph(Surv(os_time, os_event) ~ strata(cohort) + ECM_signature + CCDC80_ECM_program * cohort, data = pooled)
lr <- anova(f_add, f_int, test = "LRT")
write.csv(as.data.frame(lr), file.path(out, "CCDC80_cohort_interaction_LRT.csv"), row.names = FALSE)
write.csv(as.data.frame(summary(f_int)$coefficients), file.path(out, "CCDC80_cohort_interaction_coefficients.csv"))
writeLines(c("Completed cohort-specific and interaction Cox analysis.", "Scores standardized within cohort.", "Interaction test evaluates heterogeneity of CCDC80 program effect across cohorts."), file.path(out, "COHORT_HETEROGENEITY_STATUS.txt"))
print(by_cohort)
print(lr)
