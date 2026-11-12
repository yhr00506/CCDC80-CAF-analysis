suppressPackageStartupMessages(library(survival))

arguments <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(arguments) >= 1) arguments[[1]] else "."
input_file <- file.path(input_dir, "TCGA_patient_level_mycaf_scores.csv")
reference_file <- file.path(input_dir, "incremental_LRT_reference.csv")
stopifnot(file.exists(input_file), file.exists(reference_file))

data <- read.csv(input_file, check.names = FALSE)
stopifnot(nrow(data) == 400, !anyDuplicated(data$patient_id), sum(data$os_event) == 176)
stopifnot(all(is.finite(data$os_time)), all(data$os_time >= 0))
data$stage_simple <- factor(data$stage_simple, levels = c("I_II", "III_IV"))
data$sex <- factor(data$sex)
data$age_z <- as.numeric(scale(data$age))
data$program_z <- as.numeric(scale(data$CCDC80_program))
data$myCAF_derived_z <- as.numeric(scale(data$myCAF_derived))

comparators <- c("CAF_9gene", "ECM_3gene", "collagen_core", "matrisome_2014", "Hallmark_TGF_Beta_2025.1")
for (name in comparators) data[[paste0(name, "_z")]] <- as.numeric(scale(data[[name]]))
families <- list(clinical = character(), myCAF = "myCAF_derived_z")
for (name in comparators) families[[name]] <- paste0(name, "_z")
families$all_existing <- paste0(comparators, "_z")
families$all_existing_plus_myCAF <- c(paste0(comparators, "_z"), "myCAF_derived_z")

results <- lapply(names(families), function(name) {
  base_formula <- reformulate(c("stage_simple", "age_z", "sex", families[[name]]), response = "Surv(os_time, os_event)")
  base <- coxph(base_formula, data = data, ties = "efron")
  full <- coxph(update(base_formula, . ~ . + program_z), data = data, ties = "efron")
  likelihood_ratio <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(base)))
  probability <- pchisq(likelihood_ratio, df = 1, lower.tail = FALSE)
  confidence <- summary(full)$conf.int["program_z", ]
  data.frame(
    benchmark = name, n = full$n, events = full$nevent,
    HR = unname(confidence["exp(coef)"]),
    CI_low = unname(confidence["lower .95"]),
    CI_high = unname(confidence["upper .95"]),
    Wald_P = unname(summary(full)$coefficients["program_z", "Pr(>|z|)"]),
    LR = likelihood_ratio, LRT_P = probability
  )
})
results <- do.call(rbind, results)
results$FDR_BH <- p.adjust(results$LRT_P, method = "BH")
reference <- read.csv(reference_file, check.names = FALSE)
reference <- reference[match(results$benchmark, reference$benchmark), ]
stopifnot(identical(results$benchmark, reference$benchmark))
stopifnot(max(abs(results$LR - reference$LR)) < 1e-6)
stopifnot(max(abs(results$LRT_P - reference$p)) < 1e-8)
stopifnot(max(abs(results$FDR_BH - reference$FDR_BH)) < 1e-8)

write.csv(results, file.path(input_dir, "tcga_reproduced_models.csv"), row.names = FALSE)
cat("Verified", nrow(results), "model comparisons for", nrow(data), "unique patients and", sum(data$os_event), "events.\n")
