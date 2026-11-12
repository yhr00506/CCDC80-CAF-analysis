#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(dplyr); library(flextable); library(officer)})

root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
final <- file.path(root, "final_analysis_20260907")
out <- file.path(root, "analysis_revisions_20260909", "supplementary_tables")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

fmt_p <- function(x) ifelse(is.na(x), "NA", ifelse(x < .001, "<0.001", sprintf("%.3f", x)))
fmt_num <- function(x, d=3) ifelse(is.na(x), "NA", sprintf(paste0("%.",d,"f"), x))
style_tbl <- function(x) {
  ft <- flextable(x)
  ft <- theme_booktabs(ft)
  ft <- fontsize(ft, size=8, part="all")
  ft <- font(ft, fontname="Arial", part="all")
  ft <- bold(ft, part="header")
  ft <- align(ft, align="center", part="all")
  ft <- align(ft, j=1, align="left", part="all")
  ft <- autofit(ft)
  ft <- set_table_properties(ft, layout="autofit", width=1)
  ft
}
write_doc <- function(ft, title, note, file) {
  doc <- read_docx()
  doc <- body_add_par(doc, title, style="heading 1")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, note, style="Normal")
  print(doc, target=file)
}

# ST1: cohort and CAF-state composition.
cohort <- read_csv(file.path(final,"Figure1_CAF_atlas/Figure1E_cohort_subtype_composition_source.csv"), show_col_types=FALSE) %>%
  transmute(Cohort=dataset, `CAF state`=CAF_subtype_global, `CAF cells`=n, `Fraction`=sprintf("%.1f%%", 100*fraction))
st1 <- style_tbl(cohort)
write_doc(st1, "Supplementary Table 1. Integrated BLCA CAF cohort and state composition",
          "Counts and within-cohort fractions are derived from the globally annotated integrated CAF object (four cohorts; 5,589 CAF cells).",
          file.path(out,"Supplementary_Table_1_CAF_cohort_composition.docx"))

# ST2: program gene stability.
stab <- read_csv(file.path(final,"SupplementaryFigure3_program_gene_stability/SuppFig3_source_gene_stability.csv"), show_col_types=FALSE)
loo <- read_csv(file.path(final,"SupplementaryFigure3_program_gene_stability/SuppFig3_source_leave_one_gene_out.csv"), show_col_types=FALSE)
st2d <- left_join(stab, loo, by="gene") %>% transmute(Gene=gene, `Repeat runs`=n_runs,
  `Mean virtual-KO FC`=fmt_num(mean_FC), `SD of FC`=fmt_num(sd_FC), `Mean |Z|`=fmt_num(mean_abs_Z),
  `Leave-one-gene-out Spearman r`=fmt_num(correlation))
st2 <- style_tbl(st2d)
write_doc(st2, "Supplementary Table 2. Stability of the 19-gene CCDC80-associated ECM program",
          "Mean FC, SD, and mean absolute Z are calculated across three virtual-knockout repeats. The final column is the Spearman correlation between the full score and the score after removing the indicated gene across 5,589 integrated CAF cells.",
          file.path(out,"Supplementary_Table_2_program_gene_stability.docx"))

# ST3: external validation.
forest <- read_csv(file.path(final,"Figure5_external_validation/Figure5_source_forest.csv"), show_col_types=FALSE) %>%
  transmute(Evidence="Program survival association", Cohort=cohort, `N`=NA_character_, Estimate="HR",
            `Estimate (95% CI)`=sprintf("%.2f (%.2f–%.2f)", HR, lo, hi), `P value`=NA_character_)
corr <- read_csv(file.path(final,"Figure5_external_validation/Figure5_source_correlations.csv"), show_col_types=FALSE) %>%
  transmute(Evidence=metric, Cohort=cohort, N=as.character(n), Estimate="Spearman ρ",
            `Estimate (95% CI)`=sprintf("%.2f", rho), `P value`=fmt_p(p))
st3d <- bind_rows(forest, corr)
st3 <- style_tbl(st3d)
write_doc(st3, "Supplementary Table 3. External validation of program associations",
          "HRs summarize external program survival associations; Spearman correlations summarize associations with CAF and myeloid scores. Confidence intervals were not available in the source forest table and are therefore not imputed.",
          file.path(out,"Supplementary_Table_3_external_validation.docx"))

# ST4: clinical candidate and mechanism-hypothesis evidence.
cox <- read_csv(file.path(final,"Figure7_clinical_prioritization/Figure7_source_Cox.csv"), show_col_types=FALSE) %>%
  transmute(Evidence="Program Cox model", Candidate=term, `N`=NA_character_, `Events`=NA_character_,
            `HR (95% CI)`=sprintf("%.2f (%.2f–%.2f)", HR, lo, hi), `P value`=fmt_p(p), FDR=NA_character_)
cand <- read_csv(file.path(final,"Figure7_clinical_prioritization/Figure7_source_candidate_priority.csv"), show_col_types=FALSE) %>%
  transmute(Evidence="Candidate gene", Candidate=gene, N=as.character(n), Events=as.character(events),
            `HR (95% CI)`=sprintf("%.2f (%.2f–%.2f)", HR, CI_lower, CI_upper), `P value`=fmt_p(p_value), FDR=fmt_p(FDR))
lr <- read_csv(file.path(final,"Figure6_spatial_mechanistic/Figure6_source_LR_top8.csv"), show_col_types=FALSE) %>%
  transmute(Evidence="LR hypothesis", Candidate=pair, N=as.character(n_patients), Events=NA_character_,
            `HR (95% CI)`=sprintf("evidence score %.0f", evidence_score), `P value`=NA_character_, FDR=fmt_p(FDR_LR_myeloid_score))
st4 <- style_tbl(bind_rows(cox,cand,lr))
write_doc(st4, "Supplementary Table 4. Clinical prioritization and mechanism-hypothesis evidence",
          "Candidate-gene rows summarize TCGA clinical prioritization. LR rows are hypothesis-generating rankings and should not be interpreted as causal communication evidence. Program Cox models use the final analysis definitions.",
          file.path(out,"Supplementary_Table_4_clinical_mechanistic_evidence.docx"))

# Combined Word file with each table on a separate page.
doc <- read_docx()
for (i in 1:4) {
  dat <- list(st1,st2,st3,st4)[[i]]
  ttl <- c("Supplementary Table 1. Integrated BLCA CAF cohort and state composition", "Supplementary Table 2. Stability of the 19-gene CCDC80-associated ECM program", "Supplementary Table 3. External validation of program associations", "Supplementary Table 4. Clinical prioritization and mechanism-hypothesis evidence")[[i]]
  note <- c("Counts and within-cohort fractions are derived from the globally annotated integrated CAF object.", "Results are based on three virtual-knockout repeats and leave-one-gene-out score stability.", "HRs and correlations are reproduced from the external-validation source tables.", "LR results are hypothesis-generating and not causal evidence.")[[i]]
  if (i > 1) doc <- body_add_break(doc, pos="after")
  doc <- body_add_par(doc, ttl, style="heading 1")
  doc <- body_add_flextable(doc, dat)
  doc <- body_add_par(doc, note, style="Normal")
}
print(doc, target=file.path(out,"Supplementary_Tables_1-4.docx"))
writeLines(c("Four supplementary tables exported to Word.", "ST1 cohort/state composition; ST2 program stability; ST3 external validation; ST4 clinical/mechanistic evidence."), file.path(out,"STATUS.txt"))
