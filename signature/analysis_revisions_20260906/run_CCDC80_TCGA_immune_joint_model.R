library(SummarizedExperiment)
library(survival)

base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
out <- file.path(base, "analysis_revisions_20260906", "CCDC80_TCGA_immune_joint_validation")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

se <- readRDS(file.path(base, "00_data_audit", "TCGA_expression.rds"))
expr <- assay(se)
rd <- as.data.frame(rowData(se))
gene <- as.character(rd$gene_name)
keep <- !is.na(gene) & gene != ""
expr <- expr[keep,,drop=FALSE]
gene <- gene[keep]
expr <- expr[!duplicated(gene),,drop=FALSE]
rownames(expr) <- gene[!duplicated(gene)]
expr <- log2(as.matrix(expr) + 1)

clinical <- read.csv(file.path(base, "00_data_audit", "TCGA_clinical.csv"), check.names=FALSE, stringsAsFactors=FALSE)
sample_id <- substr(colnames(expr), 1, 12)
clinical <- clinical[match(sample_id, clinical$submitter_id),,drop=FALSE]

get_gene <- function(g) {
  if (!g %in% rownames(expr)) return(rep(NA_real_, ncol(expr)))
  as.numeric(expr[g,])
}
score <- function(genes) {
  genes <- intersect(genes, rownames(expr))
  if (!length(genes)) return(rep(NA_real_, ncol(expr)))
  z <- t(scale(t(expr[genes,,drop=FALSE])))
  colMeans(z, na.rm=TRUE)
}

dat <- data.frame(
  sample_id=sample_id,
  CCDC80=get_gene("CCDC80"),
  myeloid_score=score(c("CD68","CD163","C1QA","C1QB","C1QC","LILRB1","FCER1G","TYROBP","APOC1","SPP1","MSR1")),
  cytotoxicity=score(c("CD8A","CD8B","NKG7","GNLY","PRF1","GZMB","GZMK")),
  T_exhaustion=score(c("PDCD1","CTLA4","TIGIT","LAG3","HAVCR2","TOX","CXCL13")),
  CAF_score=score(c("COL1A1","COL1A2","DCN","LUM","COL3A1","COL5A1","SPARC","FAP","PDGFRA")),
  stage=clinical$ajcc_pathologic_stage,
  vital_status=clinical$vital_status,
  days_to_death=clinical$days_to_death,
  days_to_last_follow_up=clinical$days_to_last_follow_up,
  stringsAsFactors=FALSE
)
dat$os_time <- suppressWarnings(as.numeric(ifelse(!is.na(dat$days_to_death) & dat$days_to_death != "", dat$days_to_death, dat$days_to_last_follow_up)))
dat$os_event <- as.integer(tolower(dat$vital_status) == "dead")
dat$stage_simple <- sub("^Stage ", "", dat$stage)
dat$stage_simple[dat$stage_simple %in% c("I","II")] <- "I_II"
dat$stage_simple[dat$stage_simple %in% c("III","IV")] <- "III_IV"
dat <- dat[complete.cases(dat[,c("CCDC80","myeloid_score","cytotoxicity","T_exhaustion","CAF_score","os_time","os_event","stage_simple")]),,drop=FALSE]
dat$CCDC80_z <- as.numeric(scale(dat$CCDC80))

write.csv(dat, file.path(out, "TCGA_CCDC80_immune_joint_input.csv"), row.names=FALSE)

cor_genes <- c("myeloid_score","cytotoxicity","T_exhaustion","CAF_score")
cor_out <- do.call(rbind, lapply(cor_genes, function(v) {
  tt <- suppressWarnings(cor.test(dat$CCDC80, dat[[v]], method="spearman", exact=FALSE))
  data.frame(variable=v, n=nrow(dat),rho=unname(tt$estimate),p=tt$p.value)
}))
cor_out$FDR <- p.adjust(cor_out$p, method="BH")
write.csv(cor_out, file.path(out, "TCGA_CCDC80_immune_correlations.csv"), row.names=FALSE)

fits <- list(
  univariable = c("CCDC80_z"),
  immune_adjusted = c("CCDC80_z","myeloid_score","cytotoxicity","T_exhaustion","stage_simple"),
  CAF_immune_adjusted = c("CCDC80_z","CAF_score","myeloid_score","cytotoxicity","T_exhaustion","stage_simple")
)
cox_rows <- list()
for (nm in names(fits)) {
  f <- as.formula(paste0("Surv(os_time,os_event) ~ ", paste(fits[[nm]], collapse=" + ")))
  fit <- coxph(f, data=dat)
  sm <- summary(fit)
  for (i in seq_along(coef(fit))) {
    cox_rows[[length(cox_rows)+1]] <- data.frame(model=nm, variable=names(coef(fit))[i], n=nrow(dat), events=sum(dat$os_event), HR=unname(exp(coef(fit)[i])), CI_low=unname(sm$conf.int[i,"lower .95"]), CI_high=unname(sm$conf.int[i,"upper .95"]), p=unname(sm$coefficients[i,"Pr(>|z|)"]))
  }
}
cox_out <- do.call(rbind, cox_rows)
cox_out$FDR <- ave(cox_out$p, cox_out$model, FUN=function(x) p.adjust(x, "BH"))
write.csv(cox_out, file.path(out, "TCGA_CCDC80_immune_joint_Cox.csv"), row.names=FALSE)

writeLines(c(
  "Analysis: CCDC80 patient-level immune joint validation",
  paste0("N=", nrow(dat), "; events=", sum(dat$os_event)),
  "Scores are exploratory z-score averages from predefined immune/CAF marker panels.",
  "The single-sample scRNA-seq TME object was not used for patient-level inference."
), file.path(out, "RUN_STATUS.txt"))
