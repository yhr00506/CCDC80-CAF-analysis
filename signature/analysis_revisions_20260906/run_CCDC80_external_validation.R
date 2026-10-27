library(Biobase)
library(survival)

base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
out <- file.path(base, "analysis_revisions_20260906", "CCDC80_external_validation")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

cohorts <- list(
  GSE32894=list(file="00_data_audit/GSE32894_expression.rds", event="dod_event_(yes/no):ch1", time="time_to_dod_(months):ch1"),
  GSE13507=list(file="00_data_audit/GSE13507_expression.rds", event="overall survival:ch1", time="survival month:ch1")
)

get_probe_expr <- function(eset) {
  x <- exprs(eset)
  fd <- as.data.frame(fData(eset))
  sym_col <- intersect(c("Symbol","ILMN_Gene"), names(fd))[1]
  sym <- toupper(trimws(as.character(fd[[sym_col]])))
  ok <- !is.na(sym) & sym != "" & sym != "---"
  x <- x[ok,,drop=FALSE]; sym <- sym[ok]
  # Average duplicate probes after mapping to gene symbols.
  summed <- rowsum(x, group=sym, reorder=FALSE)
  out <- summed / as.numeric(table(sym)[rownames(summed)])
  out
}
get_score <- function(x, genes) {
  genes <- intersect(genes, rownames(x))
  if (!length(genes)) return(rep(NA_real_, ncol(x)))
  z <- t(scale(t(x[genes,,drop=FALSE])))
  colMeans(z, na.rm=TRUE)
}
parse_event <- function(x) {
  y <- tolower(trimws(as.character(x)))
  as.integer(y %in% c("yes","dead","death","deceased","1","true","progression"))
}
parse_time <- function(x) {
  y <- suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", as.character(x))))
  y[y <= 0] <- NA_real_
  y
}

caf_genes <- c("COL1A1","COL1A2","DCN","LUM","COL3A1","COL5A1","SPARC","FAP","PDGFRA")
myeloid_genes <- c("CD68","CD163","C1QA","C1QB","C1QC","LILRB1","FCER1G","TYROBP","APOC1","SPP1","MSR1")
cyto_genes <- c("CD8A","CD8B","NKG7","GNLY","PRF1","GZMB","GZMK")
exhaust_genes <- c("PDCD1","CTLA4","TIGIT","LAG3","HAVCR2","TOX","CXCL13")

all_cor <- list(); all_cox <- list(); all_input <- list()
for (nm in names(cohorts)) {
  cfg <- cohorts[[nm]]; eset <- readRDS(file.path(base, cfg$file)); x <- get_probe_expr(eset); pd <- pData(eset)
  dat <- data.frame(sample_id=rownames(pd), CCDC80=as.numeric(x["CCDC80",]), CAF_score=get_score(x,caf_genes), myeloid_score=get_score(x,myeloid_genes), cytotoxicity=get_score(x,cyto_genes), T_exhaustion=get_score(x,exhaust_genes), stringsAsFactors=FALSE)
  event_col <- if (cfg$event %in% names(pd)) cfg$event else grep("survival|dod_event|death", names(pd), ignore.case=TRUE, value=TRUE)[1]
  time_col <- if (cfg$time %in% names(pd)) cfg$time else grep("survival month|time_to_dod|survival", names(pd), ignore.case=TRUE, value=TRUE)[1]
  dat$os_event <- parse_event(pd[[event_col]])
  dat$os_time <- parse_time(pd[[time_col]])
  dat <- dat[complete.cases(dat),,drop=FALSE]
  dat$CCDC80_z <- as.numeric(scale(dat$CCDC80))
  dat$cohort <- nm
  write.csv(dat, file.path(out, paste0(nm,"_CCDC80_immune_input.csv")), row.names=FALSE)
  all_input[[nm]] <- dat
  for (v in c("CAF_score","myeloid_score","cytotoxicity","T_exhaustion")) {
    tt <- suppressWarnings(cor.test(dat$CCDC80, dat[[v]], method="spearman", exact=FALSE))
    all_cor[[length(all_cor)+1]] <- data.frame(cohort=nm, variable=v, n=nrow(dat), rho=unname(tt$estimate), p=tt$p.value)
  }
  for (model_name in c("univariable","immune_adjusted")) {
    rhs <- if (model_name == "univariable") "CCDC80_z" else "CCDC80_z + myeloid_score + cytotoxicity + T_exhaustion"
    fit <- coxph(as.formula(paste0("Surv(os_time,os_event) ~ ",rhs)), data=dat)
    sm <- summary(fit)
    for (i in seq_along(coef(fit))) all_cox[[length(all_cox)+1]] <- data.frame(cohort=nm,model=model_name,variable=names(coef(fit))[i],n=nrow(dat),events=sum(dat$os_event),HR=exp(coef(fit)[i]),CI_low=sm$conf.int[i,"lower .95"],CI_high=sm$conf.int[i,"upper .95"],p=sm$coefficients[i,"Pr(>|z|)"])
  }
}

cor_out <- do.call(rbind, all_cor); cor_out$FDR <- ave(cor_out$p, cor_out$variable, FUN=function(x) p.adjust(x,"BH")); write.csv(cor_out,file.path(out,"CCDC80_external_immune_correlations.csv"),row.names=FALSE)
cox_out <- do.call(rbind, all_cox); cox_out$FDR <- ave(cox_out$p, cox_out$model, FUN=function(x) p.adjust(x,"BH")); write.csv(cox_out,file.path(out,"CCDC80_external_Cox.csv"),row.names=FALSE)
writeLines(c("External validation: GSE32894 and GSE13507",paste0("Cohorts=",paste(names(cohorts),collapse=", ")),"Microenvironment scores are exploratory marker-panel z-score averages; interpretation is at patient level."),file.path(out,"RUN_STATUS.txt"))
