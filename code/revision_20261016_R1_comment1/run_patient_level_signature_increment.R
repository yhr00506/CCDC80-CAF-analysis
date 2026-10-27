#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(SummarizedExperiment); library(survival)})
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
out <- file.path(root, "revision_20261016_R1_comment1")
assets <- file.path(out, "assets")
resdir <- file.path(out, "patient_level_analysis")
dir.create(resdir, recursive=TRUE, showWarnings=FALSE)
expr_file <- file.path(root, "00_data_audit/TCGA_expression.rds")
clin_file <- file.path(root, "00_data_audit/TCGA_clinical.csv")
mat_file <- file.path(assets, "MatrisomePortal_Homo_sapiens_masterlist_rev2014.csv")
msig_file <- file.path(assets, "msigdb.2025.1.Hs.rds")
program <- read.csv(file.path(assets,"CCDC80_associated_19_gene_program.csv"),check.names=FALSE)$gene
caf <- c("COL1A1","COL1A2","DCN","LUM","COL3A1","COL5A1","SPARC","FAP","PDGFRA")
ecm3 <- c("COL1A1","COL3A1","FN1")
stopifnot(length(program)==19, !anyDuplicated(program), file.exists(msig_file))

# Gene-set source extraction from the package-pinned MSigDB 2025.1 human RDS.
ms <- readRDS(msig_file)
stopifnot(all(c("gs_name","db_gene_symbol") %in% names(ms)))
tgf <- sort(unique(toupper(ms$db_gene_symbol[ms$gs_name=="HALLMARK_TGF_BETA_SIGNALING"])))
stopifnot(length(tgf)>0)
writeLines(tgf,file.path(assets,"MSigDB_2025.1_HALLMARK_TGF_BETA_SIGNALING_symbols.txt"))
rm(ms); gc()
mat <- read.csv(mat_file,skip=1,check.names=FALSE,stringsAsFactors=FALSE)
mat <- mat[!is.na(mat$`Gene Symbol`),]
mat$`Gene Symbol` <- toupper(trimws(mat$`Gene Symbol`))
matrisome <- sort(unique(mat$`Gene Symbol`))
collagen <- sort(unique(mat$`Gene Symbol`[mat$`Matrisome Category`=="Collagens"]))
sets <- list(CCDC80_program=program, CAF_9gene=caf, ECM_3gene=ecm3,
             collagen_core=collagen, matrisome_2014=matrisome,
             Hallmark_TGF_Beta_2025.1=tgf)

# Read the project TCGA object; retain primary tumor aliquots (TCGA sample type 01).
se <- readRDS(expr_file)
assay_names <- assayNames(se)
if (!"unstranded" %in% assay_names) stop("Expected unstranded count assay is absent")
barcodes <- colnames(se)
pat <- substr(barcodes,1,12)
sample_type <- substr(barcodes,14,15)
keep <- sample_type=="01" & !is.na(pat) & pat!=""
if (!any(keep)) stop("No TCGA primary tumor aliquots (sample type 01) found")
rd <- as.data.frame(rowData(se))
if (!"gene_name" %in% names(rd)) stop("rowData gene_name field missing")
symbol <- toupper(trimws(as.character(rd$gene_name)))
valid <- !is.na(symbol) & symbol!="" & !duplicated(symbol)
all_genes <- unique(unlist(sets,use.names=FALSE))
sel <- which(valid & symbol %in% all_genes)
counts <- as.matrix(assay(se,"unstranded")[sel,keep,drop=FALSE])
rownames(counts) <- symbol[sel]
# Match the existing project transform, then average aliquot-level log2(count+1) values within patient.
log_expr <- log2(counts+1)
patient_ids <- pat[keep]
patient_order <- unique(patient_ids)
expr_patient <- vapply(patient_order,function(id) rowMeans(log_expr[,patient_ids==id,drop=FALSE]),numeric(nrow(log_expr)))
if (is.null(dim(expr_patient))) expr_patient <- matrix(expr_patient,ncol=1,dimnames=list(rownames(log_expr),patient_order))
rownames(expr_patient) <- rownames(log_expr); colnames(expr_patient) <- patient_order

# Per-signature score = mean gene-wise z score over patients; missing genes are listed and scores require >=80% coverage.
score_rows <- list(); coverage <- list()
for (nm in names(sets)) {
  g <- intersect(sets[[nm]],rownames(expr_patient))
  coverage[[nm]] <- data.frame(signature=nm,n_defined=length(unique(sets[[nm]])),n_detected=length(g),fraction_detected=length(g)/length(unique(sets[[nm]])))
  if (length(g)<2 || length(g)/length(unique(sets[[nm]]))<0.8) stop(paste("Insufficient coverage:",nm,length(g),"of",length(sets[[nm]])))
  z <- t(scale(t(expr_patient[g,,drop=FALSE])))
  score_rows[[nm]] <- colMeans(z,na.rm=TRUE)
}
scores <- as.data.frame(score_rows,check.names=FALSE)
scores$patient_id <- rownames(scores)
coverage <- do.call(rbind,coverage)
write.csv(coverage,file.path(resdir,"signature_gene_coverage.csv"),row.names=FALSE)

clin <- read.csv(clin_file,check.names=FALSE,stringsAsFactors=FALSE)
clin$patient_id <- substr(as.character(clin$submitter_id),1,12)
clin <- clin[!duplicated(clin$patient_id),]
clin$os_event <- as.integer(tolower(trimws(clin$vital_status)) %in% c("dead","deceased"))
days_death <- suppressWarnings(as.numeric(clin$days_to_death))
days_follow <- suppressWarnings(as.numeric(clin$days_to_last_follow_up))
clin$os_time <- ifelse(clin$os_event==1,days_death,days_follow)
clin$age <- suppressWarnings(as.numeric(clin$age_at_index))
clin$sex <- factor(tolower(trimws(clin$sex_at_birth)))
st <- trimws(as.character(clin$ajcc_pathologic_stage))
st <- sub("^Stage\\s+","",st,ignore.case=TRUE)
clin$stage_simple <- ifelse(st %in% c("I","II"),"I_II",ifelse(st %in% c("III","IV"),"III_IV",NA_character_))
clin$stage_simple <- factor(clin$stage_simple,levels=c("I_II","III_IV"))
d <- merge(scores,clin[,c("patient_id","os_time","os_event","age","sex","stage_simple","ajcc_pathologic_stage")],by="patient_id")
d$os_time <- d$os_time/365.25
# Apply a single complete-case set for all primary nested model comparisons.
model_vars <- c("os_time","os_event","age","sex","stage_simple",names(sets))
d <- d[complete.cases(d[,model_vars]),]
if(nrow(d)<100 || sum(d$os_event)<30) stop("Unexpectedly small complete-case analysis set")
d$age_z <- as.numeric(scale(d$age))
for (nm in names(sets)) d[[paste0(nm,"_z")]] <- as.numeric(scale(d[[nm]]))
write.csv(d,file.path(resdir,"TCGA_patient_level_scores_and_clinical.csv"),row.names=FALSE)

base_formula <- as.formula("survival::Surv(os_time,os_event) ~ stage_simple + age_z + sex")
fit_safe <- function(f,d) coxph(f,data=d,ties="efron",x=TRUE,y=TRUE,model=TRUE)
base_fit <- fit_safe(base_formula,d)
score_terms <- paste0(names(sets),"_z")
extract_fit <- function(f,name,comparison,fitdata=d) {
  s <- summary(f); co <- as.data.frame(s$coefficients); ci <- as.data.frame(s$conf.int)
  terms <- rownames(co)
  data.frame(model=name,comparison=comparison,term=terms,HR=ci[,"exp(coef)"],
             CI_low=ci[,"lower .95"],CI_high=ci[,"upper .95"],p=co[,"Pr(>|z|)"],
             n=nrow(fitdata),events=sum(fitdata$os_event),AIC=AIC(f),
             logLik=as.numeric(logLik(f)),C_index=unname(s$concordance[1]),
             C_index_se=unname(s$concordance[2]),stringsAsFactors=FALSE)
}
model_records <- list(); lrt_records <- list(); ph_records <- list()
model_records$clinical_only <- extract_fit(base_fit,"clinical_only","clinical_covariates",d)
f_program <- fit_safe(update(base_formula,paste("~.+", "CCDC80_program_z")),d)
model_records$clinical_plus_program <- extract_fit(f_program,"clinical_plus_CCDC80_program","program_added_to_clinical",d)
lrt_program <- anova(base_fit,f_program,test="LRT")
lrt_records$CCDC80_vs_clinical <- data.frame(benchmark="clinical_covariates",n=nrow(d),events=sum(d$os_event),LR=as.numeric(lrt_program$Chisq[2]),df=as.numeric(lrt_program$Df[2]),p=as.numeric(lrt_program$`Pr(>|Chi|)`[2]))
zph <- cox.zph(f_program,transform="km")
ph_records$program_only <- data.frame(model="clinical_plus_CCDC80_program",term=rownames(zph$table),chisq=zph$table[,"chisq"],df=zph$table[,"df"],p=zph$table[,"p"])
for (nm in setdiff(names(sets),"CCDC80_program")) {
  term <- paste0(nm,"_z")
  f0 <- fit_safe(update(base_formula,paste(".~.+",term)),d)
  f1 <- fit_safe(update(formula(f0),paste(".~.+CCDC80_program_z")),d)
  model_records[[paste0(nm,"_base")]] <- extract_fit(f0,paste0("clinical_plus_",nm),"benchmark",d)
  model_records[[paste0(nm,"_plus_program")]] <- extract_fit(f1,paste0("clinical_plus_",nm,"_plus_program"),"benchmark_plus_program",d)
  lrt <- anova(f0,f1,test="LRT")
  lrt_records[[nm]] <- data.frame(benchmark=nm,n=nrow(d),events=sum(d$os_event),
      LR=as.numeric(lrt$Chisq[2]),df=as.numeric(lrt$Df[2]),p=as.numeric(lrt$`Pr(>|Chi|)`[2]))
  zph <- cox.zph(f1,transform="km")
  ph_records[[nm]] <- data.frame(model=paste0(nm,"_plus_program"),term=rownames(zph$table),
                                  chisq=zph$table[,"chisq"],df=zph$table[,"df"],p=zph$table[,"p"])
}
# Combined model tests added information conditional on all comparator scores.
full_formula <- as.formula(paste("survival::Surv(os_time,os_event) ~ stage_simple + age_z + sex +",
                                  paste(c(score_terms,"CCDC80_program_z"),collapse=" + ")))
full_fit <- fit_safe(full_formula,d)
without_program_formula <- update(full_formula,.~.-CCDC80_program_z)
without_program_fit <- fit_safe(without_program_formula,d)
model_records$all_benchmarks <- extract_fit(without_program_fit,"clinical_plus_all_benchmarks","all_benchmarks",d)
model_records$all_benchmarks_plus_program <- extract_fit(full_fit,"clinical_plus_all_benchmarks_plus_program","all_benchmarks_plus_program",d)
lrt <- anova(without_program_fit,full_fit,test="LRT")
lrt_records$all_benchmarks <- data.frame(benchmark="all_benchmarks",n=nrow(d),events=sum(d$os_event),LR=as.numeric(lrt$Chisq[2]),df=as.numeric(lrt$Df[2]),p=as.numeric(lrt$`Pr(>|Chi|)`[2]))
zph <- cox.zph(full_fit,transform="km")
ph_records$all_benchmarks <- data.frame(model="all_benchmarks_plus_program",term=rownames(zph$table),chisq=zph$table[,"chisq"],df=zph$table[,"df"],p=zph$table[,"p"])

# Score correlations, variance inflation diagnostics, and influence/deviance summaries.
score_names <- names(sets)
cor_mat <- cor(d[,score_names],method="spearman",use="pairwise.complete.obs")
write.csv(cor_mat,file.path(resdir,"signature_score_spearman_matrix.csv"))
vifs <- data.frame(signature=score_names,VIF=NA_real_)
for (i in seq_along(score_names)) {
  y <- d[[score_names[i]]]; xx <- d[,setdiff(score_names,score_names[i]),drop=FALSE]
  vifs$VIF[i] <- 1/(1-summary(lm(y~.,data=xx))$r.squared)
}
write.csv(vifs,file.path(resdir,"signature_score_VIF.csv"),row.names=FALSE)

a <- anova(full_fit)
write.csv(do.call(rbind,model_records),file.path(resdir,"Cox_model_coefficients_per_SD.csv"),row.names=FALSE)
lrts <- do.call(rbind,lrt_records); lrts$FDR_BH <- p.adjust(lrts$p,method="BH")
write.csv(lrts,file.path(resdir,"CCDC80_program_incremental_LRT.csv"),row.names=FALSE)
write.csv(do.call(rbind,ph_records),file.path(resdir,"PH_assumption_checks_cox_zph.csv"),row.names=FALSE)
write.csv(data.frame(term=rownames(a),chisq=a$Chisq,df=a$Df,p=a$`Pr(>|Chi|)`),file.path(resdir,"full_model_term_tests.csv"),row.names=FALSE)
write.csv(data.frame(patient_n=nrow(d),events=sum(d$os_event),aliquots_used=sum(keep),unique_patients_expression=length(unique(patient_ids)),clinical_rows=nrow(clin),
  molecular_subtype_available=FALSE,tumor_purity_available=FALSE,analysis_unit="unique TCGA patient",expression="log2(unstranded counts + 1); mean across primary tumor aliquots per patient",score="mean gene-wise z-score; z-scores calculated across unique patients",survival_time="years",stage="pathologic stage I-II vs III-IV",covariates="stage, age, sex",ph_checks="cox.zph KM transform",R_version=R.version.string,survival_version=as.character(packageVersion("survival")),msigdbr_version=as.character(packageVersion("msigdbr")),MSigDB="2025.1.Hs",stringsAsFactors=FALSE),file.path(resdir,"analysis_metadata.csv"),row.names=FALSE)
writeLines(c("Patient-level incremental signature analysis completed.",paste0("N patients=",nrow(d),"; events=",sum(d$os_event)),"Primary tumor aliquots (sample type 01); duplicate aliquots averaged after log2(count+1).","Scores are equal-weight mean gene-wise z scores, standardized across unique patients.","Nested Cox models adjusted for stage, age, and sex; reported increment is an association, not causality.","Tumor purity and TCGA molecular subtype fields were not present in the audited local TCGA clinical inputs.","See CSV tables for coverage, coefficients, LRT, PH checks, C-index/AIC, score correlation, and VIF."),file.path(resdir,"RUN_STATUS.txt"))
cat(readLines(file.path(resdir,"RUN_STATUS.txt")),sep="\n")
cat("Signature gene coverage:\n");print(coverage)
cat("Incremental LRT:\n");print(lrts)
