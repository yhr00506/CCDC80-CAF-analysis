library(survival)
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/analysis_revisions_20260906"
out <- file.path(root, "ECM_CAF_signature_validation")
d <- read.csv(file.path(out, "TCGA_BLCA_ECM_signature_input.csv"), stringsAsFactors=FALSE)
d$stage_simple <- sub("^Stage ", "", d$stage)
d$stage_simple[d$stage_simple %in% c("I","II")] <- "I_II"
d$stage_simple[d$stage_simple %in% c("III","IV")] <- "III_IV"
clin <- read.csv("/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/00_data_audit/TCGA_clinical.csv", check.names=FALSE, stringsAsFactors=FALSE)
m <- match(substr(d$sample_id,1,12), clin$submitter_id)
d$age <- suppressWarnings(as.numeric(clin$age_at_diagnosis[m])) / 365.25
d <- d[complete.cases(d[,c("ECM_z","CAF_score","myeloid_score","os_time","os_event","stage_simple","age")]),,drop=FALSE]
models <- list(univariable="ECM_z", stage_adjusted="ECM_z + stage_simple", stage_age_adjusted="ECM_z + stage_simple + age", stage_CAF_adjusted="ECM_z + stage_simple + CAF_score", full="ECM_z + stage_simple + age + CAF_score + myeloid_score")
rows <- list(); ph <- list()
for (nm in names(models)) {
 fit <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ ", models[[nm]])), data=d)
 sm <- summary(fit)
 for (i in seq_along(coef(fit))) rows[[length(rows)+1]] <- data.frame(model=nm, term=names(coef(fit))[i], n=nrow(d), events=sum(d$os_event), HR=exp(coef(fit)[i]), CI_low=sm$conf.int[i,"lower .95"], CI_high=sm$conf.int[i,"upper .95"], p=sm$coefficients[i,"Pr(>|z|)"])
 z <- cox.zph(fit); ph[[length(ph)+1]] <- data.frame(model=nm, term=rownames(z$table), p=z$table[,"p"])
}
res <- do.call(rbind, rows); res$FDR <- ave(res$p, res$model, FUN=function(x) p.adjust(x,"BH")); write.csv(res,file.path(out,"TCGA_ECM_signature_multivariable_Cox.csv"),row.names=FALSE)
write.csv(do.call(rbind,ph),file.path(out,"TCGA_ECM_signature_proportional_hazards.csv"),row.names=FALSE)
write.csv(d,file.path(out,"TCGA_ECM_signature_multivariable_input.csv"),row.names=FALSE)
writeLines(c(paste0("N=",nrow(d),"; events=",sum(d$os_event)),"Models: univariable, stage-adjusted, stage+age, stage+CAF, full.","ECM signature is the mean z-score of COL1A1, COL3A1, FN1."),file.path(out,"MULTIVARIABLE_RUN_STATUS.txt"))
print(res[res$term=="ECM_z",])
