library(data.table)
base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
out <- file.path(base,"analysis_revisions_20260906/CCDC80_meta_GSE293189_GSE135337")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
a <- fread(file.path(base,"analysis_revisions_20260906/GSE293189_patient_level_CCDC80_TME/GSE293189_patient_level_CCDC80_TME.csv"))
b <- fread(file.path(base,"analysis_revisions_20260906/GSE135337_patient_level_CCDC80_TME/GSE135337_patient_level_CCDC80_TME.csv"))
a[,cohort:="GSE293189"]; b[,cohort:="GSE135337"]
keep <- c("patient_id","n_cells","n_CAF","n_Myeloid","n_Tcell","CAF_fraction","Myeloid_fraction","CAF_CCDC80_mean","CAF_CCDC80_detected","CAF_Myeloid_score","Myeloid_score","Exhaustion_score","Cytotoxicity_score","cohort")
z <- rbindlist(list(a[,intersect(keep,names(a)),with=FALSE],b[,intersect(keep,names(b)),with=FALSE]),fill=TRUE)
z <- z[is.finite(CAF_CCDC80_mean) & n_CAF>0]
z[,CCDC80_CAF_z:=as.numeric(scale(CAF_CCDC80_mean)),by=cohort]
z[,CCDC80_high_CAF:=CCDC80_CAF_z>=0]
fwrite(z,file.path(out,"combined_patient_level_CCDC80_TME.csv"))
vars <- c("Myeloid_fraction","Myeloid_score","CAF_Myeloid_score","Exhaustion_score","Cytotoxicity_score")
within <- rbindlist(lapply(vars,function(v){rbindlist(lapply(unique(z$cohort),function(cc){d<-z[cohort==cc];ok<-is.finite(d$CAF_CCDC80_mean)&is.finite(d[[v]]);if(sum(ok)>=3){tt<-suppressWarnings(cor.test(d$CAF_CCDC80_mean[ok],d[[v]][ok],method="spearman",exact=FALSE));data.table(cohort=cc,variable=v,n=sum(ok),rho=unname(tt$estimate),p=tt$p.value)}else data.table(cohort=cc,variable=v,n=sum(ok),rho=NA_real_,p=NA_real_)}))}));within[,FDR:=p.adjust(p,"BH")]
pooled <- rbindlist(lapply(vars,function(v){ok<-is.finite(z$CCDC80_CAF_z)&is.finite(z[[v]]);d<-z[ok];fit<-lm(as.formula(paste(v,"~CCDC80_CAF_z+cohort")),data=d);data.table(variable=v,n=nrow(d),beta=coef(fit)["CCDC80_CAF_z"],p=coef(summary(fit))["CCDC80_CAF_z","Pr(>|t|)"],spearman_rho=cor(d$CCDC80_CAF_z,d[[v]],method="spearman"))}));pooled[,FDR:=p.adjust(p,"BH")]
meta <- rbindlist(lapply(vars,function(v){rr<-within[variable==v & is.finite(rho)];if(nrow(rr)==2){w<-rr$n-3;zz<-atanh(pmax(pmin(rr$rho,.999999),-.999999));est<-sum(w*zz)/sum(w);se<-sqrt(1/sum(w));p<-2*pnorm(-abs(est/se));data.table(variable=v,n_cohorts=2,rho_fisher=tanh(est),p_fisher=p)}else data.table(variable=v,n_cohorts=nrow(rr),rho_fisher=NA_real_,p_fisher=NA_real_)}));meta[,FDR:=p.adjust(p_fisher,"BH")]
highlow <- rbindlist(lapply(vars,function(v){ok<-is.finite(z$CCDC80_high_CAF)&is.finite(z[[v]]);d<-z[ok];tt<-wilcox.test(d[[v]]~d$CCDC80_high_CAF,exact=FALSE);data.table(variable=v,n=nrow(d),high_mean=mean(d[CCDC80_high_CAF==TRUE][[v]],na.rm=TRUE),low_mean=mean(d[CCDC80_high_CAF==FALSE][[v]],na.rm=TRUE),p=tt$p.value)}));highlow[,FDR:=p.adjust(p,"BH")]
fwrite(within,file.path(out,"within_cohort_CCDC80_TME_correlations.csv"));fwrite(pooled,file.path(out,"pooled_cohort_adjusted_regression.csv"));fwrite(meta,file.path(out,"fisher_effect_size_meta_analysis.csv"));fwrite(highlow,file.path(out,"CCDC80_high_low_TME_comparison.csv"));fwrite(z[,.(n_patients=.N,n_CAF=sum(n_CAF),n_Myeloid=sum(n_Myeloid)),by=cohort],file.path(out,"meta_audit.csv"))
writeLines(c("Direct CCDC80 CAF meta-analysis: GSE293189 + GSE135337","Only cohorts with direct CCDC80 expression in CAF were included.","GSE211388 excluded because CCDC80 is absent from its raw matrix.","Pooled regression adjusts for cohort; Fisher analysis combines within-cohort Spearman effects."),file.path(out,"RUN_STATUS.txt"))
