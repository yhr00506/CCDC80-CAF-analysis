library(data.table)
base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
lrdir <- file.path(base,"analysis_revisions_20260906/GSE135337_CAF_myeloid_LR")
out <- file.path(base,"analysis_revisions_20260906/GSE135337_CCDC80_LR_association")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
lr <- fread(file.path(lrdir,"GSE135337_patient_LR_potential.csv"))
p <- fread(file.path(base,"analysis_revisions_20260906/GSE135337_patient_level_CCDC80_TME/GSE135337_patient_level_CCDC80_TME.csv"),select=c("patient_id","CAF_CCDC80_mean","CCDC80_high_CAF"))
lr <- merge(lr,p,by="patient_id",all=FALSE); lr[,high:=CAF_CCDC80_mean>=median(CAF_CCDC80_mean,na.rm=TRUE)]
corres <- lr[, {ok<-is.finite(LR_potential)&is.finite(CAF_CCDC80_mean); if(sum(ok)>=4){tt<-suppressWarnings(cor.test(LR_potential[ok],CAF_CCDC80_mean[ok],method="spearman",exact=FALSE)); list(n=sum(ok),rho=unname(tt$estimate),p=tt$p.value)}else list(n=sum(ok),rho=NA_real_,p=NA_real_)}, by=.(pair,ligand,receptor)]
corres[,FDR:=p.adjust(p,"BH")]
groups <- lr[, {x<-LR_potential[high==TRUE & is.finite(LR_potential)];y<-LR_potential[high==FALSE & is.finite(LR_potential)];pp<-if(length(x)>=2 & length(y)>=2) suppressWarnings(wilcox.test(x,y,exact=FALSE)$p.value) else NA_real_;list(n_high=length(x),n_low=length(y),high_mean=mean(x,na.rm=TRUE),low_mean=mean(y,na.rm=TRUE),p=pp)},by=.(pair,ligand,receptor)]
groups[,FDR:=p.adjust(p,"BH")]
ranked <- merge(corres,groups,by=c("pair","ligand","receptor")); ranked[,FDR:=FDR.x]; ranked[,priority_score:= -log10(pmax(FDR,1e-300))*abs(rho)]; ranked<-ranked[order(-priority_score)]
fwrite(lr,file.path(out,"patient_level_CCDC80_LR_merged.csv")); fwrite(ranked,file.path(out,"CCDC80_associated_LR_ranked.csv"))
writeLines(c("CCDC80-associated CAF-to-myeloid ligand-receptor analysis","LR potential = CAF ligand mean x myeloid receptor mean.","Association tested with continuous CAF CCDC80 and high/low comparison.","Scores are prioritization evidence, not causal proof."),file.path(out,"RUN_STATUS.txt"))
