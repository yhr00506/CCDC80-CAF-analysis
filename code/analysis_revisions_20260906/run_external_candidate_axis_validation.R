library(Biobase)
library(survival)
base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
out <- file.path(base,"analysis_revisions_20260906","external_candidate_axis_validation")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cohorts <- list(
 GSE32894=list(file="00_data_audit/GSE32894_expression.rds", event="dod_event_(yes/no):ch1", time="time_to_dod_(months):ch1"),
 GSE13507=list(file="00_data_audit/GSE13507_expression.rds", event="overall survival:ch1", time="survival month:ch1"))
get_expr <- function(eset){
 x<-exprs(eset); fd<-as.data.frame(fData(eset)); sc<-intersect(c("Symbol","ILMN_Gene"),names(fd))[1]
 s<-toupper(trimws(as.character(fd[[sc]]))); ok<-!is.na(s)&s!=""&s!="---"; x<-x[ok,,drop=FALSE]; s<-s[ok]
 sm<-rowsum(x,s,reorder=FALSE); sm/as.numeric(table(s)[rownames(sm)])
}
score <- function(x,g){g<-intersect(g,rownames(x)); if(!length(g)) return(rep(NA_real_,ncol(x))); colMeans(t(scale(t(x[g,,drop=FALSE]))),na.rm=TRUE)}
pe <- function(x){y<-tolower(trimws(as.character(x))); as.integer(y%in%c("yes","dead","death","deceased","1","true","progression"))}
pt <- function(x){y<-suppressWarnings(as.numeric(gsub("[^0-9.\\-]","",as.character(x)))); y[y<=0]<-NA_real_; y}
caf<-c("COL1A1","COL1A2","DCN","LUM","COL3A1","COL5A1","SPARC","FAP","PDGFRA")
myeloid<-c("CD68","CD163","C1QA","C1QB","C1QC","LILRB1","FCER1G","TYROBP","APOC1","SPP1","MSR1")
genes<-c("CSF1","COL1A1","IL6","CCL5","TGFB1","CXCL12","PDGFA","COL3A1","FN1")
all<-list(); cors<-list(); coxs<-list()
for(nm in names(cohorts)){
 cfg<-cohorts[[nm]]; e<-readRDS(file.path(base,cfg$file)); x<-get_expr(e); pd<-pData(e)
 d<-data.frame(sample_id=colnames(x),CAF_score=score(x,caf),myeloid_score=score(x,myeloid),stringsAsFactors=FALSE)
 for(g in genes) d[[g]]<-if(g%in%rownames(x)) as.numeric(x[g,]) else NA_real_
 ec<-if(cfg$event%in%names(pd))cfg$event else grep("survival|dod_event|death",names(pd),ignore.case=TRUE,value=TRUE)[1]
 tc<-if(cfg$time%in%names(pd))cfg$time else grep("survival month|time_to_dod|survival",names(pd),ignore.case=TRUE,value=TRUE)[1]
 d$os_event<-pe(pd[[ec]]); d$os_time<-pt(pd[[tc]]); d<-d[complete.cases(d[,c("os_event","os_time","CAF_score","myeloid_score")]),,drop=FALSE]
 for(g in genes){ if(all(is.na(d[[g]])))next; d[[paste0(g,"_z")]]<-as.numeric(scale(d[[g]]));
  for(v in c("CAF_score","myeloid_score")){z<-complete.cases(d[,c(g,v)]); tt<-suppressWarnings(cor.test(d[[g]][z],d[[v]][z],method="spearman",exact=FALSE)); cors[[length(cors)+1]]<-data.frame(cohort=nm,gene=g,variable=v,n=sum(z),rho=unname(tt$estimate),p=tt$p.value)}
  z<-complete.cases(d[,c(paste0(g,"_z"),"os_time","os_event")]); if(sum(z)>20&&sum(d$os_event[z])>5){fit<-coxph(as.formula(paste0("Surv(os_time,os_event)~",g,"_z")),data=d[z,]); sm<-summary(fit); coxs[[length(coxs)+1]]<-data.frame(cohort=nm,gene=g,n=sum(z),events=sum(d$os_event[z]),HR=exp(coef(fit)),lo=sm$conf.int[1,"lower .95"],hi=sm$conf.int[1,"upper .95"],p=sm$coefficients[1,"Pr(>|z|)"])}
 }
 d$cohort<-nm; all[[nm]]<-d; write.csv(d,file.path(out,paste0(nm,"_candidate_axis_input.csv")),row.names=FALSE)
}
co<-do.call(rbind,cors); co$FDR<-ave(co$p,co$variable,FUN=function(x)p.adjust(x,"BH")); write.csv(co,file.path(out,"candidate_gene_microenvironment_correlations.csv"),row.names=FALSE)
cx<-do.call(rbind,coxs); cx$FDR<-p.adjust(cx$p,"BH"); write.csv(cx,file.path(out,"candidate_gene_OS_Cox.csv"),row.names=FALSE)
writeLines(c("External candidate-axis validation completed","Cohorts: GSE32894 and GSE13507","Candidate genes: CSF1, COL1A1, IL6, CCL5, TGFB1, CXCL12, PDGFA, COL3A1, FN1","All analyses are patient-level; exploratory external validation."),file.path(out,"RUN_STATUS.txt"))
print(co); print(cx)
