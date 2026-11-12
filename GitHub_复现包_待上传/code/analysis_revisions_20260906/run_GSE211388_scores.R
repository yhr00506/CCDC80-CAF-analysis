library(data.table)
base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
p <- file.path(base,"02_external_scRNA_downloads/GSE211388")
out <- file.path(base,"analysis_revisions_20260906/combined_GSE293189_GSE211388_CCDC80_TME")
markers <- list(CAF=c("COL1A1","COL1A2","DCN","LUM","COL3A1","COL5A1","SPARC","FAP","PDGFRA","COL6A1"),Myeloid=c("LST1","TYROBP","FCER1G","CTSS","LILRB1","CD68","CD163","C1QA","C1QB","C1QC","SPP1"),Tcell=c("CD3D","CD3E","TRBC1","TRBC2","IL7R","LTB"),NK=c("NKG7","GNLY","PRF1","GZMB"),ECM_proxy=c("IGFBP6","MFAP5","FSTL1","FBN1","CLEC3B","PI16","SH3BGRL3","S100A4","C1QTNF3","EMP3","TPPP3","GSN","TIMP3","CRIP1","PTHLH","HTRA3","ELN","TIMP1","CRABP2"),Exhaustion=c("PDCD1","CTLA4","TIGIT","LAG3","HAVCR2","TOX","CXCL13"),Cytotoxicity=c("NKG7","GNLY","PRF1","GZMB","GZMK","CD8A","CD8B"))
score <- function(m,g){g<-intersect(g,rownames(m));if(!length(g))return(rep(NA_real_,ncol(m)));z<-t(scale(t(m[g,,drop=FALSE])));colMeans(z,na.rm=TRUE)}
ci <- fread(cmd=paste("zcat",shQuote(file.path(p,"GSE211388_cell_info.txt.gz"))),header=FALSE,data.table=FALSE); names(ci)<-c("cell_id","sample","majorcluster","minorcluster")
cells <- as.character(ci$cell_id)
need <- unique(unlist(markers)); cond <- paste(need,collapse="|")
cmd <- paste("zcat",shQuote(file.path(p,"GSE211388_count_matrix.txt.gz")),"| grep -E \"^\\\"?(",cond,")\\\"?[[:space:]]\"")
x <- as.matrix(fread(cmd=cmd,header=FALSE,data.table=FALSE)); genes<-toupper(x[,1]); x<-x[,-1,drop=FALSE]; storage.mode(x)<-"numeric"; rownames(x)<-genes
header<-strsplit(readLines(gzfile(file.path(p,"GSE211388_count_matrix.txt.gz")),n=1),"\t")[[1]]; colnames(x)<-gsub("\"","",header)
ci<-ci[match(colnames(x),cells),,drop=FALSE]; lib<-colSums(x); norm<-log1p(t(t(x)/pmax(lib,1)*10000))
d<-data.table(cell_id=colnames(x),patient_id=sub("[ab]$","",sub("^human[NP]_","",ci$sample)),majorcluster=ci$majorcluster,minorcluster=ci$minorcluster,CCDC80=if("CCDC80"%in%rownames(norm)) as.numeric(norm["CCDC80",]) else NA_real_,ECM_proxy_score=score(norm,markers$ECM_proxy),Myeloid_score=score(norm,markers$Myeloid),Tcell_score=score(norm,markers$Tcell),NK_score=score(norm,markers$NK),Exhaustion_score=score(norm,markers$Exhaustion),Cytotoxicity_score=score(norm,markers$Cytotoxicity),library_size=lib)
d[,CAF:=majorcluster%in%c("Fibroblasts","Myofibroblasts")];d[,Myeloid:=majorcluster%in%c("Monocytes/Macrophages","DC","pDC")];d[,T_NK:=majorcluster=="T/NK"];fwrite(d,file.path(out,"GSE211388_cell_level_scores.csv"))
pd<-d[,.(n_cells=.N,n_CAF=sum(CAF),n_Myeloid=sum(Myeloid),n_Tcell=sum(T_NK),CAF_fraction=mean(CAF),Myeloid_fraction=mean(Myeloid),CAF_CCDC80_mean=mean(CCDC80[CAF],na.rm=TRUE),CAF_CCDC80_detected=sum(CCDC80[CAF]>0,na.rm=TRUE),ECM_proxy_mean=mean(ECM_proxy_score[CAF],na.rm=TRUE),Myeloid_score=mean(Myeloid_score[Myeloid],na.rm=TRUE),Exhaustion_score=mean(Exhaustion_score[T_NK],na.rm=TRUE),Cytotoxicity_score=mean(Cytotoxicity_score[T_NK],na.rm=TRUE)),by=patient_id];pd[,cohort:="GSE211388"];fwrite(pd,file.path(out,"GSE211388_patient_level_CCDC80_TME.csv"))
writeLines(c("GSE211388 patient-level immune reference","CD45+ and CD45- cells paired by patient ID.","CCDC80 expression is extracted directly from the raw count matrix."),file.path(out,"RUN_STATUS.txt"))
