library(data.table)
base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
tar <- file.path(base,"02_external_scRNA_downloads/GSE135337/GSE135337_RAW.tar")
out <- file.path(base,"analysis_revisions_20260906/GSE135337_patient_level_CCDC80_TME")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

panels <- list(
  CAF=c("COL1A1","COL1A2","DCN","LUM","COL3A1","COL5A1","SPARC","FAP","PDGFRA","COL6A1"),
  Myeloid=c("LYZ","LST1","TYROBP","FCER1G","CTSS","LILRB1","CD68","CD163","C1QA","C1QB","C1QC","SPP1"),
  Monocyte=c("LYZ","S100A8","S100A9","FCN1","CTSS","LGALS3"),
  Macrophage=c("CD68","CD163","C1QA","C1QB","C1QC","APOE","MSR1"),
  Tcell=c("CD3D","CD3E","TRBC1","TRBC2","IL7R","LTB"),
  NK=c("NKG7","GNLY","PRF1","GZMB"),
  Exhaustion=c("PDCD1","CTLA4","TIGIT","LAG3","HAVCR2","TOX","CXCL13"),
  Cytotoxicity=c("NKG7","GNLY","PRF1","GZMB","GZMK","CD8A","CD8B")
)
score <- function(m, genes) {
  genes <- intersect(genes, rownames(m))
  if (!length(genes)) return(rep(NA_real_, ncol(m)))
  z <- t(scale(t(m[genes,,drop=FALSE])))
  colMeans(z, na.rm=TRUE)
}
files <- system2("tar", c("-tf", shQuote(tar)), stdout=TRUE)
files <- files[grepl("gene_cell_exprs_table\\.(txt|xls)\\.gz$", files)]
files <- files[grepl("_BC[1-7]_gene_cell_exprs_table", basename(files))]
all_cells <- list(); patient_rows <- list()
for (f in files) {
  message("Reading ", f)
  cmd <- paste("tar -xOf", shQuote(tar), shQuote(f), "| gzip -dc")
  x <- fread(cmd=cmd, data.table=FALSE, check.names=FALSE)
  if (ncol(x)<4) next
  sample_id <- sub("_gene_cell_exprs_table.*$", "", basename(f))
  patient_id <- sub("^GSM[0-9]+_", "", sample_id)
  patient_id <- sub("_gene.*$", "", patient_id)
  gene <- toupper(trimws(x[[2]])); keep <- !duplicated(gene) & nzchar(gene)
  x <- x[keep,,drop=FALSE]; gene <- gene[keep]
  mat <- as.matrix(x[,-c(1,2),drop=FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- gene
  lib <- colSums(mat); norm <- log1p(t(t(mat)/pmax(lib,1)*10000))
  sc <- lapply(panels, function(g) score(norm,g))
  cell <- data.table(cell_id=colnames(mat), sample_id=sample_id, patient_id=patient_id,
    CCDC80=if("CCDC80"%in%rownames(norm)) as.numeric(norm["CCDC80",]) else NA_real_,
    CAF_score=sc$CAF, Myeloid_score=sc$Myeloid, Monocyte_score=sc$Monocyte,
    Macrophage_score=sc$Macrophage, Tcell_score=sc$Tcell, NK_score=sc$NK,
    Exhaustion_score=sc$Exhaustion, Cytotoxicity_score=sc$Cytotoxicity,
    library_size=lib)
  cell[,CAF:=CAF_score>0 & CAF_score>Myeloid_score & CAF_score>Tcell_score]
  cell[,Myeloid:=Myeloid_score>0 & Myeloid_score>CAF_score & Myeloid_score>Tcell_score]
  cell[,T_NK:=Tcell_score>0 & Tcell_score>=CAF_score & Tcell_score>=Myeloid_score]
  cell[,Macrophage:=Myeloid & Macrophage_score>=Monocyte_score]
  cell[,Monocyte:=Myeloid & Monocyte_score>Macrophage_score]
  all_cells[[sample_id]] <- cell
  patient_rows[[sample_id]] <- cell[,.(n_cells=.N,n_CAF=sum(CAF),n_Myeloid=sum(Myeloid),n_Monocyte=sum(Monocyte),n_Macrophage=sum(Macrophage),n_TNK=sum(T_NK),CAF_fraction=mean(CAF),Myeloid_fraction=mean(Myeloid),CAF_CCDC80_mean=mean(CCDC80[CAF],na.rm=TRUE),CAF_CCDC80_detected=mean(CCDC80[CAF]>0,na.rm=TRUE),CAF_Myeloid_score=mean(Myeloid_score[CAF],na.rm=TRUE),Myeloid_score=mean(Myeloid_score[Myeloid],na.rm=TRUE),Exhaustion_score=mean(Exhaustion_score[T_NK],na.rm=TRUE),Cytotoxicity_score=mean(Cytotoxicity_score[T_NK],na.rm=TRUE)),by=.(patient_id,sample_id)]
}
cells <- rbindlist(all_cells, fill=TRUE); patients <- rbindlist(patient_rows, fill=TRUE)
patients[,CCDC80_high_CAF:=CAF_CCDC80_mean>=median(CAF_CCDC80_mean,na.rm=TRUE)]
fwrite(cells,file.path(out,"GSE135337_cell_level_scores.csv")); fwrite(patients,file.path(out,"GSE135337_patient_level_CCDC80_TME.csv"))
vars <- c("Myeloid_fraction","Myeloid_score","n_Monocyte","n_Macrophage","Exhaustion_score","Cytotoxicity_score")
cor_out <- rbindlist(lapply(vars,function(v){ok=is.finite(patients$CAF_CCDC80_mean)&is.finite(patients[[v]]); if(sum(ok)>=3){tt=suppressWarnings(cor.test(patients$CAF_CCDC80_mean[ok],patients[[v]][ok],method="spearman",exact=FALSE));data.table(variable=v,n=sum(ok),rho=unname(tt$estimate),p=tt$p.value)}else data.table(variable=v,n=sum(ok),rho=NA_real_,p=NA_real_)})); cor_out[,FDR:=p.adjust(p,"BH")]
fwrite(cor_out,file.path(out,"GSE135337_CCDC80_TME_correlations.csv")); fwrite(data.table(n_patients=uniqueN(patients$patient_id),n_samples=nrow(patients),n_cells=nrow(cells),n_CAF=sum(cells$CAF),n_Myeloid=sum(cells$Myeloid),CCDC80_detected_cells=sum(cells$CCDC80>0,na.rm=TRUE)),file.path(out,"GSE135337_audit.csv"))
writeLines(c("GSE135337 patient-level CCDC80 CAF-TME analysis","Cell labels are exploratory marker-score labels, not submitter annotations.","Patient-level statistics; CAF CCDC80 high/low split at the patient median.","Cancer-adjacent sample BCN is excluded by sample list (BC1-BC7 only)."),file.path(out,"RUN_STATUS.txt"))
