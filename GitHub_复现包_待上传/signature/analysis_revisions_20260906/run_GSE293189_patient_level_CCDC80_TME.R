library(data.table)

base <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
rawdir <- file.path(base, "02_external_scRNA_downloads/GSE293189/raw_files")
metafile <- file.path(base, "02_external_scRNA_downloads/metadata/GSE293189_family.soft.gz")
out <- file.path(base, "analysis_revisions_20260906/GSE293189_patient_level_CCDC80_TME")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

parse_soft <- function(path) {
  z <- readLines(gzfile(path), warn=FALSE)
  ix <- which(grepl("^\\^SAMPLE = ", z))
  ans <- vector("list", length(ix))
  for (j in seq_along(ix)) {
    a <- ix[j]; b <- if (j < length(ix)) ix[j+1]-1 else length(z)
    block <- z[a:b]
    gsm <- sub("^\\^SAMPLE = ", "", block[1])
    ttl <- sub("^!Sample_title = ", "", block[grep("^!Sample_title = ", block)[1]])
    patient <- sub("^.*Patient ([^ ]+).*$", "\\1", ttl)
    tissue <- if (grepl("Tumor", ttl, ignore.case=TRUE)) "Tumor" else if (grepl("Normal", ttl, ignore.case=TRUE)) "Normal" else "Other"
    ans[[j]] <- data.frame(gsm=gsm,title=ttl,patient_id=patient,tissue=tissue,stringsAsFactors=FALSE)
  }
  do.call(rbind, ans)
}

meta <- parse_soft(metafile)
files <- list.files(rawdir, pattern="_dge\\.txt\\.gz$", full.names=TRUE)
meta$file <- files[match(meta$gsm, sub("_.*$", "", basename(files)))]
meta <- meta[!is.na(meta$file),,drop=FALSE]

panels <- list(
  CAF=c("COL1A1","COL1A2","DCN","LUM","COL3A1","COL5A1","SPARC","FAP","PDGFRA","COL6A1"),
  Myeloid=c("LST1","TYROBP","FCER1G","CTSS","LILRB1","CD68","CD163","C1QA","C1QB","C1QC","SPP1"),
  Tcell=c("CD3D","CD3E","TRBC1","TRBC2","IL7R","LTB","MALAT1"),
  NK=c("NKG7","GNLY","PRF1","GZMB"),
  Exhaustion=c("PDCD1","CTLA4","TIGIT","LAG3","HAVCR2","TOX","CXCL13"),
  Cytotoxicity=c("NKG7","GNLY","PRF1","GZMB","GZMK","CD8A","CD8B")
)
panels$Tcell <- setdiff(panels$Tcell,"MALAT1")

score_panel <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (!length(genes)) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[genes,,drop=FALSE])))
  colMeans(z, na.rm=TRUE)
}

one_sample <- function(row) {
  f <- row$file
  needed <- unique(c("CCDC80", unlist(panels)))
  cond <- paste(sprintf("$1==\"%s\"", needed), collapse=" || ")
  cmd <- paste("zcat", shQuote(f), "| awk 'NR==1 ||", cond, "'")
  x <- as.matrix(fread(cmd=cmd, data.table=FALSE, check.names=FALSE))
  genes <- toupper(x[,1]); x <- x[,-1,drop=FALSE]
  keep <- !duplicated(genes) & genes != ""
  x <- x[keep,,drop=FALSE]; genes <- genes[keep]
  storage.mode(x) <- "numeric"; rownames(x) <- genes
  lib <- colSums(x)
  norm <- log1p(t(t(x) / pmax(lib,1) * 10000))
  vals <- list(CCDC80=if ("CCDC80" %in% rownames(norm)) as.numeric(norm["CCDC80",]) else rep(NA_real_,ncol(norm)))
  for (nm in names(panels)) vals[[paste0(nm,"_score")]] <- score_panel(norm,panels[[nm]])
  lineage <- cbind(CAF=vals$CAF_score, Myeloid=vals$Myeloid_score, Tcell=vals$Tcell_score, NK=vals$NK_score)
  lineage[!is.finite(lineage)] <- -Inf
  labels <- colnames(lineage)[max.col(lineage, ties.method="first")]
  labels[apply(lineage,1,max) <= 0] <- "Other"
  # Keep a minimum lineage signal and assign each cell to one dominant compartment.
  labels[!is.finite(apply(lineage,1,max))] <- "Other"
  d <- data.frame(cell_id=colnames(norm),label=labels,CCDC80=vals$CCDC80,CAF_score=vals$CAF_score,Myeloid_score=vals$Myeloid_score,Tcell_score=vals$Tcell_score,NK_score=vals$NK_score,Exhaustion_score=vals$Exhaustion_score,Cytotoxicity_score=vals$Cytotoxicity_score,library_size=lib,stringsAsFactors=FALSE)
  d$gsm <- row$gsm; d$title <- row$title; d$patient_id <- row$patient_id; d$tissue <- row$tissue
  d
}

cell_results <- rbindlist(lapply(seq_len(nrow(meta)), function(i) { message(i,"/",nrow(meta)," ",meta$gsm[i]); one_sample(meta[i,]) }), fill=TRUE)
fwrite(cell_results, file.path(out,"GSE293189_cell_level_scores.csv"))

sample_results <- cell_results[, .(
  n_cells=.N,
  n_CAF=sum(label=="CAF"), n_Myeloid=sum(label=="Myeloid"), n_Tcell=sum(label=="Tcell"), n_NK=sum(label=="NK"),
  CAF_fraction=mean(label=="CAF"), Myeloid_fraction=mean(label=="Myeloid"),
  CAF_CCDC80_mean=mean(CCDC80[label=="CAF"],na.rm=TRUE), CAF_CCDC80_detected=mean(CCDC80[label=="CAF"]>0,na.rm=TRUE),
  CAF_Myeloid_score=mean(Myeloid_score[label=="CAF"],na.rm=TRUE), Myeloid_score=mean(Myeloid_score[label=="Myeloid"],na.rm=TRUE),
  Exhaustion_score=mean(Exhaustion_score[label %in% c("Tcell","NK")],na.rm=TRUE), Cytotoxicity_score=mean(Cytotoxicity_score[label %in% c("Tcell","NK")],na.rm=TRUE)
), by=.(gsm,title,patient_id,tissue)]
fwrite(sample_results,file.path(out,"GSE293189_sample_level_CCDC80_TME.csv"))

tum <- sample_results[tissue=="Tumor"]
patient_results <- tum[, .(
  n_libraries=.N, n_cells=sum(n_cells), n_CAF=sum(n_CAF), n_Myeloid=sum(n_Myeloid), n_Tcell=sum(n_Tcell), n_NK=sum(n_NK),
  CAF_fraction=sum(n_CAF)/sum(n_cells), Myeloid_fraction=sum(n_Myeloid)/sum(n_cells),
  CAF_CCDC80_mean=weighted.mean(CAF_CCDC80_mean,pmax(n_CAF,1),na.rm=TRUE), CAF_CCDC80_detected=weighted.mean(CAF_CCDC80_detected,pmax(n_CAF,1),na.rm=TRUE),
  CAF_Myeloid_score=weighted.mean(CAF_Myeloid_score,pmax(n_CAF,1),na.rm=TRUE), Myeloid_score=weighted.mean(Myeloid_score,pmax(n_Myeloid,1),na.rm=TRUE),
  Exhaustion_score=weighted.mean(Exhaustion_score,pmax(n_Tcell+n_NK,1),na.rm=TRUE), Cytotoxicity_score=weighted.mean(Cytotoxicity_score,pmax(n_Tcell+n_NK,1),na.rm=TRUE)
), by=patient_id]
patient_results[, CCDC80_high_CAF := CAF_CCDC80_mean >= median(CAF_CCDC80_mean,na.rm=TRUE)]
fwrite(patient_results,file.path(out,"GSE293189_patient_level_CCDC80_TME.csv"))

cor_vars <- c("Myeloid_fraction","Myeloid_score","Exhaustion_score","Cytotoxicity_score")
cor_out <- rbindlist(lapply(cor_vars, function(v) {tt<-suppressWarnings(cor.test(patient_results$CAF_CCDC80_mean,patient_results[[v]],method="spearman",exact=FALSE)); data.table(variable=v,n=nrow(patient_results),rho=unname(tt$estimate),p=tt$p.value)}))
cor_out[,FDR:=p.adjust(p,"BH")]
fwrite(cor_out,file.path(out,"GSE293189_patient_level_CCDC80_TME_correlations.csv"))
fwrite(data.table(n_patients=nrow(patient_results),n_tumor_libraries=nrow(tum),n_tumor_cells=sum(tum$n_cells)),file.path(out,"GSE293189_patient_level_audit.csv"))
writeLines(c("GSE293189 patient-level CCDC80-TME analysis","Technical replicates were collapsed by patient.","Tumor libraries only were used for patient-level association.","Cell labels are marker-score based and exploratory; all inference uses patient-level summaries."),file.path(out,"RUN_STATUS.txt"))
