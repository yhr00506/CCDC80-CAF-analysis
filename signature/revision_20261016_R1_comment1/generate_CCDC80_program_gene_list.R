#!/usr/bin/env Rscript
summary_file <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/analysis_revisions_20260906/virtual_knockout_repeat/CCDC80/CCDC80_repeat_summary.csv"
asset_dir <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/revision_20261016_R1_comment1/assets"
dir.create(asset_dir,recursive=TRUE,showWarnings=FALSE)
x <- read.csv(summary_file,check.names=FALSE,stringsAsFactors=FALSE)
x <- x[x$gene!="CCDC80" & x$n_runs==3,,drop=FALSE]
x <- x[order(-x$mean_FC,x$gene),,drop=FALSE]
x <- head(x,19)
stopifnot(nrow(x)==19, !anyDuplicated(x$gene))
x$rank <- seq_len(nrow(x))
x <- x[,c("rank","gene","n_runs","n_FDR05","min_FDR","max_abs_Z","mean_FC")]
write.csv(x,file.path(asset_dir,"CCDC80_associated_19_gene_program.csv"),row.names=FALSE)
writeLines(x$gene,file.path(asset_dir,"CCDC80_associated_19_gene_symbols.txt"))
if(any(x$n_FDR05>0)) warning("At least one selected gene has FDR < 0.05 in a run.")
cat("Selected",nrow(x),"genes present in all 3 runs; genes with any run FDR<0.05:",sum(x$n_FDR05>0),"\n")
