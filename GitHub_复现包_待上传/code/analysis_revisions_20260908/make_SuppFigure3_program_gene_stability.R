#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(readr);library(svglite);library(scales);library(ggrepel);library(Seurat)})

# Supplementary Figure 3: 19-gene program composition and stability. Quantitative panels only.
text_size_pt <- 8; line_width <- .45; data_line_width <- .75; point_size <- 1.5
width_mm <- 300; height_mm <- 280
wong <- c("#E69F00","#56B4E9","#009E73","#CC79A7","#D55E00","#F0E442","#0072B2","#000000")
continuous_map <- "viridis" # continuous values use equal-count quantile mapping
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
rev <- file.path(root,"analysis_revisions_20260906"); final <- file.path(root,"final_analysis_20260907")
out <- file.path(final,"SupplementaryFigure3_program_gene_stability"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_pub <- theme_classic(base_family="Times",base_size=text_size_pt)+theme(text=element_text(family="Times",size=text_size_pt),plot.title=element_blank(),plot.subtitle=element_blank(),axis.title=element_text(size=text_size_pt),axis.text=element_text(size=text_size_pt),legend.text=element_text(size=text_size_pt),legend.title=element_text(size=text_size_pt),panel.grid.major=element_line(linewidth=line_width,colour="#E6E6E6"),panel.grid.minor=element_blank(),plot.margin=margin(5,6,5,6),legend.position="bottom")

repeat_all <- read_csv(file.path(rev,"virtual_knockout_repeat/CCDC80/CCDC80_all_repeat_results.csv"),show_col_types=FALSE)
repeat_sum <- read_csv(file.path(rev,"virtual_knockout_repeat/CCDC80/CCDC80_repeat_summary.csv"),show_col_types=FALSE)
program_genes <- repeat_sum %>% filter(gene!="CCDC80",n_runs==3) %>% slice_head(n=19) %>% pull(gene)
rr <- repeat_all %>% filter(gene %in% program_genes) %>% mutate(seed=gsub("seed_","",seed))

# A: repeat-level FC heatmap.
pA <- ggplot(rr,aes(seed,factor(gene,levels=rev(program_genes)),fill=FC))+geom_tile(colour="white",linewidth=.25)+geom_text(aes(label=sprintf("%.2f",FC)),size=text_size_pt/3.1,family="Times")+scale_fill_viridis_c(option="viridis",name="Virtual-KO FC")+labs(x="Repeat seed",y=NULL)+theme_pub

# B: gene-level consistency across repeats.
stab <- rr %>% group_by(gene) %>% summarise(mean_FC=mean(FC,na.rm=TRUE),sd_FC=sd(FC,na.rm=TRUE),mean_abs_Z=mean(abs(Z),na.rm=TRUE),n_runs=n(),.groups="drop")
pB <- ggplot(stab,aes(mean_FC,mean_abs_Z,size=n_runs,colour=sd_FC,label=gene))+geom_point()+geom_text_repel(size=text_size_pt/3.2,family="Times",box.padding=.4,point.padding=.25,seed=20260908,max.overlaps=Inf)+scale_colour_viridis_c(option="viridis",name="SD of FC")+labs(x="Mean virtual-KO FC",y="Mean |Z|",size="Runs")+theme_pub

# C: actual leave-one-gene-out score stability from the integrated CAF object.
obj <- readRDS(file.path(final,"Figure1_CAF_atlas/GSE135337_GSE190888_GSE130001_GSE129845_CAF_global_annotated.rds"))
present <- intersect(program_genes,rownames(obj))
expr <- as.matrix(GetAssayData(obj,assay=DefaultAssay(obj),slot="data")[present,,drop=FALSE])
zs <- t(scale(t(expr))); full <- colMeans(zs,na.rm=TRUE)
loo <- bind_rows(lapply(present,function(g) data.frame(gene=g,correlation=cor(full,colMeans(zs[setdiff(present,g),,drop=FALSE],na.rm=TRUE),method="spearman"),n_genes=length(present)-1))) %>% mutate(gene=factor(gene,levels=rev(program_genes)))
pC <- ggplot(loo,aes(correlation,gene))+geom_vline(xintercept=.995,linetype=2,linewidth=line_width)+geom_segment(aes(x=.995,xend=correlation,y=gene,yend=gene),linewidth=data_line_width,colour=wong[7])+geom_point(size=point_size,colour=wong[7])+geom_text(aes(x=.999,label=sprintf("r=%.3f",correlation)),hjust=1,size=text_size_pt/3.2,family="Times")+scale_x_continuous(limits=c(min(loo$correlation)-.01,1.005),expand=expansion(mult=c(.02,.2)))+labs(x="Spearman correlation with full 19-gene score",y=NULL)+theme_pub

# D: gene-level evidence rank.
rank_df <- stab %>% mutate(rank_score=rank(-mean_abs_Z)+rank(-mean_FC),gene=factor(gene,levels=rev(gene)))
pD <- ggplot(rank_df,aes(mean_abs_Z,gene,fill=mean_FC))+geom_col(width=.68)+geom_text(aes(label=sprintf("FC %.2f; |Z| %.2f",mean_FC,mean_abs_Z)),hjust=1.05,size=text_size_pt/3.2,family="Times",colour="white")+scale_fill_viridis_c(option="viridis",name="Mean FC")+scale_x_continuous(expand=expansion(mult=c(0,.35)))+labs(x="Mean absolute Z",y=NULL)+theme_pub

# E: enrichment dot plot from GO, Reactome and Hallmark outputs.
go <- read_csv(file.path(rev,"CCDC80_ECM_program_validation/CCDC80_ECM_program_GO_ORA.csv"),show_col_types=FALSE) %>% transmute(source="GO",term=Description,score=-log10(p.adjust),size=Count)
re <- read_csv(file.path(rev,"CCDC80_ECM_program_validation/CCDC80_ECM_program_Reactome_ORA.csv"),show_col_types=FALSE) %>% transmute(source="Reactome",term=Description,score=-log10(p.adjust),size=Count)
ha <- read_csv(file.path(rev,"CCDC80_ECM_program_validation/TCGA_program_GSEA_Hallmark_fgsea.csv"),show_col_types=FALSE) %>% transmute(source="Hallmark",term=pathway,score=-log10(padj),size=size,leadingEdge=leadingEdge)
enr <- bind_rows(go,re,ha) %>% mutate(term=substr(term,1,48)) %>% group_by(source) %>% slice_max(score,n=6,with_ties=FALSE) %>% ungroup() %>% mutate(term=factor(term,levels=rev(term)))
pE <- ggplot(enr,aes(score,term,size=size,colour=source))+geom_point()+facet_wrap(~source,scales="free_y",ncol=1)+scale_colour_manual(values=wong[1:3])+labs(x="−log10 adjusted P",y=NULL,size="Gene count",colour=NULL)+theme_pub+theme(legend.position="bottom")

# F: leading-edge overlap count across Hallmark pathways.
lead <- ha %>% mutate(pathway=substr(term,1,42),term=substr(term,1,42),leading_n=sapply(strsplit(gsub(";","|",ifelse(is.na(leadingEdge),"",leadingEdge)),"\\|"),function(x) sum(x!=""))) %>% arrange(desc(score)) %>% slice_head(n=10) %>% mutate(pathway=factor(pathway,levels=rev(pathway)))
pF <- ggplot(lead,aes(leading_n,pathway,fill=score))+geom_col(width=.68)+geom_text(aes(label=sprintf("n=%d; FDR %.2g",leading_n,10^(-score))),hjust=1.05,size=text_size_pt/3.2,family="Times",colour="white")+scale_fill_viridis_c(option="viridis",name="−log10 FDR")+scale_x_continuous(expand=expansion(mult=c(0,.50)))+labs(x="Leading-edge gene count",y=NULL) +theme_pub

pEF <- (pE|pF)+plot_layout(widths=c(1,1.45))
fig <- (free(pA)|pB)/(free(pC)|pD)/pEF+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(family="Times",size=text_size_pt,face="bold"))
ggsave(file.path(out,"SupplementaryFigure3_program_gene_stability.pdf"),fig,width=width_mm,height=height_mm,units="mm",device=cairo_pdf,bg="transparent")
ggsave(file.path(out,"SupplementaryFigure3_program_gene_stability.svg"),fig,width=width_mm,height=height_mm,units="mm",device=svglite,bg="transparent")
ggsave(file.path(out,"SupplementaryFigure3_program_gene_stability.tif"),fig,width=width_mm,height=height_mm,units="mm",dpi=600,compression="lzw",device="tiff",bg="white")
write_csv(rr,file.path(out,"SuppFig3_source_repeat_results.csv"));write_csv(stab,file.path(out,"SuppFig3_source_gene_stability.csv"));write_csv(loo,file.path(out,"SuppFig3_source_leave_one_gene_out.csv"));write_csv(rank_df,file.path(out,"SuppFig3_source_gene_rank.csv"));write_csv(enr,file.path(out,"SuppFig3_source_enrichment.csv"));write_csv(lead,file.path(out,"SuppFig3_source_leading_edge.csv"))
writeLines(c("Supplementary Figure 3 19-gene program composition and stability completed.",paste0("Program genes=",length(program_genes),"; integrated cells=",ncol(obj)),"No schematic or workflow panel was generated."),file.path(out,"SuppFig3_STATUS.txt"))
