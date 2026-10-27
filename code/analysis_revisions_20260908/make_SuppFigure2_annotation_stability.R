#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(readr);library(svglite);library(scales)})

# Supplementary Figure 2: CAF annotation stability and state-call audit only.
text_size_pt <- 8; line_width <- .45; data_line_width <- .75; point_size <- 1.4
width_mm <- 260; height_mm <- 320
wong <- c("#E69F00","#56B4E9","#009E73","#CC79A7","#D55E00","#F0E442","#0072B2","#000000")
continuous_map <- "viridis" # continuous values use equal-count quantile mapping
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
final <- file.path(root,"final_analysis_20260907"); out <- file.path(final,"SupplementaryFigure2_annotation_stability"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_pub <- theme_classic(base_family="Times",base_size=text_size_pt)+theme(text=element_text(family="Times",size=text_size_pt),plot.title=element_blank(),plot.subtitle=element_blank(),axis.title=element_text(size=text_size_pt),axis.text=element_text(size=text_size_pt),legend.text=element_text(size=text_size_pt),legend.title=element_text(size=text_size_pt),panel.grid.major=element_line(linewidth=line_width,colour="#E6E6E6"),panel.grid.minor=element_blank(),plot.margin=margin(5,6,5,6),legend.position="bottom")

meta <- read_csv(file.path(root,"external_CAF_data/integrated/CAF_specific_integrated_metadata.csv"),show_col_types=FALSE)
conf <- read_csv(file.path(final,"Figure1_CAF_atlas/global_annotation_old_new_confusion.csv"),show_col_types=FALSE)
bycohort <- read_csv(file.path(final,"Figure1_CAF_atlas/global_annotation_by_cohort.csv"),show_col_types=FALSE)
marker <- read_csv(file.path(final,"Figure1_CAF_atlas/Figure1_advanced_marker_dotplot_source.csv"),show_col_types=FALSE)
sample_scores <- read_csv(file.path(final,"Figure1_CAF_atlas/Figure1_advanced_sample_scores_source.csv"),show_col_types=FALSE)
cell_scores <- read_csv(file.path(final,"Figure1_CAF_atlas/global_annotation_cell_scores.csv"),show_col_types=FALSE)

# A: old/new annotation confusion matrix.
conf$old_subtype <- factor(conf$old_subtype); conf$new_global_subtype <- factor(conf$new_global_subtype)
pA <- ggplot(conf,aes(old_subtype,new_global_subtype,fill=Freq))+geom_tile(colour="white",linewidth=.35)+geom_text(aes(label=Freq),family="Times",size=text_size_pt/3.1)+scale_fill_viridis_c(option="viridis",name="Cells")+labs(x="Previous state",y="Global state")+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1))

# B: cohort composition with exact n and percentage labels.
comp <- bycohort %>% group_by(dataset) %>% mutate(total=sum(Freq),fraction=Freq/total) %>% ungroup() %>% filter(Freq>0) %>% mutate(label=sprintf("%d (%.1f%%)",Freq,100*fraction))
pB <- ggplot(comp,aes(dataset,fraction,fill=CAF_subtype_global))+geom_col(width=.72)+geom_text(aes(label=label),position=position_stack(vjust=.5),size=text_size_pt/3.2,family="Times",colour="white")+scale_y_continuous(labels=percent_format())+scale_fill_manual(values=wong[1:6])+labs(x=NULL,y="CAF state fraction",fill="CAF state")+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1))

# C: marker specificity dot plot, using mean expression and detection.
mk <- marker %>% group_by(CAF_subtype_global) %>% mutate(specificity=mean_expression-max(mean_expression,na.rm=TRUE)+mean_expression) %>% ungroup() %>% filter(pct_detected>0)
pC <- ggplot(mk,aes(gene,CAF_subtype_global,size=pct_detected,colour=mean_expression))+geom_point()+scale_colour_viridis_c(option="viridis",name="Mean expression")+scale_size(range=c(1,5),name="Detected (%)")+labs(x=NULL,y=NULL)+theme_pub+theme(axis.text.x=element_text(angle=60,hjust=1),legend.position="right")

# D: sample-level confidence margin: first minus second state score.
score_cols <- c("CAF_score_global_ECM_CAF","CAF_score_global_TGFb_CAF","CAF_score_global_iCAF","CAF_score_global_Vascular_CAF","CAF_score_global_Hypoxia_CAF","CAF_score_global_apCAF")
score_mat <- as.matrix(sample_scores[,intersect(score_cols,names(sample_scores))]); score_mat <- apply(score_mat,2,as.numeric)
top2 <- t(apply(score_mat,1,function(x) sort(x,decreasing=TRUE)[1:2]))
margin <- data.frame(dataset=sample_scores$dataset,top_score=top2[,1],second_score=top2[,2],margin=top2[,1]-top2[,2])
pD <- ggplot(margin,aes(dataset,margin,fill=dataset))+geom_violin(scale="width",trim=TRUE,colour=NA,alpha=.6)+geom_boxplot(width=.12,outlier.size=.2)+geom_hline(yintercept=0,linetype=2,linewidth=line_width)+scale_fill_manual(values=wong[1:4])+labs(x=NULL,y="Top score − second score")+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")

# E: sample-level state composition; counts are explicit.
cell_comp <- cell_scores %>% select(cell_id,state_for_plot=global_subtype) %>% left_join(meta %>% select(cell_id,dataset,sample_id),by="cell_id") %>% filter(!is.na(dataset),!is.na(sample_id),!is.na(state_for_plot),state_for_plot!="") %>% count(dataset,sample_id,state_for_plot) %>% group_by(dataset,sample_id) %>% mutate(fraction=n/sum(n),label=sprintf("%s\n%d",state_for_plot,n)) %>% ungroup()
pE <- ggplot(cell_comp,aes(sample_id,fraction,fill=state_for_plot))+geom_col(width=.72)+scale_y_continuous(labels=percent_format())+scale_fill_manual(values=wong[1:6])+labs(x="Sample",y="CAF state fraction",fill="CAF state")+theme_pub+theme(axis.text.x=element_text(angle=55,hjust=1,size=text_size_pt-1),legend.position="right")

# F: state-score correlation matrix at sample level.
score_long <- sample_scores %>% select(dataset,sample_id,starts_with("CAF_score_global_"))
score_cor <- cor(score_long %>% select(-dataset,-sample_id),use="pairwise.complete.obs",method="spearman") %>% as.data.frame() %>% tibble::rownames_to_column("x") %>% pivot_longer(-x,names_to="y",values_to="rho")
score_cor$x <- gsub("CAF_score_global_","",score_cor$x); score_cor$y <- gsub("CAF_score_global_","",score_cor$y)
pF <- ggplot(score_cor,aes(x,y,fill=rho))+geom_tile(colour="white",linewidth=.25)+geom_text(aes(label=sprintf("%.2f",rho)),size=text_size_pt/3.0,family="Times")+scale_fill_viridis_c(option="viridis",limits=c(-1,1),name="Spearman rho")+labs(x=NULL,y=NULL)+theme_pub+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="right")

fig <- (pA|pB)/(pC|pD)/(pE|pF)+plot_layout(widths=c(1,1),heights=c(1,1.25,1.25))+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(family="Times",size=text_size_pt,face="bold"))
ggsave(file.path(out,"SupplementaryFigure2_annotation_stability.pdf"),fig,width=width_mm,height=height_mm,units="mm",device=cairo_pdf,bg="transparent")
ggsave(file.path(out,"SupplementaryFigure2_annotation_stability.svg"),fig,width=width_mm,height=height_mm,units="mm",device=svglite,bg="transparent")
ggsave(file.path(out,"SupplementaryFigure2_annotation_stability.tif"),fig,width=width_mm,height=height_mm,units="mm",dpi=600,compression="lzw",device="tiff",bg="white")
write_csv(conf,file.path(out,"SuppFig2_source_confusion_matrix.csv"));write_csv(comp,file.path(out,"SuppFig2_source_cohort_composition.csv"));write_csv(mk,file.path(out,"SuppFig2_source_marker_specificity.csv"));write_csv(margin,file.path(out,"SuppFig2_source_score_margin.csv"));write_csv(cell_comp,file.path(out,"SuppFig2_source_sample_composition.csv"));write_csv(score_cor,file.path(out,"SuppFig2_source_score_correlation.csv"))
writeLines(c("Supplementary Figure 2 CAF annotation stability and state-call audit completed.",paste0("Integrated CAF cells=",nrow(meta),"; datasets=",n_distinct(meta$dataset),"; sample IDs=",n_distinct(meta$sample_id)),"No UMAP, schematic, or workflow panel was generated."),file.path(out,"SuppFig2_STATUS.txt"))
