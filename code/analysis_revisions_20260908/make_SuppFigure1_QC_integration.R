#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(readr);library(svglite);library(scales)})

# Supplementary Figure 1: quantitative QC and integration audit only.
text_size_pt <- 8; line_width <- .45; data_line_width <- .75; point_size <- 1.2
width_mm <- 240; height_mm <- 280
wong <- c("#E69F00","#56B4E9","#009E73","#CC79A7","#D55E00","#F0E442","#0072B2","#000000")
continuous_map <- "viridis" # finite continuous values use equal-count quantile mapping
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
final <- file.path(root,"final_analysis_20260907"); out <- file.path(final,"SupplementaryFigure1_QC_integration"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_pub <- theme_classic(base_family="Times",base_size=text_size_pt)+theme(text=element_text(family="Times",size=text_size_pt),plot.title=element_blank(),plot.subtitle=element_blank(),axis.title=element_text(size=text_size_pt),axis.text=element_text(size=text_size_pt),legend.text=element_text(size=text_size_pt),legend.title=element_text(size=text_size_pt),panel.grid.major=element_line(linewidth=line_width,colour="#E6E6E6"),panel.grid.minor=element_blank(),plot.margin=margin(5,6,5,6),legend.position="bottom")

meta <- read_csv(file.path(root,"external_CAF_data/integrated/CAF_specific_integrated_metadata.csv"),show_col_types=FALSE)
audit <- read_csv(file.path(root,"external_CAF_data/integrated/CAF_specific_candidate_audit.csv"),show_col_types=FALSE)
axes <- read_csv(file.path(root,"analysis_revisions_20260906/NNMT_data_driven_target_screen/sample_level_CAF_state_axes.csv"),show_col_types=FALSE)

# A: cohort scale.
cohort <- meta %>% count(dataset,name="CAF_cells") %>% left_join(audit %>% group_by(dataset) %>% summarise(CAF_samples=n_distinct(sample),.groups="drop"),by="dataset") %>% mutate(label=sprintf("%s\n%d cells; %d samples",dataset,CAF_cells,CAF_samples))
pA <- ggplot(cohort,aes(CAF_cells,reorder(dataset,CAF_cells)))+geom_col(fill=wong[7],width=.65)+geom_text(aes(label=label),hjust=-.05,size=text_size_pt/3.2,family="Times",lineheight=.9)+scale_x_continuous(expand=expansion(mult=c(0,.25)),labels=comma)+labs(x="Retained CAF cells",y=NULL)+theme_pub

# B: QC distributions.
qc <- meta %>% select(dataset,nCount_RNA,nFeature_RNA,percent.mt) %>% pivot_longer(-dataset,names_to="metric",values_to="value") %>% mutate(metric=recode(metric,nCount_RNA="RNA counts",nFeature_RNA="Detected genes",percent.mt="Mitochondrial %")) %>% filter(is.finite(value))
pB <- ggplot(qc,aes(dataset,value,fill=dataset))+geom_violin(scale="width",trim=TRUE,colour=NA,alpha=.65)+geom_boxplot(width=.13,outlier.size=.25,colour="#333333")+facet_wrap(~metric,scales="free_y",ncol=1)+scale_fill_manual(values=wong[1:4])+labs(x=NULL,y=NULL,fill=NULL)+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")

# C: CAF marker score and detection quality.
marker <- meta %>% group_by(dataset) %>% summarise(CAF_score=median(CAF_score,na.rm=TRUE),CAF_score_q25=quantile(CAF_score,.25,na.rm=TRUE),CAF_score_q75=quantile(CAF_score,.75,na.rm=TRUE),epithelial=median(epithelial_score,na.rm=TRUE),immune=median(immune_score,na.rm=TRUE),endothelial=median(endothelial_score,na.rm=TRUE),.groups="drop") %>% pivot_longer(-dataset,names_to="metric",values_to="value")
pC <- ggplot(marker,aes(value,reorder(dataset,value),colour=metric))+geom_point(size=point_size)+facet_wrap(~metric,ncol=2,scales="free_x")+scale_colour_manual(values=wong[1:6])+labs(x="Median score",y=NULL,colour=NULL)+theme_pub+theme(legend.position="none")

# D: contamination score audit at cell level.
cont <- meta %>% select(dataset,epithelial_score,immune_score,endothelial_score) %>% pivot_longer(-dataset,names_to="contamination",values_to="score") %>% mutate(contamination=recode(contamination,epithelial_score="Epithelial",immune_score="Immune",endothelial_score="Endothelial")) %>% filter(is.finite(score))
pD <- ggplot(cont,aes(contamination,score,fill=contamination))+geom_violin(scale="width",trim=TRUE,colour=NA,alpha=.6)+geom_boxplot(width=.12,outlier.size=.15)+facet_wrap(~dataset,ncol=2)+scale_fill_manual(values=wong[1:3])+labs(x=NULL,y="Contamination module score",fill=NULL)+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")

# E: sample-level CAF retention audit.
sample_audit <- meta %>% count(dataset,sample,name="CAF_cells") %>% mutate(included=CAF_cells>=20, status=ifelse(included,"Included","Below 20 cells"))
pE <- ggplot(sample_audit,aes(CAF_cells,reorder(sample,paste(dataset,CAF_cells)),colour=status))+geom_point(size=point_size)+geom_vline(xintercept=20,linetype=2,linewidth=line_width)+facet_wrap(~dataset,scales="free_y",ncol=1)+scale_x_log10(labels=comma)+scale_colour_manual(values=c("Included"=wong[7],"Below 20 cells"=wong[5]))+labs(x="CAF cells per sample (log10)",y=NULL,colour=NULL)+theme_pub+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

# F: sample-level state-axis matrix; quantitative heatmap, not a conceptual diagram.
heat <- axes %>% pivot_longer(c(ECM_CAF_module1,Activated_myCAF_module1,iCAF_module1,Hypoxia_CAF_module1,TGFb_CAF_module1,Vascular_CAF_module1),names_to="axis",values_to="score") %>% mutate(axis=recode(axis,ECM_CAF_module1="ECM",Activated_myCAF_module1="Activated myCAF",iCAF_module1="iCAF",Hypoxia_CAF_module1="Hypoxia",TGFb_CAF_module1="TGFβ",Vascular_CAF_module1="Vascular"))
pF <- ggplot(heat,aes(axis,reorder(sample_id,score),fill=score))+geom_tile(colour="white",linewidth=.25)+scale_fill_viridis_c(option="viridis",name="Score")+labs(x=NULL,y="Sample")+theme_pub+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(angle=35,hjust=1),legend.position="right")

fig <- (pA|pB)/(pC|pD)/(pE|pF)+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(family="Times",size=text_size_pt,face="bold"))
ggsave(file.path(out,"SupplementaryFigure1_QC_integration.pdf"),fig,width=width_mm,height=height_mm,units="mm",device=cairo_pdf,bg="transparent")
ggsave(file.path(out,"SupplementaryFigure1_QC_integration.svg"),fig,width=width_mm,height=height_mm,units="mm",device=svglite,bg="transparent")
ggsave(file.path(out,"SupplementaryFigure1_QC_integration.tif"),fig,width=width_mm,height=height_mm,units="mm",dpi=600,compression="lzw",device="tiff",bg="white")
write_csv(cohort,file.path(out,"SuppFig1_source_cohort_scale.csv"));write_csv(qc,file.path(out,"SuppFig1_source_QC_distributions.csv"));write_csv(marker,file.path(out,"SuppFig1_source_marker_contamination_summary.csv"));write_csv(cont,file.path(out,"SuppFig1_source_contamination_scores.csv"));write_csv(sample_audit,file.path(out,"SuppFig1_source_sample_retention.csv"));write_csv(heat,file.path(out,"SuppFig1_source_state_axes.csv"))
writeLines(c("Supplementary Figure 1 quantitative QC and integration audit completed.",paste0("Integrated CAF cells=",nrow(meta),"; datasets=",n_distinct(meta$dataset),"; sample IDs=",n_distinct(meta$sample_id)),"No UMAP, schematic, or workflow panel was generated."),file.path(out,"SuppFig1_STATUS.txt"))
