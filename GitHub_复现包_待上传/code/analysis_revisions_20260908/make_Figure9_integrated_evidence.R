#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(readr);library(survival);library(svglite);library(scales);library(ggrepel)})

# Quantitative evidence integration only: no schematic or workflow panel.
text_size_pt <- 8
line_width <- 0.45
data_line_width <- 0.75
point_size <- 1.8
width_mm <- 240
height_mm <- 280
wong <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#000000")
continuous_map <- "viridis" # finite continuous values use equal-count quantile mapping
root <- "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
final <- file.path(root,"final_analysis_20260907"); rev <- file.path(root,"analysis_revisions_20260906")
out <- file.path(final,"Figure9_integrated_evidence"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_pub <- theme_classic(base_family="Times",base_size=text_size_pt)+theme(text=element_text(family="Times",size=text_size_pt),plot.title=element_blank(),plot.subtitle=element_blank(),axis.title=element_text(size=text_size_pt),axis.text=element_text(size=text_size_pt),legend.text=element_text(size=text_size_pt),legend.title=element_text(size=text_size_pt),panel.grid.major=element_line(linewidth=line_width,colour="#E6E6E6"),panel.grid.minor=element_blank(),plot.margin=margin(5,6,5,6),legend.position="bottom")

# A: signed effect estimates across independent evidence layers.
f4 <- read_csv(file.path(final,"Figure4_CAF_myeloid/Figure4_source_meta_effects.csv"),show_col_types=FALSE)
f5 <- read_csv(file.path(final,"Figure5_external_validation/Figure5_source_forest.csv"),show_col_types=FALSE)
f6 <- read_csv(file.path(final,"Figure6_spatial_mechanistic/Figure6_source_spatial_summary.csv"),show_col_types=FALSE)
ev <- bind_rows(f5 %>% transmute(label=label,estimate=log(HR),lo=log(lo),hi=log(hi),type="log HR",text=sprintf("HR %.2f",HR)),f4 %>% transmute(label=metric,estimate=rho_fisher,lo=NA_real_,hi=NA_real_,type="rho",text=sprintf("rho %.2f; FDR %.3g",rho_fisher,FDR)),f6 %>% summarise(label="Spatial ECM–myeloid median",estimate=median(rho_ECM_myeloid,na.rm=TRUE),lo=quantile(rho_ECM_myeloid,.25,na.rm=TRUE),hi=quantile(rho_ECM_myeloid,.75,na.rm=TRUE),type="rho",text=sprintf("rho %.2f",median(rho_ECM_myeloid,na.rm=TRUE)))) %>% mutate(label=make.unique(label),label=factor(label,levels=rev(label)))
pA <- ggplot(ev,aes(estimate,label,colour=type))+geom_vline(xintercept=0,linetype=2,linewidth=line_width)+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.18,linewidth=data_line_width,na.rm=TRUE)+geom_point(size=point_size)+geom_text_repel(aes(label=text),direction="y",hjust=0,nudge_x=.12,box.padding=.35,point.padding=.25,min.segment.length=0,seed=20260908,size=text_size_pt/3.2,family="Times",na.rm=TRUE)+scale_colour_manual(values=c("log HR"=wong[7],"rho"=wong[1]))+scale_x_continuous(expand=expansion(mult=c(.04,.35)))+labs(x="Signed effect estimate",y=NULL,colour=NULL)+theme_pub

# B/C: tertile gradients, distinct from the Figure 7 median split.
tc <- read_csv(file.path(rev,"independent_ECM_comparator/TCGA_incremental_value_input.csv"),show_col_types=FALSE) %>% filter(complete.cases(CCDC80_ECM_program,os_time,os_event)) %>% mutate(tertile=factor(ntile(CCDC80_ECM_program,3),levels=1:3,labels=c("Low","Middle","High")))
gr <- tc %>% group_by(tertile) %>% summarise(n=n(),events=sum(os_event),rate=events/n,lo=binom.test(events,n)$conf.int[1],hi=binom.test(events,n)$conf.int[2],.groups="drop")
pB <- ggplot(gr,aes(tertile,rate,colour=tertile))+geom_linerange(aes(ymin=lo,ymax=hi),linewidth=data_line_width)+geom_point(size=point_size)+geom_text(aes(label=sprintf("%d/%d; %.1f%%",events,n,100*rate)),vjust=-1,size=text_size_pt/3.2,family="Times")+scale_colour_manual(values=wong[1:3])+scale_y_continuous(labels=percent_format(),limits=c(0,max(gr$hi)*1.2))+labs(x="CCDC80 program tertile",y="Observed OS event rate",colour=NULL)+theme_pub
oe <- gr %>% mutate(expected=mean(tc$os_event),ratio=rate/expected,lo_ratio=lo/expected,hi_ratio=hi/expected)
pC <- ggplot(oe,aes(tertile,ratio,colour=tertile))+geom_hline(yintercept=1,linetype=2,linewidth=line_width)+geom_linerange(aes(ymin=lo_ratio,ymax=hi_ratio),linewidth=data_line_width)+geom_point(size=point_size)+geom_text(aes(label=sprintf("O/E %.2f",ratio)),vjust=-1,size=text_size_pt/3.2,family="Times")+scale_colour_manual(values=wong[1:3])+labs(x="CCDC80 program tertile",y="Observed / expected event rate",colour=NULL)+theme_pub

# D: likelihood-ratio increment beyond the independent ECM signature.
l1 <- read_csv(file.path(rev,"independent_ECM_comparator/TCGA_incremental_value_LRT_ECM_plus_CCDC80.csv"),show_col_types=FALSE)
l2 <- read_csv(file.path(rev,"independent_ECM_comparator/TCGA_incremental_value_LRT_CCDC80_plus_ECM.csv"),show_col_types=FALSE)
inc <- data.frame(addition=c("CCDC80 added to ECM","ECM added to CCDC80"),chisq=c(l1$Chisq[2],l2$Chisq[2]),p=c(l1$`Pr(>|Chi|)`[2],l2$`Pr(>|Chi|)`[2]))
pD <- ggplot(inc,aes(chisq,reorder(addition,chisq)))+geom_vline(xintercept=3.84,linetype=2,linewidth=line_width)+geom_segment(aes(x=0,xend=chisq,y=addition,yend=addition),linewidth=data_line_width,colour=wong[7])+geom_point(size=point_size,colour=wong[7])+geom_text(aes(label=sprintf("χ² %.2f; P=%.3g",chisq,p)),hjust=-.05,size=text_size_pt/3.2,family="Times")+scale_x_continuous(expand=expansion(mult=c(.02,.35)))+labs(x="Likelihood-ratio χ² (1 df)",y=NULL)+theme_pub

# E: patient-level quantitative evidence matrix.
pat <- read_csv(file.path(final,"Figure4_CAF_myeloid/Figure4_source_patient_level.csv"),show_col_types=FALSE) %>% mutate(patient=paste0(cohort," / ",patient_id)) %>% select(patient,CCDC80_CAF_z,CAF_Myeloid_score,Myeloid_score,Exhaustion_score,Cytotoxicity_score) %>% pivot_longer(-patient,names_to="metric",values_to="value") %>% group_by(metric) %>% mutate(value=as.numeric(scale(value))) %>% ungroup()
pE <- ggplot(pat,aes(metric,patient,fill=value))+geom_tile(colour="white",linewidth=.25)+scale_fill_viridis_c(option="viridis",name="Within-metric z")+labs(x=NULL,y="Direct CAF patient")+theme_pub+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(angle=35,hjust=1),legend.position="right")

# F: compact quantitative evidence-domain summary, not a conceptual diagram.
dom <- data.frame(domain=c("CAF–myeloid association","Spatial ECM–myeloid","External OS replication","Bootstrap HR > 1","CCDC80 incremental LRT"),value=c(f4$p_fisher[f4$metric=="CAF–myeloid score"],median(f6$FDR_ECM_myeloid,na.rm=TRUE),min(f5$HR[f5$cohort!="TCGA"],na.rm=TRUE),mean(read_csv(file.path(final,"Figure8_robustness_incremental/Figure8_source_bootstrap.csv"),show_col_types=FALSE)$HR>1),inc$chisq[1]))
dom$display <- c(sprintf("FDR %.3g",dom$value[1]),sprintf("median FDR %.3g",dom$value[2]),sprintf("min HR %.2f",dom$value[3]),sprintf("HR>1 %.1f%%",100*dom$value[4]),sprintf("χ² %.2f",dom$value[5]))
pF <- ggplot(dom,aes(1,reorder(domain,value),fill=value))+geom_tile(width=1.1,height=.72)+geom_text(aes(label=display),colour="white",fontface="bold",size=text_size_pt/3.0,family="Times")+scale_fill_viridis_c(option="viridis",guide="none")+scale_x_continuous(breaks=NULL)+labs(x=NULL,y=NULL)+theme_pub+theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())

fig <- (pA|pB)/(pC|pD)/(pE|pF)+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(family="Times",size=text_size_pt,face="bold"))
ggsave(file.path(out,"Figure9_integrated_evidence.pdf"),fig,width=width_mm,height=height_mm,units="mm",device=cairo_pdf,bg="transparent")
ggsave(file.path(out,"Figure9_integrated_evidence.svg"),fig,width=width_mm,height=height_mm,units="mm",device=svglite,bg="transparent")
ggsave(file.path(out,"Figure9_integrated_evidence.tif"),fig,width=width_mm,height=height_mm,units="mm",dpi=600,compression="lzw",device="tiff",bg="white")
write_csv(ev,file.path(out,"Figure9_source_cross_evidence.csv"));write_csv(gr,file.path(out,"Figure9_source_tertile_event_gradient.csv"));write_csv(oe,file.path(out,"Figure9_source_observed_expected.csv"));write_csv(inc,file.path(out,"Figure9_source_incremental_LRT.csv"));write_csv(pat,file.path(out,"Figure9_source_patient_matrix.csv"));write_csv(dom,file.path(out,"Figure9_source_evidence_domains.csv"))
writeLines(c("Figure 9 quantitative cross-evidence consistency and translational value completed.",paste0("Direct CAF patients=",n_distinct(pat$patient),"; TCGA N=",nrow(tc),"; events=",sum(tc$os_event)),"No schematic or workflow panel was generated."),file.path(out,"Figure9_STATUS.txt"))
