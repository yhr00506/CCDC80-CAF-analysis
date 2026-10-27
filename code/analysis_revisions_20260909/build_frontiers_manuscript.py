from pathlib import Path
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
import csv

ROOT=Path('/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project')
BASE=ROOT/'analysis_revisions_20260909'
TEMPLATE=BASE/'frontiers_template/Frontiers_Word_Templates/Frontiers_Template.docx'
SOURCE=BASE/'manuscript_outputs/CCDC80_CAF_ECM_myeloid_manuscript_clean.docx'
FINAL=ROOT/'final_analysis_20260907'
OUT=BASE/'frontiers_manuscripts'; OUT.mkdir(exist_ok=True)

def clear_body(doc):
    body=doc._element.body
    for child in list(body):
        if child.tag != qn('w:sectPr'):
            body.remove(child)

def addp(doc,text='',style='Normal'):
    p=doc.add_paragraph(style=style)
    p.add_run(text)
    return p

def heading(doc,text,level=1): return addp(doc,text,f'Heading {level}')

def set_cell_text(cell,text,bold=False):
    cell.text=''
    p=cell.paragraphs[0]
    r=p.add_run(str(text)); r.bold=bold; r.font.name='Arial'; r.font.size=Pt(8)

def add_table(doc,title,rows,headers,note):
    addp(doc,title,'Heading 2')
    t=doc.add_table(rows=1, cols=len(headers))
    t.alignment=WD_TABLE_ALIGNMENT.CENTER
    t.style='Table Grid'
    for c,h in zip(t.rows[0].cells,headers): set_cell_text(c,h,True)
    for row in rows:
        cells=t.add_row().cells
        for c,v in zip(cells,row): set_cell_text(c,v)
    addp(doc,note,'Normal')
    addp(doc,'','Normal')

def readcsv(path):
    with open(path,newline='',encoding='utf-8-sig') as f: return list(csv.DictReader(f))

src=Document(SOURCE)
allp=[p.text.strip() for p in src.paragraphs if p.text.strip()]
try:
    ref_start=allp.index('References')+1
except ValueError:
    ref_start=0
refs=[x for x in allp[ref_start:] if x.split('.',1)[0].isdigit()]
refs=refs[:55]

doc=Document(TEMPLATE); clear_body(doc)
addp(doc,'CCDC80-associated CAF extracellular-matrix programs and myeloid microenvironment remodeling in bladder cancer','Title')
addp(doc,'Penghao Li1, Zhihong Zhang1','Author List')
addp(doc,'1Tianjin Institute of Urology, The Second Hospital of Tianjin Medical University, Tianjin, China','Normal')
addp(doc,'* Correspondence: Penghao Li; lipenghao@tmu.edu.cn','Normal')
addp(doc,'Keywords: bladder cancer; cancer-associated fibroblast; CCDC80; extracellular matrix; tumor microenvironment; single-cell RNA sequencing; myeloid cells','Author List')
addp(doc,'Abstract','Author List')
addp(doc,'Background: Cancer-associated fibroblasts (CAFs) are heterogeneous stromal cells that organize extracellular matrix (ECM) and influence immune contexture in bladder cancer, but reproducible CAF programs that connect these processes remain incompletely defined. Methods: We integrated four publicly available bladder-cancer single-cell cohorts to define CAF states, performed stringent lineage-contamination audits, and constructed a CCDC80-associated ECM program from replicated virtual-knockout analyses. Patient-level CAF–myeloid associations, independent bulk-cohort validation, spatial transcriptomic analyses, ligand–receptor prioritization, CellOracle perturbation analyses, and TCGA clinical models were evaluated in a discovery–validation framework. Results: The integrated atlas contained 5,589 CAFs from four cohorts. Six CAF states were resolved, including ECM, hypoxia, inflammatory, antigen-presenting, TGFβ, and vascular states. A 19-gene CCDC80-associated ECM program was recovered across three virtual-knockout repeats and showed stable leave-one-gene-out scores. In 23 directly profiled patients, the program was associated with myeloid-state variation, with cohort-dependent effect sizes. In TCGA-BLCA, the program was associated with overall survival in a univariable model (HR 1.57, P=0.0005), but not after adjustment for CAF score and stage (HR 1.24, P=0.357). Conclusions: CCDC80 is best interpreted as a representative marker of an ECM-associated CAF program rather than as an independently validated driver or therapeutic target. The integrated evidence supports a matrix-remodeling CAF state linked to myeloid microenvironment variation and provides a testable framework for spatial and experimental validation.','Normal')

heading(doc,'Introduction',1)
addp(doc,'Bladder cancer comprises biologically diverse diseases with distinct risks of recurrence, progression, and treatment resistance [1,2]. Molecular studies show that urothelial tumors are shaped not only by malignant epithelial cells but also by stromal and immune compartments that influence tissue architecture, inflammation, and treatment response [3–5]. Molecular subtype and stromal organization are particularly relevant in muscle-invasive disease and in the interpretation of immune checkpoint therapy [6–8]. Current clinical guidance provides the treatment context in which these biological distinctions are evaluated [9–11].')
addp(doc,'CAFs are heterogeneous mesenchymal cells that remodel ECM, regulate tissue mechanics, produce paracrine signals, and influence leukocyte recruitment and function [12,13]. Myofibroblastic, inflammatory, antigen-presenting, vascular, and matrix-associated programs can coexist within the same tumor and may vary with tissue context [14–16]. A consensus framework therefore emphasizes CAF plasticity and the need to distinguish cellular identity from transient functional state [18].')
addp(doc,'Single-cell and spatial transcriptomics provide complementary strategies for resolving stromal heterogeneity and linking CAF states to malignant and immune compartments [19–21]. Computational tools for integration, pathway scoring, and cell–cell communication make such analyses practical, but their outputs require independent validation and conservative causal interpretation [23–25].')
addp(doc,'CCDC80, also known as DRO1, has been studied as a stromal- and tumor-associated extracellular protein with context-dependent links to cancer biology. Prior work has described CCDC80 loss or altered expression in experimental tumor models and prognostic gene signatures, including muscle-invasive bladder cancer [29–31]. These studies motivate evaluation of CCDC80 in a cell-state framework, while not establishing that CCDC80 itself drives CAF function.')
addp(doc,'We reanalyzed four public bladder-cancer single-cell cohorts to define reproducible CAF states and identify a CCDC80-associated ECM program. We then assessed patient-level CAF–myeloid associations, external transcriptomic reproducibility, spatial co-localization, ligand–receptor and regulatory hypotheses, and clinical relevance in TCGA-BLCA and IMvigor210. The analysis was designed to distinguish discovery, external validation, and mechanism-hypothesis generation.')

heading(doc,'Materials and Methods',1)
heading(doc,'Study design and public datasets',2)
addp(doc,'This retrospective study used public, de-identified datasets and involved no new human recruitment, intervention, or identifiable information. The integrated single-cell discovery atlas included GSE135337, GSE190888, GSE130001, and GSE129845. External analyses used TCGA-BLCA, GSE32894, GSE13507, GSE211388, GSE293189, GSE171351 spatial data, and IMvigor210. The source studies and accession records are publicly available [40,54,55]. Analyses were organized into discovery, external validation, and mechanism-hypothesis layers.')
heading(doc,'Single-cell preprocessing, integration, and contamination audit',2)
addp(doc,'Count matrices were processed with Seurat using normalization, variable-feature selection, dimensional reduction, neighborhood construction, clustering, and marker identification [23]. Dataset effects were assessed with Harmony [24]. CAFs were identified using standardized fibroblast marker modules and subjected to epithelial, immune, and endothelial contamination audits. The audit was performed before state assignment and required that candidate CAF populations retain fibroblast identity without excessive lineage-marker signal.')
addp(doc,'The integrated object was evaluated at cell, sample, and cohort levels. We recorded detected genes, transcript counts where available, mitochondrial content, cell yield, and sample identifiers. Samples with no retained CAFs or incomplete source packages were excluded before integration. The final object contained 5,589 CAFs across four cohorts and 13 sample identifiers. These checks were used to prevent cell-level abundance from being mistaken for independent patient replication.')
heading(doc,'CAF state annotation and functional scoring',2)
addp(doc,'CAF states were assigned using standardized marker-module scores, marker specificity, cohort composition, and functional program coherence. The final state labels were ECM CAF, hypoxia CAF, inflammatory CAF, antigen-presenting CAF, TGFβ CAF, and vascular CAF. Functional scores included ECM, activated myCAF, iCAF, apCAF, hypoxia, TGFβ, and vascular modules.')
addp(doc,'Annotation stability was assessed by comparing score margins, marker specificity, cohort composition, and state calls across the integrated object. No single marker was used as a sufficient state definition. This approach follows the principle that CAF identity and CAF function should be analyzed as related but non-identical dimensions [12,13,18].')
heading(doc,'Definition and stability of the CCDC80-associated ECM program',2)
addp(doc,'CCDC80-associated genes were prioritized from replicated virtual-knockout analyses. Genes recovered in all three repeat runs were retained to define a 19-gene program. Program stability was evaluated using repeat-level fold changes, mean absolute Z scores, leave-one-gene-out score correlations, and enrichment consistency. The program was treated as a state-associated score and not as a causal perturbation signature.')
addp(doc,'For each sample or cell, the program score was calculated from the normalized expression of the retained genes after harmonized gene-symbol matching. Leave-one-gene-out analyses recalculated the score after removing each gene and compared it with the full score using Spearman correlation. Enrichment results were interpreted as biological context for the program rather than independent evidence that any pathway is mechanistically activated.')
heading(doc,'CAF–myeloid association and cell–cell communication',2)
addp(doc,'Patient-level CAF and myeloid scores were calculated from annotated single-cell data. Associations were evaluated within directly profiled cohorts and combined using Fisher’s method where appropriate. Ligand–receptor analyses used CellChat-derived communication summaries and focused on ranked hypotheses including CSF1–CSF1R, COL1A1–ITGA1, and IL6–IL6R [27,42].')
addp(doc,'Communication scores represent expression-based communication potential and depend on ligand and receptor abundance, database content, and aggregation. Accordingly, ligand–receptor results were used to prioritize experiments and were not interpreted as evidence of ligand binding, receptor activation, or directional recruitment.')
heading(doc,'External validation and clinical analysis',2)
addp(doc,'The 19-gene program was projected into independent bulk cohorts using sample-level expression scores. Associations with CAF, ECM, immune, and myeloid modules were assessed by Spearman correlation. TCGA-BLCA survival was analyzed with Cox proportional-hazards models; complete-case analyses included 480 patients and 231 events. The univariable model was compared with a model adjusted for CAF score and AJCC stage. IMvigor210 analyses were treated as external context rather than as a prespecified predictive validation.')
addp(doc,'Hazard ratios are reported with 95% confidence intervals. Multiple testing was controlled using the Benjamini–Hochberg procedure [32]. Robustness analyses included bootstrap estimation, leave-one-cohort-out comparisons, calibration summaries, and incremental model comparisons. The analysis was performed with two-sided tests where applicable, and the observational nature of all associations was retained in interpretation.')
heading(doc,'Spatial and regulatory hypothesis generation',2)
addp(doc,'Spatial associations were evaluated using patient-level spatial summaries and permutation-based comparisons in GSE171351. CellOracle analyses and virtual-knockout results were used to prioritize transcriptional regulators and perturbation-sensitive genes. Spatial and regulatory outputs were considered complementary evidence domains and were not used to claim causal regulation.')
heading(doc,'Data and code availability',2)
addp(doc,'The public datasets used in this study are available through GEO, TCGA-linked resources, and the published IMvigor210 research resource. Source tables, analysis scripts, final figures, and supplementary tables are retained in the project archive and can be made available by the corresponding author.')

heading(doc,'Results',1)
heading(doc,'An integrated bladder-cancer CAF atlas contains 5,589 cells across four cohorts',2)
addp(doc,'After standardized CAF extraction and contamination auditing, the integrated object contained 5,589 CAFs from GSE135337, GSE190888, GSE130001, and GSE129845. The cohort-level composition and sample retention were documented before downstream state scoring. QC analyses supported retention of the four cohorts and identified no basis for treating the integrated object as a single biological sample (Figure 1; Supplementary Figure 1; Supplementary Table 1).')
addp(doc,'The integrated atlas was used as a discovery resource rather than a clinical cohort. Cells were retained as the unit of state discovery, whereas patient-level summaries were used for association analyses. This separation reduced the risk of pseudoreplication and allowed the analysis to distinguish state structure from cohort-specific sampling density.')
heading(doc,'Six CAF states define partially overlapping matrix and immune-associated axes',2)
addp(doc,'Global annotation resolved ECM, hypoxia, inflammatory, antigen-presenting, TGFβ, and vascular CAF states. State calls were supported by marker specificity, cohort composition, and functional score separation. The states were related rather than mutually exclusive biological lineages, consistent with a model of context-dependent CAF plasticity (Figure 2; Supplementary Figure 2; Supplementary Table 2).')
addp(doc,'ECM and activated myCAF scores were concentrated in matrix-rich populations, whereas TGFβ, hypoxia, inflammatory, antigen-presenting, and vascular scores contributed distinct but partially overlapping axes. The annotation therefore supports functional specialization without implying that each state is a fixed lineage or that one state is uniformly present across all tumors.')
heading(doc,'A replicated CCDC80-associated ECM program contains 19 genes',2)
addp(doc,'Replicated virtual-knockout analyses identified a 19-gene program associated with the CCDC80 state. The program showed consistent repeat-level fold changes and mean absolute Z scores, while leave-one-gene-out analyses indicated that no single gene fully determined the score. Enrichment results were consistent with ECM organization and stromal remodeling, supporting interpretation of the program as a coordinated state-associated axis (Figure 3; Supplementary Figure 3; Supplementary Table 2).')
addp(doc,'The program was intentionally defined at the program level. CCDC80 was retained as a representative candidate, but the downstream analyses used the multigene score to reduce dependence on dropout, platform-specific detection, and expression of a single transcript. This distinction is central to interpreting the program as a CAF state marker rather than a demonstrated molecular driver.')
heading(doc,'CCDC80-associated CAF scores are linked to patient-level myeloid variation',2)
addp(doc,'Across two directly profiled CCDC80 cohorts comprising 23 patients, the CAF program showed patient-level associations with myeloid scores and related immune-state measures. The combined CAF–myeloid score association yielded Fisher meta P=0.00185 and FDR=0.00925, while individual effect sizes varied by cohort. These results support association between the matrix-associated CAF state and myeloid context but do not establish recruitment or directional signaling (Figure 4; Supplementary Table 4).')
addp(doc,'The cohort-dependent effect sizes are important for interpretation. They indicate that the program–myeloid relationship is not a universal fixed effect that can be separated from tissue composition, clinical context, or cohort design using the present sample size. The patient-level analysis therefore serves as an external association check rather than a definitive estimate of a clinical effect.')
heading(doc,'Independent cohorts support the program–stromal and immune associations',2)
addp(doc,'The 19-gene program was evaluated in independent bulk cohorts and IMvigor210. Directionally concordant associations were observed for CAF, ECM, and immune-related scores, with cohort-specific differences in magnitude. The external analyses support portability of the program beyond the discovery atlas, while platform and cohort heterogeneity limit direct comparison of effect sizes (Figure 5; Supplementary Table 3).')
addp(doc,'The validation strategy focused on direction and domain consistency rather than a single universal cutoff. Because bulk profiles combine malignant, stromal, and immune transcripts, a positive association with the program cannot by itself establish that the corresponding genes originate from CAFs in every sample. The single-cell atlas and the external analyses are therefore interpreted jointly.')
heading(doc,'Spatial, ligand–receptor, and regulatory analyses generate mechanistic hypotheses',2)
addp(doc,'Spatial analyses showed ECM–myeloid co-variation in subsets of patients, whereas CCDC80 itself was less stable across spatial samples. Ranked ligand–receptor candidates included CSF1–CSF1R, COL1A1–ITGA1, and IL6–IL6R. CellOracle and virtual-knockout analyses provided complementary transcription-factor and perturbation hypotheses. Together, these results prioritize experimentally testable axes but do not demonstrate physical communication or causal regulation (Figure 6; Supplementary Table 4).')
addp(doc,'The convergence of spatial, ligand–receptor, and regulatory evidence was used to rank hypotheses, not to create a causal chain. A candidate axis will require validation by spatial protein measurements, perturbation in CAF–myeloid co-culture, receptor blockade, and ideally patient-matched tissue analysis. The present analysis identifies which experiments are most informative.')
heading(doc,'Clinical prioritization indicates association but limited independence',2)
addp(doc,'In TCGA-BLCA, the program was associated with overall survival in a univariable Cox model (HR 1.57, 95% CI 1.22–2.02, P=0.0005). After adjustment for CAF score and AJCC stage, the association was attenuated and not statistically significant (HR 1.24, 95% CI 0.79–1.95, P=0.357). Bootstrap analysis gave a median HR of 1.228 with a 95% interval of 0.996–1.497. Candidate-gene prioritization therefore supports further study rather than an independent prognostic or therapeutic claim (Figure 7; Supplementary Table 4).')
addp(doc,'The attenuation after adjustment indicates that the program may capture stromal burden or disease extent that overlaps with established covariates. Accordingly, CCDC80 is not presented as an independent prognostic biomarker. The clinical analysis instead provides a quantitative boundary for how strongly the program can be interpreted using the available retrospective datasets.')
heading(doc,'Robustness and integrated evidence define the interpretation boundary',2)
addp(doc,'Leave-one-cohort-out analyses, bootstrap performance assessment, calibration summaries, and cross-evidence matrices were used to assess robustness. These analyses support reproducibility of the state-associated program and its relationship with stromal and myeloid context, but the small number of directly profiled patients and the observational design constrain clinical inference. The integrated evidence is therefore classified as discovery, external association, and mechanism-hypothesis generation (Figures 8 and 9).')

heading(doc,'Discussion',1)
addp(doc,'This study defines a CCDC80-associated ECM program within a multi-cohort bladder-cancer CAF atlas and evaluates its relationship with the myeloid microenvironment. The central result is not that CCDC80 alone drives immune remodeling, but that CCDC80 marks a reproducible matrix-associated CAF state whose 19-gene program is detectable across independent analytical layers. The program was supported by single-cell state annotation, patient-level CAF–myeloid associations, external transcriptomic analyses, spatial summaries, and complementary computational hypotheses.')
addp(doc,'The results extend the view that CAFs are functionally heterogeneous and that matrix organization, inflammatory signaling, antigen presentation, and vascular features can occupy partially overlapping state axes [12–14] and [15–17]. The integrated atlas provides a framework for separating CAF abundance from CAF state. This distinction is important because a bulk score can reflect both the number of fibroblasts and the activity of a matrix-associated program.')
addp(doc,'CCDC80 should be interpreted within this state framework. Prior studies support biological relevance of CCDC80/DRO1 in tumor and stromal contexts, but they do not establish that the gene is a universal CAF driver. The present results add a bladder-cancer CAF-specific association and show why multigene state scores are preferable to single-marker claims when expression is sparse or context-dependent [29–31].')
addp(doc,'The association with myeloid scores is biologically plausible but should not be converted into a recruitment mechanism. Fibroblast–macrophage interactions can be reciprocal and context-dependent [18,43]. Expression-based ligand–receptor inference depends on transcript abundance, database annotation, and aggregation assumptions [27,42]. Spatial co-variation and ranked LR candidates therefore identify priorities for multiplex imaging, organoid co-culture, and perturbation experiments rather than causal pathways.')
addp(doc,'The clinical analysis also requires restraint. The univariable TCGA association was attenuated after adjustment for CAF score and stage, and the directly profiled CCDC80 cohorts included only 23 patients. External validation supports directionality across datasets but cannot replace prospective validation. Future studies should test the program with prespecified cutoffs, treatment covariates, molecular subtype, PD-L1, tumor mutational burden, and established immune-response biomarkers [49–51] and [52,53].')
addp(doc,'Several limitations apply. The study is observational and cannot establish that CCDC80 or any program gene causes ECM remodeling or myeloid variation. Bulk validation cannot fully resolve cellular sources. Spatial analyses were available for limited cohorts and showed patient-dependent stability. Virtual-knockout and CellOracle outputs are model-dependent. Finally, cohort and platform heterogeneity, small subgroup sizes, and incomplete clinical covariates limit generalizability.')
addp(doc,'The most direct next step is experimental validation of the prioritized CAF–myeloid axes. CAF-enriched cultures or organoid co-cultures can test whether perturbing CCDC80-associated program activity changes myeloid chemotaxis, macrophage state, or matrix-dependent migration. Multiplex imaging can test whether the spatial associations occur at protein level. Prospective cohorts can then determine whether the program adds information beyond CAF abundance, stage, molecular subtype, and established immune biomarkers.')
addp(doc,'In conclusion, the study supports a CCDC80-associated CAF–ECM program as a reproducible stromal state linked to myeloid microenvironment variation in bladder cancer. CCDC80 should be used as a representative state marker, while causal and therapeutic interpretations require experimental and prospective validation.')

heading(doc,'Ethics statement',1); addp(doc,'This study used publicly available, de-identified datasets and did not involve new human participant recruitment or identifiable information. Ethics committee approval was therefore not required.')
heading(doc,'Author contributions',1); addp(doc,'Penghao Li: Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Visualization, and Writing – original draft. Zhihong Zhang: Conceptualization, Supervision, Validation, Resources, Project administration, and Writing – review and editing. Both authors read and approved the final manuscript and accept responsibility for the integrity of the work.')
heading(doc,'Funding',1); addp(doc,'This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.')
heading(doc,'Acknowledgments',1); addp(doc,'The authors thank the investigators and participants who generated and deposited the public datasets used in this study.')
heading(doc,'Conflict of interest',1); addp(doc,'The authors declare that the research was conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest.')
heading(doc,'Data availability statement',1); addp(doc,'The datasets analyzed in this study are available through GEO, TCGA-linked resources, and the published IMvigor210 research resource. Processed tables, analysis scripts, figures, and supplementary materials are available from the corresponding author on reasonable request.')

heading(doc,'References',1)
for r in refs: addp(doc,r)

heading(doc,'Tables',1)
c=readcsv(FINAL/'Figure1_CAF_atlas/Figure1E_cohort_subtype_composition_source.csv')
add_table(doc,'Supplementary Table 1. Integrated BLCA CAF cohort and state composition',[[x['dataset'],x['CAF_subtype_global'],x['n'],f"{100*float(x['fraction']):.1f}%"] for x in c],['Cohort','CAF state','CAF cells','Fraction'],'Counts and within-cohort fractions are derived from the globally annotated integrated CAF object (four cohorts; 5,589 CAF cells).')
s=readcsv(FINAL/'SupplementaryFigure3_program_gene_stability/SuppFig3_source_gene_stability.csv'); l=readcsv(FINAL/'SupplementaryFigure3_program_gene_stability/SuppFig3_source_leave_one_gene_out.csv'); lm={x['gene']:x for x in l}
add_table(doc,'Supplementary Table 2. Stability of the 19-gene CCDC80-associated ECM program',[[x['gene'],x['n_runs'],f"{float(x['mean_FC']):.3f}",f"{float(x['sd_FC']):.3f}",f"{float(x['mean_abs_Z']):.3f}",f"{float(lm[x['gene']]['correlation']):.3f}"] for x in s],['Gene','Repeat runs','Mean FC','SD of FC','Mean |Z|','Leave-one-out Spearman r'],'Mean FC, SD, and mean absolute Z are calculated across three virtual-knockout repeats. The final column is the correlation between the full score and the score after removing the indicated gene.')
f=readcsv(FINAL/'Figure5_external_validation/Figure5_source_forest.csv'); co=readcsv(FINAL/'Figure5_external_validation/Figure5_source_correlations.csv')
rows=[[x['label'],x['cohort'],'','HR',f"{float(x['HR']):.2f} ({float(x['lo']):.2f}–{float(x['hi']):.2f})",''] for x in f]
rows += [[x['metric'],x['cohort'],x['n'],'Spearman ρ',f"{float(x['rho']):.2f}",'<0.001' if float(x['p'])<.001 else f"{float(x['p']):.3f}"] for x in co]
add_table(doc,'Supplementary Table 3. External validation of program associations',rows,['Evidence','Cohort','N','Estimate','Value (95% CI or ρ)','P value'],'HRs summarize external program survival associations; Spearman correlations summarize associations with CAF and myeloid scores. Confidence intervals were not imputed when unavailable.')
cox=readcsv(FINAL/'Figure7_clinical_prioritization/Figure7_source_Cox.csv'); ca=readcsv(FINAL/'Figure7_clinical_prioritization/Figure7_source_candidate_priority.csv'); lr=readcsv(FINAL/'Figure6_spatial_mechanistic/Figure6_source_LR_top8.csv')
rows=[["Program Cox model",x['term'],'','','%.2f (%.2f–%.2f)'%(float(x['HR']),float(x['lo']),float(x['hi'])),x['p'],''] for x in cox]
rows += [["Candidate gene",x['gene'],x['n'],x['events'],'%.2f (%.2f–%.2f)'%(float(x['HR']),float(x['CI_lower']),float(x['CI_upper'])),x['p_value'],x['FDR']] for x in ca]
rows += [["LR hypothesis",x['pair'],x['n_patients'],'',f"evidence score {float(x['evidence_score']):.0f}",'',x['FDR_LR_myeloid_score']] for x in lr]
add_table(doc,'Supplementary Table 4. Clinical prioritization and mechanism-hypothesis evidence',rows,['Evidence','Candidate/term','N','Events','HR (95% CI) or score','P value','FDR'],'Candidate-gene rows summarize TCGA clinical prioritization. LR rows are hypothesis-generating rankings and should not be interpreted as causal communication evidence.')

heading(doc,'Figures',1)
figs=[(1,'Figure1_CAF_atlas','Figure 1. Integrated bladder-cancer CAF atlas.'),(2,'Figure2_CAF_programs','Figure 2. CAF state programs.'),(3,'Figure3_CCDC80_ECM','Figure 3. CCDC80-associated ECM program definition.'),(4,'Figure4_CAF_myeloid','Figure 4. CAF–myeloid association.'),(5,'Figure5_external_validation','Figure 5. External validation.'),(6,'Figure6_spatial_mechanistic','Figure 6. Spatial and regulatory mechanism hypotheses.'),(7,'Figure7_clinical_prioritization','Figure 7. Clinical prioritization.'),(8,'Figure8_robustness_incremental','Figure 8. Robustness and incremental value.'),(9,'Figure9_integrated_evidence','Figure 9. Integrated evidence.')]
figs += [(10,'SupplementaryFigure1_QC_integration','Supplementary Figure 1. QC and integration audit.'),(11,'SupplementaryFigure2_annotation_stability','Supplementary Figure 2. CAF annotation stability.'),(12,'SupplementaryFigure3_program_gene_stability','Supplementary Figure 3. 19-gene program stability.')]
for i,dname,cap in figs:
    addp(doc,cap,'Normal')
    if i>1: doc.add_page_break()
    fp=FINAL/dname/(('Figure'+str(i)+'.tif') if False else '')
    candidates=list((FINAL/dname).glob('*.tif'))
    if not candidates: candidates=list((FINAL/dname).glob('*.tiff'))
    if candidates: doc.add_picture(str(candidates[0]),width=Inches(6.4))

for name in ['CCDC80_CAF_ECM_myeloid_Frontiers_clean.docx','CCDC80_CAF_ECM_myeloid_Frontiers_redline.docx']:
    doc.save(OUT/name)
print(OUT)
