import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

ROOT=os.path.dirname(__file__); IN='/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/external_CAF_data/spatial_validation_20260906/GSE171351/GSE171351_combined_visium.h5ad'; OUT=os.path.join(ROOT,'GSE171351_spatial_validation'); os.makedirs(OUT,exist_ok=True)
a=ad.read_h5ad(IN); a=a[a.obs['in_tissue'].astype(bool)].copy(); genes=list(a.var_names); ix={g:i for i,g in enumerate(genes)}
def score(glist):
 ids=[ix[g] for g in glist if g in ix]
 if not ids:return np.full(a.n_obs,np.nan)
 x=a.X[:,ids]; x=x.toarray() if hasattr(x,'toarray') else np.asarray(x); lib=np.asarray(a.obs['total_counts']).reshape(-1,1); x=np.log1p(x/np.maximum(lib,1)*1e4); z=(x-x.mean(0))/np.where(x.std(0)>0,x.std(0),1); return np.nanmean(z,axis=1)
for nm,g in {'CCDC80':['CCDC80'],'ECM_signature':['COL1A1','COL3A1','FN1'],'CAF_score':['COL1A1','COL1A2','DCN','LUM','COL3A1','COL5A1','SPARC','FAP','PDGFRA'],'Myeloid_score':['CD68','CD163','C1QA','C1QB','C1QC','LILRB1','FCER1G','TYROBP','APOC1','SPP1','MSR1'],'CSF1':['CSF1'],'IL6':['IL6'],'CXCL12':['CXCL12']}.items(): a.obs[nm]=score(g)
spot=a.obs[['sampleID','Patient','array_row','array_col','CCDC80','ECM_signature','CAF_score','Myeloid_score','CSF1','IL6','CXCL12']].copy(); spot.insert(0,'spot_id',spot.index); spot.to_csv(os.path.join(OUT,'GSE171351_spot_scores.csv'),index=False)
rows=[]
for sid,d in spot.groupby('sampleID'):
 row={'sampleID':sid,'Patient':d.Patient.iloc[0],'n_spots':len(d)}
 for x,y,nm in [('ECM_signature','Myeloid_score','ECM_myeloid'),('CCDC80','Myeloid_score','CCDC80_myeloid'),('CSF1','Myeloid_score','CSF1_myeloid'),('IL6','Myeloid_score','IL6_myeloid'),('CXCL12','Myeloid_score','CXCL12_myeloid')]:
  z=d[[x,y]].dropna(); r,p=spearmanr(z[x],z[y]); row[f'rho_{nm}']=r; row[f'p_{nm}']=p
 row['high_high_ECM_myeloid_fraction']=((d.ECM_signature>=d.ECM_signature.quantile(.75))&(d.Myeloid_score>=d.Myeloid_score.quantile(.75))).mean()
 # nearest-neighbor spatial autocorrelation proxy: correlate each spot's ECM with neighbor-average myeloid score
 xy=d[['array_row','array_col']].to_numpy(); tree=cKDTree(xy); _,nn=tree.query(xy,k=min(7,len(d))); neigh=d.Myeloid_score.to_numpy()[nn[:,1:]].mean(1); row['rho_ECM_neighbor_myeloid']=spearmanr(d.ECM_signature,neigh,nan_policy='omit').statistic; row['p_ECM_neighbor_myeloid']=spearmanr(d.ECM_signature,neigh,nan_policy='omit').pvalue
 rows.append(row)
summary=pd.DataFrame(rows); pcols=[c for c in summary if c.startswith('p_')];
for c in pcols: summary['FDR_'+c[2:]]=multipletests(summary[c].fillna(1),method='fdr_bh')[1]
summary.to_csv(os.path.join(OUT,'GSE171351_spatial_patient_summary.csv'),index=False)
with open(os.path.join(OUT,'README.txt'),'w') as f:f.write('GSE171351 human bladder Visium spatial validation. Scores are spot-level log-normalized z-score averages. Neighbor analysis uses six nearest spatial spots.\n')
print(summary.to_string(index=False))
