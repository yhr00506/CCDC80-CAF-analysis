import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

ROOT=os.path.dirname(__file__); OUT=os.path.join(ROOT,'patient_internal_colocalization_proxy'); os.makedirs(OUT,exist_ok=True)
files=[('GSE135337',os.path.join(ROOT,'GSE135337_patient_level_CCDC80_TME','GSE135337_cell_level_scores.csv')),('GSE293189',os.path.join(ROOT,'GSE293189_patient_level_CCDC80_TME','GSE293189_cell_level_scores.csv'))]
allrows=[]
for cohort,fn in files:
 d=pd.read_csv(fn); d['patient_id']=d['patient_id'].astype(str)
 if 'label' in d.columns:
  d['CAF']=d['label'].eq('CAF'); d['Myeloid']=d['label'].eq('Myeloid')
 else:
  for c in ['CAF','Myeloid']: d[c]=d[c].astype(str).str.lower().eq('true')
 rows=[]
 for pid,x in d.groupby('patient_id'):
  caf=x[x.CAF]; my=x[x.Myeloid]
  rows.append({'cohort':cohort,'patient_id':pid,'n_cells':len(x),'n_CAF':len(caf),'n_Myeloid':len(my),'CAF_fraction':x.CAF.mean(),'Myeloid_fraction':x.Myeloid.mean(),'CAF_score_mean':caf.CAF_score.mean() if len(caf) else np.nan,'Myeloid_score_mean':my.Myeloid_score.mean() if len(my) else np.nan,'CAF_Myeloid_fraction_product':x.CAF.mean()*x.Myeloid.mean(),'CAF_Myeloid_score_product':caf.CAF_score.mean()*my.Myeloid_score.mean() if len(caf) and len(my) else np.nan})
 allrows.append(pd.DataFrame(rows))
res=pd.concat(allrows,ignore_index=True)
labels=pd.concat([pd.read_csv(os.path.join(ROOT,'GSE135337_patient_level_CCDC80_TME','GSE135337_patient_level_CCDC80_TME.csv')).assign(cohort='GSE135337'),pd.read_csv(os.path.join(ROOT,'GSE293189_patient_level_CCDC80_TME','GSE293189_patient_level_CCDC80_TME.csv')).assign(cohort='GSE293189')],ignore_index=True)
labels['patient_id']=labels['patient_id'].astype(str); labels['CCDC80_high_CAF']=labels['CCDC80_high_CAF'].astype(str).str.lower().eq('true').astype(int)
res=res.merge(labels[['cohort','patient_id','CAF_CCDC80_mean','CCDC80_high_CAF']],on=['cohort','patient_id'],how='left'); res.to_csv(os.path.join(OUT,'patient_internal_colocalization_proxy.csv'),index=False)
rows=[]
for cohort,x in res.groupby('cohort'):
 for y in ['CAF_Myeloid_fraction_product','CAF_Myeloid_score_product','Myeloid_score_mean']:
  z=x[['CAF_CCDC80_mean',y]].dropna()
  if len(z)>=4 and z.iloc[:,0].nunique()>1 and z.iloc[:,1].nunique()>1:
   r,p=spearmanr(z.iloc[:,0],z.iloc[:,1]); rows.append({'cohort':cohort,'n':len(z),'proxy':y,'rho_CCDC80':r,'p':p})
cor=pd.DataFrame(rows); cor['FDR']=multipletests(cor.p,method='fdr_bh')[1] if len(cor) else []
cor.to_csv(os.path.join(OUT,'proxy_CCDC80_correlations.csv'),index=False)
with open(os.path.join(OUT,'README.txt'),'w') as f:f.write('Patient-internal co-occurrence proxy, not true spatial transcriptomics. Products combine patient-level CAF and myeloid abundance/activity.\n')
print(cor.to_string(index=False)); print(res.groupby('cohort').size())
