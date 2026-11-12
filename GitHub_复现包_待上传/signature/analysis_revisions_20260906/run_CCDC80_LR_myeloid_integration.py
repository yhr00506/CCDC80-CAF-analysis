import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = os.path.dirname(__file__)
MET = pd.read_csv(os.path.join(ROOT, 'GSE135337_patient_level_CCDC80_TME', 'GSE135337_patient_level_CCDC80_TME.csv'))
LR = pd.read_csv(os.path.join(ROOT, 'GSE135337_CAF_myeloid_LR', 'GSE135337_patient_LR_potential.csv'))
OUT = os.path.join(ROOT, 'CCDC80_LR_myeloid_integration')
os.makedirs(OUT, exist_ok=True)

MET['patient_id'] = MET['patient_id'].astype(str)
LR['patient_id'] = LR['patient_id'].astype(str)
for c in ['CAF_ligand', 'Myeloid_receptor', 'LR_potential']:
    LR[c] = pd.to_numeric(LR[c], errors='coerce').fillna(0.0)
dat = LR.merge(MET[['patient_id','CAF_CCDC80_mean','Myeloid_score','Myeloid_fraction','CCDC80_high_CAF']], on='patient_id', how='inner')
dat['CCDC80_high_CAF'] = dat['CCDC80_high_CAF'].astype(str).str.lower().eq('true').astype(int)
rows = []
for pair, d in dat.groupby('pair'):
    ligand = d['CAF_ligand']; receptor = d['Myeloid_receptor']; pot = d['LR_potential']
    def sp(x, y):
        if x.nunique() < 2 or y.nunique() < 2: return np.nan, np.nan
        return spearmanr(x, y)
    r1,p1 = sp(ligand, d['CAF_CCDC80_mean'])
    r2,p2 = sp(ligand, d['Myeloid_score'])
    r3,p3 = sp(receptor, d['Myeloid_score'])
    r4,p4 = sp(pot, d['Myeloid_score'])
    hi, lo = d.loc[d.CCDC80_high_CAF==1, 'LR_potential'], d.loc[d.CCDC80_high_CAF==0, 'LR_potential']
    if len(hi) and len(lo):
        try: u,p5 = mannwhitneyu(hi, lo, alternative='two-sided')
        except ValueError: p5 = np.nan
    else: p5 = np.nan
    rows.append({'pair':pair, 'ligand':d.ligand.iloc[0], 'receptor':d.receptor.iloc[0], 'n_patients':len(d),
                 'mean_LR_potential':pot.mean(), 'rho_ligand_CCDC80':r1, 'p_ligand_CCDC80':p1,
                 'rho_ligand_myeloid_score':r2, 'p_ligand_myeloid_score':p2,
                 'rho_receptor_myeloid_score':r3, 'p_receptor_myeloid_score':p3,
                 'rho_LR_myeloid_score':r4, 'p_LR_myeloid_score':p4,
                 'median_LR_high':hi.median() if len(hi) else np.nan, 'median_LR_low':lo.median() if len(lo) else np.nan,
                 'p_LR_high_vs_low':p5})
res = pd.DataFrame(rows)
for col in [c for c in res.columns if c.startswith('p_')]:
    res['FDR_'+col[2:]] = multipletests(res[col].fillna(1), method='fdr_bh')[1]
res['evidence_score'] = ((res['rho_ligand_CCDC80'] > 0).astype(int) +
                         (res['rho_LR_myeloid_score'] > 0).astype(int) +
                         (res['p_ligand_CCDC80'] < .1).astype(int) +
                         (res['p_LR_myeloid_score'] < .1).astype(int))
res = res.sort_values(['evidence_score','mean_LR_potential'], ascending=[False,False])
res.to_csv(os.path.join(OUT, 'CCDC80_LR_myeloid_integrated_ranking.csv'), index=False)
dat.to_csv(os.path.join(OUT, 'CCDC80_LR_myeloid_patient_matrix.csv'), index=False)
with open(os.path.join(OUT, 'README.txt'),'w') as f:
    f.write('Integrates GSE135337 patient-level CCDC80 CAF state, CAF ligand, myeloid receptor, LR potential, and myeloid score.\n')
    f.write('Missing ligand/receptor expression is treated as zero, consistent with the upstream LR calculation.\n')
    f.write('This is exploratory (7 independent patients); p-values are descriptive and require external validation.\n')
print(res[['pair','n_patients','mean_LR_potential','rho_ligand_CCDC80','p_ligand_CCDC80','rho_LR_myeloid_score','p_LR_myeloid_score','evidence_score']].head(15).to_string(index=False))
