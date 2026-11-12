import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

ROOT = os.path.dirname(__file__)
INFILE = os.path.join(ROOT, 'CCDC80_robustness_audit', 'robust_patient_level_metrics.csv')
OUTDIR = os.path.join(ROOT, 'patient_level_CCDC80_myeloid_association')
os.makedirs(OUTDIR, exist_ok=True)
df = pd.read_csv(INFILE)
df['CCDC80_high'] = df['CCDC80_high'].astype(str).str.lower().eq('true').astype(int)
df['cohort'] = df['cohort'].astype(str)
predictors = ['CCDC80_z', 'CAF_CCDC80_mean', 'CCDC80_high']
responses = ['Myeloid_fraction', 'n_Myeloid', 'Myeloid_score_mean']
rows = []
for cohort, d in df.groupby('cohort'):
    for x in predictors:
        for y in responses:
            z = d[[x, y]].dropna()
            if len(z) < 4 or z[x].nunique() < 2 or z[y].nunique() < 2:
                continue
            rho, p = spearmanr(z[x], z[y])
            rows.append({'cohort': cohort, 'n': len(z), 'predictor': x, 'response': y, 'rho': rho, 'p': p})
for x in predictors:
    for y in responses:
        z = df[[x, y, 'cohort', 'CAF_fraction']].dropna()
        if len(z) < 8 or z[x].nunique() < 2 or z[y].nunique() < 2:
            continue
        rho, p = spearmanr(z[x], z[y])
        rows.append({'cohort': 'pooled', 'n': len(z), 'predictor': x, 'response': y, 'rho': rho, 'p': p})
corr = pd.DataFrame(rows)
if len(corr):
    corr['FDR'] = multipletests(corr['p'], method='fdr_bh')[1]
corr.to_csv(os.path.join(OUTDIR, 'patient_level_spearman_correlations.csv'), index=False)
models = []
for y in responses:
    z = df[['CCDC80_z', y, 'cohort', 'CAF_fraction']].dropna().copy()
    if len(z) >= 10 and z['CCDC80_z'].nunique() > 1 and z[y].nunique() > 1:
        fit = smf.ols(f'Q("{y}") ~ CCDC80_z + CAF_fraction + C(cohort)', data=z).fit(cov_type='HC3')
        models.append({'response': y, 'n': len(z), 'beta_CCDC80_z': fit.params.get('CCDC80_z', np.nan),
                       'p_CCDC80_z': fit.pvalues.get('CCDC80_z', np.nan), 'r2': fit.rsquared})
model_df = pd.DataFrame(models)
if len(model_df):
    model_df['FDR'] = multipletests(model_df['p_CCDC80_z'], method='fdr_bh')[1]
model_df.to_csv(os.path.join(OUTDIR, 'pooled_adjusted_models.csv'), index=False)
df.to_csv(os.path.join(OUTDIR, 'analysis_input_patient_metrics.csv'), index=False)
with open(os.path.join(OUTDIR, 'README.txt'), 'w') as f:
    f.write('Patient-level CCDC80 CAF to myeloid association analysis.\n')
    f.write('Input is the robustness-audit table; CAF state assignment is independent of myeloid outcomes.\n')
    f.write('Spearman analyses are cohort-stratified and pooled; adjusted OLS includes CAF_fraction and cohort.\n')
    f.write('Exploratory because independent patient count is small (GSE293189=13, GSE135337=7).\n')
print(corr.to_string(index=False))
print('\nAdjusted models')
print(model_df.to_string(index=False))
