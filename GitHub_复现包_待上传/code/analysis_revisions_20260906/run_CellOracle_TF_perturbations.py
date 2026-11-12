import os
import numpy as np
import pandas as pd
import celloracle as co

BASE = '/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project/analysis_revisions_20260906/CellOracle_GSE135337'
ORACLE_PATH = os.path.join(BASE, 'GSE135337_CAF_oracle_postGRN.celloracle.oracle')
LINKS_PATH = os.path.join(BASE, 'GSE135337_CAF_CCDC80_state_links.celloracle.links')
OUT = os.path.join(BASE, 'CellOracle_TF_perturbation_summary.csv')

# Candidates are selected from state-differential contextual GRNs, using
# multi-target support and effect size; they are upstream TF perturbations,
# not direct CCDC80 perturbation.
candidates = ['STAT1', 'RXRA', 'PRDM1', 'IRF1', 'ATF3', 'KLF5']
targets = ['CSF1','PDGFA','IL6','CXCL12','COL1A1','COL3A1','FN1','TGFB1','TGFB2',
           'CCDC80','FSTL1','MFAP5','FBN1','COL5A1']

oracle = co.load_hdf5(ORACLE_PATH)
links = co.load_hdf5(LINKS_PATH)
links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True, verbose_level=0)

rows = []
for tf in candidates:
    if tf not in oracle.adata.var_names:
        rows.append({'TF': tf, 'status': 'not_in_input'})
        continue
    oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3,
                          ignore_warning=True, clip_delta_X=True)
    delta = oracle.adata.to_df(layer='delta_X')
    row = {'TF': tf, 'status': 'ok', 'n_cells': delta.shape[0]}
    for g in targets:
        if g in delta.columns:
            row[f'delta_mean_{g}'] = float(delta[g].mean())
            row[f'delta_median_{g}'] = float(delta[g].median())
            row[f'frac_decrease_{g}'] = float((delta[g] < 0).mean())
    for state in oracle.adata.obs['CCDC80_state'].unique():
        ix = oracle.adata.obs['CCDC80_state'].values == state
        row[f'n_{state}'] = int(ix.sum())
        for g in targets:
            if g in delta.columns:
                row[f'delta_mean_{state}_{g}'] = float(delta.loc[ix, g].mean())
    rows.append(row)

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(res[['TF','status','n_cells'] + [c for c in res.columns if c.startswith('delta_mean_') and '_' not in c[len('delta_mean_'):]]].to_string(index=False))
print('saved', OUT)
