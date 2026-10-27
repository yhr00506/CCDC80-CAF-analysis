"""Permutation-based spatial neighborhood validation for GSE171351.

This implements the core logic of spatialCooccur without requiring a package
installation: the observed spatial graph is fixed, while myeloid scores are
permuted within each patient/sample. It tests high-high co-occurrence and
ECM-to-neighbor myeloid enrichment.
"""
import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(__file__)
INFILE = os.path.join(HERE, "GSE171351_spatial_validation", "GSE171351_spot_scores.csv")
OUTDIR = os.path.join(HERE, "GSE171351_spatial_validation")
N_PERM = 2000
SEED = 20260907

df = pd.read_csv(INFILE)
rng = np.random.default_rng(SEED)
rows = []

def corr_from_ranks(x, y):
    x = rankdata(x)
    y = rankdata(y)
    x = x - x.mean()
    y = y - y.mean()
    den = np.sqrt(np.sum(x*x) * np.sum(y*y))
    return float(np.sum(x*y) / den) if den else np.nan

for sample_id, d in df.groupby("sampleID", sort=True):
    d = d.reset_index(drop=True)
    keep = d[["ECM_signature", "Myeloid_score", "array_row", "array_col"]].notna().all(axis=1)
    d = d.loc[keep].reset_index(drop=True)
    n = len(d)
    if n < 20:
        continue

    e = d["ECM_signature"].to_numpy(float)
    m = d["Myeloid_score"].to_numpy(float)
    xy = d[["array_row", "array_col"]].to_numpy(float)
    k = min(6, n - 1)
    _, nn = cKDTree(xy).query(xy, k=k + 1)
    nn = nn[:, 1:]

    e_hi = e >= np.quantile(e, 0.75)
    m_hi = m >= np.quantile(m, 0.75)
    obs_high_high = float(np.mean(e_hi & m_hi))
    obs_edge_rate = float(np.mean(m_hi[nn][e_hi])) if e_hi.any() else np.nan
    m_rank = rankdata(m)
    obs_neighbor_rho = corr_from_ranks(e, m_rank[nn].mean(axis=1))

    null_high_high = np.empty(N_PERM)
    null_edge_rate = np.empty(N_PERM)
    null_neighbor_rho = np.empty(N_PERM)
    for i in range(N_PERM):
        perm = rng.permutation(n)
        mp = m[perm]
        mp_hi = mp >= np.quantile(mp, 0.75)
        null_high_high[i] = np.mean(e_hi & mp_hi)
        null_edge_rate[i] = np.mean(mp_hi[nn][e_hi]) if e_hi.any() else np.nan
        null_neighbor_rho[i] = corr_from_ranks(e, rankdata(mp)[nn].mean(axis=1))

    p_high_high = (1 + np.sum(null_high_high >= obs_high_high)) / (N_PERM + 1)
    p_edge = (1 + np.sum(null_edge_rate >= obs_edge_rate)) / (N_PERM + 1)
    p_neighbor = (1 + np.sum(null_neighbor_rho >= obs_neighbor_rho)) / (N_PERM + 1)
    rows.append({
        "sampleID": sample_id,
        "Patient": d.get("Patient", pd.Series(["NA"])).iloc[0] if "Patient" in d else "NA",
        "n_spots": n,
        "n_neighbors": k,
        "ECM_high_definition": "top_quartile",
        "Myeloid_high_definition": "top_quartile",
        "observed_high_high_fraction": obs_high_high,
        "null_high_high_mean": null_high_high.mean(),
        "null_high_high_sd": null_high_high.std(ddof=1),
        "high_high_perm_p": p_high_high,
        "observed_ECM_to_myeloid_neighbor_rate": obs_edge_rate,
        "null_neighbor_rate_mean": null_edge_rate.mean(),
        "neighbor_enrichment_ratio": obs_edge_rate / null_edge_rate.mean() if null_edge_rate.mean() else np.nan,
        "neighbor_enrichment_perm_p": p_edge,
        "observed_ECM_neighbor_myeloid_spearman": obs_neighbor_rho,
        "null_neighbor_spearman_mean": null_neighbor_rho.mean(),
        "neighbor_spearman_perm_p": p_neighbor,
    })

res = pd.DataFrame(rows)
pcols = ["high_high_perm_p", "neighbor_enrichment_perm_p", "neighbor_spearman_perm_p"]
for c in pcols:
    res["FDR_" + c.replace("_perm_p", "")] = multipletests(res[c], method="fdr_bh")[1]
res.to_csv(os.path.join(OUTDIR, "GSE171351_spatial_permutation_summary.csv"), index=False)

with open(os.path.join(OUTDIR, "SPATIAL_PERMUTATION_STATUS.txt"), "w") as f:
    f.write("Completed patient-internal spatial permutation analysis.\n")
    f.write(f"N_PERM={N_PERM}; k-nearest-neighbor graph; top-quartile labels.\n")
    f.write("Spatial graph fixed; myeloid scores permuted within sample.\n")
    f.write("This is permutation-based association evidence, not causal proof.\n")

print(res.to_string(index=False))
