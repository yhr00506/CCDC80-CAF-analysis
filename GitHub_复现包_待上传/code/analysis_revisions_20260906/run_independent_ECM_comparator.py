"""Compare the CCDC80-associated 19-gene program with an independent ECM score."""
import os
import pandas as pd
from scipy.stats import spearmanr

BASE = "/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project"
REV = os.path.join(BASE, "analysis_revisions_20260906")
OUT = os.path.join(REV, "independent_ECM_comparator")
os.makedirs(OUT, exist_ok=True)

cohorts = {
    "TCGA": (os.path.join(REV, "CCDC80_ECM_program_validation", "TCGA_ECM_program_scores.csv"), os.path.join(REV, "ECM_CAF_signature_validation", "TCGA_BLCA_ECM_signature_input.csv")),
    "GSE13507": (os.path.join(REV, "CCDC80_ECM_program_validation", "GSE13507_ECM_program_scores.csv"), os.path.join(REV, "ECM_CAF_signature_validation", "GSE13507_ECM_signature_input.csv")),
    "GSE32894": (os.path.join(REV, "CCDC80_ECM_program_validation", "GSE32894_ECM_program_scores.csv"), os.path.join(REV, "ECM_CAF_signature_validation", "GSE32894_ECM_signature_input.csv")),
}

rows = []
merged_all = []
for cohort, (program_file, signature_file) in cohorts.items():
    p = pd.read_csv(program_file).rename(columns={"program_score": "CCDC80_ECM_program"})
    s = pd.read_csv(signature_file).rename(columns={"ECM_signature": "independent_ECM_signature"})
    m = p.merge(s, on="sample_id", how="inner")
    m.insert(0, "cohort", cohort)
    merged_all.append(m)
    x = m[["CCDC80_ECM_program", "independent_ECM_signature"]].dropna()
    rho, pval = spearmanr(x.iloc[:, 0], x.iloc[:, 1])
    rows.append({"cohort": cohort, "n_merged": len(x), "spearman_rho": rho, "p_value": pval})

merged = pd.concat(merged_all, ignore_index=True)
merged.to_csv(os.path.join(OUT, "CCDC80_program_vs_independent_ECM_sample_scores.csv"), index=False)
pooled = merged[["CCDC80_ECM_program", "independent_ECM_signature"]].dropna()
rho, pval = spearmanr(pooled.iloc[:, 0], pooled.iloc[:, 1])
rows.append({"cohort": "pooled_unadjusted", "n_merged": len(pooled), "spearman_rho": rho, "p_value": pval})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "program_signature_correlations.csv"), index=False)

current = pd.read_csv(os.path.join(REV, "CCDC80_ECM_program_validation", "external_GEO_ECM_program_survival.csv"))
independent = pd.read_csv(os.path.join(REV, "ECM_CAF_signature_validation", "ECM_signature_OS_Cox.csv"))
current["score_type"] = "CCDC80_ECM_program"
independent["score_type"] = "independent_ECM_signature"
pd.concat([current, independent], ignore_index=True, sort=False).to_csv(os.path.join(OUT, "program_signature_survival_comparison.csv"), index=False)

with open(os.path.join(OUT, "RUN_STATUS.txt"), "w") as f:
    f.write("Independent ECM comparator completed.\n")
    f.write("Comparator: COL1A1/COL3A1/FN1 ECM signature, independent of CCDC80 virtual knockout.\n")
    f.write("Primary analysis: sample-level Spearman correlation with the CCDC80-associated 19-gene program.\n")
    f.write("Existing survival summaries were placed side-by-side; this does not establish independence or causality.\n")

print(pd.DataFrame(rows).to_string(index=False))
