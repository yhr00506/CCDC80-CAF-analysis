import os, re, gzip, tarfile
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

BASE='/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project'
TAR=os.path.join(BASE,'02_external_scRNA_downloads/GSE135337/GSE135337_RAW.tar')
SCORE=os.path.join(BASE,'analysis_revisions_20260906/GSE135337_patient_level_CCDC80_TME/GSE135337_cell_level_scores.csv')
OUT=os.path.join(BASE,'analysis_revisions_20260906/CellOracle_GSE135337')
os.makedirs(OUT,exist_ok=True)
scores=pd.read_csv(SCORE,usecols=['cell_id','sample_id','patient_id','CCDC80','CAF'])
scores=scores[scores['CAF'].astype(bool)].copy()
chunks=[]; obs=[]
with tarfile.open(TAR,'r') as tf:
    members=[m for m in tf.getmembers() if re.search(r'_BC[1-7]_gene_cell_exprs_table\.(txt|xls)\.gz$',m.name)]
    for m in members:
        sid=re.sub(r'_gene_cell_exprs_table.*$','',os.path.basename(m.name))
        ids=scores.loc[scores.sample_id.eq(sid),'cell_id'].tolist()
        if not ids: continue
        f=tf.extractfile(m)
        with gzip.open(f,'rt') as gz:
            header=pd.read_csv(gz,sep='\t',nrows=0).columns.tolist()
        wanted=['Gene_ID','Symbol']+[x for x in ids if x in header]
        f=tf.extractfile(m)
        with gzip.open(f,'rt') as gz:
            x=pd.read_csv(gz,sep='\t',usecols=wanted)
        x['Symbol']=x['Symbol'].astype(str).str.upper()
        x=x[x['Symbol'].ne('') & ~x['Symbol'].duplicated()]
        genes=x['Symbol'].tolist(); mat=x.drop(columns=['Gene_ID','Symbol']).to_numpy(dtype=np.float32)
        chunks.append(pd.DataFrame(mat,index=genes,columns=x.columns[2:])); chunks[-1].columns=[sid+"|"+c for c in chunks[-1].columns]
        oo=scores[scores.sample_id.eq(sid)].set_index("cell_id").reindex(x.columns[2:]); oo.index.name="cell_id"; oo=oo.reset_index()
        oo["cell_id"]=sid+"|"+oo["cell_id"].astype(str); obs.append(oo)
expr=pd.concat(chunks,axis=1,join='outer').fillna(0).astype(np.float32)
ob=pd.concat(obs,ignore_index=True).set_index('cell_id').loc[expr.columns]
varmean=expr.mean(axis=1); var=expr.var(axis=1); keep=var.sort_values(ascending=False).head(min(3000,len(var))).index
expr=expr.loc[keep]
adata=ad.AnnData(X=sparse.csr_matrix(expr.T.to_numpy(dtype=np.float32)),obs=ob,var=pd.DataFrame(index=expr.index))
adata.layers['raw_count']=adata.X.copy()
adata.write_h5ad(os.path.join(OUT,'GSE135337_CAF_CellOracle_input.h5ad',),compression='gzip')
audit=pd.DataFrame([{'n_cells':adata.n_obs,'n_genes':adata.n_vars,'n_patients':adata.obs.patient_id.nunique(),'n_samples':adata.obs.sample_id.nunique(),'CCDC80_detected_cells':int((adata.obs.CCDC80>0).sum())}])
audit.to_csv(os.path.join(OUT,'input_audit.csv'),index=False)
print(audit.to_string(index=False))
