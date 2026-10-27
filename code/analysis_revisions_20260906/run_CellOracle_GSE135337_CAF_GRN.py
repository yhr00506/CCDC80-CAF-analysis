import os
import numpy as np
import scanpy as sc
import celloracle as co
from sklearn.decomposition import TruncatedSVD

base='/data/y01014/diliyaer/BLCA_CAF_project/GDF15_CAF_project'
inp=os.path.join(base,'analysis_revisions_20260906/CellOracle_GSE135337/GSE135337_CAF_CellOracle_input.h5ad')
out=os.path.join(base,'analysis_revisions_20260906/CellOracle_GSE135337')
adata=sc.read_h5ad(inp)
pmap=__import__('pandas').read_csv(os.path.join(base,'analysis_revisions_20260906/GSE135337_patient_level_CCDC80_TME/GSE135337_patient_level_CCDC80_TME.csv')).set_index('patient_id')['CCDC80_high_CAF'].to_dict()
patient_mean=adata.obs['patient_id'].astype(str).map(pmap).astype(bool)
adata.obs['CCDC80_state']=np.where(patient_mean,'CCDC80_high','CCDC80_low').astype(str)
X=adata.layers['raw_count']
adata.obsm['X_pca']=TruncatedSVD(n_components=30,random_state=123).fit_transform(X)
oracle=co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata,cluster_column_name='CCDC80_state',embedding_name='X_pca')
base_grn=co.data.load_human_promoter_base_GRN()
oracle.import_TF_data(TF_info_matrix=base_grn)
oracle.perform_PCA(n_components=30)
oracle.knn_imputation(n_pca_dims=30,k=30,balanced=True,b_sight=240,b_maxl=120,n_jobs=4)
oracle.to_hdf5(os.path.join(out,'GSE135337_CAF_oracle_preGRN.celloracle.oracle'))
links=oracle.get_links(cluster_name_for_GRN_unit='CCDC80_state',alpha=10,bagging_number=10,verbose_level=1,test_mode=False,n_jobs=4)
links.to_hdf5(os.path.join(out,'GSE135337_CAF_CCDC80_state_links.celloracle.links'))
oracle.to_hdf5(os.path.join(out,'GSE135337_CAF_oracle_postGRN.celloracle.oracle'))
for k,v in links.links_dict.items():
    v.to_csv(os.path.join(out,'GRN_'+str(k)+'.csv'),index=False)
print('cells',adata.n_obs,'genes',adata.n_vars,'states',adata.obs['CCDC80_state'].value_counts().to_dict())
print('GRN states',list(links.links_dict.keys()))
