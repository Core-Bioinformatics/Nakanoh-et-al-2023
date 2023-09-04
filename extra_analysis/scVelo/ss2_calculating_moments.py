import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import cellrank as cr

adata_subset = sc.read('saved_object/raw_scvelo_adata_29_33_34_35_zhai2022.h5ad')

##################

fil_norm_adata = scv.pp.filter_and_normalize(adata_subset, 
                            min_counts=2, min_counts_u=2,
                           min_cells=10, min_cells_u=10,
                           min_shared_counts= 2,min_shared_cells=10,
                           subset_highly_variable=False, copy= True)

################

df_umap = fil_norm_adata.obs.loc[:,['UMAP_1','UMAP_2']].copy()
df_umap.shape


df_pca = pd.read_csv('../scVelo/saved_object/pca_56636cells_50dim_zhai2022_fixed_index.csv',index_col=0)

df_pca = df_pca.loc[fil_norm_adata.obs_names,:]

Counter(fil_norm_adata.obs_names == df_umap.index)

################

fil_norm_adata.obsm['X_pca'] = df_pca.values
fil_norm_adata.obsm['X_umap'] = df_umap.values

#################

mo_fil_norm_adata = scv.pp.moments(fil_norm_adata, 
                                       n_neighbors=30, n_pcs=30,
                                       method="umap", use_rep='X_pca',
                                       use_highly_variable=False, copy= True)




sc.write('saved_object/ss2_mo_scvelo_adata_29_33_34_35_zhai2022.h5ad',mo_fil_norm_adata)


