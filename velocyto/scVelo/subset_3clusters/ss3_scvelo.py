import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import cellrank as cr

adata_4_1_10 = sc.read('write/raw_scvelo_adata_4_1_10.h5ad')


fil_norm_adata = scv.pp.filter_and_normalize(adata_4_1_10, 
                            min_counts=2, min_counts_u=2,
                           min_cells=10, min_cells_u=10,
                           min_shared_counts= 2,min_shared_cells=10,
                           subset_highly_variable=False, copy= True)

hvg_df = pd.read_csv('hvg3000_seurat_subset_4_1_10.csv',index_col=0)

hvg_fil_adata = fil_norm_adata[:,fil_norm_adata.var_names.isin(hvg_df['x'])]


pca_df = pd.read_csv('pca_df__subset_4_1_10.csv',index_col=0)

umap_df = pd.read_csv('umap_df__subset_4_1_10.csv',index_col=0)


hvg_fil_adata.obsm['X_pca'] = pca_df.values
hvg_fil_adata.obsm['X_umap'] = umap_df.values

mo_hvg_fil_norm_adata = scv.pp.moments(hvg_fil_adata, 
                                       n_neighbors=30, n_pcs=30,
                                       method="umap", use_rep='X_pca',
                                       use_highly_variable=False, copy= True)

scv.tl.recover_dynamics(mo_hvg_fil_norm_adata, copy=False, n_jobs=8)


scv.tl.velocity(mo_hvg_fil_norm_adata,
                mode='dynamical', copy=False, use_highly_variable=False)
scv.tl.velocity_graph(mo_hvg_fil_norm_adata, n_jobs= 8)


cr.tl.terminal_states(mo_hvg_fil_norm_adata, 
                     cluster_key="clusters_medres",
                    n_jobs= 1)

cr.tl.initial_states(mo_hvg_fil_norm_adata, 
                     cluster_key="clusters_medres",
                    n_jobs= 1)

scv.tl.latent_time(mo_hvg_fil_norm_adata,
                   root_key="initial_states_probs", end_key="terminal_states_probs")


scv.tl.paga(
    mo_hvg_fil_norm_adata,
    groups="clusters_medres",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="latent_time",
)


root_idx = np.where(mo_hvg_fil_norm_adata.obs["latent_time"] == 0)[0][0]


mo_hvg_fil_norm_adata.uns["iroot"] = root_idx

sc.tl.dpt(mo_hvg_fil_norm_adata)

with open('write/adata_latent_time_paga_dpt_3clusters.pkl', 'wb') as file:
      
    # A new file will be created
    pickle.dump(mo_hvg_fil_norm_adata, file)



