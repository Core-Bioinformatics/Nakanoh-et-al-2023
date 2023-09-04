import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import cellrank as cr
import pickle
from cellrank.tl.kernels import VelocityKernel

with open('saved_object/ss5_adata_latent_time_paga_dpt_noEPI.pkl', 'rb') as file:
      
    # Call load method to deserialze
    dynamic_adata = pickle.load(file)
    
root_idx = dynamic_adata.uns["iroot"] 
root_idx


fig, ax = plt.subplots(3, 2,figsize=(11, 11))
ax = ax.flatten()

scv.pl.scatter(dynamic_adata,
    color=root_idx,    fontsize=16, basis='umap',
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title="root cell diffusion pseudotime",show=False,ax=ax[0]
)

scv.pl.scatter(
    dynamic_adata,
    color="latent_time",  basis='umap',
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title="latent time",show=False,ax=ax[1]
)

scv.pl.scatter(
    dynamic_adata,
    color="dpt_pseudotime", basis='umap',
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title="dpt pseudotime",show=False,ax=ax[2]
)

scv.pl.scatter(
    dynamic_adata,
    color="velocity_pseudotime", basis='umap',
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title="velocity_pseudotime",show=False,ax=ax[3]
)

sc.pl.umap(dynamic_adata, color='cell_cluster', show=False,ax=ax[4], legend_loc='on data')
ax[5].axis('off')

plt.savefig('saved_plot/UMAP_ss6_different_pseudotime_noEPI.pdf',bbox_inches='tight')

plt.show()


