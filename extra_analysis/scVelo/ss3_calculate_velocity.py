import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import cellrank as cr

from cellrank.tl.kernels import VelocityKernel


adata = sc.read('saved_object/ss2_mo_scvelo_adata_29_33_34_35_zhai2022.h5ad')

scv.tl.recover_dynamics(adata, copy=False, n_jobs=10)

scv.tl.velocity(adata,
                mode='dynamical', copy=False, use_highly_variable=False)

scv.tl.velocity_graph(adata, n_jobs= 10)

fig, ax = plt.subplots(1, 1,figsize=(8, 6))
scv.pl.velocity_embedding_stream(adata, 
                                 basis='umap', color='cell_cluster',
                                show=False, ax=ax)
plt.xlabel('UMAP1')  
plt.ylabel('UMAP2')  

plt.savefig('saved_plot/umap_velocity_stream__Zhai2022_subset_noEPI.png',bbox_inches='tight')

plt.show()

fig, ax = plt.subplots(1, 1,figsize=(8, 6))
scv.pl.velocity_embedding_grid(adata, basis='umap', color='cell_cluster',ax=ax,show=False)

plt.xlabel('UMAP1')  
plt.ylabel('UMAP2')  

plt.savefig('saved_plot/umap_velocity_grid__Zhai2022_subset_noEPI.pdf',bbox_inches='tight')

plt.show()


sc.write('saved_object/ss3_dynamic_adata_29_33_34_35_zhai2022.h5ad',adata)



