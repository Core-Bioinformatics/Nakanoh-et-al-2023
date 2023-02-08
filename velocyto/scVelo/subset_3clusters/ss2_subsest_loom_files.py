import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc


root_dir='/path/to/data/'

path_looms = ["velocyto/SLX-21814_SITTA3_H7T7CDRX2.loom", 
              "velocyto/SLX-21814_SITTB3_H7T7CDRX2.loom",
              "velocyto/SLX-21814_SITTC3_H7T7CDRX2.loom", 
              "velocyto/SLX-21814_SITTD3_H7T7CDRX2.loom"]
sample_ids = ["1", "2", "3", "4"]


adata_rna_A3 = scv.read(f'{root_dir}{path_looms[0]}')

adata_rna_B3 = scv.read(f'{root_dir}{path_looms[1]}')

adata_rna_C3 = scv.read(f'{root_dir}{path_looms[2]}')

adata_rna_D3 = scv.read(f'{root_dir}{path_looms[3]}')
list_adata = [adata_rna_A3,
              adata_rna_B3,
              adata_rna_C3,
              adata_rna_D3]
for i, temp_adata in enumerate(list_adata):
    temp_adata.obs_names = [x.split(':')[1][:-1] + "-" + sample_ids[i] for x in temp_adata.obs_names]
    temp_adata.var_names_make_unique()
    
    print(temp_adata.obs_names[:3])
    
    
com_adata = list_adata[0].concatenate(list_adata[1],list_adata[2],list_adata[3],index_unique=None)

meta = pd.read_csv('meta__subset_4_1_10.csv',index_col=0)


adata_4_1_10 = com_adata[meta.index,:]

adata_4_1_10.obs['sample'] = meta['sample'].values
adata_4_1_10.obs['clusters_medres'] = meta['clusters_medres'].values

adata_4_1_10.obs['clusters_medres']= adata_4_1_10.obs['clusters_medres'].astype("category")
adata_4_1_10.obs['sample']= adata_4_1_10.obs['sample'].astype("category")

sc.write('write/raw_scvelo_adata_4_1_10.h5ad',adata_4_1_10)



