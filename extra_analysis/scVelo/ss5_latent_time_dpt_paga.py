import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import cellrank as cr
import pickle
from cellrank.tl.kernels import VelocityKernel
import seaborn as sns

with open('saved_object/ss4_dynamic_adata_macrostates.pkl', 'rb') as file:
    # A new file will be created
    dynamic_adata = pickle.load(file)
    
    
umap_33_cluster = dynamic_adata.obs[dynamic_adata.obs['terminal_states']=='33']
    
    
    
root_list=[]

end_list=[]

for index, row in dynamic_adata.obs.iterrows():
#     print(index)
    if row['initial_states'] =='29':
        root_list.append(True)
    else:
        root_list.append(False)
        
    if row['terminal_states'] =='33' and row['UMAP_1'] > 2:
        end_list.append(True)
    else:
        end_list.append(False)
    
    
dynamic_adata.obs['root_key'] = root_list
dynamic_adata.obs['end_key'] = end_list

#########################

scv.tl.latent_time(dynamic_adata,
                   root_key='root_key', end_key='end_key')
    
    
########################
    
scv.tl.paga(
    dynamic_adata,
    groups="cell_cluster",
#     root_key="initial_states_probs",
#     end_key="end_key",
    use_time_prior="latent_time",
)
    
    
#######################


root_idx = np.where(dynamic_adata.obs["latent_time"] == 0)[0][0]
dynamic_adata.uns["iroot"] = root_idx
sc.tl.dpt(dynamic_adata)
    
#########################

with open('saved_object/ss5_adata_latent_time_paga_dpt_noEPI.pkl', 'wb') as file:
      
    # A new file will be created
    pickle.dump(dynamic_adata, file)

meta = dynamic_adata.obs.to_csv('saved_object/meta_ss5_latent_time_paga_dpt_noEPI.csv')
    
    
    
    
    
    
    
    