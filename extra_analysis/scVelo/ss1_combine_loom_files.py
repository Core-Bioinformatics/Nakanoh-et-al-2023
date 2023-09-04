import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import glob

dir_meta='/path/to/data'

meta = pd.read_csv('meta_new_index.csv',index_col=0)

root_dir='/path/to/data'

list_loom = glob.glob(f'{root_dir}/**/velocyto/*.loom')

list_adata = []

for filename in list_loom:
    temp_adata = scv.read(filename, cache=True)

    list_adata.append(temp_adata)

#####################
    
dict_sample_embryo={'SRR17489634':'CS8-e1','SRR17489635':'CS8-e2',
                    'SRR17489636':'CS9-e1','SRR17489637':'CS9-e2',
                    'SRR17489638':'CS11-e1','SRR17489639':'CS11-ys1',
                    'SRR17489640':'CS11-e2'}

dict_sample_embryo_2={'SRR17489634':'E1','SRR17489635':'E2',
                    'SRR17489636':'E1','SRR17489637':'E2',
                    'SRR17489638':'E1','SRR17489639':'Y1',
                    'SRR17489640':'E1'}

dict_sample_stage={'SRR17489634':'E20','SRR17489635':'E20',
                    'SRR17489636':'E23','SRR17489637':'E23',
                    'SRR17489638':'E26','SRR17489639':'E26',
                    'SRR17489640':'E29'}

#####################

dict_adata=dict()

for i , temp_adata in enumerate(list_adata):

    
    temp_index_list=[]
    
    for index in temp_adata.obs_names:
        new_index = index.split(':')[-1]
        stage_name = dict_sample_stage[index.split(':')[0]]
        embryo_name = dict_sample_embryo_2[index.split(':')[0]]
        temp_index_list.append(f'{stage_name}-{embryo_name}_{new_index[:-1]}')
        
    sample_list=[dict_sample_embryo[index.split(':')[0]]]*len(temp_index_list)

    stage_list=[stage_name]*len(temp_index_list)

    temp_adata.obs_names = temp_index_list
    
    temp_adata.obs['sample'] = sample_list

    temp_adata.obs['stage'] = stage_list
    
    dict_adata[index.split(':')[0]] = temp_adata
    
    print(len(temp_index_list))
    print(len(set(temp_index_list)))
    print('==========')
    
for key, value in dict_adata.items():
    
    print(key)
    
    value.var_names_make_unique()   
    
    
####################
    
com_adata = dict_adata['SRR17489634'].concatenate(dict_adata['SRR17489635'],
                                                  dict_adata['SRR17489636'],
                                                  dict_adata['SRR17489637'],
                                                  dict_adata['SRR17489638'],
                                                  dict_adata['SRR17489639'],
                                                  dict_adata['SRR17489640'],
                                                  index_unique=None)

####################
    
subset_meta=meta[meta['cell_type'].isin(['ECT','SE1','SE2','AM'])]
    
adata_subset = com_adata[com_adata.obs_names.isin(subset_meta['new_index']),:]

    
subset_meta.index = subset_meta['new_index'].to_list()
    
indata_meta = subset_meta.loc[adata_subset.obs_names,:]
    
adata_subset.obs = indata_meta
    
    
adata_subset.obs['cell_type'] = adata_subset.obs['cell_type'].astype('category')

adata_subset.obs['cell_cluster'] = adata_subset.obs['cell_cluster'].astype('category')
adata_subset.obs['stage'] = adata_subset.obs['stage'].astype('category')

adata_subset.obs['sample'] = adata_subset.obs['sample'].astype('category')
    
    
    
sc.write('saved_object/raw_scvelo_adata_29_33_34_35_zhai2022.h5ad',adata)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    