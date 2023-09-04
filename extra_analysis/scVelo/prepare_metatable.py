import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import glob

dir_meta='/path/to/data'

meta = pd.read_csv(f'{dir_meta}/MFE56636-meta_fixed_index.csv',index_col=0)


orig_meta = pd.read_csv(f'{dir_meta}/MFE56636-meta.csv',index_col=0)
orig_meta


embryo_list=[]

stage_list=[]

new_index=[]

for index in orig_meta.index:
    
    first_part = index.split('_')[0]
    
    temp_stage = first_part.split('-')[0]
    
    temp_embryo = first_part.split('-')[-1]
    
    embryo_list.append(temp_embryo)
    
    stage_list.append(temp_stage)
    
    new_index.append(index[:-2])


meta['new_index'] = new_index
meta.to_csv('meta_new_index.csv')






