import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import cellrank as cr

meta = pd.read_csv('meta_new_index.csv',index_col=0)


df_pca = pd.read_csv('saved_object/pca_56636cells_50dim_zhai2022.csv',index_col=0)


Counter(meta.index == df_pca.index)


df_pca.index = meta['new_index'].to_list()


df_pca.to_csv('saved_object/pca_56636cells_50dim_zhai2022_fixed_index.csv')












