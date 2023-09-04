import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import scanpy as sc
import cellrank as cr
import pickle
from cellrank.tl.kernels import VelocityKernel


adata = sc.read('saved_object/ss3_dynamic_adata_29_33_34_35_zhai2022.h5ad')


cr.tl.terminal_states(adata, 
                     cluster_key="cell_cluster",n_states=1,
                    n_jobs= 1)

cr.tl.initial_states(adata, 
                     cluster_key="cell_cluster",n_states=1,
                    n_jobs= 1)


with open('saved_object/ss4_dynamic_adata_macrostates.pkl', 'wb') as file:
    # A new file will be created
    pickle.dump(adata, file)





