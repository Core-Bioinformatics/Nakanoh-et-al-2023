import pandas as pd

meta = pd.read_csv('MFE56636-meta.csv',index_col=0)

new_index = []

for cell_id in meta.index:
    new_id = cell_id.split('_')[-1]
    new_index.append(new_id)
    
    
    
meta.index = new_index
meta.to_csv('MFE56636-meta_fixed_index.csv')
    