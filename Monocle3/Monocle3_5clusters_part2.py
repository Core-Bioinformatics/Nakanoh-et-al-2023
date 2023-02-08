import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import rankdata


meta = pd.read_csv('meta_5clusters_monocle_pseudotime.csv',index_col=0)

rank_numbers = rankdata(meta.monocle_pseudotime, method='dense')
rank_numbers[:4]

meta['pseudotime_rank'] = rank_numbers

pseudo_order = meta.sort_values(by=['monocle_pseudotime'],ascending = True)

difference = pseudo_order['monocle_pseudotime'].diff()

reverse_list=[]

start_value= np.max(pseudo_order['monocle_pseudotime'])

continue_value = 0

for i, value in enumerate(difference):
    if i == 0:
        reverse_list.append(start_value)
        
        continue_value = start_value
    else:   
        continue_value = continue_value - value

        reverse_list.append(continue_value)
        
pseudo_order['python_reverse_pseudotime'] = reverse_list

pseudo_order['clusters_medres'] = pseudo_order['clusters_medres'].astype(str)


reverse_order = pseudo_order.sort_values(by=['python_reverse_pseudotime'],ascending = True)

reverse_rank_numbers = rankdata(reverse_order.python_reverse_pseudotime, method='dense')


reverse_order['python_reverse_rank'] = reverse_rank_numbers

reverse_order.to_csv('meta_5clusters_monocle_pseudotime_reverse.csv')


