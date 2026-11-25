# Combine the results into a single dataframe for each slide
import numpy as np
import pandas as pd
import scanpy as sc

# Slides and methods
slides = ['V1', 'V2', 'V3', 'V4','Vf1', 'Vf2', 'Vf3', 'Vf4', 'Vf9', 'Vf10', 'Vf11', 'Vf12']
svg_methods = ['SPARKX', 'normhsic']


results_dir = 'results/output'
data_dir = 'results/data'

for slide in slides:
    adata = sc.read_h5ad(f'{data_dir}/{slide}.h5ad')
    adata.var['sparsity'] = 1 - np.array((np.sum(adata.X != 0, axis = 0)/adata.X.shape[0])).reshape(-1) # number of zeros values
    adata.var['total_counts'] = np.array(np.sum(adata.X, axis = 0), int).reshape(-1)
    df = adata.var[['n_cells', 'total_counts', 'sparsity']]
    df.index.name = 'gene'
    
    
    for svg_method in svg_methods:
        if svg_method == 'SPARKX':
            df_2 = pd.read_csv(f'{results_dir}/{svg_method}/{slide}.csv').sort_values(by = 'adjustedPval', ascending = True)
            df_2.index = df_2['Row.names']
        else:
            df_2 = pd.read_csv(f'{results_dir}/{svg_method}/{slide}.csv').sort_values(by = 'test_stat', ascending = False)
            df_2.index = df_2['gene']
        df_2[f'rank_{svg_method}'] = np.arange(len(df_2)) + 1
        df_2.index.name = 'gene'
        if svg_method == 'SPARKX':
            df_2 = df_2[[f'rank_{svg_method}', 'adjustedPval']]
            df_2.columns = [f'rank_{svg_method}', f'pval_{svg_method}']
        else:
            df_2 = df_2[[f'rank_{svg_method}', 'pval', 'test_stat']]
            df_2.columns = [f'rank_{svg_method}', f'pval_{svg_method}', f'test_stat_{svg_method}']

        df = pd.merge(df, df_2, left_index=True, right_index=True)
    df.to_csv(f'results/summary/{slide}.csv')
