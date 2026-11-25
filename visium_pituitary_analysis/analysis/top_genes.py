import pandas as pd

input_dir = 'results/output'
out_dir = 'results/top_N_genes'

def get_genes(df, method, n_svgs):
    if method == 'normhsic':
        df = df.sort_values(['test_stat'], ascending=False)
        df = df.head(n_svgs).reset_index()
        df = df[['gene']]

    # For SPARKX
    if method == 'SPARKX':
        df = df.sort_values(['combinedPval'], ascending=False)
        df = df.head(n_svgs).reset_index()
        df = df[['Row.names']]
    return df


method_list = ['normhsic', 'SPARKX']
N_genes_list = [200, 500, 1000, 2000]
dataset_list = ['V1', 'V2', 'V3', 'V4', 'Vf1', 'Vf2', 'Vf3', 'Vf4', 'Vf9', 'Vf10', 'Vf11', 'Vf12']

for method in method_list:
    for dataset in dataset_list:
        df = pd.read_csv(f'{input_dir}/{method}/{dataset}.csv')
        for N_genes in N_genes_list:
            df_genes = get_genes(df, method, N_genes)
            df_genes.to_csv(f'{out_dir}/{method}/{dataset}_{N_genes}.csv', index=False)