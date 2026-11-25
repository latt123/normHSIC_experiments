
import numpy as np 
import pandas as pd

results_dir = '../results/summary'

slides = ['V1', 'V2', 'V3', 'V4', 'Vf1', 'Vf2', 'Vf3', 'Vf4', 'Vf9', 'Vf10', 'Vf11', 'Vf12']


slide = slides[0]
columns = ['gene']
df = pd.read_csv(f'{results_dir}/{slide}.csv')[['gene', 'rank_normhsic']]
columns.append(f'rank_{slide}')
df.columns = columns


for slide in slides[1:]:
    df_2 = pd.read_csv(f'{results_dir}/{slide}.csv')[['gene', 'rank_normhsic']]
    df = pd.merge(df, df_2, how = 'outer', on = 'gene')
    columns.append(f'rank_{slide}')
    df.columns = columns

df_normhsic = df

N_genes = 500
df_normhsic[f'Slides in Top {N_genes} Genes'] = (df_normhsic.values[:, 1:(len(slides)+1)] <= N_genes).sum(axis = 1) 
df_normhsic['Gene in N Slides'] = len(slides) - df_normhsic.iloc[:,1:13].isna().sum(axis = 1)
genes = list(df_normhsic[(df_normhsic['Slides in Top 500 Genes'] == 12)]['gene'].values)


slide = slides[0]
columns = ['gene']
df = pd.read_csv(f'{results_dir}/{slide}.csv')[['gene', 'rank_SPARKX']]
columns.append(f'rank_{slide}')
df.columns = columns
for slide in slides[1:]:
    df_2 = pd.read_csv(f'{results_dir}/{slide}.csv')[['gene', 'rank_SPARKX']]
    df = pd.merge(df, df_2, how = 'outer', on = 'gene')
    columns.append(f'rank_{slide}')
    df.columns = columns

df_SPARKX = df

N_genes = 500
df_SPARKX[f'Slides in Top {N_genes} Genes'] = (df_SPARKX.values[:, 1:(len(slides)+1)] <= N_genes).sum(axis = 1) 
df_SPARKX['Gene in N Slides'] = len(slides) - df_SPARKX.iloc[:,1:13].isna().sum(axis = 1)


genes = list(df_normhsic[(df_normhsic['Slides in Top 500 Genes'] == 12)]['gene'].values)
df_rank_normhsic = pd.DataFrame(genes, columns = ['gene'])
df_pval_normhsic = pd.DataFrame(genes, columns = ['gene'])
df_rank_SPARKX = pd.DataFrame(genes, columns = ['gene'])
df_pval_SPARKX = pd.DataFrame(genes, columns = ['gene'])
input_dir = '../results/output'

for slide in slides:
    method = 'normhsic'
    df_2 = pd.read_csv(f'{input_dir}/{method}/{slide}.csv', index_col  = 0)
    df_2[slide] = len(df_2) -  np.argsort(np.argsort(df_2['test_stat']) + 1)
    df_2 = df_2[df_2['gene'].isin(genes)]
    df_rank_normhsic = df_rank_normhsic.merge(df_2[['gene', slide]], on = 'gene')
    df_2[slide] = df_2['pval']
    df_pval_normhsic = df_pval_normhsic.merge(df_2[['gene', slide]], on =  'gene')

    method = 'SPARKX'
    df_2 = pd.read_csv(f'{input_dir}/{method}/{slide}.csv', index_col  = 0)
    df_2['gene'] = df_2['Row.names']
    df_2[f'{slide}'] = np.argsort(np.argsort(df_2['adjustedPval']) + 1)
    df_2 = df_2[df_2['gene'].isin(genes)]
    df_rank_SPARKX = df_rank_SPARKX.merge(df_2[['gene', slide]], on = 'gene')
    df_2[slide] = df_2['adjustedPval']
    df_pval_SPARKX = df_pval_SPARKX.merge(df_2[['gene', slide]], on = 'gene')


# df
df_rank_normhsic['method'] = 'normHSIC'
df_rank_SPARKX['method'] = 'SPARKX'
df_ranks = pd.concat((df_rank_normhsic, df_rank_SPARKX))

df_pval_normhsic['method'] = 'normHSIC'
df_pval_SPARKX['method'] = 'SPARKX'
df_pvals = pd.concat((df_pval_normhsic, df_pval_SPARKX))

df_ranks.to_csv('rank_top_genes_all_slides.csv')
df_pvals.to_csv('pval_top_genes_all_slides.csv')
