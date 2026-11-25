import warnings
warnings.filterwarnings(action='ignore')
import os
import scanpy as sc
import squidpy as sq
import numpy as np

output_file = 'results/data'
os.makedirs(output_file, exist_ok=True)

slide_list = ['V1', 'V2', 'V3', 'V4', 'Vf1', 'Vf2', 'Vf3', 'Vf4', 'Vf9', 'Vf10', 'Vf11', 'Vf12']

for slide in slide_list:
    input_file = f'../../data/visium_crcl/{slide}/outs/'
    adata = sq.read.visium(input_file)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=200)
    sc.pp.filter_cells(adata, min_genes=200)
    adata.var['total_counts'] = np.array(np.sum(adata.X, axis = 0), int).reshape(-1)
    adata.obs['total_counts'] = np.array(np.sum(adata.X, axis = 1), int).reshape(-1)
    adata.write(f"{output_file}/{slide}.h5ad")