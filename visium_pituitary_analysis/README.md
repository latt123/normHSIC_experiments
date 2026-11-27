Code to reproduce Section 5.2, 5.3 of "A Normalised Kernel Independence Test". 

Experiments comparing the output of spatially variable gene (SVG) selection methods normHSIC ("A Normalised Kernel Independence Test") and SPARKX (Zhu, J., S. Sun, and X. Zhou. 2021. “SPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies.”).

Analysis of 12 spatial single cell RNAseq visium datasets of pituitary tumors. 


# Structure
- `simulate_data.py' Generating simulate data sets following experiment design in Gretton et al. 2007, A Kernel Statistical Test of Independence. Datasets simulated varying sample size, dimension, source distribution and parameter for strength of dependence signal. 
- `snakemake/` snakemake pipeline for running both methods (normHSIC and SPARKX) on all 12 datasets. 
used for kernels. 
- `analysis/` 
    - `analysis.qmd` Quarto file summarising analysis
    - `top_genes.py` Python script to generate csv files of lists of the highest ranked (SVG) for each dataset and method.
    - `pval_top_genes_all.py` Python script to generate csv files of the p-values and relative ranks assigned by both methods for the 91 genes that are consistently ranked top 500 across all slides by SVG by normHSIC. 
- `plots_for_paper.ipynb` Jupyter notebook to create Figures used in Section 5.2, 5.3 of "A Normalised Kernel Independence Test". 
- `preparing_data.py` Preprocessing step applied to anndata files. Filtering to remove genes expressed in less than 200 spots and spots with less than 200 unique genes expressed.
- `summarise_results.py` Combine results for both methods and add rank of p-value of each gene among datasets. 


# Using snakemake

1. **Simulate data:**
```bash 
python simulate_data.py
```

2. Ensure **snakemake** and **conda** are installed

3. **Create `knit_env` and `SPARKX_env` conda environment:**
```bash
conda env create -f knit_env.yml -n knit_env
```
```bash
conda env create -f SPARKX_env.yml -n SPARKX_env
```
Environment requirements are in `config` folder.

4. **Configure paths:** Edit paths in Snakefile and ``.smk' rules files.

5. **Run snakemake**
