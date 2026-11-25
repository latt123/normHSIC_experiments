Code to reproduce Section 4 of "A Normalised Kernel Independence Test". 

Experiments measuring the empirical power and level of kernel independence tests normHSIC ("A Normalised Kernel Independence Test") and HSIC (Gretton et al. 2007, A Kernel Statistical Test of Independence). 

# Structure
- `simulate_data.py' Generating simulate data sets following experiment design in Gretton et al. 2007, A Kernel Statistical Test of Independence. Datasets simulated varying sample size, dimension, source distribution and parameter for strength of dependence signal. 

- `snakemake/` Snakemake pipeline for running both independence tests (normHSIC and HSIC) on all datasets. 
- `snakemake_proportion_train` Snakemake pipeline for running normHSIC for various train set sizes (normHSIC requires splitting dataset)
- `snakemake_quantiles` Snakemake pipeline for running both independence tests (normHSIC and HSIC) on all datasets, varying the bandwidth parameters used for kernels. 

# Usage
1. **Simulate data:**
```bash 
python simulate_data.py
```

2. Ensure **Snakemake** and **conda** are installed

3. **Create `knit_env` conda environement:**
```bash
conda env create -f knit_env.yml -n knit_env
```

4. **Configure paths:** Edit paths in Snakefile and ``.smk' rules files.

5. **Run Snakemake**
