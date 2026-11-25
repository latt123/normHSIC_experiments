Code to reproduce the analysis in "A Normalised Kernel Test for Independence". 

# Data Availability
- Experimental Visium Spatial scRNAseq datasets. 

# Folders
- `{un_test` Environment files to recreate conda environements to run normHSIC (`knit_env`) and SPARKX (`SPARKX_env`).
- `hsic_comparision` Anaylsis comparing normHSIC and HSIC on simulated datasets. Corresponding to Section 4 of "A Normalised Kernel Test for Independence".
- `visium_pituitary_analysis` Analysis of normHSIC and SPARKX on Visium Spatial scRNAseq samples of pituitary tumours. Correspoding to Sections 5.2 and 5.3 of "A Normalised Kernel Test for Independence".

# knit package 
NormHSIC implemented with `knit` (Kernel Normalised Independence Test) package.
```bash
    python -m pip install git+https://gitbio.ens-lyon.fr/lattwood/knit
```
