Code to reproduce the analysis in "A Normalised Kernel Test for Independence". 

# Data Availability
- Experimental Visium Spatial scRNAseq datasets. 

# Folders
- `config`: Environment files to recreate conda environments to run normHSIC (`knit_env`) and SPARKX (`SPARKX_env`).
- `hsic_comparison`: Analysis comparing normHSIC and HSIC on simulated datasets. Corresponding to Section 4 of "A Normalised Kernel Test for Independence".
- `visium_pituitary_analysis`: Analysis of normHSIC and SPARKX on Visium Spatial scRNAseq samples of pituitary tumours. Corresponding to Sections 5.2 and 5.3 of "A Normalised Kernel Test for Independence".

# knit package 
NormHSIC implemented with [`knit`](https://github.com/latt123/knit) (Kernel Normalised Independence Test) package.

Install from GitHub:
```console
$ pip install git+https://github.com/latt123/knit
```

# License
See [License](License).