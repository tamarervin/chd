# Coronal Hole Detection and Modeling

Methods based on [CHMAP](https://zenodo.org/record/5039440) from Predictive Science Inc. (PSI) 

## Environment Setup
This requires a custom environment.
Check for conda-forge:
```
conda config --show channels
```
If needed, add conda-forge:
```
conda config --append channels conda-forge
```
Now we are ready to build the custom 'chd' conda environment. Navigate to the folder containing the configuration file and run the below snippet.
```
conda env create --file conda_env.yml
```
## [Carrington Map CHD]:
- download files from time period of interest
- prepare images for algorithm
- loop through maps to create minimum merge composite map
  - apply PSI's EZSEG algorithm for coronal hole detection
  - combine maps using PSI's minimum merge mapping techniques
