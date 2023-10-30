Example notebooks
=================

Below are descriptions regarding example Jupyter notebookes using the subsettools API. Currently, only there are only examples for subsetting the CONUS1 domain but notebooks supporting the subsetting of the CONUS2 domain will soon be supported. 

### Information Regarding the CONUS1 domain

The example notebooks in the following section subset datasets from the baseline CONUS1 model (conus1_baseline_mod and conus1_domain) stored on the database Hydrodata. 
Please refer to [Dataset Information](https://hydroframe-ml.github.io/readthedocs/tables/dataset.html) for details about CONUS datasets contained on Hydrodata.

CONUS1 is a box domain over the contiguous United States. It does not cover areas near the coasts. This model has 4 soil layers (top 2 meters) and 1 geologic layer that is 100 m thick. The cell size is 1 kilometer. 

Please refer to and cite the following DOIs if you want more details or to release work with the datasets used in the example notebooks:
1. [10.5194/gmd-14-7223-2021](https://gmd.copernicus.org/articles/14/7223/2021/)
2. [10.1016/j.advwatres.2015.04.008](https://www.sciencedirect.com/science/article/pii/S0309170815000822)
3. [10.1002/2014WR016774](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014WR016774)
4. [10.1002/2015GL066916](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066916)  

### Example Notebook Descriptions

#### Doing a transient simulation with ParFlow-CLM

This notebook walks through an example of subsetting a HUC8 from the CONUS1 domain. This example will subset everything needed to do a transient run with ParFlow-CLM. 
This is includes all hydrogeologic datasets and climate forcing data from NLDAS2 (see additional DOI: [10.1002/2016GL069964](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL069964)). All of the data is written to a folder for the specified days to run. This example uses the template runscript conus1_pfclm_transient_solid.yaml and edits it to correspond with the domain subset. It also sets-up and performs the designed simulation. The result will be model output pressure and saturation pfbs according to the days specified.

```{toctree}
:maxdepth: 1
conus1_subsetting_transient.ipynb
```

#### Performing a model initialization (spin up) with ParFlow

This notebook will subset a HUC8 from the CONUS1 domain as in the conus1_subsetting.ipynb example. However, it does not subset climate forcing data as it is an example of model initialization or spin up. So, it only runs with ParFlow. The template runscript used is conus1_pf_spinup_solid.yaml which has the necessary timing keys to run with a growth timestep and a longterm recharge mask (PmE) as opposed to a constant time step and time varying climate forcing data. The example notebook is not set to run for the entire period it would normally take to spin-up a model, you are encouraged to do your own evaluation to determine if your model is at steady state. 

```{toctree}
:maxdepth: 1
conus1_subsetting_spinup.ipynb
```

