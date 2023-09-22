Example notebooks
=================

Below are descriptions regarding example Jupyter notebookes using the subsettools API. Currently, only there are only examples for subsetting the CONUS1 domain but notebooks supporting the subsetting of the CONUS2 domain will soon be supported. 

### Information Regarding the CONUS1 domain

The example notebooks in the following section subset datasets from the baseline CONUS1 model (conus1_baseline_mod and conus1_domain) stored on the database Hydrodata. 
Please refer to [Dataset Information](https://hydroframe-ml.github.io/readthedocs/tables/dataset.html) for more information about CONUS datasets contained on hydrodata.

CONUS1 is a box domain over the contiguous United States. It does not cover areas near the coasts. This model has 4 soil layers (top 2 meters) and 1 geologic layer that is 100 m thick. The cell size is 1 kilometer. Please refer to and cite the following DOIs if you want more details or to release work with the datasets used in the example notebooks: 10.5194/gmd-14-7223-2021; 10.1016/j.advwatres.2015.04.008; 10.1002/2014WR016774; 10.1002/2015GL066916.  

### Example Notebook Descriptions

#### Doing a transient simulation with ParFlow-CLM - conus1_subsetting.ipynb 

This notebook walks through an example of subsetting a HUC6 from the CONUS1 domain. This example will subset everything needed to do a transient run with ParFlow-CLM. 
This is includes all hydrogeologic datasets and climate forcing data from NLDAS2 (see DOIs: 10.5194/gmd-14-7223-2021; 10.1002/2016GL069964). All of the data is written to a folder for the specified days to run. 
Additionally, this example uses the template runscript conus1_pfclm_transient_solid.yaml stored in this repository and edits it to correspond with the domain subset. Next, the example will walk through how to use additional functions to set-up and alter runs and at the end of the notebook will perform the designed simulation. The result will be model output pressure and saturation pfbs according to the timeseries specified.

#### Performing a model initialization (spin up) with ParFlow - conus1_subsetting_spinup.ipynb

This notebook will subset the same HUC6 from the CONUS1 domain as in the conus1_subsetting.ipynb example. However, it does not subset climate forcing data as it is an example of spin-up or model initialization. So, it only runs with ParFlow. The template runscript used (conus1_pf_spinup_solid.yaml) has the necessary timing keys to run with a growth timestep and a longterm recharge mask (PmE) as opposed to climate forcing data. The example notebook is not set to run for the entire period it would normally take to spin-up a model, you are encouraged to do your own evaluation if your model is at steady state especially if you make any changes to your modeling domain. 

```{toctree}
:maxdepth: 1
conus1_subsetting.ipynb
conus1_subsetting_spinup.ipynb
```

