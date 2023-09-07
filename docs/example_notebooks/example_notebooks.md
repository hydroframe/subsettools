Example notebooks
=================

Here is a collection of examples using the subsettools API:

CONUS1 Example Notebooks:
The example notebooks in the following section subset datasets from the baseline CONUS1 model (conus1_baseline_mod and conus1_domain) stored on Hydrodata.
It is a box domain over the contiguous United States. It does not cover areas near the coasts. This model has 4 soil layers (top 2 meters) and 1 geologic layer that is 100 m thick. The cell size is 1 kilometer. For more information about the CONUS1 datasets being subset in these example notebooks please refer to the following DOIs: 10.5194/gmd-14-7223-2021; 10.1016/j.advwatres.2015.04.008; 10.1002/2014WR016774; 10.1002/2015GL066916.  

1. conus1_subsetting.ipynb 
This notebook walks through an example of subsetting a HUC6 from the CONUS1 domain described above. This example will subset everything needed to do a transient run with ParFlow-CLM. 
This is includes all hydrogeologic datasets and climate forcing data from NLDAS2 (see DOIs: 10.5194/gmd-14-7223-2021; 10.1002/2016GL069964). All of the data is written to a folder for the specified days to run. 
Additionally, this example uses the template runscript conus1_pfclm_transient_solid.yaml stored in this repository and edits it to correspond with the domain subset. Next, the example will walk through how to use additional functions to set-up and alter runs and at the end of the notebook will perform the designed simulation. The result will be model output pressure and saturation pfbs according to the timeseries specified. 
2. conus1_subsetting_spinup.ipynb
This notebook will subset the same HUC6 from the CONUS1 domain as in the conus1_subsetting.ipynb example. However, it does not subset climate forcing data as it is an example of spin-up or model initialization. So, it only runs with ParFlow. The template runscript used (conus1_pf_spinup_solid.yaml) has the necessary timing keys to run with a growth timestep and a longterm recharge mask (PmE) as opposed to climate forcing data. The example notebook is not set to run for the entire period it would normally take to spin-up a model, you are encouraged to do your own evaluation if your model is at steady state especially if you make any changes to your modeling domain. 

```{toctree}
:maxdepth: 1
conus1_subsetting.ipynb
conus1_subsetting_spinup.ipynb
```

CONUS2 Example Notebooks:
The following example notebooks subset datasets from the CONUS2 model stored on Hydrodata.The CONUS2 modeling domain covers all areas of the contiguous United States. The model has 4 soil and 6 geologic layers. The cell size is 1 kilometer. For more information about the CONUS2 domain and simulations, refer to (PUBLICATIONS).(NEED TO REDO THIS ONCE IT IS IN  HYDRODATA). The example notebooks perform exactly the same functions as the correspondign CONUS1 notebook above, except for the CONUS2 domain, so they will not be described in as much detail. 

1. conus2_subsetting.ipynb
2. conus2_subsetting_spinup.ipynb
