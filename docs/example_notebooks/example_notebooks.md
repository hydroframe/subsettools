Example notebooks
=================

Below are descriptions regarding example Jupyter notebookes using the subsettools API. Currently, only there are only examples for subsetting the CONUS1 domain but notebooks supporting the subsetting of the CONUS2 domain will soon be supported. 

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

#### Short examples

```{toctree}
:maxdepth: 1
short_examples.ipynb
```
