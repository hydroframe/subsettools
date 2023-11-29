Example notebooks
=================

Below are descriptions regarding example Jupyter notebooks using the subsettools API. Currently, only there are only examples for subsetting the CONUS1 domain but notebooks supporting the subsetting of the CONUS2 domain will soon be supported. 

#### Doing a transient simulation with ParFlow-CLM

This notebook walks through an example of subsetting a HUC8 from the CONUS1 domain. This example will subset everything needed to do a transient run with ParFlow-CLM. 

```{toctree}
:maxdepth: 1
conus1_subsetting_transient.ipynb
```

#### Performing a model initialization (spin up) with ParFlow

This notebook will subset a HUC8 from the CONUS1 domain as in the above example. However, it does not subset climate forcing data as it is an example of model initialization or spin up. 

```{toctree}
:maxdepth: 1
conus1_subsetting_spinup.ipynb
```
