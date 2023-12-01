Example notebooks
=================

In addition to the tutorials provided in the `How To` section we provide several Jupyter notebooks that demonstrate complete end-to-end workflows for common tasks.  

**1. Performing a model initialization (spin up) with ParFlow**: This notebook will subset a HUC8 from the CONUS1 domain and will create a ParFlow run script to run a spinup simulation uses a constant recharge forcing to achieve a steady state water table depth (*NOTE*: This run does not include land surface simulation).

```{toctree}
:maxdepth: 1
conus1_subsetting_spinup.ipynb
```


**2. Transient simulation with ParFlow-CLM**: This notebook walks through an example of subsetting a HUC8 from the CONUS1 domain. This example will subset everything needed to do a transient run with ParFlow-CLM and create a run script. This simulation is designed to start from a steady state groundwater condition.  Here we subset the initial conditions from the pre-generated national steady state water table depth simulations, but the starting point could also be generated from the spin up workbook.  

```{toctree}
:maxdepth: 1
conus1_subsetting_transient.ipynb
```

