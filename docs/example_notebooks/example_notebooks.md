Example notebooks
=================

In addition to the tutorials provided in the `How To` section we provide several Jupyter notebooks that demonstrate complete end-to-end workflows for common tasks. 

You can get the tutorials by cloning the subsettools repo and navigating to the "example_notebooks" folder, or you can download individual notebooks from the GitHub repo.

```bash
git clone https://github.com/hydroframe/subsettools.git
cd subsettools/docs/example_notebooks
```

There are three ways to run the notebooks, as explained in more detail in the [Getting started](https://hydroframesubsettools.readthedocs.io/en/latest/getting_started.html) section.

Each notebook can be run in `Binder`. Click the `launch binder` button at the top of each example notebook.

The notebooks can also be run in `Docker`. Follow the instructions in [Getting started](https://hydroframesubsettools.readthedocs.io/en/latest/getting_started.html) to get the `subsettools` Docker container and run the notebooks.

Finally, the notebooks can be run on your local machine. Clone the `subsettools` repo or download an individual notebook from GitHub, and run the notebooks with a kernel that has the package pip-installed. Note that if you're running locally, in order to complete some workflows, you will also need to have ParFlow installed. Please refer to the [HydroFrame](https://hydroframe.org/) for help getting started with ParFlow.


**1. [Performing a model initialization (spin up) with ParFlow on CONUS1](https://hydroframesubsettools.readthedocs.io/en/latest/example_notebooks/conus1_subsetting_spinup.html)**: This notebook will subset a HUC8 from the CONUS1 domain and will create a ParFlow run script to run a spinup simulation uses a constant recharge forcing to achieve a steady state water table depth (*NOTE*: This run does not include land surface simulation).

**2. [Transient simulation with ParFlow-CLM on CONUS1](https://hydroframesubsettools.readthedocs.io/en/latest/example_notebooks/conus1_subsetting_transient.html)**: This notebook walks through an example of subsetting a HUC8 from the CONUS1 domain. This example will subset everything needed to do a transient run with ParFlow-CLM and create a run script. This simulation is designed to start from a steady state groundwater condition.  Here we subset the initial conditions from the pre-generated national steady state water table depth simulations, but the starting point could also be generated from the spin up workbook.

**3. [Transient simulation with ParFlow-CLM on CONUS2](https://hydroframesubsettools.readthedocs.io/en/latest/example_notebooks/conus2_subsetting_transient.html)**: This notebook walks through an example of subsetting a HUC8 from the CONUS2 domain. This example will subset everything needed to do a transient run with ParFlow-CLM and create a run script. This simulation is designed to start from a steady state groundwater condition.  Here we subset the initial conditions from the pre-generated national steady state water table depth simulations, but the starting point could also be generated from the spin up workbook.
