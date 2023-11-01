# Getting started

## Installation

Subsettools can be installed in a python virtual environment with pip. You can get the latest stable release from PyPI with:

```bash
$ pip install subsettools
```

Or you can get the latest development version directly from GitHub with:

```bash
$ pip install git+https://github.com/hydroframe/subsettools.git
```

In addition, we provide a reproducible computational environment using [MyBinder](https://mybinder.org/v2/gh/hydroframe/subsettools-binder/HEAD), where you can execute the example notebooks without the need to install the subsettools package or ParFlow. Please note that the MyBinder project has limited resources, so this is not an appropriate place to run large simulations - but it's a great environment to get started with the subsettools API.

If you prefer using Docker, you can get an image with JupyterLab, subsettools and ParFlow installed from DockerHub with:

```bash
$ docker pull george135/parflow:latest
```

## Creating a HydroGEN API Account

Users must create a HydroGEN API account and register their PIN before using the 
``subsettools`` package.

First, please visit our [HydroGEN PIN Page](https://hydrogen.princeton.edu/pin) to 
sign up for an account and create a 4-digit PIN.

After creating your PIN, you must register that PIN on the machine that you intend
to use. You can run the following code one time to register your PIN.::  

```python
from hf_hydrodata.gridded import register_api_pin
register_api_pin("<your_email>", "<your_pin>")
```

Your PIN will expire after 2 days of non-use. If your PIN expires, you must return to
the [HydroGEN PIN Page](https://hydrogen.princeton.edu/pin) and create a new PIN. 
You only need to re-register this PIN with the `register_api_pin` method if the 
new 4-digit PIN is different from your previous 4-digit PIN (the PIN is allowed
to stay the same).


## Example Notebooks

In this section you can find a collection of example workflows using the subsettools API. 

The current list of example notebooks is given below. A more detailed explanation of each notebook and how it should be used can be found in the Example notebooks tab. 

*Note: Only CONUS1 examples are currently being supported, but CONUS2 examples will be added in the future.* 

1. conus1_subsetting.ipynb - Subsets a CONUS1 domain and runs a transient simulation with ParFlow-CLM.
2. conus1_subsetting_spinup.ipynb - Subsets a CONUS1 domain and performs a model initialization (spin up) with ParFlow.

## Template Model Runscripts

In addition to example notebooks, several reference .yaml files are provided with the subsettools package. 
You may use these as a template for a ParFlow run that most closely meets the specifications of the model 
you are trying to build. Additional details can be found under the Template runscripts tab. 
