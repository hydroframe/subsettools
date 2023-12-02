# Getting started

## Installation

The best way to install `subsettools` is using pip. This installs our latest stable release with fully-supported features.

```bash
$ pip install subsettools
```

You can also install the latest development version by cloning the GitHub repository and using pip to install from the local directory:

```bash
$ pip install git+https://github.com/hydroframe/subsettools.git
```

In addition, we provide a reproducible computational environment using [MyBinder](https://mybinder.org/v2/gh/hydroframe/subsettools-binder/HEAD), where you can execute the example notebooks without the need to install the subsettools package or ParFlow. Please note that the MyBinder project has limited resources, so this is not an appropriate place to run large simulations - but it's a great environment to get started with the subsettools API.

If you prefer using Docker, you can get an image with JupyterLab, subsettools and ParFlow installed from DockerHub with:

```bash
$ docker pull george135/parflow:latest
$ docker run -dp 8888:8888 --platform=linux/amd64 george135/parflow:latest
```
You should now be able use the container if you open a browser at http://localhost:8888/lab.

## Creating a HydroGEN, HydroFrame, HydroData account and registering a PIN

Users must create a [HydroGEN](https://hydro-generation.org) or [HydroFrame](https://hydroframe.org) account and register their PIN before using the ``subsettools`` package.

First, please visit our [Signup page](https://hydrogen.princeton.edu/signup) to sign up for a HydroGEN/HydroFrame/HydroData account.  This single sign up process is free and allows for user tracking and security.

Second, please visit the [HydroGEN PIN Page](https://hydrogen.princeton.edu/pin) to log in and create a 4-digit PIN.

After creating your PIN, you must register that PIN on the machine that you intend
to use to access our platforms and APIs. You can run the following code one time to register your PIN.::  

```python
from hf_hydrodata.gridded import register_api_pin
register_api_pin("<your_email>", "<your_pin>")
```

Your PIN will expire after 2 days of non-use. If your PIN expires, you must return to
the [HydroGEN PIN Page](https://hydrogen.princeton.edu/pin) and create a new PIN. 
You only need to re-register this PIN with the `register_api_pin` method if the 
new 4-digit PIN is different from your previous 4-digit PIN (the PIN is allowed
to stay the same).

## Building a ParFlow run from a template runscript

In order to run ParFlow, you must have a model object in the form of a .pfidb or .yaml file. 
These files can be generated initially from a python script using pftools keys.
The template runscripts provided with the subsettools package are all in .yaml format and provide examples of common use configurations for ParFlow. You may use these as a template for a ParFlow run that most closely meets the specifications of the model 
you are trying to build. 

We recommend a user select a template runscript corresponding to following three guidelines:
 
1. Select a source data domain you wish to subset from:
   a. CONUS1
   b. CONUS2
   
2. Select an input file type:
   a. Solid file (if you are providing a HUC or HUC list)
   b. A box domain (you are providing lat/lon coordinates)
   
3. Select a template corresponding to your planned way to run ParFlow:
   a. Transient run, fully coupled to the climate model CLM
   b. Model Initialization running ParFlow on its own with a longterm recharge (PmE) mask.
   
We provide 8 template runscripts which correspond to all unique combinations of the above three guidelines.
These template runscripts have values set to the values used to run simulations of CONUS1 or CONUS2 in the past.

To get a template runscript provided with the package, you should use the `get_template_runscript` function (API reference [here](https://hydroframesubsettools.readthedocs.io/en/edit-docs/autoapi/subsettools/datasets/index.html#subsettools.datasets.get_template_runscript)). The function will get the correct template runscript based on your choices and write it to your chosen directory. For example, to get the template for a ParFlow-CLM coupled run on the CONUS1 grid with a solid input file, you can do:

```python
from subsettools.datasets import get_template_runscript
reference_run = get_template_runscript(grid="conus1", mode="transient", input_file_type="solid",
                                          write_dir="/path/to/your/chosen/directory"
                                      )
```

*If you want to use your own runscript:*

You can provide your own .pfidb or .yaml file to subsettools and are not required to use these templates to use the subsettools functions. 
However, we encourage starting from one of these templates before making other changes to your model unless you are an experienced ParFlow user as it is possible the settings in your runscript will be incompatible with running data from the CONUS1 and 2 domains. 

