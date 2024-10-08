# Getting started

## Installation

`subsettools` can be installed using pip. This installs our latest stable release with fully-supported features. `subsettools` currently supports Python versions 3.10, 3.11 and 3.12.

```bash
$ pip install subsettools
```

You can also install the latest development version by cloning the GitHub repository and using pip to install from the local directory:

```bash
$ pip install git+https://github.com/hydroframe/subsettools.git
```

If you encounter issues with installation, check out our FAQ page. If you still can't find an answer to your problem, please open an issue at the [subsettools GitHUB page](https://github.com/hydroframe/subsettools).

In order to use a small subset of the `subsettools` functions, you also need to have ParFlow installed. Please refer to the [HydroFrame](https://hydroframe.org/) for help getting started with ParFlow.

In addition, we provide a reproducible computational environment using [Binder](https://mybinder.org/v2/gh/hydroframe/subsettools-binder/HEAD), where you can execute the example notebooks without the need to install the subsettools package or ParFlow. Please note that the Binder project has limited resources, so this is not an appropriate place to run large simulations - but it's a great environment to get started with the subsettools API. Note that Binder may take several minutes to launch.

### Using Docker

If you prefer using Docker, you can get an image with JupyterLab, subsettools and ParFlow installed from DockerHub. Make sure you have an up-to-date version of Docker. Follow the instructions to set up your ParFlow Docker environment:

1. Go to [Docker](https://www.docker.com/products/docker-desktop/) and download Docker Desktop. Make sure that you download the correct version based on your computer’s operating system and architecture. Follow the instructions to install Docker Desktop. ![alt text](Docker1.png)
2. Launch Docker Desktop. The application should be running while you do the next steps. (You might need to create a Docker account and sign in if you don’t already have one.)
3. Start a new terminal session (outside Docker) to type the commands in the next two steps:
4. Pull the parflow/subsettools image from DockerHub. Choose the correct version based on your computer’s architecture.
**NOTE:** *The Docker image might take a while to download. Also, if your Docker is running out of space, you might need to use [docker system prune](https://docs.docker.com/engine/reference/commandline/system_prune/) with the appropriate options to clear out old containers and make space for the new one.*

- For the x86_64/amd64 (Intel Chip for Mac or Windows) architecture:
```bash
docker pull george135/subsettools_amd64
```
- For the arm64 (M1/M2 chip for Mac) architecture:
```bash
docker pull george135/subsettools_arm64
```
5. Once the image has finished downloading, you can run the container with:
- For the x86_64/amd64 architecture:
```bash
docker run -dp 8888:8888 george135/subsettools_amd64:latest start-notebook.sh --NotebookApp.token=''
```
- For the arm64 architecture:
```bash
docker run -dp 8888:8888 george135/subsettools_arm64:latest start-notebook.sh --NotebookApp.token=''
```
6. Use a browser to navigate to your [JupyterLab container](http://localhost:8888/lab?) or use the link that will appear next to your container on the Docker Desktop application: ![alt text](Docker2.png)
7. You should see a JupyterLab environment like this: ![alt text](Docker3.png)
8. Click on the Terminal application to start a terminal session *inside* the container.
9. Clone the GitHub repository that contains the example workflows for `subsettools` and `hf_hydrodata`:
```bash
git clone https://github.com/hydroframe/subsettools-binder.git
```
10. Navigate to `subsettools-binder` -> `subsettools` folder and click on the `definte_subset.ipynb` notebook. ![alt text](Docker4.png)
11. Make sure the notebook runs successfully without errors. **You will need to provide your Hydrogen email and pin in the first code cell of the notebook.**
12. Congratulations, your ParFlow Docker environment is ready to use!

## Creating a HydroGEN, HydroFrame, HydroData account and registering a PIN

Users must create a HydroGEN API account and register their PIN before using the `hf_hydrodata` package.

1. If its your first time using this package you will need to sign up for a HydroFrame account on the [HydroGEN Signup Page](https://hydrogen.princeton.edu/signup) (*Note: This only needs to be done once*)

2. Visit our [HydroGEN PIN Page](https://hydrogen.princeton.edu/pin) to create a 4-digit PIN.

3. After creating your PIN, you must register that PIN on the machine that you intend to use. You can run the following code one time to register your PIN::  


```python
from hf_hydrodata.gridded import register_api_pin
register_api_pin("<your_email>", "<your_pin>")
```

Your PIN will expire after 2 days of non-use. If your PIN expires, you must return to
the [HydroGEN PIN Page](https://hydrogen.princeton.edu/pin) and create a new PIN. 
You only need to re-register this PIN with the `register_api_pin` method if the 
new 4-digit PIN is different from your previous 4-digit PIN (the PIN is allowed
to stay the same).

## Things you can do with SubsetTools

SubsetTools is a package that leverages the [HydroData](https://www.hydroframe.org/hydrodata) platform and API to automate a number of data tasks.  These currently include subsetting (clipping) large data products in space or in time.  For example, static products like a Digital Elevation Model (DEM) or time-varying products like meterological forcing data such as precipitation or temperature can be clipped and formatted using the example workflows.  SubsetTools also rapidly accelerates model development and can clip and formulate all the inputs to run a hydrologic model, for example over a HUC in the US for a water year.

## Building a ParFlow run from a template runscript

One of the primary workflows for SubsetTools is to automate building an integrated hydrologic simulation over the Continental US (CONUS).  In addition to static model inputs, the HydroData platform also has four decades of hourly forcing at 1km resolution.  Currently, we support the integrated hydrology model [ParFlow](https://www.parflow.org) coupled to CLM, [PFCLM](https://www.hydroframe.org/parflow-resources). 

Using SubsetTools you can automatically generate ParFlow inputs (e.g. a `.pfidb` or `.yaml` file) from a Python script using pftools keys.  The template runscripts provided with the SubsetTools package are all in .yaml format and provide examples of common use configurations for ParFlow. You may use these as a template for a ParFlow run that most closely meets the specifications of the model 
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

To get a template runscript provided with the package, you should use the `get_template_runscript` function (API reference [here](https://hydroframesubsettools.readthedocs.io/en/latest/autoapi/subsettools/index.html#subsettools.get_template_runscript)). The function will get the correct template runscript based on your choices and write it to your chosen directory. For example, to get the template for a ParFlow-CLM coupled run on the CONUS1 grid with a solid input file, you can do:

```python
from subsettools.datasets import get_template_runscript
reference_run = get_template_runscript(grid="conus1", mode="transient", input_file_type="solid",
                                          write_dir="/path/to/your/chosen/directory"
                                      )
```

*If you want to use your own runscript:*

You can provide your own .pfidb or .yaml file to SubsetTools and are not required to use these templates to use the subsettools functions. 
However, we encourage starting from one of these templates before making other changes to your model unless you are an experienced ParFlow user as it is possible the settings in your runscript will be incompatible with running data from the CONUS1 and 2 domains. 

