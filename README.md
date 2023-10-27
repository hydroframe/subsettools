# subsettools

## Usage

The subsettools package provides functions that simplify the process of setting up a modeling domain in the contiguous United States (CONUS) and running the hydrologic model ParFlow. The package allows the user to subset all required hydrogeologic and climate forcing datasets from the database Hydrodata. It also contains various functions to edit and run your subset model. These capabilities allow for more rapid and replicable application of ParFlow for hydrologic simulations. 

## Installation

Subsettools can be installed in a python virtual environment with pip: 

```bash
$ pip install git+https://github.com/hydroframe/subsettools.git
```

In addition, we provide a reproducible computational environment using [MyBinder](https://mybinder.org/v2/gh/hydroframe/subsettools-binder/HEAD), where you can execute the example notebooks without the need to install the subsettools package or ParFlow. Please note that the MyBinder project has limited resources, so this is not an appropriate place to run large simulations - but it's a great environment to get started with the subsettools API.

If you prefer using Docker, you can get an image with JupyterLab, subsettools and ParFlow installed from DockerHub with:

```bash
$ docker pull george135/parflow:latest
```

## Getting started

Detailed documentation can be found at [Read the Docs](https://hydroframesubsettools.readthedocs.io/en/latest/).

You will need to register on the [hydrogen website](https://hydrogen.princeton.edu/pin) before you can use the hydrodata utilities.

### Example Notebooks

In this section you can find a collection of example jupyter notebooks using the subsettools API. 

The current list of example notebooks is given below. A more detailed explanation of each notebook and how it should be used can be found in the Example notebooks tab at the Read the Docs link above. 

*Note: Only CONUS1 examples are currently being supported, but CONUS2 examples will be added in the future.* 

1. conus1_subsetting.ipynb - Subsets a CONUS1 domain and runs a transient simulation with ParFlow-CLM.
2. conus1_subsetting_spinup.ipynb - Subsets a CONUS1 domain and performs a model initialization (spin up) with ParFlow.

### Template Model Runscripts

In addition to example notebooks, several reference .yaml files are provided with the subsettools package. You may use these as a template for a ParFlow run that most closely meets the specifications of the model you are trying to build. Additional details can be found in the package documentation (linked above) under the `Template runscripts` tab. 

A list of currently provided run templates is provided below:

1. conus1_pf_spinup_box.yaml
2. conus1_pf_spinup_solid.yaml
3. conus1_pfclm_transient_box.yaml
4. conus1_pfclm_transient_solid.yaml
5. conus2_pf_spinup_box.yaml
6. conus2_pf_spinup_solid.yaml
7. conus2_pfclm_transient_box.yaml
8. conus2_pfclm_transient_solid.yaml

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`subsettools` was created by George Artavanis, Amanda Triplett. It is licensed under the terms of the MIT license.

## Credits

`subsettools` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
