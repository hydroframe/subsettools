# subsettools

## Usage

The subsettools package provides functions to simplify the process of setting up a modeling domain to run the hydrologic model ParFlow. The package allows the user to subset all required hydrogeologic and climate forcing datasets from Hydrodata. It also contains various functions to edit and run models. These capabilities allow for more rapid application of ParFlow to hydrologic simulation. 

## Installation

Subsettools can be installed in a python virtual environment with pip: 

```bash
$ pip install git+ssh://git@github.com/hydroframe/subsettools
```

Note for verde users: the subsettools package is part of the parflow-ml module, so you can access it via

```bash
$ module load parflow-ml
```

## Getting started

A collection of examples using the subsettools API is provided in the example jupyter notebooks.

To get the entire folder of example notebooks, start a terminal session, navigate to your chosen download location and do:

```bash
$ svn checkout https://github.com/hydroframe/subsettools/trunk/docs/example_notebooks
```

If you want to download a single notebook, copy the notebook title from the github URL and do:

```bash
$ curl https://raw.githubusercontent.com/hydroframe/subsettools/main/docs/example_notebooks/<notebook_title>.ipynb -o <notebook_title>.ipynb
```
The current list of example notebooks with brief explanations is given below. If you would like more information or to explore the notebooks within readthedocs please go to the following:

1. conus1_subsetting.ipynb - Subsets data from the CONUS2 modeling domain stored on Hydrodata as well as climate forcing to perform a simulation with ParFlow-CLM, that is ParFlow fully coupled to the climate model
. 
2. conus1_subsetting_spinup.ipynb - Subsets data from the CONUS1 modeling domain stored on Hydrodata to perform a spin-up (model initialization) to steady state with only ParFlow and a longterm recharge mask (PmE).
3. conus2_subsetting.ipynb - Subsets data from the CONUS2 modeling domain stored on Hydrodata as well as climate forcing to perform a simulation with ParFlow-CLM, that is ParFlow fully coupled to the climate model.
4. conus2_subsetting_spinup.ipynb - Subsets data from the CONUS2 modeling domain stored on Hydrodata to perform a spin-up (model initialization) to steady state with only ParFlow and a longterm recharge mask (PmE).

In addition to example notebooks, several reference .yaml files are provided at subsettools/src/subsettools/ref_yamls. You may use these as a template for a ParFlow (or ParFlow-CLM) run that most closely meets the specifications of the model you are trying to build. A list of currently provided run templates is provided below. The first word in the filename refers to the larger modeling domain to be subset from (either conus1 or conus2).The second part refers to pf_spinup, which is a model file set up for initializing a modeling domain with only parflow or pfclm_transient which is a model run with parflow fully linked to the climate model clm. The last specifies if the model file is for a box or solid file domain. You may always use your own .yaml or .pfidb that defines your model with the subsettools package.

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
