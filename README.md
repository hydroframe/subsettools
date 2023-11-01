# subsettools

## Overview

The subsettools package provides functions that simplify the process of setting up a modeling domain in the contiguous United States (CONUS) and running the hydrologic model ParFlow. The package allows the user to subset all required hydrogeologic and climate forcing datasets from the database Hydrodata. It also contains various functions to edit and run your subset model. These capabilities allow for more rapid and replicable application of ParFlow for hydrologic simulations. 

Detailed documentation can be found at [Read the Docs](https://hydroframesubsettools.readthedocs.io/en/latest/).

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

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`subsettools` was created by George Artavanis, Amanda Triplett. It is licensed under the terms of the MIT license.

## Credits

`subsettools` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
