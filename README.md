# subsettools

Subsetting tools and utilities for ParFlow

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

## Usage

- TODO

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`subsettools` was created by George Artavanis, Amanda Triplett. It is licensed under the terms of the MIT license.

## Credits

`subsettools` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
