# subsettools

new lline
second line

Subsetting tools and utilities for ParFlow

There are two options to use subsettools: use it as part of the parflow-ml module on verde, or install it in a virtual environment of your own.
If you want to use it within parflow-ml follow the 'Getting started' instructions, otherwise go to 'Installation'.

## Getting started

subsettools is already installed in the parflow-ml module on verde. You only need to get the example notebook from the repo and run it.

To get the example notebook, start a terminal session and clone the repository:

```bash
$ cd ~
$ git clone https://github.com/hydroframe/subsettools
$ cp subsettools/docs/conus1_subsetting.ipynb .
```

The notebook should now be in your home directory. To run the notebook, go to verde.princeton.edu. From the interactive apps tab, choose "Jupyter". \
"Number of hours": choose how many hours you want to run the notebook session. \
"Anaconda3 version used for starting up jupyter interface": choose "parflow-ml" from the list of environments. \
"Extra slurm options": --ntasks=16 (for the example notebook). In general, set ntasks to the number of processes you will need when you distribute the parflow run. \
Finally, launch the notebook. Feel free to run through the whole notebook as it is, or play with the functions and change their arguments!


## Installation

```bash
$ pip install git+ssh://git@github.com/hydroframe/subsettools
```

## Usage

- TODO

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`subsettools` was created by George Artavanis, Amanda Triplett. It is licensed under the terms of the MIT license.

## Credits

`subsettools` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
