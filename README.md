# subsettools

Subsetting tools and utilities for ParFlow

## Installation

```bash
$ pip install git+ssh://git@github.com/hydroframe/subsettools
```

### Detailed instructions to install on verde and run the example notebook:

Create a virtual environment ENV_NAME:

```bash
$ module load anaconda3/2021.11
$ conda create -n ENV_NAME python=3.9
$ conda activate ENV_NAME
$ pip install notebook ipykernel matplotlib
$ pip install git+ssh://git@github.com/hydroframe/subsettools
$ conda deactivate
```

To get the example notebook, clone the repository:

```bash
$ cd ~
$ git clone https://github.com/hydroframe/subsettools .
$ cp subsettools/docs/conus1_subsetting.ipynb .
```

To run the notebook, go to the [verde on demand](verde.princeton.edu).
From the interactive apps tab, choose "Jupyter".
"Number of hours": choose how many hours you want to run the notebook session.
"Anaconda3 version used for starting up jupyter interface": choose the virtual environment you just created.
"Extra slurm options": --ntasks=16 (for the example notebook). In general, set ntasks to the number of processes you will need when you distribute the parflow run.
Finally, launch the notebook. Feel free to run through the whole notebook as it is, or play with the functions and change their arguments!


## Usage

- TODO

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`subsettools` was created by George Artavanis, Amanda Triplett. It is licensed under the terms of the MIT license.

## Credits

`subsettools` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
