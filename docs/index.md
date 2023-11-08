# subsettools

## Introduction

The subsettools package provides functions that simplify the process of setting up a modeling domain in the contiguous United States (CONUS) and running the hydrologic model ParFlow. The package allows the user to subset all required hydrogeologic and climate forcing datasets from the database Hydrodata. It also contains various functions to edit and run your subset model. These capabilities allow for more rapid and replicable application of ParFlow for hydrologic simulations.

Checkout our Getting Started guide for installation instructions and information on the tutorials, example workflows and templates provided with the package!

## Hydrodata overview

>>> TODO: put information on Hydrodata and link to the docs: [Hydrodata documentation](https://maurice.princeton.edu/hydroframe/docs/index.html#). 

## CONUS1 domain overview

CONUS1 is a box domain over the contiguous United States. It does not cover areas near the coasts. This model has 4 soil layers (top 2 meters) and 1 geologic layer that is 100 m thick. The cell size is 1 kilometer. 

Please refer to our [grid documentation](https://maurice.princeton.edu/hydroframe/docs/gridded_data.html#dataset-type-parameters) for more details about the CONUS datasets contained on Hydrodata.

Please refer to and cite the following DOIs if you want more details or to release work with the datasets used in the example notebooks:
1. [10.5194/gmd-14-7223-2021](https://gmd.copernicus.org/articles/14/7223/2021/)
2. [10.1016/j.advwatres.2015.04.008](https://www.sciencedirect.com/science/article/pii/S0309170815000822)
3. [10.1002/2014WR016774](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014WR016774)
4. [10.1002/2015GL066916](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066916)

## CONUS2 domain overview

>>> TODO: put information on CONUS2 with references

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`subsettools` is part of the [HydroFrame project](https://hydroframe.org/). It is licensed under the terms of the MIT license.


```{toctree}
:maxdepth: 2
:hidden:

getting_started.md
examples/gallery.md
example_notebooks/example_notebooks.md
template_runscripts.md
autoapi/index
changelog.md
contributing.md
conduct.md
```
