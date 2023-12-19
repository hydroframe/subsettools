# SubsetTools

## Introduction

The [SubsetTools](https://hydroframe.org/subsettools) package is developed and maintained by the [HydroFrame](https://hydroframe.org) project. It is designed to subset model inputs and outputs from the national ParFlow modeling framework. We provide tools to subset all required hydrogeologic and climate forcing datasets for a ParFlow simulation as well as obtaining pre-configured run scripts for your desired domain based on the most common use cases.  

Example workflows are provided for working with both the first ([PFCONUS1](https://hydroframe.org/parflow-conus1)) and second ([PFCONUS2](https://hydroframe.org/parflow-conus2)) generation of the national model. Refer to the HydroFrame website for more information on these domains. 

Subset tools has been configured to work with the [HydroData](https://hydroframe.org/hydrodata) data catalog which houses national ParFlow inputs and simulation results as well as a broad array of other hydrologic variables.  Refer to the [HydroData documentation](https://hf-hydrodata.readthedocs.io/en/latest/index.html) for information on how to access additional datasets. 

SubsetTools will provide you with ParFlow scripts that can run locally, but it should be noted that depending on the size of your domain you may want to deploy your runs on HPC resources. ParFlow is designed to run efficiently in parallel and all simulations can easily be distributed across multiple processors.

Checkout our Getting Started guide for installation instructions, information on creating a Hydrogen account, setting up a ParFlow run from templates and more! The HowTo section contains short examples, while the Example notebooks sections contains longer workflows that setup a ParFlow run at the end. The API reference contains the full list of available functions.

## Citing SubsetTools

If you use our tools please cite this package in your work. In addition please make sure to cite all of the datasets that you subset. Examples for obtaining the DOIs for any dataset you use are provided in the examples. 

## Contributing

Interested in contributing? Check out the contributing guidelines on the [SubsetTools GitHub repo](https://github.com/hydroframe/subsettools). Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License
`subsettools` is part of the [HydroFrame project](https://hydroframe.org/). It is licensed under the terms of the MIT license.
