---
title: 'SubsetTools: A Python package to subset data to build and run ParFlow hydrologic models'
tags:
  - Python
  - hydrology
  - modeling
  - simulation
authors:
  - name: Amanda K. Triplett
	orcid: 0009-0009-8085-3938
	equal-contrib: true
	affiliation: 2
  - name: Georgios Artavanis
	equal-contrib: true
	affiliation: 1
  - name: William M. Hasling
	affiliation: 1
  - name: Amy C. Defnet
	affiliation: 1
  - name: Amy Johnson
	affiliation: "2, 3"
  - name: Will E. Lytle
	affiliation: 2
  - name: Elena Leonarduzzi
	affiliation: 1
  - name: Andrew Bennett
	affiliation: 2
  - name: Laura E. Condon
	affiliation: 2
  - name: Reed M. Maxwell
	affiliation: 1
affiliations:
 - name: Princeton University, USA
   index: 1
 - name: Department of Hydrology and Atmospheric Sciences, University of Arizona, USA
   index: 2
 - name: CyVerse, USA
   index: 3

date: 31 January 2024
bibliography: paper.bib

---

# Summary
Hydrologic models are an integral part of understanding and managing water supply. There are countless hydrologic models available that differ in their complexity, scale and focus on different parts of the hydrologic cycle. ParFlow is a fully integrated, physics-based model that simulates surface and subsurface flow simultaneously.  ParFlow is also coupled with a land surface model which allows it to simulate the full terrestrial hydrologic cycle from the bedrock to the top of the treetops. ParFlow has been applied to a myriad of watersheds across the US and around the world to answer questions of water supply and groundwater–surface water interactions (citations). 

ParFlow is a scientifically rigorous hydrologic model, however its application by the broader community has been limited to a degree by its technical complexity which creates a high barrier to entry for new users. Intensive training and hydrologic expertise is required to appropriately build a ParFlow model from scratch.

`SubsetTools` is a Python package that seeks to lower the barrier to entry by allowing a user to subset published and verified ParFlow inputs and model configurations to build their own watershed models.  These tools allow a user to set up and run a model in a matter of minutes, rather than weeks or months. `SubsetTools` is designed to interface with the first and second generation ParFlow configurations which provide model inputs for the contiguous United States (CONUS). 


# Statement of need
There are several big barriers to building a hydrologic model from scratch: 
Assuring that model inputs have the correct format, units, spatial resolution, and orientation is a time consuming process, but this is just the beginning. It takes significant time and expertise to assemble and process all of the input datasets that the model will require. It has taken years of development to build a national geofabric for the ParFlow CONUS simulations. Our team also conducted large data assembly and analysis projects to develop hydrologically consistent topographic datasets and spatially consistent and continuous hydrostratigraphy. Rather than repeating this effort, `SubsetTools` users can start from all of the input datasets that have already been developed and tested for hydrologic consistency. 
It requires modeling expertise to set up a ParFlow run script. A run script often includes more than a hundred input keys and parameters that need to be set and tuned for the simulation to run smoothly.  We have multiple working model configurations already developed for our national platform and can easily adapt these scripts for watershed simulation. 
Groundwater models require a very long initialization known as ‘spinup’ to develop a steady state groundwater configuration. This has to be completed before any transient simulations are run. Because we have already developed steady state conditions at the national level users of subset tools can skip this step and go straight to running their model. 

The `SubsetTools` package provides functions that simplify the process of setting up and running a ParFlow model in CONUS. It allows the user to subset all required hydrogeologic and climate forcing datasets from the `HydroData` database (reference or link for hydrodata?). It also provides template model runscripts which are designed to link seamlessly with functions that edit the model keys corresponding to the domain and model configuration specified by the user. These features enable a more rapid and replicable application of ParFlow for hydrologic simulations.

`SubsetTools` is designed to be used by both hydrology students and researchers. For students, the functions and examples provided in the package can be run with little programming or hydrologic knowledge to start teaching concepts. However, the functions have been thoughtfully designed to be flexible and transparent so that more advanced users can develop customized models that meet their needs.

The source code for `SubsetTools` is available on [GitHub](https://github.com/hydroframe/subsettools). The documentation for the package is hosted on [ReadTheDocs](https://hydroframesubsettools.readthedocs.io/en/latest/) and includes installation instructions, short tutorials, example notebooks, a complete API reference as well as contributing guidelines.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

This research has been supported by the U.S. Department of Energy Office of Science (DE-AC02-05CH11231) and the US National Science Foundation Office of Advanced Cyberinfrastructure (OAC- 2054506 and OAC-1835855).

# References

