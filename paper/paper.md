---
title: 'SubsetTools: A Python package to subset data to build and run ParFlow hydrologic models'
tags:
  - Python
  - hydrology
  - modeling
  - simulation
authors:
  - name: Amanda K. Triplett
    orcid: 0000-0000-0000-0000
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
    affiliation: "2, 3"
  - name: Reed M. Maxwell
    corresponding: true 
    affiliation: 1
affiliations:
 - name: Princeton University, USA
   index: 1
 - name: University of Arizona, USA
   index: 2
 - name: CyVerse, USA
   index: 3

date: 31 January 2024
bibliography: paper.bib

---

# Summary
Hydrologic models are an integral part of understanding and managing water supplies in the contiguous United States (CONUS). There are countless models available to use that differ in their complexity, scale and focus on different parts of the hydrologic cycle. ParFlow is a fully integrated, physics-based model that simulates surface and subsurface flow. ParFlow has been used in various publications to answer questions of water supply, quality, and use in natural and managed states as well as to examine climate change impacts on hydrologic systems (references). It is a robust tool, however its application by the broader community has been limited to a degree by the high barrier to entry. Intensive training and guidance is required to appropriately build a ParFlow model from scratch. `SubsetTools` is a Python package that seeks to lower the barrier to entry by allowing a user to subset published and verified ParFlow inputs and set up a script to run the model for their domain of interest. These tools allow a user to set up and run a model in a matter of minutes, rather than weeks or months.


# Statement of need
Assuring that model inputs have the correct format, units, spatial resolution, and orientation is a time consuming process. The majority of inputs to the ParFlow model are in formats that are not commonly used such as a ParFlow Binary or .pfb. This may cause difficulties in learning to manipulate the files. Further, the inputs needed vary depending on what configuration you want to run the model in. For instance, running base ParFlow or ParFlow linked to the climate model CLM. 

The ParFlow runscript must also be customized to accept these modified files and the specified model configuration. This is done with a model runscript which is defined by hundreds of keys or parameters. For a beginner or even experienced user, changing a key can introduce significant errors which can be difficult to trace. 

The `SubsetTools` package provides functions that simplify the process of setting up and running a ParFlow model in CONUS. It allows the user to subset all required hydrogeologic and climate forcing datasets from the `HydroData` database (reference or link for hydrodata?). It also provides template model runscripts which are designed to link seamlessly with functions that edit the model keys corresponding to the domain and model configuration specified by the user. These features enable a more rapid and replicable application of ParFlow for hydrologic simulations.

`SubsetTools` is designed to be used by both hydrology students and researchers. For students, the functions and examples provided in the package can be run with little programming or hydrologic knowledge to start teaching concepts. However, the functions have been thoughtfully designed to be flexible and transparent so that more advanced users can develop customized models that meet their needs. `SubsetTools` has already been used in a graduate hydrology course at Princeton, as well as at a DOE IDEAS-watersheds workshop at Stanford. 

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

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

This project was funded by the US National Science Foundation Convergence Accelerator
Program, Grant No. CA-2040542 as well as US National Science Foundation Grant No. 1835855.

# References
