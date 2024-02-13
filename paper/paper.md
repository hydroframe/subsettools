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
    corresponding: true
    affiliation: 3
  - name: Georgios Artavanis
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2"
  - name: William M. Hasling
    affiliation: "1, 2"
  - name: Reed M. Maxwell
    affiliation: "2, 5, 6"
  - name: Amy C. Defnet
    affiliation: "1, 2"
  - name: Amy Johnson
    affiliation: "3, 4"
  - name: Will E. Lytle
    affiliation: 3
  - name: Andrew Bennett
    affiliation: 3
  - name: Elena Leonarduzzi
    affiliation: 2
  - name: Laura E. Condon
    affiliation: 3

affiliations:
 - name: Research Software Engineering, Princeton University, USA
   index: 1
 - name: Integrated GroundWater Modeling Center, Princeton University, USA
   index: 2
 - name: Department of Hydrology and Atmospheric Sciences, University of Arizona, USA
   index: 3
 - name: CyVerse, USA
   index: 4
 - name: Department of Civil and Environmental Engineering, Princeton University, USA
   index: 5
 - name: High Meadows Environmental Institute, Princeton University, USA
   index: 6


date: 31 January 2024
bibliography: paper.bib

---

# Summary
Hydrologic models are an integral part of understanding and managing water supply. There are countless hydrologic models available that differ in their complexity, scale and focus on different parts of the hydrologic cycle. ParFlow is a fully integrated, physics-based model that simulates surface and subsurface flow simultaneously [@RN351; @RN316; @RN255; @RN320]. ParFlow is also coupled with a land surface model which allows it to simulate the full terrestrial hydrologic cycle from the bedrock to the top of the treetops [@RN322; @RN321]. ParFlow has been applied to a myriad of watersheds across the US and around the world to answer questions of water supply and groundwater–surface water interactions.

ParFlow is a scientifically rigorous hydrologic model, however its application by the broader community has been limited to a degree by its technical complexity which creates a high barrier to entry for new users. Intensive training and hydrologic expertise is required to appropriately build a ParFlow model from scratch.

`SubsetTools` is a Python package that seeks to lower the barrier to entry by allowing a user to subset published and verified ParFlow inputs and model configurations to build their own watershed models. These tools allow a user to set up and run a model in a matter of minutes, rather than weeks or months. `SubsetTools` is designed to interface with the [first](https://hydroframe.org/parflow-conus1) [@RN257; @RN257; @RN353]and [second](https://hydroframe.org/parflow-conus2) [@RN352] generation ParFlow configurations which provide model inputs for the contiguous United States (CONUS). 


# Statement of need
There are several big barriers to building a hydrologic model from scratch. `SubsetTools` helps to resolve the issues caused by three primary ones:  

Finding quality data and then using it within a model. It requires significant time and expertise to assemble and process all of the input datasets that the model will require. It has taken years of development to build a national geofabric for the ParFlow CONUS simulations [@RN257; @RN352]. Our team also conducted large data assembly and analysis projects to develop hydrologically consistent topographic datasets [@RN354] and spatially consistent and continuous hydrostratigraphy [@RN355]. Rather than repeating this effort, `SubsetTools` users can start from all of the input datasets that have already been developed and tested for hydrologic consistency. This assures that model inputs have the correct format, units, spatial resolution, and orientation to run a new subset model.
It requires modeling expertise to set up a ParFlow run script. A run script often includes more than a hundred input keys and parameters that need to be set and tuned for the simulation to run smoothly.  We have multiple working model configurations already developed for our national platform and can easily adapt these scripts for watershed simulation. 
Groundwater models require a very long initialization known as ‘spinup’ to develop a steady state groundwater configuration. This has to be completed before any transient simulations are run. Because we have already developed steady state conditions at the national level, users of `SubsetTools` can skip this step and go straight to running their model. 

The `SubsetTools` package provides functions that simplify the process of setting up and running a ParFlow model in CONUS. It allows the user to subset all required hydrogeologic and climate forcing datasets from [`HydroData`](https://hf-hydrodata.readthedocs.io/en/latest/). It also provides template model runscripts which are designed to link seamlessly with functions that edit the model keys corresponding to the domain and model configuration specified by the user. These features enable a more rapid and replicable application of ParFlow for hydrologic simulations.

`SubsetTools` is designed to be used by both hydrology students and researchers. For students, the functions and examples provided in the package can be run with little programming or hydrologic knowledge to start teaching concepts. However, the functions have been thoughtfully designed to be flexible and transparent so that more advanced users can develop customized models that meet their needs.

# Functionality 

The source code for `SubsetTools` is available on [GitHub](https://github.com/hydroframe/subsettools). The documentation for the package is hosted on [ReadTheDocs](https://hydroframesubsettools.readthedocs.io/en/latest/) and includes installation instructions, short tutorials, example notebooks, a complete API reference, as well as contributing guidelines. In the section below we will go through an abbreviated outline of how a user can interact with `SubsetTools` to build and run their own ParFlow model. 

First, the user should consider what domain they want to use for source data, CONUS1 or CONUS2. We then require information regarding their geographic area of interest, which may be given as a hydrologic unit code (HUC) or list of HUCs, a bounding box or a single point of latitude and longitude. The user should also specify timing information such as a date range for their simulation. 

With the information provided above, users can subset all required model input files. An example is shown below for subsetting static and climate forcing data on the CONUS2 grid for the Upper Verde region.

```python
import subsettools as st

static_paths = st.subset_static(
    ij_bounds=[815, 1200, 924, 1290], # CONUS2 grid bounds for the Upper Verde region
    dataset=”conus2_domain”,
    write_dir=”/path/to/your/output/directory”,
)

forcing_paths = st.subset_forcing(
    ij_bounds=[815, 1200, 924, 1290],
    grid=”conus2”,
    start=”2012-10-11”,
    end=”2013-10–11”,
    dataset=”CW3E”,
    write_dir=”/path/to/your/output/directory”,
    forcing_vars=(‘precipitation’, ‘air_temp’,),
)
```

![Two example subset outputs for HUC 15060202, the Upper Verde Watershed in Arizona. (a) shows the subset Gauckler-Manning coefficient n for this domain as a result of the function st.subset_static(). (b) is the air temperature forcing, one of the variables output by the function st.subset_forcing().](fig1.png)


An appropriate run script must also be selected based on the kind of ParFlow simulation to be performed. The `SubsetTools` package provides eight different templates, which can be used as a starting point for building a ParFlow model. The example function call shown below specifies a transient run using ParFlow-CLM over a solid file domain on the CONUS2 grid. 

```python
import subsettools as st

reference_run = st.get_template_runscript(
    grid=”conus2”,
    mode=”transient”,
    input_file_type=”solid”,
    write_dir=”/path/to/your/output/directory”
)
runscript_path = st.edit_runscript_for_subset(
    ij_bounds=[815, 1200, 924, 1290], # CONUS2 grid bounds for the Upper Verde region
    runscript_path=reference_run,
    runname=”your_runname”,
    forcing_dir=”/path/to/your/forcing/directory”,
)
```

The `SubsetTools` package also provides functions to customize a template runscript, for example by specifying the desired subset domain to match the subset inputs, modifying the file paths of the model input files, and changing the processor topology for the ParFlow run. Once the customized Parflow runscript is ready, the user can launch a ParFlow simulation using the [pftools](https://pypi.org/project/pftools/) package utilities. 

# Acknowledgements

This research has been supported by the U.S. Department of Energy Office of Science (DE-AC02-05CH11231) and the US National Science Foundation Office of Advanced Cyberinfrastructure (OAC- 2054506 and OAC-1835855).

# References

