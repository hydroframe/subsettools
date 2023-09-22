Template runscripts
===================

In order to run ParFlow, you must have a model object in the form of a .pfidb or .yaml file. 
These files can be generated initially from a python script using pftools keys or a .tcl using tcl keys. 
The template runscripts contained in this repo are all in .yaml format and provide examples of common use configurations for ParFlow. 

We recommend a user select a template runscript corresponding to following three guidelines:
 
1. Select a source data domain you wish to subset from:
   a. CONUS1
   b. CONUS2
   
2. Select a subset type:
   a. Solid file (if you are providing a HUC or HUC list)
   b. A box domain (you are providing lat/lon coordinates)
   
3. Select a template corresponding to your planned way to run ParFlow:
   a. Transient run, fully coupled to the climate model CLM
   b. Model Initialization running ParFlow on its own with a longterm recharge (PmE) mask.
   
We provide 8 template runscripts which correspond to all unique combinations of the above three guidelines:

1. Conus1_pf_spinup_box.yaml
2. Conus1_pf_spinup_solid.yaml
3. Conus1_pfclm_transient_box.yaml
4. Conus1_pfclm_transient_solid.yaml
5. Conus2_pf_spinup_box.yaml
6. Conus2_pf_spinup_solid.yaml
7. Conus2_pfclm_transient_box.yaml
8. Conus2_pfclm_transient_solid.yaml

These template runscripts have values set to the values used to run simulations of CONUS1 or CONUS2 in the past.

*If you want to use your own runscript:*

You can provide your own .pfidb or .yaml file to subsettools and are not required to use these templates to use the subsettools functions. 
However, we encourage starting from one of these templates before making other changes to your model unless you are an experienced ParFlow user as it is possible the settings in your runscript will be incompatible with running data from the CONUS1 and 2 domains. 
