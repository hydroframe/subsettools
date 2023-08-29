Template runscripts
===================

In order to run ParFlow, you must have a model object in the form of a .pfidb or .yaml file. 
These files can be generated initially from a python script using pftools keys or a .tcl using tcl keys. 
The template runscripts contained in this repo are provided to give examples of common use configurations for both the CONUS1 and 2 domains. 

We recommend a user select a template runscript corresponding to following three guidelines:
 
1. Select a template corresponding to the domain you wish to subset from (CONUS1 or CONUS2) 
2. Select a template corresponding to your planned way to run ParFlow, either a transient run fully coupled to the climate model, or spin up running ParFlow on its own with a longterm recharge (PmE) mask. 
3. Select a template which uses a solid file (if you are providing a HUC to subsettools for your subset) or a box domain (you are providing lat/lon coordinates)

We provide 8 template runscripts which correspond to all unique combinations of the above three guidelines:

1. Conus1_pf_spinup_box.yaml
2. Conus1_pf_spinup_solid.yaml
3. Conus1_pfclm_transient_box.yaml
4. Conus1_pfclm_transient_solid.yaml
5. Conus2_pf_spinup_box.yaml
6. Conus2_pf_spinup_solid.yaml
7. Conus2_pfclm_transient_box.yaml
8. Conus2_pfclm_transient_solid.yaml

These template runscripts have values set to the values used to run simulations of CONUS1 and CONUS2 in the past.

You can provide your own .pfidb or .yaml file to subsettools and are not required to use these templates to use the subsettools functions. 
However, we encourage starting from one of these templates before making other changes to your model unless you are an experienced ParFlow user. 
