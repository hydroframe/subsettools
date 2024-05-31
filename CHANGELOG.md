# Changelog

<!--next-version-placeholder-->

## 2.0.0 (30/05/2024)
- subset_static will raise a ValueError if a variable is requested that does not
  exist in the dataset instead of just printing an error message.
- subset_forcing has been refactored to improve cooperation between threads.

## 1.0.5 (17/05/2024)
- huc_to_ij, latlon_to_ij, upstream_area_to_ij and create_mask_solid have been deprecated.
  They have been replaced by define_huc_domain, define_latlon_domain, define_upstream_domain
  and write_mask_solid.

## 1.0.4 (03/05/2024)
- Added upstream_area_to_ij function to subset based on the upstream area of a point or collection of points in the grid

## 1.0.0 (07/12/2023)
- Update versions for pftools, numpy, etc
- Rename variables to keep consinsent with changes in the datacatalog
-  Add CONUS2 examples in the documentation

## 0.4.0 (29/11/2023)
- All functions that write files now return a dictionary where the keys are variable names and the values are paths to the files written.
- subset_press_init gets data at midnight on the date specified (before it was getting data an hour before midnight on the date specified)
- Short tutorials were added and documentation reorganized
- Improved error messages and argument-checking code was added

## 0.3.0 (31/10/2023)

- Created multi-threaded version of subset_forcing to optimize performance.
- Added optional forcing_vars argument to subset_forcing in case users want to request only some forcing variables.

## 0.2.5 (31/10/2023)

- First release of `subsettools` in PyPI!
- Added timezone handling functionality in subset_press_init, subset_forcing, config_clm

## 0.1.0 (22/06/2023)

- First release of `subsettools`!
