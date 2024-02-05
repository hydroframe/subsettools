"""This module contains subsetting functions for Parflow runs.
"""

import os
import shutil
import pathlib
import subprocess
from datetime import datetime, timedelta
import threading 
import numpy as np
import hf_hydrodata
import re
from parflow import Run
from parflow.tools.io import read_pfb, write_pfb
from .subset_utils import (
    write_land_cover,
    edit_drvclmin,
    get_UTC_time,
)
from ._error_checking import (
    _validate_huc_list,
    _validate_grid,
    _validate_dir,
    _validate_grid_bounds,
    _validate_date,
)


def huc_to_ij(huc_list, grid):
    """Get the grid ij bounds of a bounding box that encompasses HUC IDs provided in huc_list.

       All HUC IDs in huc_list must be the same length (HUCs of the same level). All HUCs 
       should be adjacent. If a HUC is only partially covered by the provided grid, the grid 
       bounds for the covered area will be returned.

    Args:
        huc_list (list[str]): a list of USGS HUC IDs
        grid (str): The spatial grid that the ij indices are calculated relative to and that the subset data will be returned on. Possible values: “conus1” or “conus2”

    Returns:
        tuple[int]: A tuple of the form (imin, jmin, imax, jmax) representing the bounds in the conus grid of the area defined by the huc IDs in huc_list. imin, jmin, imax, jmax are the west, south, east and north sides of the box respectively and all i,j indices are calculated relative to the lower southwest corner of the domain.

    Raises:
        ValueError: If all HUC IDs are not the same length or if the area defined by the provided HUCs 
            is not part of the given grid.

    Example:

    .. code-block:: python

        grid_bounds = huc_to_ij(huc_list=["14080201", "14080202", "14080203"], grid="conus1")
    """
    _validate_huc_list(huc_list)
    _validate_grid(grid)
    conus_hucs, _, indices_j, indices_i = _get_conus_hucs_indices(huc_list, grid.lower())
    if indices_i.size == 0 or indices_j.size == 0:
        raise ValueError(f"The area defined by the provided HUCs is not part of the {grid} grid.")  
    return _indices_to_ij(conus_hucs, indices_j, indices_i)


def _get_conus_hucs_indices(huc_list, grid):
    """Get the huc datafile as an ndarray and three mask arrays representing the selected hucs.
    
    Args:
        huc_list (list[str]): a list of huc IDs 
        grid (str): The spatial grid that the ij indices are calculated relative to and that the subset 
            data will be returned on. Possible values: “conus1” or “conus2”

    Returns:   
        A tuple (conus_hucs, sel_hucs, indices_j, indices_i) where                                               
        conus_hucs is an ndarray of the huc datafile, sel_hucs is
        a mask array for the selected hucs, and indices_i and
        indices_j mask arrays in the j and i directions.
    """
    huc_len = len(huc_list[0])
    huc_list = [int(huc) for huc in huc_list]
    entry = hf_hydrodata.get_catalog_entry(
        dataset="huc_mapping", grid=grid, file_type="tiff"
    )
    if entry is None:
        raise ValueError(f"There is no HUC mapping entry for grid {grid}.")
    conus_hucs = hf_hydrodata.gridded.get_ndarray(entry, level=str(huc_len))
    sel_hucs = np.isin(conus_hucs, huc_list).squeeze()
    indices_j, indices_i = np.where(sel_hucs > 0)
    return conus_hucs, sel_hucs, indices_j, indices_i


def _indices_to_ij(conus_hucs, indices_j, indices_i):
    """Get the conus ij-bounds for the conus_hucs boundary defined by indices_j and indices_i.                                                  

    Args:                                                                                                                           
        conus_hucs (numpy.ndarray): conus huc data                                                                                             
        indices_j (numpy.ndarray): mask in the j direction for selected hucs                                                                    
        indices_i (numpy.ndarray): mask in the i direction for selected hucs                                                                   
 
    Returns:                                                                                                         
        A tuple of the form (imin, jmin, imax, jmax) representing the bounds                                                                    
        in conus_hucs defined by the two mask arrays indices_j and indices_i.                                                                   
    """
    imin = np.min(indices_i)
    imax = np.max(indices_i) + 1  # right bound inclusive
    jmin = np.min(indices_j)
    jmax = np.max(indices_j) + 1  # right bound inclusive                                                                                   
    return (int(imin), int(jmin), int(imax), int(jmax))


def latlon_to_ij(latlon_bounds, grid):
    """Get the ij bounds of the area defined by the the latitute/longitude bounds (latlon_bounds) relative to a selected conus grid.

    Args:
        latlon_bounds (List[List[float]]): list of the form [[lat1, lon1], [lat2, lon2]].
            [lat1, lon1] and [lat2, lon2] define the northwest and southeast corners of the desired box respectively.
        grid (str):  The spatial grid that the ij indices are calculated relative to and that 
            the subset data will be returned on. Possible values: “conus1” or “conus2”

    Returns:
        tuple[int]: A tuple of the form (imin, jmin, imax, jmax) representing the bounds in the conus grid of the area defined by latlon_bounds.  imin, jmin, imax, jmax are the west, south, east and north sides of the box respectively and all i,j indices are calculated relative to the lower southwest corner of the domain.

    Example:

    .. code-block:: python

        grid_bounds = latlon_to_ij(latlon_bounds=[[37.91, -91.43], [37.34, -90.63]], grid="conus2")
    """
    _validate_grid(grid)    
    if not isinstance(latlon_bounds, list):
        raise TypeError("latlon_bounds must be a list")
    if len(latlon_bounds) != 2:
        raise ValueError("latlon_bounds must contain exactly two lat-lon points: [[lat1, lon1], [lat2, lon2]]")
    for point in latlon_bounds:
        if not (isinstance(point, list) and
                len(point) == 2 and
                all(isinstance(value, (int, float)) for value in point)):
            raise ValueError("latlon_bounds must contain exactly two lat-lon points: [[lat1, lon1], [lat2, lon2]]")
    grid = grid.lower()
    point0 = hf_hydrodata.to_ij(grid, latlon_bounds[0][0], latlon_bounds[0][1])
    point1 = hf_hydrodata.to_ij(grid, latlon_bounds[1][0], latlon_bounds[1][1])
    imin, imax = [
        min(point0[0], point1[0]),
        max(point0[0], point1[0]),
    ]
    jmin, jmax = [
        min(point0[1], point1[1]),
        max(point0[1], point1[1]),
    ]
    return (imin, jmin, imax, jmax)


def create_mask_solid(huc_list, grid, write_dir):
    """Create mask and solid files for a ParFlow simulation and write them in a user specified directory. 

    Given a list of HUC IDs, this function will define a bounding box that encompasses all of the selected HUCs and will create:
        1. a 2D mask file that indicates which cells inside the box domain are part of the selected HUCS.
        2. a solid file that defines a 3D domain extending to the depth of whichever grid has been selected and tracing the boundaries of the selected HUCS. 

    Args:  
        huc_list (list[str]): a list of HUC IDs   
        grid (str): The spatial grid that the ij indices are calculated relative to and that 
            the subset data will be returned on. Possible values: “conus1” or “conus2”
        write_dir (str): directory path where the mask and solid files will be written  

    Returns:
        dict: A dictionary of paths with keys ("mask", "mask_vtk", "solid") and values filepaths to the created files.

    Raises:  
        FileNotFoundError: If write_dir is not a valid directory.  
        ValueError: If the area defined by the provided HUCs is not part of the given grid.

    Example:

    .. code-block:: python

        filepaths = create_mask_solid(
            huc_list=["14080201", "14080202", "14080203"],
            grid="conus2",
            write_dir="/path/to/your/chosen/directory"
        )
    """
    _validate_huc_list(huc_list)
    _validate_grid(grid)
    _validate_dir(write_dir)
    grid = grid.lower()
    conus_hucs, sel_hucs, indices_j, indices_i = _get_conus_hucs_indices(huc_list, grid)
    if indices_i.size == 0 or indices_j.size == 0:
        raise ValueError(f"The area defined by the provided HUCs is not part of the {grid} grid.")      
    arr_jmin = np.min(indices_j)
    arr_jmax = np.max(indices_j) + 1  # right bound inclusive
    ij_bounds = _indices_to_ij(conus_hucs, indices_j, indices_i)

    nj = ij_bounds[3] - ij_bounds[1]
    ni = ij_bounds[2] - ij_bounds[0]

    # checks conus1 / 2 grid and assigns appripriate dz and z_total for making the mask and solid file
    if grid  == "conus1":
        print("grid is conus1")
        layz = 100
        z_total = str(500)
    else:
        print("grid is conus2")
        layz = 200
        z_total = str(2000)
    # could look this up in dC and get the information

    # create and write the pfb mask
    mask_clip = np.zeros((1, nj, ni))
    mask_clip[0, :, :] = sel_hucs[
        arr_jmin:arr_jmax, ij_bounds[0] : ij_bounds[2]
    ]  # we need to use numpy iymin / iymax because we are subsetting the tif file
#    mask_clip = np.flip(
#        mask_clip, 1
#    )  # This flip tooks the section we just subset and puts it in the appropriate parflow orientation
    mask_clip = mask_clip.astype(float)
    mask_file_path = os.path.join(write_dir, "mask.pfb")
    write_pfb(mask_file_path, mask_clip, dx=1000, dy=1000, dz=layz, dist=False)
    print("Wrote mask.pfb")
    mask_vtk_path = os.path.join(write_dir, "mask_vtk.vtk")
    solid_file_path = os.path.join(write_dir, "solidfile.pfsol")

    try:
        parflow_dir = os.environ["PARFLOW_DIR"]
    except KeyError:
        raise KeyError('The environment variable PARFLOW_DIR has not been defined. Please make sure you have ParFlow installed ' \
                       'and os.environ["PARFLOW_DIR"] points to that installation.')
    file_path = os.path.join(parflow_dir, "bin", "pfmask-to-pfsol")
    if not os.path.exists(file_path):
        raise FileNotFoundError('pfmask-to-pfsol file not found. Please make sure you have ParFlow installed ' \
                       'and os.environ["PARFLOW_DIR"] points to that installation.')
    try:
        subprocess.run(
            [
                file_path,
                "--mask",
                mask_file_path,
                "--pfsol",
                solid_file_path,
                "--vtk",
                mask_vtk_path,
                "--z-bottom",
                "0.0",
                "--z-top",
                z_total,
            ],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError("pfmask-to-pfsol error:", e.stderr)

    print(f"Wrote solidfile.pfsol and mask_vtk.vtk with total z of {z_total} meters")
    file_paths = {"mask": mask_file_path, "mask_vtk": mask_vtk_path, "solid": solid_file_path}
    return file_paths 


def subset_static(
    ij_bounds,
    dataset,
    write_dir,
    var_list=(
        "slope_x",
        "slope_y",
        "pf_indicator",
        "mannings",
        "pf_flowbarrier",
        "pme",
        "ss_pressure_head",
    ),
):
    """Subset static input files from national datasets and write the subset values out as ParFlow binary files (pfbs) in a user specified directory.

    By default the following variables will be subset:
        - Slope in the east/west direction (slope_x)
        - Slope in the north/south direction (slope_y)
        - Subsurface units indicator file (pf_indicator)
        - Mannings roughness coefficients (mannings)
        - Depth to bedrock (pf_flowbarrier)
        - Long term average preciptation minus evaporation (i.e. recharge) (pme)
        - Steady state pressure head used to initialize transient simulations (ss_pressure_head)

    Note that some datasets might not contain all 7 static input variables. In that case, the subset_static function is going to 
    get and save the data for the variables supported by the dataset and print out a message for those that are not. 

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as i,j index values where 0,0 is the lower left hand 
            corner of a domain. ij_bounds are given relative to whatever grid is being used for the subset. Use the latlon_to_ij function to 
            determine ij indices from lat long values.  
        dataset (str): static inputs dataset name from the HydroData catalog e.g. "conus1_domain"
        write_dir (str): directory where the subset files will be written
        var_list (tuple[str]): tuple of variables to subset from the dataset. By default all 7 variables above will be subset. The user can specify as subset of these variables or list additional variables that are available in their dataset of choice.

    Returns:
        dict: A dictionary in which the keys are the static variable names and the values are file paths where the subset data were written. 
        
    Raises:
        FileNotFoundError: If write_dir is not a valid directory.

    Example:

    .. code-block:: python

        filepaths = subset_static(
            ij_bounds=(375, 239, 487, 329),
            dataset="conus1_domain", 
            write_dir="/path/to/your/chosen/directory",
            var_list=("slope_x", "slope_y")
        )
    """
    _validate_grid_bounds(ij_bounds)
    if not isinstance(dataset, str):
        raise TypeError("dataset name must be a string.")
    _validate_dir(write_dir)
    if not all(isinstance(var, str) for var in var_list):
        raise TypeError("All variable names should be strings.")
    file_paths = {}
    for var in var_list:
        entry = hf_hydrodata.get_catalog_entry(
            dataset=dataset, file_type="pfb", period="static", variable=var
        )
        if entry is not None:
            subset_data = hf_hydrodata.gridded.get_ndarray(entry, grid_bounds=ij_bounds)
            file_path = os.path.join(write_dir, f"{var}.pfb")
            write_pfb(file_path, subset_data, dist=False)
            file_paths[var] = file_path
            print(f"Wrote {var}.pfb in specified directory.")
        else:
            print(f"{var} not found in dataset {dataset}")
    return file_paths


def subset_press_init(ij_bounds, dataset, date, write_dir, time_zone="UTC"):
    """Subset a pressure file and write it as a pfb file to a user specified directory to be used as an initial pressure for a ParFlow simulation.

    This function will select the pressure file for midnight on the date provided and subset the selected pressure file to the ij_bounds provided. The subset data will be written out as a ParFlow binary file (pfb) to be used as an initial pressure file for a ParFlow simulation.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as i,j index values where 0,0 is 
            the lower left hand corner of a domain. ij_bounds are given relative to whatever grid is being used for the 
            subset. Use the latlon_to_ij function to determine ij indices from lat long values.  
        dataset (str): dataset name from the HydroData catalog that the pressure file will be subset from e.g. "conus1_baseline_mod"
        date (str): The date of the pressure file that you would like to subset, in the form 'yyyy-mm-dd'
        write_dir (str): directory where the subset file will be written
        time_zone (str): timezone information for subset date. Data will be subset at midnight in the specified timezone.
            Defaults to "UTC".

    Returns:
        str: The filepath of the subset file, which includes datetime information, so that it can be used by later functions (e.g. edit_runscript_for_subset).

    Raises:
        FileNotFoundError: If write_dir is not a valid directory.

    Example:

    .. code-block:: python
    
        filepath = subset_press_init(
            ij_bounds=(375, 239, 487, 329),
            dataset="conus1_baseline_mod",
            date="2005-12-15", 
            write_dir="/path/to/your/chosen/directory",
            time_zone="EST"
        )
    """
    _validate_grid_bounds(ij_bounds)
    if not isinstance(dataset, str):
        raise TypeError("dataset name must be a string.")
    _validate_date(date)
    _validate_dir(write_dir)
    if not isinstance(time_zone, str):
        raise TypeError("time_zone must be a string.")
    
    entry = hf_hydrodata.get_catalog_entry(
        dataset=dataset, file_type="pfb", variable="pressure_head", period="hourly"
    )
    if entry is None:
        print(f"No pressure file found for in dataset {dataset}")
        return None
    
    new_date = get_UTC_time(date, time_zone)
    print(f"UTC Date: {new_date}")
        
    date_string = new_date.strftime("%Y.%m.%d:%H.%M.%S_UTC0")

    subset_data = hf_hydrodata.gridded.get_ndarray(
        entry, grid_bounds=ij_bounds, start_time=new_date
    )

    file_path = os.path.join(write_dir, f"{dataset}_{date_string}_press.pfb")
    write_pfb(file_path, subset_data[0, :, :, :], dist=False)
    print(f"Wrote {file_path} in specified directory.")
    return file_path


def config_clm(ij_bounds, start, end, dataset, write_dir, time_zone="UTC"):
    """Modify template clm driver files for a desired subdomain and run duration and write them to a user specified directory.

    This function will obtain template clm driver files (specifically vegm, vep and drv_clmin) from the existing national simulations on HydroData and modify them to reflect the desired subdomain (indicated by the ij_bounds) and run duration (indicated by the start and end dates). The modified files will be written out to a user specified directory. These files are required if you are going to run a ParFlow-CLM simulation. 

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as i,j index values where 0,0 is 
            the lower left hand corner of a domain. ij_bounds are given relative to whatever grid is being used for the 
            subset. Use the latlon_to_ij function to determine ij indices from lat long values.  
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        dataset (str): the dataset that the files should be obtained from name e.g. "conus1_baseline_mod"
        write_dir (str): directory where the subset files will be written
        timezone (str): timezone information for start and end dates. Defaults to "UTC".

    Returns:
        dict: A dictionary in which the keys are ("vegp", "vegm", "drv_clm") and the values are file paths where the CLM files were written.

    Raises:
        FileNotFoundError: If write_dir is not a valid directory.     

    Example:

    .. code-block:: python

        filepaths = config_clm(
            ij_bounds=(375, 239, 487, 329),
            start="2005-10-01",
            end="2006-10-01",
            dataset="conus1_baseline_mod",
            write_dir="/path/to/your/chosen/directory"
        )
    """
    _validate_grid_bounds(ij_bounds)
    _validate_date(start)
    _validate_date(end)
    if not isinstance(dataset, str):
        raise TypeError("dataset name must be a string.")
    _validate_dir(write_dir)
    if not isinstance(time_zone, str):
        raise TypeError("time_zone must be a string.")

    file_type_list = ["vegp", "pfb", "drv_clm"] # get the pfb version of the vegm file
    file_paths = {}
    for file_type in file_type_list:
        entry = hf_hydrodata.get_catalog_entry(dataset=dataset,
                                                       file_type=file_type,
                                                       variable="clm_run",
                                                       period="static"
        )
        if entry is None:
            raise ValueError("No {file_type} entry matches your request.")
        print(f"processing {file_type}")
        if file_type == "vegp":
            file_path = os.path.join(write_dir, "drv_vegp.dat")
            hf_hydrodata.get_raw_file(
                file_path,
                dataset=dataset,
                file_type=file_type,
                variable="clm_run",
                period="static"
            )
            file_paths[file_type] = file_path
            print("copied vegp")
        elif file_type == "pfb":
            subset_data = hf_hydrodata.gridded.get_ndarray(entry, grid_bounds=ij_bounds)
            land_cover_data = _reshape_ndarray_to_vegm_format(subset_data)
            file_path = write_land_cover(land_cover_data, write_dir)
            file_paths[file_type] = file_path
            print("subset vegm")
        elif file_type == "drv_clm":
            file_path = os.path.join(write_dir, "drv_clmin.dat")
            hf_hydrodata.get_raw_file(file_path,
                            dataset=dataset,
                            file_type=file_type,
                            variable="clm_run",
                            period="static"
            )
            print("copied drv_clmin")
            edit_drvclmin(
                file_path=file_path,
                start=start,
                end=end,
                time_zone=time_zone
            )
            file_paths[file_type] = file_path
            print("edited drv_clmin")
    return file_paths


def _reshape_ndarray_to_vegm_format(data):
    """Reshape ndarray returned by datacatalog to vegm format.

    Parameters:
        data (ndarray) – raw subset vegm data (2d array)

    Returns:
        Ndarray reshaped to vegm format.
    """       
    _, nj, ni = data.shape
    indices = np.indices((nj, ni)) + 1
    indices = indices[::-1, :, :]
    data = np.vstack([indices, data])  # stack x,y indices on vegm                                                                              
    # transpose and reshape back into expected 2D vegm file format for the subset                                                               
    return data.transpose(1, 2, 0).reshape(-1, 25)


def subset_forcing(
        ij_bounds,
        grid,
        start,
        end,
        dataset,
        write_dir,
        time_zone="UTC",
        forcing_vars = (
            "precipitation",
            "downward_shortwave",
            "downward_longwave",
            "specific_humidity",
            "air_temp",
            "atmospheric_pressure",
            "east_windspeed",
            "north_windspeed",
        )
):
    """Subset forcing files for a given box or point of interest and write them as pfb files to a user specified directory.

    Subset forcing data will be written out as pfb files formatted for a ParFlow run with 24 hours per forcing file. Per ParFlow-CLM 
    convention separate files will be written for each variable following the standard clm variable naming convention. 

    Forcing file outputs will be numbered starting with 0000 and data will start at midnight local time for the timezone that has been 
    provided. If no timezone is provided it will default to midnight UTC.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as i,j index values where 0,0 is the lower left hand corner 
            of a domain. ij_bounds are given relative to whatever grid is being used for the subset. Use the latlon_to_ij function to determine ij 
            indices from lat long values.  
        grid (str): The spatial grid that the ij indices are calculated relative to and that the subset data will be returned on. Possible 
            values: "conus1" or "conus2"
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        dataset (str): forcing dataset name from the HydroData catalog that the forcing files will be subset from e.g. "NLDAS2". 
        write_dir (str): directory where the subset files will be written
        time_zone (str): timezone information for start and end dates. Data will be subset starting at midnight in the specified timezone. 
            Defaults to "UTC".
        forcing_vars (tuple[str]): tuple of forcing variables to subset. By default all 8 variables needed to run ParFlow-CLM will be subset. 

    Returns:
        dict: A dictionary in which the keys are the forcing variable names and the values are lists of file paths where the subset data were written.

    Raises:
        FileNotFoundError: If write_dir is not a valid directory.

    Example:

    .. code-block:: python

        filepaths = subset_forcing(
            ij_bounds=(1225, 1738, 1347, 1811), 
            grid="conus2", 
            start="2005-11-01", 
            end="2005-12-01", 
            dataset="CW3E", 
            write_dir="/path/to/your/chosen/directory",
            forcing_vars=("precipitation", "air_temp")
        )
    """
    _validate_grid_bounds(ij_bounds)
    _validate_grid(grid)
    _validate_date(start)
    _validate_date(end)
    if not isinstance(dataset, str):
        raise TypeError("dataset name must be a string.")
    _validate_dir(write_dir)
    if not isinstance(time_zone, str):
        raise TypeError("time_zone must be a string.")
    if not all(isinstance(var, str) for var in forcing_vars):
        raise TypeError("All variable names should be strings.")    
    outputs = {}
    start_date = get_UTC_time(start, time_zone)
    end_date = get_UTC_time(end, time_zone)
    exit_event = threading.Event()    
    threads = [threading.Thread(target=_subset_forcing_variable,
                                args=(variable, ij_bounds, grid, start_date, end_date, dataset, write_dir, time_zone, outputs, exit_event))
        for variable in forcing_vars]
    
    for thread in threads:
        thread.start()

    try:
        for thread in threads:
            thread.join()
    except KeyboardInterrupt:
        print("Interrupted. Stopping threads...")
        exit_event.set()
        
        for thread in threads:
            thread.join()

        print("All threads stopped.")
    finally:
        return outputs


def _subset_forcing_variable(variable, ij_bounds, grid, start_date, end_date, dataset, write_dir, time_zone, outputs, exit_event):
    """Helper function for subset_forcing that subsets data for one variable for the specified dates and dataset."""

    entry = hf_hydrodata.get_catalog_entry(
        dataset=dataset, variable=variable, grid=grid, file_type="pfb", period="hourly"
    )
    if entry is None:
        raise ValueError(f"No {variable} entry matches your request.")
    
    day = 1
    date = start_date
    delta = timedelta(days=1)
    outputs[variable] = []
    print(f"Reading {variable} pfb sequence")
    
    while date < end_date and not exit_event.is_set():
        start_time = date
        end_time = date + delta
        # we need to distinguish between UTC and non-UTC as the datacatalog returns the wrong answer
        # for requests that start reading from the middle of a file and span multiple files

        if time_zone == "UTC":
            subset_data = hf_hydrodata.gridded.get_ndarray(
                entry,
                start_time=start_time,
                end_time=end_time,
                grid_bounds=ij_bounds,
            )
        else:
            next_day_midnight = datetime(end_time.year, end_time.month, end_time.day)
            data1 = hf_hydrodata.gridded.get_ndarray(
                entry,
                start_time=start_time,
                end_time=next_day_midnight,
                grid_bounds=ij_bounds,
            )
            data2 = hf_hydrodata.gridded.get_ndarray(
                entry,
                start_time=next_day_midnight,
                end_time=end_time,
                grid_bounds=ij_bounds,
            )                
            subset_data = np.concatenate((data1, data2), axis=0)
            
        assert subset_data.shape[0] == 24, "attempted to write more than 24 hours of data to a pfb file"

        paths = hf_hydrodata.gridded.get_file_paths(
            entry, start_time=start_time, end_time=end_time
        )
        
        write_path = os.path.join(write_dir,
                                  _adjust_filename_hours(os.path.basename(paths[0]),day)
        )
        outputs[variable].append(write_path)
        write_pfb(write_path, subset_data[:, :, :], dist=False)
        date = date + delta
        day = day + 1
    if not exit_event.is_set():
        print(f"Finished writing {variable} to folder")
    

def _adjust_filename_hours(filename, day):
    """Adjust the forcing filename hours so that they match what a parflow simulation expects on each day of the simulation.

    The first day of the simulation the hours will be “.000001_to_000024.”, the second day the hours will be “.000025_to_000048.” and so on. 
    This is in case the first day of simulation does not coincide with the first day of the water year (Oct 1st), as the dataset 
    filenames assume day 1 is Oct 1st. The input and output filenames must match the regular expression “..[0-9]{6}_to_[0-9]{6}.*”

    Parameters:
        filename (str) – original forcing filename
        day (int) – day relative to the start date of forcing file subsetting

    Returns:
        The forcing filename with adjusted hours.

    Raises:
        AssertionError – If the input or output filename string do not match the above regex.
    """
    assert day >= 1
    s1, s2, s3, s4 = filename.split(".")
    assert s1 != "" and s2 != "" and s4 != "", "invalid forcing filename"
    pattern = re.compile("[0-9]{6}_to_[0-9]{6}")
    assert pattern.fullmatch(s3) is not None, "invalid forcing filename"
    start = str(24 * (day - 1) + 1).rjust(6, "0")
    end = str(24 * day).rjust(6, "0")
    s3 = start + "_to_" + end
    assert pattern.fullmatch(s3) is not None, "invalid adjusted forcing filename"
    return ".".join([s1, s2, s3, s4])
    
    
def edit_runscript_for_subset(
    ij_bounds, runscript_path, write_dir=None, runname=None, forcing_dir=None
):
    """Modify a ParFlow run script for a new subdomain run

    This function is designed to start from a national ParFlow run script template and:  
        1. Modify the geometry to reflect the bounds of the desired ij_bounds (i.e. the number of grid cells in the x and y direction and the upper bounds of the geometry)
        2. update the runname to for the desired new run
        3. update the location of the climate forcings for the new run

    If the runname is None and write_dir is the directory containing
    the runscript file, the runscript file will be overwritten.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as i,j index values where 0,0 is 
            the lower left hand corner of a domain. ij_bounds are given relative to whatever grid is being used for the 
            subset. Use the latlon_to_ij function to determine ij indices from lat long values.
        runscript_path (str): absolute path to the template parflow runscript file
        write_dir (str): directory where the new template file will be written.
            If it is None, defaults to the directory containing the runscript.
        runname (str): name for the new parflow run. If it is None, defaults to 
            the runscript's previous runname.
        forcing_dir (str): path to the directory containing the subset forcing files.
            If it is None, defaults to the runscript's previous forcing directory path.

    Returns:
        str: Path to the new runscript file that will be created.

    Raises:
        FileNotFoundError: If runscript_path is not an existing file or write_dir is not a valid directory.

    Example:

    .. code-block:: python

        runscript_path = edit_runscript_for_subset(
            ij_bounds=(375, 239, 487, 329),
            runscript_path="/path/to/your/original/runscript",
            runname="my_conus1_run",
            forcing_dir="/path/to/your/forcing/directory"
        )
    """
    _validate_grid_bounds(ij_bounds)
    if not os.path.isfile(runscript_path):
        raise FileNotFoundError("runscript_path must be a valid file path")
    if write_dir is None:
        write_dir = os.path.dirname(runscript_path)
    _validate_dir(write_dir)
    if runname is not None and not isinstance(runname, str):
        raise TypeError("runname must be a string.")
    if forcing_dir is not None:
        _validate_dir(forcing_dir)
    
    # load in the reference pfidb or yaml specified by the user
    run = Run.from_definition(runscript_path)
    _, file_extension = os.path.splitext(runscript_path)
    if runname is not None:
        run.set_name(runname)
        print(
            f"New runname: {runname} provided, a new {file_extension[1:]} file will be created"
        )
    else:
        print(f"No runname provided, old {file_extension[1:]} file will be overwritten")

    # checks if we're running with clm
    if forcing_dir is not None:
        print(
            f"Climate forcing directory has been changed to {forcing_dir} in runscript."
        )
        run.Solver.CLM.MetFilePath = forcing_dir
    else:
        print("No forcing directory provided, run.Solver.CLM.MetFilePath key not set")

    # getting the subset ni/nj to update keys
    imin, jmin, imax, jmax = ij_bounds
    ni, nj = imax - imin, jmax - jmin
    run.ComputationalGrid.NY = int(nj)
    run.ComputationalGrid.NX = int(ni)
    print(f"ComputationalGrid.NY set to {nj} and NX to {ni}")

    domain_type = run.GeomInput.domaininput.InputType
    if domain_type == "SolidFile":
        print(
            "GeomInput.domaininput.InputType detected as SolidFile, no additional keys to change for subset"
        )
    else:
        print(
            f"GeomInput.domaininput.InputType detected as Box, updating Geom.domain.Upper.X to {ni * 1000} and Y to {nj * 1000} to match subset"
        )
        run.Geom.domain.Upper.X = ni * 1000
        run.Geom.domain.Upper.Y = nj * 1000

    print(f"Updated runscript written to {write_dir}")
    file_path, _ = run.write(
        working_directory=write_dir, file_format=f"{file_extension[1:]}"
    )
    return file_path


def copy_files(read_dir, write_dir):
    """Copy all files from read_dir to write_dir.

    Args:
        read_dir (str): read-from directory path
        write_dir (str): write-to directory path

    Raises:
        FileNotFoundError: If read_dir or write_dir is not a valid directory path.

    Example:

    .. code-block:: python

        copy_files(
            read_dir="/path/to/read-from/directory",
            write_dir="/path/to/write-to/directory"
        )
    """
    _validate_dir(read_dir)
    _validate_dir(write_dir)
    for filename in os.listdir(read_dir):
        file_path = os.path.join(read_dir, filename)
        if os.path.isfile(file_path):
            shutil.copy(file_path, write_dir)


def change_filename_values(
    runscript_path,
    write_dir=None,
    runname=None,
    slopex=None,
    slopey=None,
    solidfile=None,
    init_press=None,
    indicator=None,
    depth_to_bedrock=None,
    mannings=None,
    evap_trans=None,
):
    """Change the filenames of input files.

    This function will update the paths to input files in a ParFlow run script.  The provided arguments will reset the corresponding parflow keys to match the user specified paths to input files.  File names can be specified with or without relative or absolute file paths. If no path is provided ParFlow will expect the input files to be present in the run directory at the time of simulation. 

    Note that this will only change paths for keys that already exist in the template ParFlow run script you are starting from and will not reconfigure a run to use new keys (for example if you are not starting from a run script that uses a solid file, adding a new solid file path will not configure the run to use a solid file).  

    Refer to the ParFlow manual for additional information on any of the keys listed above. 

    
    If the runname is None and write_dir is the directory containing the runscript file, the runscript file will be overwritten.

    Args:
        runscript_path (str): path to the runscript file (yaml or pfidb)
        write_dir (str): directory where the new template file will be written. If
            it is None, defaults to the directory containing the runscript file.
        runname (str): name of the new parflow run
        slopex (str): new slopex filename (and path)
        slopey (str): new slopey filename (and path)
        solidfile (str): new solidfile filename (and path)
        init_press (str): new initial pressure filename (and path)
        indicator (str): new indicator input filename (and path)
        depth_to_bedrock (str): new depth to bedrock filename (and path)
        mannings (str): new mannings filename (and path)
        evap_trans (str): new evapotranspiration filename (and path)

    Returns:
        str: Path to the new runscript file that will be created.

    Raises:
        FileNotFoundError: If runscript_path is not a valid file path or write_dir is not a valid directory.

    Example:

    .. code-block:: python

        runscript_path = change_filename_values(
            runscript_path="/path/to/your/original/runscript",
            runname="my_conus1_run",
            init_press="/filename/of/initial/pressure/pfb/file"
        )
    """
    if not os.path.isfile(runscript_path):
        raise FileNotFoundError("runscript_path must be a valid file path")
    if write_dir is None:
        write_dir = os.path.dirname(runscript_path)
    _validate_dir(write_dir)

    _, file_extension = os.path.splitext(runscript_path)
    run = Run.from_definition(runscript_path)
    if runname is not None:
        run.set_name(runname)
        print(
            f"New runname: {runname} provided, a new {file_extension[1:]} file will be created"
        )
    # check which input files are not none, to update the key
    if slopex is not None:
        run.TopoSlopesX.FileName = slopex
        print(f"X Slopes filename changed to {slopex}")
    if slopey is not None:
        run.TopoSlopesY.FileName = slopey
        print(f"Y Slopes filename changed to {slopey}")
    if solidfile is not None:
        run.GeomInput.domaininput.FileName = solidfile
        print(f"Solidfile filename changed to {solidfile}")
    if init_press is not None:
        run.Geom.domain.ICPressure.FileName = init_press
        print(f"Initial pressure filename changed to {init_press}")
    if indicator is not None:
        run.Geom.indi_input.FileName = indicator
        print(f"Indicator filename changed to {indicator}")
    if depth_to_bedrock is not None:
        run.Geom.domain.FBz.FileName = depth_to_bedrock
        print(f"Depth to bedrock filename changed to {depth_to_bedrock}")
    if mannings is not None:
        run.Mannings.FileName = mannings
        print(f"Mannings filename changed to {mannings}")
    if evap_trans is not None:
        run.Solver.EvapTrans.FileName = evap_trans
        print(f"Evaptrans filename changed to {evap_trans}")

    print(f"Updated runscript written to {write_dir}")
    file_path, _ = run.write(
        working_directory=write_dir, file_format=f"{file_extension[1:]}"
    )
    return file_path


def dist_run(P, Q, runscript_path, working_dir=None, dist_clim_forcing=True):
    """Distribute parflow input files for parallel computing.

    This function will distribute input files to P grids in the
    x direction and Q grids in the y direction. If dist_clim_forcing
    is true, forcing files will be distributed as well according to the
    same topology. If working_dir is different that the directory containing
    the runscript file, the edited runscipt file will be written to working_dir.

    Args:
        P (int): number of grids (processes) to create in the x direction
        Q (int): number of grids (processes) to create in the y direction
        runscript_path (str): path to the runscript file (yaml or pfidb)
        working_dir (str): directory containing the files to be distributed.
            If it is None, it defaults to the directory containing the runscript file.
        dist_clim_forcing (bool): if true, distribute forcing files

    Returns:
        str: Path to the edited runscript file that will be created.

    Raises:
        FileNotFoundError: If runscript_path is not a valid file path or write_dir is not a valid directory.

    Example:

    .. code-block:: python

        runscript_path = dist_run(
            P=2,
            Q=2,
            runscript_path="/path/to/your/original/runscript",
            dist_clim_forcing=False
        )
    """
    if not isinstance(P, int) or P <= 0:
        raise TypeError("P must be a positive integer.")
    if not isinstance(Q, int) or Q <= 0:
        raise TypeError("Q must be a positive integer.")
    if not os.path.isfile(runscript_path):
        raise FileNotFoundError("runscript_path must be a valid file path")    
    if working_dir is None:
        working_dir = os.path.dirname(runscript_path)
    _validate_dir(working_dir)
    
    run = Run.from_definition(runscript_path)
    run.Process.Topology.P = P
    run.Process.Topology.Q = Q

    if dist_clim_forcing:
        print("Distributing your climate forcing")
        forcing_dir = run.Solver.CLM.MetFilePath
        _validate_dir(forcing_dir)
        for filename in pathlib.Path(forcing_dir).glob("*.pfb"):
            run.dist(filename)
    else:
        print("no forcing directory provided, only distributing static inputs")

    static_input_paths = pathlib.Path(working_dir).glob("*.pfb")
    max_nz = 0
    for path in static_input_paths:
        input_array = read_pfb(path)
        nz = input_array.shape[0]
        if nz > max_nz:
            max_nz = nz
        run.dist(path)
        print(f"Distributed {os.path.basename(path)} with NZ {nz}")
    run.ComputationalGrid.NZ = max_nz

    _, file_extension = os.path.splitext(runscript_path)
    file_path, _ = run.write(working_directory=working_dir, file_format=f"{file_extension[1:]}")
    return file_path
