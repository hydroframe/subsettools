"""This module contains subsetting functions for Parflow runs.
"""

import os
import shutil
import pathlib
import subprocess
from datetime import datetime, timedelta
from threading import Thread
import numpy as np
import hf_hydrodata
import re
from parflow import Run
from parflow.tools.io import read_pfb, write_pfb
from .subset_utils import (
    get_conus_hucs_indices,
    indices_to_ij,
    write_land_cover,
    edit_drvclmin,
    get_UTC_time,
)


def huc_to_ij(huc_list, grid):
    """Get the grid ij-bounds of the area defined by the the HUC IDs in huc_list.

       All HUC IDs in huc_list must be the same length (HUCs of the same
       level). All HUCs should be adjacent. Supported grids are "conus1" 
       and "conus2". If a HUC is only partially covered by the provided grid,
       the grid bounds for the covered area will be returned.

    Args:
        huc_list (list[str]): a list of HUC IDs
        grid (str): "conus1" or "conus2"

    Returns:
        A tuple of the form (imin, jmin, imax, jmax) representing the bounds
        in the conus grid of the area defined by the huc IDs in huc_list.

    Raises:
        AssertionError: If all HUC IDs are not the same length.
        ValueError: If the area defined by the provided HUCs is not part of the given grid.
    """
    huc_len = len(huc_list[0])
    assert all(
        [len(huc) == huc_len for huc in huc_list]
    ), "All huc IDs should have the same length!"
    conus_hucs, _, indices_j, indices_i = get_conus_hucs_indices(huc_list, grid)
    if indices_i.size == 0 or indices_j.size == 0:
        raise ValueError(f"The area defined by the provided HUCs is not part of the {grid} grid.")  
    return indices_to_ij(conus_hucs, indices_j, indices_i)


def latlon_to_ij(latlon_bounds, grid):
    """Get the conus ij-bounds of the area defined by the the latitute/longitude bounds (latlon_bounds) in the conus grid.

       Supported grids are "conus1" and "conus2".

    Args:
        latlon_bounds (List[List[float]]): list of the form [[lat1, lon1], [lat2, lon2]].
            [lat1, lon1] and [lat2, lon2] are the two points defining the conus
            bounding box.
        grid (str): "conus1" or "conus2"

    Returns:
        A tuple of the form (imin, jmin, imax, jmax) representing the bounds
        in the conus grid of the area defined by latlon_bounds.
    """
    grid = grid.lower()
    assert grid in ["conus1", "conus2"], "invalid grid name"
    assert len(latlon_bounds) == 2, "please provide 2 latlon points"
    assert len(latlon_bounds[0]) == 2 and len(latlon_bounds[1]) == 2, "invalid latlon point"
    point0 = hf_hydrodata.grid.to_ij(grid, latlon_bounds[0][0], latlon_bounds[0][1])
    point1 = hf_hydrodata.grid.to_ij(grid, latlon_bounds[1][0], latlon_bounds[1][1])
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
    """Create mask and solid input files for the ParFlow simulation.  

    Args:  
        huc_list (list[str]): a list of HUC IDs  
        grid (str): "conus1" or "conus2"  
        write_dir (str): directory path where the mask and solid files will be written  

    Returns:
        A dictionary of paths with keys ("mask", "mask_vtk", "solid") and values filepaths to the created files.

    Raises:  
        AssertionError: If write_dir is not a valid directory.  
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    conus_hucs, sel_hucs, indices_j, indices_i = get_conus_hucs_indices(huc_list, grid)
    arr_jmin = np.min(indices_j)
    arr_jmax = np.max(indices_j) + 1  # right bound inclusive
    ij_bounds = indices_to_ij(conus_hucs, indices_j, indices_i)

    nj = ij_bounds[3] - ij_bounds[1]
    ni = ij_bounds[2] - ij_bounds[0]

    # checks conus1 / 2 grid and assigns appripriate dz and z_total for making the mask and solid file
    if grid.lower() == "conus1":
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
    mask_clip = np.flip(
        mask_clip, 1
    )  # This flip tooks the section we just subset and puts it in the appropriate parflow orientation
    mask_clip = mask_clip.astype(float)
    mask_file_path = os.path.join(write_dir, "mask.pfb")
    write_pfb(mask_file_path, mask_clip, dx=1000, dy=1000, dz=layz, dist=False)
    print("Wrote mask.pfb")
    mask_vtk_path = os.path.join(write_dir, "mask_vtk.vtk")
    solid_file_path = os.path.join(write_dir, "solidfile.pfsol")
    
    subprocess.run(
        [
            os.path.join(os.environ["PARFLOW_DIR"], "bin", "pfmask-to-pfsol"),
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
    )

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
        "depth_to_bedrock",
        "pme",
        "ss_pressure_head",
    ),
):
    """Subset static inputs from dataset required to do a baseline run.

    Args:
        ij_bounds (Tuple[int]): bounding box for subset
        dataset (str): dataset name e.g. "conus1_domain"
        write_dir (str): directory where the subset files will be written
        var_list (List[str]): list of variables to subset from the dataset

    Returns:
        A dictionary in which the keys are the static variable names and the values are
        file paths where the subset data were written.
        
    Raises:
        AssertionError: If write_dir is not a valid directory.
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"
    file_paths = {}
    for var in var_list:
        entry = hf_hydrodata.gridded.get_catalog_entry(
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
    """Subset the initial pressure file.

    This represents the pressure one hour before midnight on the day before the start
    date of the simulation.

    Args:
        ij_bounds (Tuple[int]): bounding box for subset
        dataset (str): dataset name e.g. "conus1_baseline_mod"
        date (str): simulation start date, in the form 'yyyy-mm-dd'
        write_dir (str): directory where the subset file will be written
        time_zone (str): time_zone to calculate initial pressure datetime. Defaults to "UTC".

    Returns:
        The filepath of the subset file, which includes datetime information, so that it can be
        used by later functions (e.g. edit_runscript_for_subset)

    Raises:
        AssertionError: If write_dir is not a valid directory.
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    entry = hf_hydrodata.gridded.get_catalog_entry(
        dataset=dataset, file_type="pfb", variable="pressure_head", period="hourly"
    )
    
    if entry is None:
        print(f"No pressure file found for in dataset {dataset}")
        return None
    
    new_date = get_UTC_time(date, time_zone) - timedelta(hours=1)
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
    """Get and subset the clm drivers associated with the run dataset.

    vegm, vep and drv_clmin files will be written in the specified static
    input directory. For consistency, dataset must be the same dataset that
    is passed to subset_press_init().

    Args:
        ij_bounds (Tuple[int]): bounding box for subset
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        dataset (str): dataset name e.g. "conus1_baseline_mod"
        write_dir (str): directory where the subset files will be written
        timezone (str): timezone information for start and end dates 

    Returns:
        A dictionary in which the keys are ("vegp", "vegm", "drv_clm") and the values are
        file paths where the CLM files were written.

    Raises:
        AssertionError: If write_dir is not a valid directory.        
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    file_type_list = ["vegp", "vegm", "drv_clm"]
    file_paths = {}
    for file_type in file_type_list:
        print(f"processing {file_type}")
        if file_type == "vegp":
            file_path = os.path.join(write_dir, "drv_vegp.dat")
            hf_hydrodata.gridded.get_raw_file(
                file_path,
                dataset=dataset,
                file_type=file_type,
                variable="clm_run",
                period="static"
            )
            file_paths[file_type] = file_path
            print("copied vegp")
        elif file_type == "vegm":
            entry = hf_hydrodata.gridded.get_catalog_entries(
                dataset=dataset,
                file_type=file_type,
                variable="clm_run",
                period="static"
            )[0]
            subset_data = hf_hydrodata.gridded.get_ndarray(entry, grid_bounds=ij_bounds)
            land_cover_data = _reshape_ndarray_to_vegm_format(subset_data)
            file_path = write_land_cover(land_cover_data, write_dir)
            file_paths[file_type] = file_path
            print("subset vegm")
        elif file_type == "drv_clm":
            file_path = os.path.join(write_dir, "drv_clmin.dat")
            hf_hydrodata.gridded.get_raw_file(file_path,
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
    """Subset forcing files for a given box or point of interest.

    Subset forcing data will be written out as pfb files formatted for a ParFlow run with 24 hours per forcing file. Per ParFlow-CLM convention separate files will be written for each variable following the standard clm variable naming convention. 

    Forcing file outputs will be numbered starting with 0000 and data will start at midnight local time for the timezone that has been provided. If no timezone is provided it will default to midnight UTC.

    Args:
        ij_bounds (Tuple[int]): bounding box for subset. This should be given as i,j index values where 0,0 is the lower left hand corner of a domain. ij_bounds are given to whatever grid is being used for the subset. Use the latlon_to_ij function to determine ij indices from lat long values.  
        grid (str): The spatial grid that the ij indices are calculated relative to and that the subset data will be returned on. Possible values: "conus1" or "conus2"
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        dataset (str): forcing dataset name from the HydroData catalog e.g. "NLDAS2". 
        write_dir (str): directory where the subset file will be written
        timezone (str): timezone information for start and end dates. Data will be subset starting at midnight in the specified timezone.
        forcing_vars (Tuple[str]): tuple of forcing variables to subset. By default all 8 variables needed to run ParFlow-CLM will be subset. 

    Returns:
        A dictionary in which the keys are the forcing variable names and the values are lists of
        file paths where the subset data were written.

    Raises:
        AssertionError: If write_dir is not a valid directory.

    Examples:
        subset_forcing(ij_bounds=(1225, 1738, 1347, 1811), 
                        grid="conus2", 
                        start="2005-11-01", 
                        end="2005-12-01", 
                        dataset="CW3E", 
                        write_dir="/path/to/your/chosen/directory",
        )
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    outputs = {}
    start_date = get_UTC_time(start, time_zone)
    end_date = get_UTC_time(end, time_zone)
    
    threads = [Thread(target=_subset_forcing_variable, args=(variable, ij_bounds, grid, start_date, end_date, dataset, write_dir, time_zone, outputs))
        for variable in forcing_vars]
    
    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

    return outputs


def _subset_forcing_variable(variable, ij_bounds, grid, start_date, end_date, dataset, write_dir, time_zone="UTC", outputs={}):
    """Helper function for subset_forcing that subsets data for one variable for the specified dates and dataset."""
    
    entry = hf_hydrodata.gridded.get_catalog_entry(
        dataset=dataset, variable=variable, grid=grid, file_type="pfb", period="hourly"
    )
    
    day = 1
    date = start_date
    delta = timedelta(days=1)
    outputs[variable] = []
    print(f"Reading {variable} pfb sequence")
    
    while date < end_date:
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
    """Edit the parflow runscript keys.

    The geometry, runname and forcing directory keys will be reset in the new runscript. 
    If the runname is None and write_dir is the directory containing the runscript file, the 
    runscript file will be overwritten.

    Args:
        ij_bounds (Tuple[int]): bounding box for subset
        runscript_path (str): absolute path to the parflow runscript file.
        write_dir (str): directory where the new template file will be written.
            If it is None, defaults to the directory containing the runscript.
        runname (str): name for the new parflow run. If it is None, defaults to 
            the runscript's previous runname.
        forcing_dir (str): path to the directory containing the subset forcing files.
            If it is None, defaults to the runscript's previous forcing directory path.

    Returns:
        Path to the new runscript file that will be created.

    Raises:
       AssertionError: If runscript_path is not a valid file path or if forcing_dir is not a valid directory path.
    """
    assert os.path.isfile(runscript_path), "runscript_path must be a valid file path"
    if forcing_dir is not None:
        assert os.path.isdir(forcing_dir), "forcing_dir must be a valid directory path"
    
    if write_dir is None:
        write_dir = os.path.dirname(runscript_path)
    
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
        AssertionError: If read_dir or write_dir is not a valid directory path.
    """
    assert os.path.isdir(read_dir), "read_dir must be a directory"
    assert os.path.isdir(write_dir), "write_dir must be a directory"
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

    The provided arguments will reset the corresponding parflow keys in the new
    runscript. If the runname is None and write_dir is the directory containing
    the runscript file, the runscript file will be overwritten.

    Args:
        runscript_path (str): path to the runscript file (yaml or pfidb)
        write_dir (str): directory where the new template file will be written. If
            it is None, defaults to the directory containing the runscript file.
        runname (str): name of the new parflow run
        slopex (str): new slopex filename
        slopey (str): new slopey filename
        solidfile (str): new solidfile filename
        init_press (str): new initial pressure filename
        indicator (str): new indicator input filename
        depth_to_bedrock (str): new depth to bedrock filename
        mannings (str): new mannings filename
        evap_trans (str): new evapotranspiration filename

    Returns:
        Path to the new runscript file that will be created.

    Raises:
        AssertionError: If runscript_path is not a valid file path.
    """
    assert os.path.isfile(
        runscript_path
    ), "runscript_path must be a valid path to an existing file"

    if write_dir is None:
        write_dir = os.path.dirname(runscript_path)

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
    """Distribute parflow input files.

    This function will distribute static input files to P grids in the
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
        Path to the edited runscript file that will be created.

    Raises:
        AssertionError: If runscript_path is not a valid file path.
    """
    assert os.path.isfile(runscript_path), "runscript_path must be a valid file path"

    if working_dir is None:
        working_dir = os.path.dirname(runscript_path)
    
    run = Run.from_definition(runscript_path)

    run.Process.Topology.P = P
    run.Process.Topology.Q = Q

    if dist_clim_forcing:
        print("Distributing your climate forcing")
        forcing_dir = run.Solver.CLM.MetFilePath
        for filename in pathlib.Path(forcing_dir).glob("*.pfb"):
            run.dist(filename)
    else:
        print("no forcing dir provided, only distributing static inputs")

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
