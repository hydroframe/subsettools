"""This module contains subsetting functions for Parflow runs.
"""

import os
import shutil
import pathlib
import subprocess
from datetime import datetime, timedelta

import numpy as np
import pytz
from hf_hydrodata import grid, gridded
from parflow import Run
from parflow.tools.io import read_pfb, write_pfb
from .subset_utils import (
    get_conus_hucs_indices,
    indices_to_ij,
    reshape_ndarray_to_vegm_format,
    write_land_cover,
    edit_drvclmin,
    adjust_filename_hours,
)


def huc_to_ij(huc_list, grid):
    """Get the grid ij-bounds of the area defined by the the huc IDs in huc_list.

       All huc IDs in huc_list must be the same length (hucs of the same
       level). Supported grids are "conus1" and "conus2".

    Args:
        huc_list (list[str]): a list of huc IDs
        grid (str): "conus1" or "conus2"

    Returns:
        A tuple of the form (imin, jmin, imax, jmax) representing the bounds
        in the conus grid of the area defined by the huc IDs in huc_list.

    Raises:
        AssertionError: If all huc IDs are not the same length.
    """
    huc_len = len(huc_list[0])
    assert all(
        [len(huc) == huc_len for huc in huc_list]
    ), "All huc IDs should have the same length!"
    conus_hucs, _, indices_j, indices_i = get_conus_hucs_indices(huc_list, grid)
    return indices_to_ij(conus_hucs, indices_j, indices_i)


def latlon_to_ij(latlon_bounds, grid):
    """Get the conus ij-bounds of the area defined by the the latitute/longitude bounds (latlon_bounds) in the conus grid.

       Supported grids are "conus1" and "conus2".

    Args:
        latlon_bounds (list[float]): list of the form [[lat1, lon1], [lat2, lon2]].
            [lat1, lon1] and [lat2, lon2] are the two points defining the conus
            bounding box.
        grid (str): "conus1" or "conus2"

    Returns:
        A tuple of the form (imin, jmin, imax, jmax) representing the bounds
        in the conus grid of the area defined by latlon_bounds.
    """
    grid = grid.lower()
    point0 = grid.from_latlon(grid, latlon_bounds[0][1], latlon_bounds[0][0])
    point1 = grid.from_latlon(grid, latlon_bounds[1][1], latlon_bounds[1][0])
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
        huc_list (list[str]): a list of huc IDs  
        grid (str): "conus1" or "conus2"  
        write_dir (str): directory path where the mask and solid files will be written  

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

    subprocess.run(
        [
            os.path.join(os.environ["PARFLOW_DIR"], "bin", "pfmask-to-pfsol"),
            "--mask",
            mask_file_path,
            "--pfsol",
            os.path.join(write_dir, "solidfile.pfsol"),
            "--vtk",
            os.path.join(write_dir, "mask_vtk.vtk"),
            "--z-bottom",
            "0.0",
            "--z-top",
            z_total,
        ],
        check=True,
    )

    print(f"Wrote solidfile.pfsol and mask_vtk.vtk with total z of {z_total} meters")


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

    Raises:
        AssertionError: If write_dir is not a valid directory.
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"
    for var in var_list:
        entry = gridded.get_catalog_entry(
            dataset=dataset, file_type="pfb", period="static", variable=var
        )
        if entry is not None:
            subset_data = gridded.get_ndarray(entry, grid_bounds=ij_bounds)
            write_pfb(os.path.join(write_dir, f"{var}.pfb"), subset_data, dist=False)
            print(f"Wrote {var}.pfb in specified directory.")
        else:
            print(f"{var} not found in dataset {dataset}")


def subset_press_init(ij_bounds, dataset, date, write_dir, time_zone="UTC"):
    """Subset the initial pressure file.

    This represent the pressure one hour before midnight on the day before the start
    date of the simulation.

    Args:
        ij_bounds (Tuple[int]): bounding box for subset
        dataset (str): dataset name e.g. "conus1_baseline_mod"
        write_dir (str): directory where the subset file will be written
        time_zone (str): time_zone to calculate initial pressure datetime. Defaults to "UTC".

    Returns:
        The filename of the subset file, which includes datetime information, so that it can be
        used by later functions (e.g. edit_runscript_for_subset)

    Raises:
        AssertionError: If write_dir is not a valid directory.
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    entry = gridded.get_catalog_entry(
        dataset=dataset, file_type="pfb", variable="pressure_head", period="hourly"
    )

    # assumes time is UTC 0 like CONUS runs, so can remain time unaware and grab the right pressure
    new_date = datetime.strptime(date, "%Y-%m-%d") - timedelta(hours=1)
    if entry is None:
        print(f"No pressure file found for {new_date} in dataset {dataset}")
        return None

    if time_zone != "UTC":
        print(f"Converting the requested datetime from UTC0 to {time_zone}")
        new_date = new_date.replace(tzinfo=pytz.UTC).astimezone(
            pytz.timezone(time_zone)
        )
    date_string = new_date.strftime("%Y.%m.%d:%H.%M.%S_UTC0")

    subset_data = gridded.get_ndarray(
        entry, grid_bounds=ij_bounds, start_time=new_date
    )

    out_file = f"{dataset}_{date_string}_press.pfb"
    write_pfb(os.path.join(write_dir, out_file), subset_data[0, :, :, :], dist=False)
    print(f"Wrote {out_file} in specified directory.")
    return out_file


def config_clm(ij_bounds, start, end, dataset, write_dir):
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

    Returns:
        The filename of the subset file, which includes datetime information,
        so that it can be used by later functions (e.g. edit_runscript_for_subset).
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    file_type_list = ["vegp", "vegm", "drv_clm"]
    for file_type in file_type_list:
        print(f"processing {file_type}")
        if file_type == "vegp":
            gridded.get_raw_file(
                os.path.join(write_dir, "drv_vegp.dat"),
                dataset=dataset,
                file_type=file_type,
                variable="clm_run",
                period="static"
            ) 
            print("copied vegp")
        elif file_type == "vegm":
            entry = gridded.get_catalog_entries(
                dataset=dataset,
                file_type=file_type,
                variable="clm_run",
                period="static"
            )[0]
            subset_data = gridded.get_ndarray(entry, grid_bounds=ij_bounds)
            land_cover_data = reshape_ndarray_to_vegm_format(subset_data)
            write_land_cover(land_cover_data, write_dir)
            print("subset vegm")
        elif file_type == "drv_clm":
            file_path = os.path.join(write_dir, "drv_clmin.dat")
            gridded.get_raw_file(file_path,
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
            )
            print("edited drv_clmin")


def subset_forcing(ij_bounds, grid, start, end, dataset, write_dir):
    """Get and subset the forcing files filtered by grid, dataset and start/end dates.

    The forcing filenames will be adjusted if the start date does not coincide with the
    start of the water year (Oct 1st), so that they match what a parflow simulation expects.

    Args:
        ij_bounds (Tuple[int]): bounding box for subset
        grid (str): "conus1" or "conus2"
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        dataset (str): forcing dataset name e.g. "NLDAS2"
        write_dir (str): directory where the subset file will be written

    Returns:
        A dictionary in which the keys are the forcing variables and the values are lists of
        subset file paths. The return value is useful for logging purposes and can be discarded
        otherwise.

    Raises:
        AssertionError: If write_dir is not a valid directory.
    """
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    var_list = (
        "precipitation",
        "downward_shortwave",
        "downward_longwave",
        "specific_humidity",
        "air_temp",
        "atmospheric_pressure",
        "east_windspeed",
        "north_windspeed",
    )
    outputs = {}

    for var in var_list:
        entries = gridded.get_catalog_entry(
            dataset=dataset, variable=var, grid=grid, file_type="pfb", period="hourly"
        )

        day = 1
        start_date = datetime.strptime(start, "%Y-%m-%d")
        end_date = datetime.strptime(end, "%Y-%m-%d")
        delta = timedelta(days=1)
        outputs[var] = []
        print(f"Reading {var} pfb sequence")

        while start_date < end_date:
            subset_data = gridded.get_ndarray(
                entries,
                start_time=start_date,
                end_time=start_date + delta,
                grid_bounds=ij_bounds,
            )
            paths = gridded.get_file_paths(
                entries, start_time=start_date, end_time=start_date + delta
            )
            write_paths = [
                os.path.join(
                    write_dir, adjust_filename_hours(os.path.basename(path), day)
                )
                for path in paths
            ]
            outputs[var] += write_paths
            for path in write_paths:
                write_pfb(path, subset_data[:, :, :], dist=False)
                day = day + 1
            start_date = start_date + delta

        print(f"Finished writing {var} to folder")

    return outputs


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
       AssertionError: If runscript_path is not a valid file path.
    """
    assert os.path.isfile(runscript_path), "runscript_path must be a valid file path"

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
