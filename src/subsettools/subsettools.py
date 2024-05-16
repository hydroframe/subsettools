"""This module contains subsetting functions for Parflow runs."""

import os
import shutil
import subprocess
from datetime import datetime, timedelta
import threading
import re
import warnings
import numpy as np
import hf_hydrodata
from parflow import Run
from parflow.tools.io import read_pfb, write_pfb
from .subset_utils import (
    write_land_cover,
    edit_drvclmin,
    get_utc_time,
    _reshape_ndarray_to_vegm_format,
)
from ._dev_utils import replace_kwargs
from ._error_checking import (
    _validate_huc_list,
    _validate_grid,
    _validate_latlon_list,
    _validate_dir,
    _validate_grid_bounds,
    _validate_date,
    _validate_mask,
)


CONUS_DX = 1000
CONUS_DY = 1000
CONUS1_DZ = 100
CONUS1_Z_TOP = 500
CONUS2_DZ = 200
CONUS2_Z_TOP = 2000
CONUS_Z_BOTTOM = 0
_HOURS_PER_FORCING_FILE = 24


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
    """Subset static input files from national datasets in HydroData.

    The subset values will be written as ParFlow binary files (pfbs) in
    write_dir. By default the following variables will be subset.
        - Slope in the east/west direction (slope_x)
        - Slope in the north/south direction (slope_y)
        - Subsurface units indicator file (pf_indicator)
        - Mannings roughness coefficients (mannings)
        - Depth to bedrock (pf_flowbarrier)
        - Long term average preciptation minus evaporation (i.e. recharge) (pme)
        - Steady state pressure head used to initialize transient simulations
          (ss_pressure_head)

    Note that some datasets might not contain all 7 static input variables. In
    that case, the subset_static function is going to get and save the data for
    the variables supported by the dataset and print out a message for those
    that are not.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as
            i,j index values where 0,0 is the lower left hand corner of a domain.
            ij_bounds are given relative to whatever grid is being used for the
            subset.
        dataset (str): static inputs dataset name from the HydroData catalog e.g.
            "conus1_domain"
        write_dir (str): directory where the subset files will be written
        var_list (tuple[str]): tuple of variables to subset from the dataset.
            By default all 7 variables above will be subset. The user can specify
            as subset of these variables or list additional variables that are
            available in their dataset of choice.

    Returns:
        A dictionary mapping the static variable names to the corresponding file
        paths where the subset data were written.

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
        try:
            subset_data = hf_hydrodata.get_gridded_data(
                dataset=dataset,
                variable=var,
                file_type="pfb",
                temporal_resolution="static",
                grid_bounds=ij_bounds,
            )
        except ValueError as err:
            print(f"Variable '{var}' not found in dataset '{dataset}':", err)
        else:
            file_path = os.path.join(write_dir, f"{var}.pfb")
            write_pfb(file_path, subset_data, dist=False)
            file_paths[var] = file_path
            print(f"Wrote {var}.pfb in specified directory.")

    return file_paths


def subset_press_init(ij_bounds, dataset, date, write_dir, time_zone="UTC"):
    """Subset a pressure file from a national dataset in HydroData.

    This function will select the pressure file for midnight on the date provided
    and subset the selected pressure file to the ij_bounds provided. The subset
    data will be written out as a ParFlow binary file (pfb) to be used as an
    initial pressure file for a ParFlow simulation.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as
            i,j index values where 0,0 is the lower left hand corner of a domain.
            ij_bounds are given relative to whatever grid is being used for the
            subset.
        dataset (str): dataset name from the HydroData catalog that the pressure
            file will be subset from e.g. "conus1_baseline_mod"
        date (str): The date of the pressure file that you would like to subset,
            in the form 'yyyy-mm-dd'
        write_dir (str): directory where the subset file will be written
        time_zone (str): timezone information for subset date. Data will be
            subset at midnight in the specified timezone. Defaults to "UTC".

    Returns:
        The filepath of the subset file, which includes datetime information, so
        that it can be used by later functions (e.g. edit_runscript_for_subset).

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

    new_date = get_utc_time(date, time_zone)
    print(f"UTC Date: {new_date}")
    date_string = new_date.strftime("%Y.%m.%d:%H.%M.%S_UTC0")

    try:
        subset_data = hf_hydrodata.get_gridded_data(
            dataset=dataset,
            variable="pressure_head",
            file_type="pfb",
            temporal_resolution="hourly",
            grid_bounds=ij_bounds,
            start_time=new_date,
        )
    except ValueError as err:
        print(f"No pressure file found in dataset '{dataset}':", err)
        return None

    file_path = os.path.join(write_dir, f"{dataset}_{date_string}_press.pfb")
    write_pfb(file_path, subset_data[0, :, :, :], dist=False)
    print(f"Wrote {file_path} in specified directory.")
    return file_path


def config_clm(ij_bounds, start, end, dataset, write_dir, time_zone="UTC"):
    """Modify template CLM driver files for a desired subdomain and run duration.

    This function will obtain template clm driver files (specifically vegm, vep
    and drv_clmin) from the existing national simulations on HydroData and modify
    them to reflect the desired subdomain (indicated by the ij_bounds) and run
    duration (indicated by the start and end dates). The modified files will be
    written out to a user specified directory. These files are required if you
    are going to run a ParFlow-CLM simulation.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as
            i,j index values where 0,0 is the lower left hand corner of a domain.
            ij_bounds are given relative to whatever grid is being used for the
            subset.
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        dataset (str): the dataset that the files should be obtained from name
            e.g. "conus1_baseline_mod"
        write_dir (str): directory where the subset files will be written
        timezone (str): timezone information for start and end dates.
            Defaults to "UTC".

    Returns:
        A dictionary mapping the CLM file types ("vegp", "vegm", "drv_clm") to
        the corresponging filepaths where the CLM files were written.

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

    # get the pfb version of the vegm file
    file_type_list = ["vegp", "pfb", "drv_clm"]
    file_paths = {}
    for file_type in file_type_list:
        if file_type == "vegp":
            file_path = os.path.join(write_dir, "drv_vegp.dat")
            try:
                hf_hydrodata.get_raw_file(
                    file_path,
                    dataset=dataset,
                    file_type=file_type,
                    variable="clm_run",
                    temporal_resolution="static",
                )
            except ValueError as err:
                print(f"Failed to get {file_type} file for dataset '{dataset}':", err)
            else:
                file_paths[file_type] = file_path
                print("copied vegp")
        elif file_type == "pfb":
            try:
                subset_data = hf_hydrodata.get_gridded_data(
                    dataset=dataset,
                    file_type=file_type,
                    variable="clm_run",
                    temporal_resolution="static",
                    grid_bounds=ij_bounds,
                )
            except ValueError as err:
                print(f"Failed to get vegm file for dataset '{dataset}':", err)
            else:
                land_cover_data = _reshape_ndarray_to_vegm_format(subset_data)
                file_path = write_land_cover(land_cover_data, write_dir)
                file_paths[file_type] = file_path
                print("subset vegm")
        elif file_type == "drv_clm":
            file_path = os.path.join(write_dir, "drv_clmin.dat")
            try:
                hf_hydrodata.get_raw_file(
                    file_path,
                    dataset=dataset,
                    file_type=file_type,
                    variable="clm_run",
                    temporal_resolution="static",
                )
            except ValueError as err:
                print(f"Failed to get {file_type} file for dataset '{dataset}':", err)
            else:
                print("copied drv_clmin")
                edit_drvclmin(
                    file_path=file_path, start=start, end=end, time_zone=time_zone
                )
                file_paths[file_type] = file_path
                print("edited drv_clmin")
    return file_paths


def subset_forcing(
    ij_bounds,
    grid,
    start,
    end,
    dataset,
    write_dir,
    time_zone="UTC",
    forcing_vars=(
        "precipitation",
        "downward_shortwave",
        "downward_longwave",
        "specific_humidity",
        "air_temp",
        "atmospheric_pressure",
        "east_windspeed",
        "north_windspeed",
    ),
):
    """Subset forcing files from national datasets in HydroData.

    Subset forcing data will be written out as pfb files formatted for a ParFlow
    run with 24 hours per forcing file. Per ParFlow-CLM convention separate files
    will be written for each variable following the standard clm variable naming
    convention.

    Forcing file outputs will be numbered starting with 0000 and data will start
    at midnight local time for the timezone that has been provided. If no
    timezone is provided it will default to midnight UTC.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as
            i,j index values where 0,0 is the lower left hand corner of a domain.
            ij_bounds are given relative to whatever grid is being used for the
            subset.
        grid (str): The spatial grid that the ij indices are calculated relative
            to and that the subset data will be returned on. Possible values:
            "conus1" or "conus2"
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        dataset (str): forcing dataset name from the HydroData catalog that the
            forcing files will be subset from e.g. "NLDAS2".
        write_dir (str): directory where the subset files will be written
        time_zone (str): timezone information for start and end dates. Data will
            be subset starting at midnight in the specified timezone. Defaults to
            "UTC".
        forcing_vars (tuple[str]): tuple of forcing variables to subset. By
            default all 8 variables needed to run ParFlow-CLM will be subset.

    Returns:
        A dictionary mapping the forcing variable names to the corresponding file
        paths where the subset data were written.

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
    start_date = get_utc_time(start, time_zone)
    end_date = get_utc_time(end, time_zone)
    exit_event = threading.Event()
    threads = [
        threading.Thread(
            target=_subset_forcing_variable,
            args=(
                variable,
                ij_bounds,
                grid,
                start_date,
                end_date,
                dataset,
                write_dir,
                time_zone,
                outputs,
                exit_event,
            ),
        )
        for variable in forcing_vars
    ]

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


def _subset_forcing_variable(
    variable,
    ij_bounds,
    grid,
    start_date,
    end_date,
    dataset,
    write_dir,
    time_zone,
    outputs,
    exit_event,
):
    """Helper for subset_forcing that subsets data for one forcing variable."""

    day = 1
    date = start_date
    delta = timedelta(hours=_HOURS_PER_FORCING_FILE)
    outputs[variable] = []
    print(f"Reading {variable} pfb sequence")

    while date < end_date and not exit_event.is_set():
        start_time = date
        end_time = date + delta
        # we need to distinguish between UTC and non-UTC as the datacatalog
        # returns the wrong answer for requests that start reading from the
        # middle of a file and span multiple files
        try:
            if time_zone == "UTC":
                subset_data = hf_hydrodata.get_gridded_data(
                    dataset=dataset,
                    variable=variable,
                    grid=grid,
                    file_type="pfb",
                    temporal_resolution="hourly",
                    start_time=start_time,
                    end_time=end_time,
                    grid_bounds=ij_bounds,
                )
            else:
                next_day_midnight = datetime(
                    end_time.year, end_time.month, end_time.day
                )
                data1 = hf_hydrodata.get_gridded_data(
                    dataset=dataset,
                    variable=variable,
                    grid=grid,
                    file_type="pfb",
                    temporal_resolution="hourly",
                    start_time=start_time,
                    end_time=next_day_midnight,
                    grid_bounds=ij_bounds,
                )
                data2 = hf_hydrodata.get_gridded_data(
                    dataset=dataset,
                    variable=variable,
                    grid=grid,
                    file_type="pfb",
                    temporal_resolution="hourly",
                    start_time=next_day_midnight,
                    end_time=end_time,
                    grid_bounds=ij_bounds,
                )
                subset_data = np.concatenate((data1, data2), axis=0)
        except Exception as exc:
            raise ValueError(
                f"Failed to get {variable} data from {start_time} to {end_time}."
            ) from exc

        paths = hf_hydrodata.get_paths(
            dataset=dataset,
            variable=variable,
            grid=grid,
            file_type="pfb",
            temporal_resolution="hourly",
            start_time=start_time,
            end_time=end_time,
        )

        write_path = os.path.join(
            write_dir, _adjust_filename_hours(os.path.basename(paths[0]), day)
        )
        outputs[variable].append(write_path)
        write_pfb(write_path, subset_data[:, :, :], dist=False)
        date = date + delta
        day = day + 1
    if not exit_event.is_set():
        print(f"Finished writing {variable} to folder")


def _adjust_filename_hours(filename, day):
    """Adjust the part of the the forcing filename representing hours.

    The aim is to match what a parflow simulation expects on each day of the
    simulation. The first day of the simulation the hours will be
    “.000001_to_000024.”, the second day the hours will be “.000025_to_000048.”
    and so on. This is in case the first day of simulation does not coincide with
    the first day of the water year (Oct 1st), as the dataset filenames assume
    day 1 is Oct 1st. The input and output filenames must match the regular
    expression “*.*.[0-9]{6}_to_[0-9]{6}.*”

    Parameters:
        filename (str): original forcing filename
        day (int): day relative to the start date of forcing file subsetting

    Returns:
        The new forcing filename with adjusted hours.
    """
    if day < 1:
        raise ValueError("day must be >= 1")
    s1, s2, s3, s4 = filename.split(".")
    if not s1 or not s2 or not s3:
        raise ValueError("invalid forcing filename")
    pattern = re.compile("[0-9]{6}_to_[0-9]{6}")
    if pattern.fullmatch(s3) is None:
        raise ValueError("invalid forcing filename")
    start = str(_HOURS_PER_FORCING_FILE * (day - 1) + 1).rjust(6, "0")
    end = str(_HOURS_PER_FORCING_FILE * day).rjust(6, "0")
    s3 = start + "_to_" + end
    if pattern.fullmatch(s3) is None:
        raise ValueError("invalid forcing filename")
    return ".".join([s1, s2, s3, s4])


def edit_runscript_for_subset(
    ij_bounds, runscript_path, write_dir=None, runname=None, forcing_dir=None
):
    """Modify a ParFlow run script for a new subdomain run.

    This function is designed to start from a national ParFlow runscript template
    and perform the following three modifications.
        1. Modify the geometry to reflect the bounds of the desired ij_bounds
           (i.e. the number of grid cells in the x and y direction and the upper
           bounds of the geometry)
        2. Update the runname to for the desired new run.
        3. Update the location of the climate forcings for the new run.

    If the runname is None and write_dir is the directory containing the
    runscript file, the runscript file will be overwritten.

    Args:
        ij_bounds (tuple[int]): bounding box for subset. This should be given as
            i,j index values where 0,0 is the lower left hand corner of a domain.
            ij_bounds are given relative to whatever grid is being used for the
            subset.
        runscript_path (str): absolute path to the template parflow runscript
            file
        write_dir (str): directory where the new template file will be written.
            If it is None, defaults to the directory containing the runscript.
        runname (str): name for the new parflow run. If it is None, defaults to
            the runscript's previous runname.
        forcing_dir (str): path to the directory containing the subset forcing
            files. If it is None, defaults to the runscript's previous forcing
            directory path.

    Returns:
        A path to the new runscript file that will be created.

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

    run = Run.from_definition(runscript_path)
    _, file_extension = os.path.splitext(runscript_path)
    if runname is not None:
        run.set_name(runname)
        print(
            f"New runname: {runname} provided, a new {file_extension[1:]} "
            "file will be created"
        )
    else:
        print(
            f"No runname provided, old {file_extension[1:]} file will be " "overwritten"
        )

    if forcing_dir is not None:
        print(
            f"Climate forcing directory has been changed to {forcing_dir} "
            " in runscript."
        )
        run.Solver.CLM.MetFilePath = forcing_dir
    else:
        print("No forcing directory provided, run.Solver.CLM.MetFilePath key not set")

    imin, jmin, imax, jmax = ij_bounds
    ni, nj = imax - imin, jmax - jmin
    run.ComputationalGrid.NY = int(nj)
    run.ComputationalGrid.NX = int(ni)
    print(f"ComputationalGrid.NY set to {nj} and NX to {ni}")

    domain_type = run.GeomInput.domaininput.InputType
    if domain_type == "SolidFile":
        print(
            "GeomInput.domaininput.InputType detected as SolidFile, no "
            "additional keys to change for subset"
        )
    else:
        run.Geom.domain.Upper.X = ni * CONUS_DX
        run.Geom.domain.Upper.Y = nj * CONUS_DY
        print(
            f"""GeomInput.domaininput.InputType detected as Box, updating 
                Geom.domain.Upper.X to {run.Geom.domain.Upper.X} and Y to 
                {run.Geom.domain.Upper.Y} to match subset.
             """
        )

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
    """Change the filenames of input files in the ParFlow runscript.

    This function will update the paths to input files in a ParFlow runscript.
    The provided arguments will reset the corresponding parflow keys to match the
    user specified paths to input files.  File names can be specified with or
    without relative or absolute file paths. If no path is provided ParFlow will
    expect the input files to be present in the run directory at the time of
    simulation.

    Note that this will only change paths for keys that already exist in the
    template ParFlow run script you are starting from and will not reconfigure a
    run to use new keys (for example if you are not starting from a run script
    that uses a solid file, adding a new solid file path will not configure the
    run to use a solid file).

    Refer to the ParFlow manual for additional information on any of the keys
    listed above.


    If the runname is None and write_dir is the directory containing the
    runscript file, the runscript file will be overwritten.

    Args:
        runscript_path (str): path to the runscript file (yaml or pfidb)
        write_dir (str): directory where the new template file will be written.
            If it is None, defaults to the directory containing the runscript
            file.
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
        A path to the new runscript file that will be created.

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
            f"New runname: {runname} provided, a new {file_extension[1:]} "
            "file will be created"
        )
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


@replace_kwargs({"P": "topo_p", "Q": "topo_q"})
def dist_run(topo_p, topo_q, runscript_path, working_dir=None, dist_clim_forcing=True):
    """Distribute ParFlow input files for parallel computing.

    This function will distribute input files to topo_p grids in the
    x direction and topo_q grids in the y direction. If dist_clim_forcing
    is true, forcing files will be distributed as well according to the
    same topology. If working_dir is different that the directory containing
    the runscript file, the edited runscipt file will be written to working_dir.

    Args:
        topo_p (int): number of grids (processes) to create in the x direction
        topo_q (int): number of grids (processes) to create in the y direction
        runscript_path (str): path to the runscript file (yaml or pfidb)
        working_dir (str): directory containing the files to be distributed.
            If it is None, it defaults to the directory containing the runscript
            file.
        dist_clim_forcing (bool): if true, distribute forcing files

    Returns:
        str: Path to the edited runscript file that will be created.

    Example:

    .. code-block:: python

        runscript_path = dist_run(
            topo_p=2,
            topo_q=2,
            runscript_path="/path/to/your/original/runscript",
            dist_clim_forcing=False
        )
    """
    if not isinstance(topo_p, int) or topo_p <= 0:
        raise TypeError("topo_p must be a positive integer.")
    if not isinstance(topo_q, int) or topo_q <= 0:
        raise TypeError("topo_q must be a positive integer.")
    if not os.path.isfile(runscript_path):
        raise FileNotFoundError("runscript_path must be a valid file path")
    if working_dir is None:
        working_dir = os.path.dirname(runscript_path)
    _validate_dir(working_dir)

    run = Run.from_definition(runscript_path)
    run.Process.Topology.P = topo_p
    run.Process.Topology.Q = topo_q

    if dist_clim_forcing:
        print("Distributing your climate forcing")
        forcing_dir = run.Solver.CLM.MetFilePath
        _validate_dir(forcing_dir)
        for file_name in os.listdir(forcing_dir):
            if file_name.endswith(".pfb") and ".out." not in file_name:
                file_path = os.path.join(forcing_dir, file_name)
                run.dist(file_path)
    else:
        print("No forcing directory provided, only distributing static inputs")

    max_nz = 0
    for file_name in os.listdir(working_dir):
        if file_name.endswith(".pfb") and ".out." not in file_name:
            file_path = os.path.join(working_dir, file_name)
            input_array = read_pfb(file_path)
            nz = input_array.shape[0]
            max_nz = max(max_nz, nz)
            run.dist(file_path)
            print(f"Distributed {os.path.basename(file_path)} with NZ {nz}")
    run.ComputationalGrid.NZ = max_nz

    _, file_extension = os.path.splitext(runscript_path)
    file_path, _ = run.write(
        working_directory=working_dir, file_format=f"{file_extension[1:]}"
    )
    return file_path
