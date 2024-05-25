"""Functions to subset gridded files from national datasets in HydroData.

The following functions can be used to subset gridded input files to set up a
ParFlow simulation.
    - subset static model inputs
    - subset meteorological forcings
    - subset initial pressure data
    - subset gridded CLM inputs (vegm)
"""

import os
from datetime import datetime, timedelta
import threading
import re
import warnings
import numpy as np
import hf_hydrodata
from parflow.tools.io import write_pfb
from ._common import (
    get_utc_time,
    get_hf_gridded_data,
)
from ._error_checking import (
    _validate_grid,
    _validate_dir,
    _validate_grid_bounds,
    _validate_date,
)


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
        - Long term average precipitation minus evaporation (i.e. recharge) (pme)
        - Steady state pressure head used to initialize transient simulations
          (ss_pressure_head)

    Note that some datasets might not contain all 7 static input variables. In
    that case, the subset_static function is going to raise a ValueError for any
    variables that do not exist in the dataset. The default variable list
    contains the necessary static variables for the CONUS2 grid. For CONUS1-based
    datasets, "mannings" and "pf_flowbarrier" should be removed from the list.

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
            a subset of these variables or list additional variables that are
            available in their dataset of choice.

    Returns:
        A dictionary mapping the static variable names to the corresponding file
        paths where the subset data were written.

    Example:

    .. code-block:: python

        # Subsetting static variables for a CONUS1 workflow
        # We need to remove "pf_flowbarrier" and "mannings" from the list
        filepaths = subset_static(
            ij_bounds=(375, 239, 487, 329),
            dataset="conus1_domain",
            write_dir="/path/to/your/chosen/directory",
            var_list=("slope_x", "slope_y", "pf_indicator", "pme",
                      "ss_pressure_head")
        )

        # Subsetting static variables for a CONUS2 workflow
        # Note that we can use the default var_list here
        filepaths = subset_static(
            ij_bounds=(3701, 1544, 3792, 1633),
            dataset="conus2_domain",
            write_dir="/path/to/your/chosen/directory",
        )
    """
    warnings.warn(
        "Note that for subsettools versions >= 2.0.0, this function will raise "
        "a ValueError if a variable in var_list is not supported in the "
        "dataset. (In older versions, it just printed an error message and "
        "continued executing normally). You can check in the HydroData "
        "documentation which variables are contained in each dataset "
        "(https://hf-hydrodata.readthedocs.io/en/latest/available_data.html).",
        DeprecationWarning,
        stacklevel=2,
    )
    _validate_grid_bounds(ij_bounds)
    if not isinstance(dataset, str):
        raise TypeError("dataset name must be a string.")
    _validate_dir(write_dir)
    if not all(isinstance(var, str) for var in var_list):
        raise TypeError("All variable names should be strings.")
    file_paths = {}
    options = {
        "dataset": dataset,
        "file_type": "pfb",
        "temporal_resolution": "static",
        "grid_bounds": ij_bounds,
    }
    for var in var_list:
        options["variable"] = var
        subset_data = get_hf_gridded_data(options)
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

    options = {
        "dataset": dataset,
        "variable": "pressure_head",
        "file_type": "pfb",
        "temporal_resolution": "hourly",
        "grid_bounds": ij_bounds,
        "start_time": new_date,
    }
    subset_data = get_hf_gridded_data(options)
    file_path = os.path.join(write_dir, f"{dataset}_{date_string}_press.pfb")
    write_pfb(file_path, subset_data[0, :, :, :], dist=False)
    print(f"Wrote {file_path} in specified directory.")
    return file_path


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
    forcing_vars = tuple(set(forcing_vars))
    outputs = {}
    start_date = get_utc_time(start, time_zone)
    end_date = get_utc_time(end, time_zone)
    exit_event = threading.Event()
    lock = threading.Lock()
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
                lock,
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

    if exit_event.is_set():
        raise ValueError("One or more threads were interrupted.")
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
    lock,
):
    """Helper for subset_forcing that subsets data for one forcing variable."""

    day = 1
    date = start_date
    delta = timedelta(hours=_HOURS_PER_FORCING_FILE)
    write_paths = []
    print(f"Reading {variable} pfb sequence")
    options = {
        "dataset": dataset,
        "variable": variable,
        "grid": grid,
        "file_type": "pfb",
        "temporal_resolution": "hourly",
        "grid_bounds": ij_bounds,
    }

    while date < end_date and not exit_event.is_set():
        start_time = date
        end_time = date + delta
        # we need to distinguish between UTC and non-UTC as the datacatalog
        # returns the wrong answer for requests that start reading from the
        # middle of a file and span multiple files
        try:
            if time_zone == "UTC":
                options["start_time"] = start_time
                options["end_time"] = end_time
                subset_data = get_hf_gridded_data(options)
            else:
                next_day_midnight = datetime(
                    end_time.year, end_time.month, end_time.day
                )
                options["start_time"] = start_time
                options["end_time"] = next_day_midnight
                data_day1 = get_hf_gridded_data(options)
                options["start_time"] = next_day_midnight
                options["end_time"] = end_time
                data_day2 = get_hf_gridded_data(options)
                subset_data = np.concatenate((data_day1, data_day2), axis=0)
        except Exception:
            exit_event.set()
            raise
        paths = hf_hydrodata.get_paths(options)
        write_path = os.path.join(
            write_dir, _adjust_filename_hours(os.path.basename(paths[0]), day)
        )
        write_paths.append(write_path)
        write_pfb(write_path, subset_data[:, :, :], dist=False)
        date = date + delta
        day = day + 1
    if not exit_event.is_set():
        with lock:
            outputs[variable] = write_paths
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
