"""Functions to subset gridded files from national datasets in HydroData.

The following functions can be used to subset gridded input files to set up a
ParFlow simulation.
    - subset static model inputs
    - subset meteorological forcings
    - subset initial pressure data
    - subset gridded CLM inputs (vegm)
"""
#pylint: disable=C0301,R0913,R0914,R0917
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
            subset at midnight in the specified timezone. This should be a 
            zoneinfo-supported time zone. Defaults to "UTC".

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
    dataset_version=None,
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
            be subset starting at midnight in the specified timezone. This
            should be a zoneinfo-supported time zone. Defaults to "UTC".
        forcing_vars (tuple[str]): tuple of forcing variables to subset. By
            default all 8 variables needed to run ParFlow-CLM will be subset.
        dataset_version (str): version of the forcing dataset. By default the
            latest version of a dataset will be returned.

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
            forcing_vars=("precipitation", "air_temp"),
            dataset_version="0.9",
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
                dataset_version,
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
    dataset_version,
    outputs,
    exit_event,
    lock,
):
    """
        Helper to read forcing data of one variable and write it daily files in write_dir folder.
    """

    (base_file_path, hf_filter_options) = _get_forcing_file_basename(
        dataset, variable, grid, dataset_version, ij_bounds)

    # Allocate a numpy array to buffer 24 hours (1 day) of data to be written to each file
    block_day_np = np.full(
        (24, ij_bounds[3] - ij_bounds[1], ij_bounds[2] - ij_bounds[0]), np.nan
    )

    # compute timezone hours offset between UTC time and requested timezone time
    start_date_midnight = datetime(start_date.year, start_date.month, start_date.day)
    timezone_offset = (start_date - start_date_midnight).total_seconds() / 3600

    # Get optimal hours to read at time from hf_hydrodata based on grid size
    block_hours_per_read = _get_read_block_size(ij_bounds)

    # Read data in blocks of data using hf_hydrodata
    start_block_date = start_date_midnight
    day = 1
    hour = 0
    skipped_offset_in_day_1 = False
    write_paths = []

    while start_block_date < end_date and not exit_event.is_set():
        (subset_data, end_block_date) = _read_hf_hydrodata_block(
            hf_filter_options, start_block_date, end_date, block_hours_per_read)

        (block_day_np, hour, day, skipped_offset_in_day_1) = _write_block_to_files(
            subset_data, block_day_np, hour, day, skipped_offset_in_day_1, ij_bounds, timezone_offset, base_file_path, write_dir, write_paths)

        # increment start_block_date (in hours) for next hf_hydrodata block read
        start_block_date = end_block_date

    if not exit_event.is_set():
        with lock:
            outputs[variable] = write_paths
        print(f"Finished writing {variable} to folder")

def _read_hf_hydrodata_block(hf_filter_options, start_block_date, end_date, block_hours_per_read):
    """
        Read a block of data from hydrodata to support subset_forcing function.

        Parameters:
            hf_filter_options:      The hf_hydrodata filter options to be used for the read.
            start_block_date:       The value to put put in start_time options of hf_filter options.
            end_date:               The final end date for the subset_forcing function.
            block_hours_per_read:   The number of hours to read to be used to set end_time option in read.

        Returns:
            A tuple (subset_data, end_block_date)

        The subset data is numpy array containing dimensions (z, y, x) where z is the hours of read data.
        The end_block_date is the start date of the next call to this function.

        The date range read may be a range of 24 hour UTC periods or
        may be range of hours less than 24 hours. This depends on the ij_bound size.
    """

    block_delta_per_read = timedelta(hours=block_hours_per_read)
    if block_hours_per_read >= 24:
        # The block is more than a day so read the block_hours_per_read hours
        end_block_date = min(start_block_date + block_delta_per_read, end_date)
    else:
        # If the block is less than 24 hours then the start/end must be in the same UTC day
        end_block_date = min(start_block_date + block_delta_per_read, end_date)
        start_block_date_midnight = datetime(
            start_block_date.year, start_block_date.month, start_block_date.day
        )
        end_block_date_midnight = datetime(
            end_block_date.year, end_block_date.month, end_block_date.day
        )
        if start_block_date_midnight != end_block_date_midnight:
            # The end block date must in the same UTC day as the start date
            end_block_date = end_block_date_midnight

    # Read the block using get_gridded_data
    hf_filter_options["start_time"] = start_block_date
    hf_filter_options["end_time"] = end_block_date
    subset_data = get_hf_gridded_data(hf_filter_options)

    return (subset_data, end_block_date)


def _write_block_to_files(subset_data, block_day_np, hour, day, skipped_offset, ij_bounds, timezone_offset, base_filename, write_dir, write_paths):
    """
        Process the subset_data of one block and write the data to daily forcing files.

        Parameters:
            subset_data:    The numpy array of a block returned by get_gridded_data.
            block_day_np:   A numpy buffer the size of 24 hours if ij_bounds to be written to file.
            hour:           The next hour to be added to the block_day_np buffer (number from 0-23).
            day:            The next day to be used to for the next daily file to be written.
            skipped_offset: True if we already skipped the first hours of the block for timezone reads.
            ij_bounds:      The grid bounds of the requested data [min_x, min_y, max_x, max_y].
            timezone_offset:The number of hours of the timezone from UTC time (is 0 for UTC).
            base_filename:  The base file name of the daily forcing files to be written.
            write_dir:      The directory to write the daily forcing files.
            write_paths:    An array of full file paths of files written to write_dir.

        Returns: block_day_np, hour, day, skipped_offset)

        The returned block_day_np may the the same as the one passed in and filled with new hours
        or it may be a newly allocated buffer if the buffer was just written to a file.

        If the block_day_np buffer is filled to 24 hours then the buffer is written to a daily
        forcing file and the hours set back to 0 to be filled by the next hours from subset_data.
        The subset_data block may contain many days of files so this may write many daily files
        or if the bounds is large this function may write no files, but only partiall fill the buffer.

        The hour, day and skipped offset are returned as updated values after processing the
        hours data in the subset_data and updating the block_day_np.

    """
    if timezone_offset == 0:
        # There is no timezone offset, but subset_data might be 1 hour or 24 hours or > 24 hours
        # Loop over all hours in the hf_hydrodata result (z is the hour dimension of the result)
        for z in range(0, subset_data.shape[0]):
            block_day_np[hour, :, :] = subset_data[z, :, :]
            hour = hour + 1
            if hour == 24:
                _write_day_to_file(
                    block_day_np, base_filename, day, write_dir, write_paths
                )
                block_day_np = np.full(
                    (24, ij_bounds[3] - ij_bounds[1], ij_bounds[2] - ij_bounds[0]),
                    np.nan,
                )
                day = day + 1
                hour = 0
    else:
        # There is a timezone offset
        # Loop over all hours in the hf_hydrodata result (z is the hour dimension of the result)
        for z in range(0, subset_data.shape[0]):
            if day == 1 and hour < timezone_offset and not skipped_offset:
                # These are the first hours before timezone offset of the first day
                # Skip these rows, until we read all the hours before timezone offset
                hour = hour + 1
                if hour == timezone_offset:
                    skipped_offset = True
                    hour = 0
            elif hour < _HOURS_PER_FORCING_FILE - timezone_offset:
                # These are the first hours of an output day, this does not fill day
                block_day_np[hour, :, :] = subset_data[z, :, :]
                hour = hour + 1
            elif hour < 24:
                # These are hours from the next 24 hour UTC block to fill up the output day
                block_day_np[hour, :, :] = subset_data[z, :, :]
                hour = hour + 1
                if hour == 24:
                    _write_day_to_file(
                        block_day_np, base_filename, day, write_dir, write_paths
                    )
                    block_day_np = np.full(
                        (
                            24,
                            ij_bounds[3] - ij_bounds[1],
                            ij_bounds[2] - ij_bounds[0],
                        ),
                        np.nan,
                    )
                    day = day + 1
                    hour = 0
            else:
                # This should not happen because hour should be zero after writting 24 hour day
                raise ValueError("Too many hours in a day for file.")
    return (block_day_np, hour, day, skipped_offset)


def _get_forcing_file_basename(dataset, variable, grid, dataset_version, ij_bounds):
    """
        Get the base file name to be used to write daily forcing files by subset_forcing function.

        Parameters:
            dataset:        The hf_hydrodata dataset version of the forcing file.
            variable:       The hf_hydrodata variable name of the focing file.
            grid:           The grid of the requested forcing data
            dataset_version:The hf_hydrodata version of the dataset of the forcing file.
            ij_bounds:      The grid_bounds of the request to put into the return hf_filter_options.
        Returns: A tuple (base_filename, hf_filter_options).

        The base_filename is the base file name of the forcing file stored in /hydrodata.
        The hf_filter_options are the options to be passed to get_gridded_data to get forcing data.
    """
        # Get base path name of output forcing file using hf_hydrodata path
    hf_filter_options = {
        "dataset": dataset,
        "variable": variable,
        "grid": grid,
        "file_type": "pfb",
        "grid_bounds": ij_bounds,
        "mask": "false",
        "temporal_resolution": "hourly",
        "dataset_version": dataset_version
    }

    path = hf_hydrodata.get_paths(hf_filter_options)[0]
    base_filename = os.path.basename(path)
    return (base_filename, hf_filter_options)


def _get_read_block_size(ij_bounds):
    """
        Get the number of hours of data to be read by the blocking reads to hf_hydrodata to
        optimize the performance of reads for the subset_forcing function.

        Parameters:
            ij_bounds:      This is the grid bounds of the request [min_x, min_y, max_x, max_y].

        Returns:
            The number of hours of data to return each each block call for subset_forcing function.

        This may return many days of data for a small subgrid or only a few hours of data
        for a large subgrid that is too large to return 24 hours of data in a single hf_hydrodata call.
    """
    float64_byte_size = 8
    max_chunking_memory_bytes = 1500000000

    # Compute the time delta between hf_hydrodata chunking block reads of data to optimize performance
    # The block of time we can read depends on the ij_bounds of the requested subset of data

    # First try to find the number of days of data we can read using the max bytes per hf_hydrodata read
    subgrid_byte_size = (
        abs(ij_bounds[0] - ij_bounds[2])
        * abs(ij_bounds[1] - ij_bounds[3])
        * float64_byte_size
    )
    block_days_per_read = int(
        max_chunking_memory_bytes / (subgrid_byte_size * _HOURS_PER_FORCING_FILE)
    )
    # Limit number of days per read to 1 year
    block_days_per_read = min(block_days_per_read, 366)

    # Get number of hours to read per block
    block_hours_per_read = block_days_per_read * _HOURS_PER_FORCING_FILE

    if block_days_per_read == 0:
        # ij_bounds is so large we cannot download 24 hours in one get_gridded_data call
        # So use the most hours we can at a time instead of blocks of 24 hours in one call
        block_hours_per_read = int(max_chunking_memory_bytes / (subgrid_byte_size))

    return block_hours_per_read


def _write_day_to_file(data_np, base_filename, day, write_dir, write_paths):
    """
        Write a file containing 1 day of hourly data to a file.

        Parameters:
            data_np:        A number array of dimension (z, y, x) where z is 24 hours.
            base_filename:  The base name of the daily file to be written.
            write_dir:      The directory to write the daily file.
            write_paths:    An array of full file paths written. The path is appended to this.
    """
    write_path = os.path.join(write_dir, _adjust_filename_hours(base_filename, day))
    write_paths.append(write_path)
    write_pfb(write_path, data_np, dist=False)


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
