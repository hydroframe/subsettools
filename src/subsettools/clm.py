"""Functions to configure CLM for coupled ParFlow-CLM simulations."""

import os
from datetime import timedelta
import numpy as np
import hf_hydrodata
from ._common import (
    get_utc_time,
)
from ._error_checking import (
    _validate_dir,
    _validate_grid_bounds,
    _validate_date,
)


_VEGM_COLUMNS = 25


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
                file_path = _write_land_cover(land_cover_data, write_dir)
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
                _edit_drvclmin(
                    file_path=file_path, start=start, end=end, time_zone=time_zone
                )
                file_paths[file_type] = file_path
                print("edited drv_clmin")
    return file_paths


def _reshape_ndarray_to_vegm_format(data):
    """Reshape ndarray returned by datacatalog to vegm format.

    Args:
        data (ndarray): raw subset vegm data (2d array)

    Returns:
        Ndarray reshaped to vegm format.
    """
    _, nj, ni = data.shape
    indices = np.indices((nj, ni)) + 1
    indices = indices[::-1, :, :]
    data = np.vstack([indices, data])  # stack x,y indices on vegm
    # transpose and reshape back into expected 2D vegm file format for the subset
    return data.transpose(1, 2, 0).reshape(-1, _VEGM_COLUMNS)


def _write_land_cover(land_cover_data, write_dir):
    """Write the land cover ndarray in vegm format.

    Read in a gridded landcover dataset and write it out as a vegm file which is
    correctly formatted for CLM.

    Args:
        land_cover_data (ndarray): formatted vegm data (2d array)
        write_dir (str): path to output directory

    Returns:
        str:
            path to output vegm file.
    """
    heading = (
        "x y lat lon sand clay color fractional coverage of grid, by "
        "vegetation class (Must/Should Add to 1.0) "
    )
    vegm_col_names = (
        "",
        "",
        "(Deg)",
        "(Deg)",
        "(%/100)",
        "",
        "index",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
    )
    header = "\n".join([heading, " ".join(vegm_col_names)])
    file_path = os.path.join(write_dir, "drv_vegm.dat")
    np.savetxt(
        fname=file_path,
        X=land_cover_data,
        delimiter=" ",
        comments="",
        header=header,
        fmt=["%d"] * 2 + ["%.6f"] * 2 + ["%.2f"] * 2 + ["%d"] * 19,
    )
    return file_path


def _edit_drvclmin(
    file_path,
    start=None,
    end=None,
    time_zone="UTC",
    startcode=2,
    vegp_name="drv_vegp.dat",
    vegm_name="drv_vegm.dat",
):
    """Edit a template CLM driver for a new simulation.

    Update the start and end dates, timezone, restart type and vegm and vegp
    input file names for a new simulation.

    Args:
        file_path (str): clm driver file path
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        time_zone (str): time_zone used to calculate start/end dates. Defaults
            to "UTC".
        startcode (int): startcode for the parflow simulation
        vegp_name (str): vegp filename
        vegm_name (str): vegm filename

    Raises:
        AssertionError: If one of start or end is None and the other is not.
    """
    assert (start is None and end is None) or (start is not None and end is not None)

    with open(file_path, encoding="utf-8") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "vegtf" in line:
            lines[i] = (
                f"{'vegtf':<15}{vegm_name:<37} Vegetation Tile Specification File\n"
            )
        elif "vegpf" in line:
            lines[i] = f"{'vegpf':<15}{vegp_name:<37} Vegetation Type Parameter\n"
        elif "startcode" in line:
            lines[i] = (
                f"{'startcode':<15}{startcode:<37} 1=restart file, 2=defined\n"
            )
        elif "clm_ic" in line:
            lines[i] = f"{'clm_ic':<15}{startcode:<37} 1=restart file, 2=defined\n"

    if start is not None:
        start_date = get_utc_time(start, time_zone)
        end_date = get_utc_time(end, time_zone) - timedelta(hours=1)

        for i, line in enumerate(lines):
            if "shr" in line:
                lines[i] = f"{'shr':<15}{start_date.hour:<37} Starting Hour\n"
            if "sda" in line:
                lines[i] = f"{'sda':<15}{start_date.day:<37} Starting Day\n"
            elif "smo" in line:
                lines[i] = f"{'smo':<15}{start_date.month:<37} Starting Month\n"
            elif "syr" in line:
                lines[i] = f"{'syr':<15}{start_date.year:<37} Starting Year\n"
            elif "ehr" in line:
                lines[i] = f"{'ehr':<15}{end_date.hour:<37} Ending Hour\n"
            elif "eda" in line:
                lines[i] = f"{'eda':<15}{end_date.day:<37} Ending Day\n"
            elif "emo" in line:
                lines[i] = f"{'emo':<15}{end_date.month:<37} Ending Month\n"
            elif "eyr" in line:
                lines[i] = f"{'eyr':<15}{end_date.year:<37} Ending Year\n"

    with open(file_path, "w", encoding="utf-8") as f:
        f.writelines(lines)
