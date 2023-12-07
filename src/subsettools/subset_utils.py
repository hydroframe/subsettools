"""This module contains helper functions for the subsetting module.
"""

import os
import shutil
from datetime import datetime, timedelta
import pytz
import numpy as np
from parflow.tools.io import read_clm


def subset_vegm(path, ij_bounds):
    """Read in vegm file and subset it according to ij_bounds.

    Subset a national vegm file to a smaller domain based on ij_bounds that are provided relative to the grid that the national file is on. 

    Args:
        path (str): path to read vegm file
        ij_bounds (Tuple[int]): bounding box for subset. This should be given as i,j index values where 0,0 is the lower left hand 
            corner of a domain. ij_bounds are given relative to whatever grid is being used for the subset. Use the latlon_to_ij function to 
            determine ij indices from lat long values. 

    Returns:
        ndarray:
            Subset vegm data.
    """
    vegm = read_clm(path, type="vegm")  # returns (j,i,k)
    vegm = np.transpose(vegm, (2, 0, 1))  # transpose to k,j,i

    imin, jmin, imax, jmax = ij_bounds
    vegm = vegm[:, jmin:jmax, imin:imax]  # slicing on k,j,i

    _, nj, ni = vegm.shape
    indices = np.indices((nj, ni)) + 1
    indices = indices[::-1, :, :]
    vegm = np.vstack([indices, vegm])  # stack x,y indices on vegm

    # transpose and reshape back into expected 2D vegm file format for the subset
    vegm = vegm.transpose(1, 2, 0).reshape(-1, 25)

    return vegm


def write_land_cover(land_cover_data, write_dir):
    """Write the land cover ndarray in vegm format.

    Read in a gridded landcover dataset and write it out as a vegm file which is correctly formatted for CLM. 

    Args:
        land_cover_data (ndarray): formatted vegm data (2d array)
        write_dir (str): path to output directory

    Returns:
        str:
            path to output vegm file.
    """
    heading = (
        "x y lat lon sand clay color fractional coverage of grid, by vegetation class (Must/Should Add to "
        "1.0) "
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


def edit_drvclmin(
    file_path,
    start=None,
    end=None,
    time_zone="UTC",
    startcode=2,
    vegp_name="drv_vegp.dat",
    vegm_name="drv_vegm.dat",
):
    """Edit a template CLM driver for a new simulation. 

    Update the start and end dates, timezone, restart type and vegm and vegp input file names for a new simulation. 

    Args:
        file_path (str): clm driver file path
        start (str): start date (inclusive), in the form 'yyyy-mm-dd'
        end (str): end date (exlusive), in the form 'yyyy-mm-dd'
        time_zone (str): time_zone used to calculate start/end dates. Defaults to "UTC".
        startcode (int): startcode for the parflow simulation
        vegp_name (str): vegp filename
        vegm_name (str): vegm filename
    
    Raises:
        AssertionError: If one of start or end is None and the other is not.
    """
    assert (start is None and end is None) or (start is not None and end is not None)
    
    with open(file_path, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "vegtf" in line:
            lines[i] = f"{'vegtf':<15}{vegm_name:<37} Vegetation Tile Specification File\n"
        elif "vegpf" in line:
            lines[i] = f"{'vegpf':<15}{vegp_name:<37} Vegetation Type Parameter\n"
        elif "startcode" in line:
            lines[i] = f"{'startcode':<15}{startcode:<37} 1=restart file, 2=defined\n"
        elif "clm_ic" in line:
            lines[i] = f"{'clm_ic':<15}{startcode:<37} 1=restart file, 2=defined\n"

    if start is not None:
        start_date = get_UTC_time(start, time_zone)
        end_date = get_UTC_time(end, time_zone) - timedelta(hours=1)
            
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

    with open(file_path, "w") as f:
        f.writelines(lines)


def get_UTC_time(date_string, time_zone):
    """Convert the given date and time_zone to UTC time. 

    Args:
        date_string (str): date in the form 'yyyy-mm-dd'
        time_zone (str):

    Returns:
        A timezone-unaware datetime object representing the time in UTC.
    """    
    date = datetime.strptime(date_string, "%Y-%m-%d")
    if time_zone != "UTC":
        date = date.replace(tzinfo=pytz.timezone(time_zone)).astimezone(pytz.UTC).replace(tzinfo=None)
    return date
