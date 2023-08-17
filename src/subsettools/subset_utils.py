"""This module contains helper functions for the subsetting module.
"""

import os
import shutil
import re
from datetime import datetime

import numpy as np
from hydrodata.data_catalog import data_access
from parflow.tools.io import read_clm


def get_conus_hucs_indices(huc_list, grid):
    """Get the huc datafile as an ndarray and three mask arrays
    representing the selected hucs.

    Args:
        huc_list (list[str]): a list of huc IDs
        grid (str): "conus1" or "conus2"

    Returns:
        A tuple (conus_hucs, sel_hucs, indices_j, indices_i) where
        conus_hucs is an ndarray of the huc datafile, sel_hucs is
        a mask array for the selected hucs, and indices_i and
        indices_j mask arrays in the j and i directions.
    """
    huc_len = len(huc_list[0])
    huc_list = [int(huc) for huc in huc_list]
    entry = data_access.get_catalog_entry(
        dataset="huc_mapping", grid=grid.lower(), file_type="tiff"
    )
    conus_hucs = data_access.get_ndarray(entry, level=str(huc_len))
    sel_hucs = np.isin(conus_hucs, huc_list).squeeze()
    indices_j, indices_i = np.where(sel_hucs > 0)
    return conus_hucs, sel_hucs, indices_j, indices_i


def indices_to_ij(conus_hucs, indices_j, indices_i):
    """Get the conus ij-bounds for the conus_hucs boundary defined by
    indices_j and indices_i.

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
    arr_jmin = np.min(indices_j)
    arr_jmax = np.max(indices_j) + 1  # right bound inclusive

    jmin = conus_hucs.shape[1] - arr_jmax
    jmax = conus_hucs.shape[1] - arr_jmin

    return (imin, jmin, imax, jmax)


def subset_vegm(path, ij_bounds):
    """Docstring: TODO"""
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
    """Write the land cover file in vegm format

    Args:
        land_cover_data (ndarray): formatted vegm data (2d array)
        write_dir (str): path to output directory

    Returns:
        out_file (str):
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
    out_file = os.path.join(write_dir, "drv_vegm.dat")
    np.savetxt(
        fname=out_file,
        X=land_cover_data,
        delimiter=" ",
        comments="",
        header=header,
        fmt=["%d"] * 2 + ["%.6f"] * 2 + ["%.2f"] * 2 + ["%d"] * 19,
    )
    return out_file


def edit_drvclmin(
    read_path,
    write_dir,
    start=None,
    end=None,
    startcode=2,
    vegp_name="drv_vegp.dat",
    vegm_name="drv_vegm.dat",
):
    """Docstring: TODO"""
    write_path = os.path.join(write_dir, "drv_clmin.dat")
    shutil.copyfile(read_path, write_path)
    with open(write_path, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "vegtf" in line:
            lines[
                i
            ] = f"vegtf           {vegm_name}                         Vegetation Tile Specification File\n"
        elif "vegpf" in line:
            lines[
                i
            ] = f"vegpf           {vegp_name}                         Vegetation Type Parameter\n"
        elif "startcode" in line:
            lines[
                i
            ] = f"startcode       {startcode}                                    1=restart file, 2=defined\n"
        elif "clm_ic" in line:
            lines[
                i
            ] = f"clm_ic          {startcode}                                    1=restart file, 2=defined\n"

    if start is not None:
        startdt = datetime.strptime(start, "%Y-%m-%d")
        sd = startdt.strftime("%d")
        sm = startdt.strftime("%m")
        enddt = datetime.strptime(end, "%Y-%m-%d")
        ed = enddt.strftime("%d")
        em = enddt.strftime("%m")

        for i, line in enumerate(lines):
            if "sda" in line:
                lines[
                    i
                ] = f"sda            {sd}                                    Starting Day\n"
            elif "smo" in line:
                lines[
                    i
                ] = f"smo            {sm}                                    Starting Month\n"
            elif "syr" in line:
                lines[
                    i
                ] = f"syr            {startdt.year}                                  Starting Year\n"
            elif "eda" in line:
                lines[
                    i
                ] = f"eda            {ed}                                    Ending Day\n"
            elif "emo" in line:
                lines[
                    i
                ] = f"emo            {em}                                    Ending Month\n"
            elif "eyr" in line:
                lines[
                    i
                ] = f"eyr            {enddt.year}                                  Ending Year\n"

    with open(write_path, "w") as f:
        f.writelines(lines)
    return write_path


def adjust_filename_hours(filename, day):
    """Adjust the forcing filename hours so that they match what a parflow simulation expects
    on each day of the simulation.

    The first day of the simulation the hours will be "*.000001_to_000024.*", the second day
    the hours will be "*.000025_to_000048.*" and so on. This is in case the first day of simulation
    does not coincide with the first day of the water year (Oct 1st), as the dataset filenames
    assume day 1 is Oct 1st. The input and output filenames must match the regular expression
    "*.*.[0-9]{6}_to_[0-9]{6}.*"

    Args:
        filename (str): original forcing filename
        day (int): day relative to the start date of forcing file subsetting

    Returns:
        The forcing filename with adjusted hours.

    Raises:
        AssertionError: If the input or output filename string do not match the above regex.
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
