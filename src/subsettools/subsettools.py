import os
import shutil
import pathlib
import subprocess
from datetime import datetime, timedelta

import numpy as np
import pytz
import re
from hydrodata.national_mapping.map_wgs84 import ConusMap
from hydrodata.data_catalog import data_access
from parflow import Run
from parflow.tools.io import read_clm, read_pfb, write_pfb


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


def huc_to_ij(huc_list, grid):
    """Get the conus ij-bounds of the area defined by the the huc IDs
       in huc_list in the conus grid. 
       
       All huc IDs in huc_list must be the same length (hucs of the same
       level). Supported grids are "conus1" and "conus2".                                                                                                                                                        
    Args:                                                                                                                                             
        huc_list (list[str]): a list of huc IDs                                                                                                       
        grid (str): "conus1" or "conus2"                                                                                                               
    Returns:
        A tuple of the form (imin, jmin, imax, jmax) representing the bounds
        in the conus grid of the area defined by the huc IDs in huc_list.

    Raises:                                                                                                                                                   AssertionError: If all huc IDs are not the same length.                                                                                           """
    huc_len = len(huc_list[0])
    assert all(
        [len(huc) == huc_len for huc in huc_list]
    ), "All huc IDs should have the same length!"
    conus_hucs, _, indices_j, indices_i = get_conus_hucs_indices(huc_list, grid)
    return indices_to_ij(conus_hucs, indices_j, indices_i)


def latlon_to_ij(latlon_bounds, grid):
    """Get the conus ij-bounds of the area defined by the the latitute/longitude
       bounds (latlon_bounds) in the conus grid. 
       
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
    conus_map = ConusMap(grid.lower())
    point0 = conus_map.map_to_grid(
        latlon_bounds[0][1], latlon_bounds[0][0]
    )
    point1 = conus_map.map_to_grid(
        latlon_bounds[1][1], latlon_bounds[1][0]
    )
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
    """Docstring: TODO
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
        capture_output=True,
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
        entry = data_access.get_catalog_entry(
            dataset=dataset, file_type="pfb", period="static", variable=var
        )
        if entry is not None:
            subset_data = data_access.get_ndarray(entry, grid_bounds=ij_bounds)
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

    entry = data_access.get_catalog_entry(
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

    subset_data = data_access.get_ndarray(
        entry, grid_bounds=ij_bounds, start_time=new_date
    )

    out_file = f"{dataset}_{date_string}_press.pfb"
    write_pfb(os.path.join(write_dir, out_file), subset_data, dist=False)
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
        entry = data_access.get_catalog_entries(
            dataset=dataset, file_type=file_type, variable="clm_run", period="static"
        )[0]
        file_path = data_access.get_file_paths(entry)[0]
        print(file_path)
        if file_type == "vegp":
            shutil.copyfile(file_path, os.path.join(write_dir, "drv_vegp.dat"))
            print("copied vegp")
        elif file_type == "vegm":
            land_cover_data = subset_vegm(file_path, ij_bounds)
            write_land_cover(land_cover_data, write_dir)
            print("subset vegm")
        elif file_type == "drv_clm":
            edit_drvclmin(
                read_path=file_path,
                write_dir=write_dir,
                start=start,
                end=end,
            )
            print("edited drv_clmin")


def subset_vegm(path, ij_bounds):
    """ Docstring: TODO
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
    """ Write the land cover file in vegm format

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
    """Docstring: TODO
    """
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
        entries = data_access.get_catalog_entry(
            dataset=dataset, variable=var, grid=grid, file_type="pfb", period="hourly"
        )

        day = 1
        start_date = datetime.strptime(start, "%Y-%m-%d")
        end_date = datetime.strptime(end, "%Y-%m-%d")
        delta = timedelta(days=1)
        outputs[var] = []
        print(f"Reading {var} pfb sequence")

        while start_date < end_date:
            subset_data = data_access.get_ndarray(
                entries,
                start_time=start_date,
                end_time=start_date + delta,
                grid_bounds=ij_bounds,
            )
            paths = data_access.get_file_paths(
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
    assert s1 != '' and s2 != '' and s4 != '', "invalid forcing filename"
    pattern = re.compile("[0-9]{6}_to_[0-9]{6}")
    assert pattern.fullmatch(s3) is not None, "invalid forcing filename"
    
    start = str(24 * (day - 1) + 1).rjust(6, "0")
    end = str(24 * day).rjust(6, "0")
    s3 = start + "_to_" + end
    assert pattern.fullmatch(s3) is not None, "invalid adjusted forcing filename"    
    return ".".join([s1, s2, s3, s4])


def edit_runscript_for_subset(
    ij_bounds, runscript_path, write_dir, runname=None, forcing_dir=None
):
    """Set up a new parflow run from a template file.

    This function will return the correct runscript file to do your run based on the grid, 
    if you're doing spin-up and if you're using a solid file for input. The runname and 
    forcing directory keys will be reset for your new run. If the runname is None and
    write_dir is the directory containing runscript_path, the original template file will 
    be overwritten.

    Args:
        ij_bounds (Tuple[int]): bounding box for subset
        runscript_path (str): path to the original template file
        write_dir (str): directory where the new template file will be written
        runname (str): name for the new parflow run
        forcing_dir (str): path to the directory containing the subset forcing files 

    Returns:
        Path to the new runscript file that will be created. 

    Raises:
        AssertionError: If write_dir is not a valid directory or runscript_path is not
        a valid file path.
    """
    assert os.path.isdir(write_dir), "write_dir must be a director"
    assert os.path.isfile(runscript_path), "runscript_path must be a valid file path"

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


def copy_static_files(read_dir, write_dir):
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
    write_dir,
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
    runscript_path, the original template file will be overwritten.

    Args:
        runscript_path (str): path to the runscript file (yaml or pfidb)
        write_dir (str): directory where the new template file will be written
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
        AssertionError: If write_dir is not a valid directory or runscript_path is not
        a valid file path.
    """    
    assert os.path.isdir(write_dir), "write_dir must be a directory"
    assert os.path.isfile(
        runscript_path
    ), "runscript_path must be a valid path to an existing file"

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


def dist_run(P, Q, runscript_path, write_dir, dist_clim_forcing=True):
    """Distribute parflow input files.

    This function will distribute static input files to P grids in the 
    x direction and Q grids in the y direction. If dist_clim_forcing
    is true, forcing files will be distributed as well according to the
    same topology.

    Args:
        P (int): number of grids (processes) to create in the x direction
        Q (int): number of grids (processes) to create in the y direction
        runscript_path (str): path to the runscript file (yaml or pfidb)
        write_dir (str): directory where the new template file will be written
        dist_clim_forcing (bool): if true, distribute forcing files

    Raises:
        AssertionError: If write_dir is not a valid directory or runscript_path is not
        a valid file path.
    """    
    assert os.path.isdir(write_dir), "write_dir must be a directory"
    assert os.path.isfile(runscript_path), "runscript_path must be a valid file path"

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

    static_input_paths = pathlib.Path(write_dir).glob("*.pfb")
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
    run.write(working_directory=write_dir, file_format=f"{file_extension[1:]}")


def restart_run(runscript_path):
    """Docstring: TODO
    """
    pf_dir = os.path.dirname(runscript_path)
    file_name, file_extension = os.path.splitext(runscript_path)
    file_extension = file_extension[1:]
    run = Run.from_definition(runscript_path)
    runname = run.get_name()

    dump_interval = run.TimingInfo.DumpInterval

    path_to_tcl = pf_dir + "/clm_restart.tcl"
    lines = open(path_to_tcl, "r").readlines()[0]
    istep = [int(i) for i in lines.split() if i.isdigit()][0]
    print(istep)

    press_base = runname + ".out.press.{:05d}.pfb"
    print(press_base)
    restart_press = press_base.format(istep - 2)
    run.Geom.domain.ICPressure.FileName = restart_press
    print(f"Initial Press filename changed to: {run.Geom.domain.ICPressure.FileName}")

    # update pf timing (needs to be done no matter how you're running pf)
    run.TimingInfo.StartCount = istep - 1
    run.TimingInfo.StartTime = (istep - 1) * dump_interval

    print(f"start time is now: {run.TimingInfo.StartTime}")
    print(f"start count is now: {run.TimingInfo.StartCount}")

    if run.Solver.LSM == "CLM":
        edit_drvclmin(read_path=f"{pf_dir}/drv_clmin.dat", startcode=1)

        print(f"Overwrote drv_clmin.dat (changed startcode to 1 (restart file))")

    run.write(working_directory=pf_dir, file_format=f"{file_extension}")
