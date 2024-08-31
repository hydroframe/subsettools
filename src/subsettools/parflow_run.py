"""Functions to get and customize a ParFlow runscript."""

import os
import re
import shutil
from importlib import resources
from datetime import datetime, timedelta
import netCDF4 as nc
from parflow import Run
from parflow.tools.io import read_pfb, write_pfb
from ._dev_utils import replace_kwargs
from ._error_checking import (
    _validate_grid,
    _validate_dir,
    _validate_grid_bounds,
)
from ._constants import (
    CONUS_DX,
    CONUS_DY,
)


def get_template_runscript(grid, mode, input_file_type, write_dir):
    """Get a ParFlow template runscript.

    The runscript is selected based on the grid, mode and input file type and
    is copied to write_dir.

    Args:
        grid (str): The spatial grid that the ij indices are calculated relative
            to and that the subset data will be returned on. Possible values:
            “conus1” or “conus2”
        mode (str): The type of simulation you would like to do. Possible values:
            "spinup" (run ParFlow with a constant recharge forcing at the upper
            boundary) and "transient" (coupled ParFlow-CLM run)
        input_file_type (str): The type of domain you will run. Possible values:
            "box" or "solid"
        write_dir (str): directory where the template runscript file will be
            copied

    Returns:
        A path to the template runscript.

    Example:

    .. code-block:: python

        runscript_path = get_template_runscript(
            grid="conus1",
            mode="spinup",
            input_file_type="solid",
            write_dir="/path/to/your/chosen/directory"
        )
    """
    _validate_grid(grid)
    _validate_dir(write_dir)
    if not isinstance(mode, str):
        raise TypeError("mode must be a string")
    if not isinstance(input_file_type, str):
        raise TypeError("input_file_type must be a string")
    if mode not in ["transient", "spinup"]:
        raise ValueError("Supported modes are 'transient' and 'spinup'")
    if input_file_type not in ["box", "solid"]:
        raise ValueError("Supported input file types are 'box' and 'solid'")

    filename = "_".join([grid, mode, input_file_type]) + ".yaml"
    with resources.as_file(
        resources.files("subsettools.template_runscripts").joinpath(filename)
    ) as f:
        shutil.copy(f, write_dir)
    return os.path.join(write_dir, filename)


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


def restart_run(
    runscript_path,
    stop_time,
    new_dir=None,
    runname=None,
    forcing_dir=None,
    output_type="netcdf",
):
    run = Run.from_definition(runscript_path)
    if runname is not None:
        run.set_name(runname)
    if forcing_dir is not None:
        run.Solver.CLM.MetFilePath = forcing_dir

    restart_timestep = _get_restart_timestep(runscript_path, output_type)
    _set_timing_parameters(run, new_dir, restart_timestep, stop_time)
    
    if new_dir is not None:
        os.mkdir(new_dir)
        working_directory = new_dir
        _copy_static_inputs(runscript_path, new_dir)
        if run.Solver.LSM == 'CLM':
            _copy_clm_restart_files(runscript_path, new_dir)
            _copy_clm_driver_files(runscript_path, new_dir)
    else:
        working_directory = os.path.dirname(runscript_path)
        if run.Solver.LSM == 'CLM':
            _rename_clm_restart_files(runscript_path, restart_timestep)

    init_press_data = _get_ic_pressure_from_old_run(
        runscript_path, output_type, restart_timestep
    )
    filename = "initial_pressure.pfb"
    write_pfb(os.path.join(working_directory, filename), init_press_data)
    run.Geom.domain.ICPressure.FileName = filename

    if run.Solver.LSM == 'CLM':
        drvclm_path = os.path.join(working_directory, "drv_clmin.dat")
        _update_drvclm_for_restart(drvclm_path, restart_timestep, stop_time)

    runscript_path, _ = run.write(
        file_format="yaml",
        working_directory=working_directory,
    )
    return runscript_path


def _set_timing_parameters(run, new_dir, restart_timestep, stop_time):
    if new_dir is not None:
        run.TimingInfo.StartCount = 0
    else:
        run.TimingInfo.StartCount = restart_timestep
    run.TimingInfo.StartTime = run.TimingInfo.StartCount
    if run.Solver.LSM == "CLM":
        run.Solver.CLM.IstepStart = run.TimingInfo.StartCount + 1
    if stop_time < run.TimingInfo.StartTime:
        raise ValueError("Simulation stop time is smaller than start time.")
    run.TimingInfo.StopTime = stop_time


def _copy_static_inputs(runscript_path, new_dir):
    filenames = _extract_filenames_from_runscript(runscript_path)
    old_dir = os.path.dirname(runscript_path)
    for filename in filenames:
        shutil.copy(os.path.join(old_dir, filename), new_dir)
    # Remove old ic pressure file from new directory
    run = Run.from_definition(runscript_path)
    ic_pressure = run.Geom.domain.ICPressure.FileName
    os.remove(os.path.join(new_dir, ic_pressure))


def _copy_clm_restart_files(runscript_path, new_dir):
    old_dir = os.path.dirname(runscript_path)
    pattern = re.compile(r"^clm\.rst\.00000\.\d+$")
    restart_files = [f for f in os.listdir(old_dir) if pattern.match(f)]
    if not restart_files:
        raise ValueError(f"No clm restart files found in {old_dir}.")
    for file in restart_files:
        shutil.copy(os.path.join(old_dir, file), new_dir)


def _copy_clm_driver_files(runscript_path, new_dir):
    filenames = ["drv_clmin.dat", "drv_vegm.dat", "drv_vegp.dat"]
    old_dir = os.path.dirname(runscript_path)
    for filename in filenames:
        shutil.copy(os.path.join(old_dir, filename), new_dir)


def _rename_clm_restart_files(runscript_path, restart_timestep):
    run = Run.from_definition(runscript_path)
    working_directory = os.path.dirname(runscript_path)
    num_procs = run.Process.Topology.P * run.Process.Topology.Q * run.Process.Topology.R
    restart_timestep = str(restart_timestep).rjust(5, "0")
    for i in range(num_procs):
        old_file = os.path.join(working_directory, f"clm.rst.00000.{i}")
        if not os.path.exists(old_file):
            raise ValueError(f"Clm restart file {old_file} not found.")
        new_file = os.path.join(working_directory, f"clm.rst.{restart_timestep}.{i}")
        os.rename(old_file, new_file)


def _extract_filenames_from_runscript(runscript_path):
    filenames = []
    with open(runscript_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("FileName: "):
                filenames.append(line[len("FileName: ") :])
    return filenames


def _get_restart_timestep(runscript_path, output_type):
    run = Run.from_definition(runscript_path)
    working_directory = os.path.dirname(runscript_path)
    if run.Solver.LSM == "CLM":
        # Get restart timestep from clm_restart.tcl file. If no such file
        # exists, set timestep to 0. (This is a bug with ParFlow, need to
        # fix this.)
        restart_file = os.path.join(working_directory, "clm_restart.tcl")
        if os.path.exists(restart_file):
            with open(restart_file) as f:
                first_line = f.readline()
            restart_timestep = int(first_line.split()[2])
        else:
            restart_timestep = 0
    else:
        if output_type == "pfb":
            files = _get_pfb_output_files(runscript_path)
            restart_timestep = len(files) - 1
        elif output_type == "netcdf":
            files = _get_netcdf_output_files(runscript_path)
            restart_timestep = 0
            for file in files:
                dataset = nc.Dataset(os.path.join(working_directory, file))
                timesteps_in_file = len(dataset.dimensions["time"])
                restart_timestep += timesteps_in_file
                dataset.close()
            restart_timestep = restart_timestep - 1
        else:
            raise ValueError("Invalid output_type provided.")
    return restart_timestep


def _get_ic_pressure_from_old_run(runscript_path, output_type, restart_timestep):
    working_directory = os.path.dirname(runscript_path)
    if output_type == "pfb":
        files = _get_pfb_output_files(runscript_path)
        return read_pfb(os.path.join(working_directory, files[restart_timestep]))
    elif output_type == "netcdf":
        files = _get_netcdf_output_files(runscript_path)
        timesteps = 0
        for file in files:
            dataset = nc.Dataset(os.path.join(working_directory, file))
            timesteps_in_file = len(dataset.dimensions["time"])
            timesteps = timesteps + timesteps_in_file
            if timesteps > restart_timestep:
                index = restart_timestep - (timesteps - timesteps_in_file)
                press_data = dataset.variables["pressure"][index, ...]
                dataset.close()
                return press_data
            dataset.close()
        raise ValueError("Invalid restart timestep.")
    else:
        raise ValueError("Invalid output_type provided.")


def _get_pfb_output_files(runscript_path):
    run = Run.from_definition(runscript_path)
    runname = run.get_name()
    pattern = re.compile(rf"^{runname}\.out\.press\.(\d{{5}})\.pfb$")
    working_directory = os.path.dirname(runscript_path)
    files = [f for f in os.listdir(working_directory) if pattern.match(f)]
    if not files:
        raise ValueError(f"No output pfb pressure files found in {working_directory}.")
    files.sort(key=lambda x: int(pattern.search(x).group(1)))
    return files


def _get_netcdf_output_files(runscript_path):
    run = Run.from_definition(runscript_path)
    runname = run.get_name()
    pattern = re.compile(rf"^{runname}\.out\.(\d{{5}})\.nc$")
    working_directory = os.path.dirname(runscript_path)
    files = [f for f in os.listdir(working_directory) if pattern.match(f)]
    if not files:
        raise ValueError(f"No output netcdf files found in {working_directory}.")
    files.sort(key=lambda x: int(pattern.search(x).group(1)))
    return files


def _update_drvclm_for_restart(file_path, restart_timestep, stop_time):
    _set_drvclm_startcode(file_path, 1)
    start_date = _read_drvclm_date(file_path, start=True)
    new_start_date = start_date + timedelta(hours=restart_timestep)
    new_end_date = start_date + timedelta(hours=stop_time)
    _write_drvclm_date(file_path, new_start_date, start=True)
    _write_drvclm_date(file_path, new_end_date, start=False)
    return file_path


def _set_drvclm_startcode(file_path, startcode):
    with open(file_path, "r") as file:
        content = file.read()

    startcode_pattern = r"(startcode\s+)(\d{1})"
    ic_pattern = r"(clm_ic\s+)(\d{1})"
    startcode_match = re.search(startcode_pattern, content)
    ic_match = re.search(ic_pattern, content)

    if not (startcode_match and ic_match):
        raise ValueError("Failed to find startcode fields.")

    content = re.sub(startcode_pattern, rf"\g<1>{startcode}", content)
    content = re.sub(ic_pattern, rf"\g<1>{startcode}", content)
    with open(file_path, "w") as file:
        file.write(content)


def _read_drvclm_date(file_path, start):
    with open(file_path, "r") as file:
        content = file.read()

    year_pattern, month_pattern, day_pattern, hour_pattern = _get_drvclm_patterns(start)
    year_match = re.search(year_pattern, content)
    month_match = re.search(month_pattern, content)
    day_match = re.search(day_pattern, content)
    hour_match = re.search(hour_pattern, content)

    if not (year_match and month_match and day_match and hour_match):
        raise ValueError("Failed to find the date fields in the file.")

    year = int(year_match.group(2))
    month = int(month_match.group(2))
    day = int(day_match.group(2))
    hour = int(hour_match.group(2))

    date = datetime(year, month, day, hour)
    return date


def _write_drvclm_date(file_path, date, start):
    with open(file_path, "r") as file:
        content = file.read()

    year = date.year
    month = date.month
    day = date.day
    hour = date.hour

    year_pattern, month_pattern, day_pattern, hour_pattern = _get_drvclm_patterns(start)
    content = re.sub(year_pattern, rf"\g<1>{year}", content)
    content = re.sub(month_pattern, rf"\g<1>{month}", content)
    content = re.sub(day_pattern, rf"\g<1>{day}", content)
    content = re.sub(hour_pattern, rf"\g<1>{hour}", content)

    with open(file_path, "w") as file:
        file.write(content)


def _get_drvclm_patterns(start):
    if start:
        year_pattern = r"(syr\s+)(\d{4})"
        month_pattern = r"(smo\s+)(\d{1,2})"
        day_pattern = r"(sda\s+)(\d{1,2})"
        hour_pattern = r"(shr\s+)(\d{1,2})"
    else:
        year_pattern = r"(eyr\s+)(\d{4})"
        month_pattern = r"(emo\s+)(\d{1,2})"
        day_pattern = r"(eda\s+)(\d{1,2})"
        hour_pattern = r"(ehr\s+)(\d{1,2})"
    return year_pattern, month_pattern, day_pattern, hour_pattern
