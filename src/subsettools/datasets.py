"""This module contains functions for accessing the template runscripts.
"""

from importlib import resources
import os
import shutil

def get_template_runscript(grid, mode, input_file_type, write_dir):
    """Get a ParFlow template runscript based on grid, mode and input file type and write it to write_dir.

    Args:
        grid (str): "conus1" or "conus2"
        mode (str): "spinup" or "transient"
        input_file_type (str): "box" or "solid"
        write_dir (str): directory where the template runscript file will be copied

    Returns:
        str: Path to the template runscript.
    """
    if not isinstance(grid, str):
        raise TypeError("grid must be a string")
    if not isinstance(mode, str):
        raise TypeError("mode must be a string")
    if not isinstance(input_file_type, str):
        raise TypeError("input_file_type must be a string")
    if grid not in ["conus1", "conus2"]:
        raise ValueError("Supported grids are 'conus1' and 'conus2'")
    if mode not in ["transient", "spinup"]:
        raise ValueError("Supported modes are 'transient' and 'spinup'")
    if input_file_type not in ["box", "solid"]:
        raise ValueError("Supported input file types are 'box' and 'solid'")
    if not os.path.isdir(write_dir):
        raise FileNotFoundError("write_dir must be a valid directory") 

    if mode == "transient":
        mode = "_pfclm_" + mode + "_"
    elif mode == "spinup":
        mode = "_pf_" + mode + "_"

    filename = grid + mode + input_file_type + ".yaml"
    with resources.path("subsettools.ref_yamls", filename) as f:
        shutil.copy(f, write_dir)
    return os.path.join(write_dir, filename)
