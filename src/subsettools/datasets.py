"""This module contains functions for accessing the template runscripts.
"""

from importlib import resources
import os
import shutil

def get_ref_yaml_path(grid, mode, input_file_type, write_dir):
    """Get a ParFlow template runscript based on grid, mode and input file type.

    Args:
        grid (str): "conus1" or "conus2"
        mode (str): "spinup" or "transient"
        input_file_type (str): "box" or "solid"
        write_dir (str): directory where the template runscript file will be copied

    Returns:
        str: Path to the template runscript.
    """
    assert grid in ["conus1", "conus2"], "invalid grid provided"
    assert mode in ["transient", "spinup"], "invalid mode"
    assert input_file_type in ["box", "solid"], "invalid input file type"
    assert os.path.isdir(write_dir), "write_dir must be a directory"

    if mode == "transient":
        mode = "_pfclm_" + mode + "_"
    elif mode == "spinup":
        mode = "_pf_" + mode + "_"

    filename = grid + mode + input_file_type + ".yaml"
    with resources.path("subsettools.ref_yamls", filename) as f:
        shutil.copy(f, write_dir)
    return os.path.join(write_dir, filename)
