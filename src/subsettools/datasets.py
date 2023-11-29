"""This module contains functions for accessing the template runscripts.
"""

from importlib import resources
import os
import shutil
from ._error_checking import (
    _validate_grid,
    _validate_dir,
)

def get_template_runscript(grid, mode, input_file_type, write_dir):
    """Get a ParFlow template runscript based on grid, mode and input file type and write it to write_dir.

    Args:
        grid (str): The spatial grid that the ij indices are calculated relative to and that the subset data 
            will be returned on. Possible values: “conus1” or “conus2”
        mode (str): The type of simulation you would like to do. Possible values: "spinup" (run ParFlow with a
            constant recharge forcing at the upper boundary) and "transient" (coupled ParFlow-CLM run)
        input_file_type (str): The type of domain you will run. Possible values: "box" or "solid"
        write_dir (str): directory where the template runscript file will be copied

    Returns:
        str: Path to the template runscript.

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

    if mode == "transient":
        mode = "_pfclm_" + mode + "_"
    elif mode == "spinup":
        mode = "_pf_" + mode + "_"

    filename = grid + mode + input_file_type + ".yaml"
    with resources.path("subsettools.ref_yamls", filename) as f:
        shutil.copy(f, write_dir)
    return os.path.join(write_dir, filename)
