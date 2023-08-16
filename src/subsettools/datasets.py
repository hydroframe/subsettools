from importlib import resources


def get_ref_yaml_path(grid, mode, input_file_type):
    """Get the correct template (yaml) runscript path based on grid, mode
    and input file type. 
       
    Args:
        grid (str): "conus1" or "conus2"
        mode (str): "spinup" or "transient"
        input_file_type (str): "box" or "solid"

    Returns:
        pathlib.Path: Path to the template runscript.
    """
    assert grid in ["conus1", "conus2"], "invalid grid provided"
    assert mode in ["transient", "spinup"], "invalid mode"
    assert input_file_type in ["box", "solid"], "invalid input file type"

    if mode == "transient":
        mode = "_pfclm_" + mode + "_"
    elif mode == "spinup":
        mode = "_pf_" + mode + "_"

    filename = grid + mode + input_file_type + ".yaml"
    with resources.path("subsettools.ref_yamls", filename) as f:
        data_file_path = f
    return data_file_path
