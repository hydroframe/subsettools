from importlib import resources


def get_ref_yaml_path(grid, mode, input_file_type):
    assert grid in ["conus1", "conus2"], "invalid grid provided"
    assert mode in [
        "transient",
        "spinup",
    ], "invalid mode (valid modes are spinup/transient)"
    assert input_file_type in ["box", "solid"], "invalid input file type"

    if mode == "transient":
        mode = "_pfclm_" + mode + "_"
    elif mode == "spinup":
        mode = "_pf_" + mode + "_"

    filename = grid + mode + input_file_type + ".yaml"
    with resources.path("subsettools.ref_yamls", filename) as f:
        data_file_path = f
    return data_file_path
