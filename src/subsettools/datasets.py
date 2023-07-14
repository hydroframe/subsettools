from importlib import resources

def get_ref_yaml_path(grid, mode, input_file_type):    
    if mode == "transient":
        mode = "_pfclm_" + mode + "_"
    else:
        mode = "_pf_" + mode + "_"
    filename = grid + mode + input_file_type + ".yaml"
    with resources.path("subsettools.ref_yamls", filename) as f:
        data_file_path = f
    return data_file_path
