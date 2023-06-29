from importlib import resources

def get_ref_yaml_path(grid):
    if grid == "conus1":
        filename = "C1_pfclm_transient_solid.yaml"
    elif grid == "conus2":
        filename = ""
        
    with resources.path("subsettools.ref_yamls", filename) as f:
        data_file_path = f
    return data_file_path
