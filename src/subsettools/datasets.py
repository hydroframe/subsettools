from importlib import resources

def get_ref_yaml_path(grid):
    if grid == "conus1":
        filename = "C1_pfclm_transient_solid.yaml"
    elif grid == "conus2":
        filename = ""

    return resources.files("subsettools.ref_yamls").joinpath(filename)
