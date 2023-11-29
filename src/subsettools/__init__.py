# read version from installed package
from importlib.metadata import version

__version__ = version("subsettools")

from subsettools.subsettools import (
    huc_to_ij,
    latlon_to_ij,
    create_mask_solid,
    subset_static,
    subset_press_init,
    config_clm,
    subset_forcing,
    edit_runscript_for_subset,
    copy_files,
    change_filename_values,
    dist_run,
)
from subsettools.datasets import (
    get_template_runscript,
)
