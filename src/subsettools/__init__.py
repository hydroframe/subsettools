# read version from installed package
from importlib.metadata import version

__version__ = version("subsettools")

from subsettools.domain import (
    define_huc_domain,
    define_latlon_domain,
    define_upstream_domain,
    write_mask_solid,
    huc_to_ij,
    latlon_to_ij,
    create_mask_solid,
)
from subsettools.subsetting import (
    subset_static,
    subset_press_init,
    subset_forcing,
)
from subsettools.subsettools import (
    config_clm,
    edit_runscript_for_subset,
    copy_files,
    change_filename_values,
    dist_run,
)
from subsettools.datasets import (
    get_template_runscript,
)
