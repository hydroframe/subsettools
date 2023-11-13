"""This module contains helper functions for error-checking user-provided arguments.
"""

import os 

def _validate_huc_list(huc_list):
    """Error checking for a user-provided list of HUC IDs for the functions huc_to_ij() and create_mask_solid()."""
    if not all(isinstance(huc, str) for huc in huc_list):
        raise TypeError("All elements of huc_list must be strings")
    if not all(huc.isdigit() for huc in huc_list):
        raise ValueError("HUC IDs must contain only digits")
    huc_len = len(huc_list[0])
    if huc_len not in [2, 4, 6, 8, 10]:
        raise ValueError("HUC IDs are 2, 4, 6, 8, or 10-digit")
    if not all([len(huc) == huc_len for huc in huc_list]):
        raise ValueError("All HUC IDs should have the same length")    
    

def _validate_grid(grid):
    """Error checking for a user-provided grid argument."""
    if not isinstance(grid, str):
        raise TypeError("grid must be a string")
    grid = grid.lower()
    if grid not in ["conus1", "conus2"]:
        raise ValueError("Supported grids are 'conus1' and 'conus2'")


def _validate_dir(dir_name):
    if not os.path.isdir(dir_name):
        raise FileNotFoundError(f"{dir_name} is not a valid existing directory") 
