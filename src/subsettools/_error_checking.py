"""This module contains helper functions for error-checking user-provided arguments.
"""

import os 
import re

def _validate_huc_list(huc_list):
    """Check errors for a user-provided list of HUC IDs."""
    if not isinstance(huc_list, list):
        raise TypeError("huc_list must be a list.")
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
    """Check errors for a user-provided grid argument."""
    if not isinstance(grid, str):
        raise TypeError("grid must be a string")
    grid = grid.lower()
    if grid not in ["conus1", "conus2"]:
        raise ValueError("Supported grids are 'conus1' and 'conus2'")


def _validate_dir(dir_name):
    """Check that dir_name is a valid directory."""
    if not os.path.isdir(dir_name):
        raise FileNotFoundError(f"{dir_name} is not a valid existing directory") 


def _validate_grid_bounds(grid_bounds):
    if not isinstance(grid_bounds, tuple):
        raise TypeError("bounds must be a tuple.")
    if not all(isinstance(val, (int, float)) for val in grid_bounds):
        raise TypeError("bounds must contain integers or floats.")
    if not len(grid_bounds) == 4:
        raise ValueError("bounds must contain four elements (imin, imax, jmin, jmax).")

    
def _validate_date(date):
    if not isinstance(date, str):
        raise TypeError("date must be a string.")
    pattern = re.compile(r'\d{4}-\d{2}-\d{2}$')
    if not bool(pattern.match(date)):
        raise ValueError("invalid date format.")
