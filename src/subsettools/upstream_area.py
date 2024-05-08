"""Functions for subsetting the upstream area of a point in the conus grid."""

import numpy as np
import hf_hydrodata
from .subsettools import _indices_to_ij
from ._error_checking import (
    _validate_grid,
    _validate_latlon_list,
)


def upstream_area_to_ij(outlets, grid):
    """Define a domain that is the upstream area of the points in outlets.

    The domain is defined by the grid ij bounds of the bounding box that
    encompasses the upstream area of all the points in outlets and a mask and a
    mask for that bounding box indicating which cells are part of the selected
    area.

    The flow_direction files that are used to define the upstream area follow the
    convention: down: 1, left: 2, up: 3, right: 4.

    Args:
        outlets (List[List[float]]): list of lat-lon points of the form
            [[lat1, lon1], [lat2, lon2], ...]
        grid (str): The spatial grid that the ij indices are calculated relative
            to and that the subset data will be returned on. Possible values:
            “conus1” or “conus2”

    Returns:
        A tuple (bounds, mask).

        Bounds is a tuple of the form (imin, jmin, imax, jmax) representing the
        bounds in the conus grid of the upstream area of the outlets. imin, jmin,
        imax, jmax are the west, south, east and north sides of the box
        respectively and all i,j indices are calculated relative to the lower
        southwest corner of the domain.

        Mask is a 2D numpy.ndarray that indicates which cells inside the bounding
        box are part of the computed upstream area of the outlets.

    Raises:
        ValueError: If the computed upstream area of the outlets is empty.

    Example:
    .. code-block:: python
        bounds, mask = upstream_area_to_ij(
            outlets=[[44.1348, -95.5084], [44.1352, -95.4949]],
            grid="conus2"
        )
    """
    _validate_latlon_list(outlets)
    _validate_grid(grid)
    grid = grid.lower()
    try:
        flow_direction = hf_hydrodata.get_gridded_data(
            variable="flow_direction", grid=grid, file_type="tiff"
        )
    except Exception as exc:
        raise ValueError(
            f"Failed to get flow direction data for the grid {grid}."
        ) from exc
    nj, ni = flow_direction.shape
    marked = np.zeros((nj, ni), dtype=int)
    flow_values = [1, 2, 3, 4]  # D4 neighbors
    directions = np.array([[0, -1], [-1, 0], [0, 1], [1, 0]])
    queue = [hf_hydrodata.to_ij(grid, outlet[0], outlet[1]) for outlet in outlets]

    for point in queue:
        i, j = point
        if not np.isnan(flow_direction[j, i]):
            marked[j, i] = 1

    while queue:
        next_queue = []
        for point in queue:
            i, j = point
            # Look for cells that drain to this cell
            for direction, flow_value in zip(directions, flow_values):
                i_upstream = i - direction[0]
                j_upstream = j - direction[1]
                if (
                    0 <= i_upstream < ni
                    and 0 <= j_upstream < nj
                    and marked[j_upstream, i_upstream] == 0
                    and flow_direction[j_upstream, i_upstream] == flow_value
                ):
                    marked[j_upstream, i_upstream] = 1
                    next_queue.append([i_upstream, j_upstream])
        queue = next_queue

    masklist = np.argwhere(marked == 1)
    if masklist.size == 0:
        raise ValueError("Empty upstream area.")
    bounds = _indices_to_ij(masklist[:, 0], masklist[:, 1])
    imin, jmin, imax, jmax = bounds
    mask = marked[jmin:jmax, imin:imax]
    return bounds, mask
