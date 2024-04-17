import numpy as np
import hf_hydrodata
from ._error_checking import (
    _validate_grid,
    _validate_latlon_list,
    )

def delin_watershed(outlets, grid):
    """Calculate the upstream area of a point or list of points in the grid.

    The flow_direction files that are used to create the upstream area mask follow the convention: down: 1, left: 2, up: 3, right: 4.

    Args:
        outlets (List[List[float]]): list of lat-lon points of the form [[lat1, lon1], [lat2, lon2], ...]
        grid (str): The spatial grid that the upstream area will be returned on. Possible values: “conus1” or “conus2”

    Returns:
        TODO

    Raises:
        TODO

    Example:

    .. code-block:: python

        TODO
    """
    _validate_latlon_list(outlets)
    _validate_grid(grid)
    grid = grid.lower()
    try:
        flow_direction = hf_hydrodata.get_gridded_data(variable="flow_direction", grid=grid, file_type="tiff")
    except Exception as e:
        raise ValueError(f"There is no flow direction entry for the grid {grid}.")
    print("flow_direction shape: ", flow_direction.shape)
    nj, ni = flow_direction.shape
    
    # Initialize a matrix to store the mask
    marked = np.zeros((nj, ni), dtype=int)
    print("marked shape: ", marked.shape)
    
    # D4 neighbors
    d4 = [1, 2, 3, 4]
    kd = np.zeros((4, 2), dtype=int)
    kd[:, 0] = [0, -1, 0, 1]
    kd[:, 1] = [-1, 0, 1, 0]
    
    # Initialized the queue with the outlet points in ij form
    queue = [hf_hydrodata.to_ij(grid, outlet[0], outlet[1]) for outlet in outlets]

    # Add the outlets to the mask
    for point in queue:
        i, j = point
        print("point: ", point)
        print("i: ", i, "j: ", j)
        marked[j, i] = 1
    
    while queue:
        next = []
        for point in queue:
            i, j = point
            # Look for cells that drain to this cell
            for d in range(4):
                i_upstream = i - kd[d, 0]
                j_upstream = j - kd[d, 1]
                if (i_upstream > 0 and j_upstream > 0 and i_upstream < ni and j_upstream < nj
                        and not np.isnan(flow_direction[j_upstream, i_upstream]) and marked[j_upstream, i_upstream] == 0):
                    if flow_direction[j_upstream, i_upstream] == d4[d]:
                        marked[j_upstream, i_upstream] = 1  # Add the upstream cell to the mask
                        next.append([i_upstream, j_upstream])
        print("next: ", next)
        queue = next

    masklist = np.argwhere(marked == 1)
    jrange = (masklist[:, 0].min(), masklist[:, 0].max())
    irange = (masklist[:, 1].min(), masklist[:, 1].max())

    return  {'watershed': marked, 'irange': irange, 'jrange': jrange}
