"""Functions to define a domain on the CONUS grids.

A domain is defined by the grid ij bounds and a mask array. 

The grid bounds is a tuple (imin, jmin, imax, jmax) that defines a bounding box
on the CONUS grid that encompases the domain. Note that the origin (0, 0) of the
CONUS grids is always the southwest corner. Thus, (imin, jmin) is the southwest
corner of the bounding box, and (imax - 1, jmax - 1) the northeast corner (upper
bounds are exclusive in keeping with Python convention).

The mask is a 0/1 NumPy array which has the same shape as the bounding box and
indicates which cells are part of the subset domain.

SubsetTools provides functions to define domains in different ways.
    - A domain that is a collections of HUCs.
    - A box that is defined by two latitude-longitude points.
    - A domain that is the upstream area of a collection of outlets.
"""

import os
import subprocess
import warnings
import numpy as np
import hf_hydrodata
from parflow.tools.io import write_pfb
from ._error_checking import (
    _validate_huc_list,
    _validate_grid,
    _validate_latlon_list,
    _validate_dir,
    _validate_mask,
)
from ._constants import (
    CONUS_DX,
    CONUS_DY,
    CONUS1_DZ,
    CONUS1_Z_TOP,
    CONUS2_DZ,
    CONUS2_Z_TOP,
    CONUS_Z_BOTTOM,
)


def define_huc_domain(hucs, grid):
    """Define a domain by a collection of HUCs.

    The domain is defined by the grid ij bounds of a bounding box that
    encompasses the HUCs in the list and a mask for that bounding box indicating
    which cells in the bounding box are part of these HUCs.

    All HUC IDs in hucs must be the same length (HUCs of the same level).
    All HUCs should be adjacent. If a HUC is only partially covered by the
    provided grid, the grid bounds for the covered area will be returned.

    Args:
        hucs (list[str]): a list of USGS HUC IDs
        grid (str): The spatial grid that the ij indices are calculated relative
            to and that the subset data will be returned on. Possible values:
            “conus1” or “conus2”

    Returns:
        A tuple (bounds, mask).

        Bounds is a tuple of the form (imin, jmin, imax, jmax) representing the
        bounds in the conus grid of the area defined by the HUC IDs in hucs.
        imin, jmin, imax, jmax are the west, south, east and north sides of the
        box respectively and all i,j indices are calculated relative to the
        lower southwest corner of the domain.

        Mask is a 2D numpy.ndarray that indicates which cells inside the bounding
        box are part of the selected HUC(s).

    Raises:
        ValueError: If the area defined by the provided HUCs is not part of the
            given grid.

    Example:

    .. code-block:: python

        grid_bounds, mask = define_huc_domain(
            hucs=["14080201", "14080202", "14080203"], grid="conus1"
        )
    """
    _validate_huc_list(hucs)
    _validate_grid(grid)

    huc_len = len(hucs[0])
    hucs = [int(huc) for huc in hucs]
    try:
        conus_hucs = hf_hydrodata.get_gridded_data(
            dataset="huc_mapping",
            grid=grid,
            file_type="tiff",
            level=str(huc_len),
        )
    except Exception as exc:
        raise ValueError(
            f"Failed to get huc mapping data for the grid {grid}."
        ) from exc
    huc_mask = np.isin(conus_hucs, hucs).squeeze()
    indices_j, indices_i = np.where(huc_mask > 0)
    if indices_i.size == 0 or indices_j.size == 0:
        raise ValueError(
            f"The area defined by the provided HUCs is not part of the {grid} grid."
        )

    bounds = _indices_to_ij(indices_j, indices_i)
    imin, jmin, imax, jmax = bounds
    return bounds, huc_mask[jmin:jmax, imin:imax].astype(int)


def huc_to_ij(huc_list, grid):
    """This function is deprecated.

    Use define_huc_domain() instead.
    """
    warnings.warn(
        "This function is deprecated. Use define_huc_domain() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    _validate_huc_list(huc_list)
    _validate_grid(grid)
    _, _, indices_j, indices_i = _get_conus_hucs_indices(huc_list, grid.lower())
    if indices_i.size == 0 or indices_j.size == 0:
        raise ValueError(
            f"The area defined by the provided HUCs is not part of the {grid} grid."
        )
    return _indices_to_ij(indices_j, indices_i)


def _get_conus_hucs_indices(huc_list, grid):
    """Get the huc datafile as an ndarray and three mask arrays representing the selected hucs.

    Args:
        huc_list (list[str]): a list of huc IDs
        grid (str): The spatial grid that the ij indices are calculated relative to and that the subset
            data will be returned on. Possible values: “conus1” or “conus2”

    Returns:
        A tuple (conus_hucs, sel_hucs, indices_j, indices_i) where
        conus_hucs is an ndarray of the huc datafile, sel_hucs is
        a mask array for the selected hucs, and indices_i and
        indices_j mask arrays in the j and i directions.
    """
    huc_len = len(huc_list[0])
    huc_list = [int(huc) for huc in huc_list]
    entry = hf_hydrodata.get_catalog_entry(
        dataset="huc_mapping", grid=grid, file_type="tiff"
    )
    if entry is None:
        raise ValueError(f"There is no HUC mapping entry for grid {grid}.")
    conus_hucs = hf_hydrodata.gridded.get_ndarray(entry, level=str(huc_len))
    sel_hucs = np.isin(conus_hucs, huc_list).squeeze()
    indices_j, indices_i = np.where(sel_hucs > 0)
    return conus_hucs, sel_hucs, indices_j, indices_i


def _indices_to_ij(indices_j, indices_i):
    """Get the grid ij bounds for the domain defined by indices_j and indices_i.

    Args:
        indices_j (numpy.ndarray): mask in the j direction for selected hucs
        indices_i (numpy.ndarray): mask in the i direction for selected hucs

    Returns:
        A tuple of the form (imin, jmin, imax, jmax) representing the bounds
        in the grid defined by the two mask arrays indices_j and indices_i.
    """
    imin = np.min(indices_i)
    imax = np.max(indices_i) + 1
    jmin = np.min(indices_j)
    jmax = np.max(indices_j) + 1
    return (int(imin), int(jmin), int(imax), int(jmax))


def define_latlon_domain(latlon_bounds, grid):
    """Define a domain by latitude/longitude bounds.

    The domain is defined by the grid ij bounds of a bounding box formed by the
    latitude/longitude bounds (latlon_bounds) relative to the selected conus grid
    and a mask for that bounding box indicating which cells are active CONUS
    points.

    Args:
        latlon_bounds (List[List[float]]): list of the form [[lat1, lon1],
            [lat2, lon2]]. [lat1, lon1] and [lat2, lon2] define the northwest
            and southeast corners of the desired box respectively.
        grid (str):  The spatial grid that the ij indices are calculated relative
            to and that the subset data will be returned on. Possible values:
            “conus1” or “conus2”.

    Returns:
        A tuple (bounds, mask).

        Bounds is a tuple of the form (imin, jmin, imax, jmax) representing the
        bounds in the conus grid of the area defined by the latlon_bounds. imin,
        jmin, imax, jmax are the west, south, east and north sides of the box
        respectively and all i,j indices are calculated relative to the lower
        southwest corner of the domain.

        Mask is a 2D numpy.ndarray that indicates which cells inside the bounding
        box are active CONUS points (for example, if ocean is part of the bounding
        box the corresponding cells will not be part of the mask).

    Example:

    .. code-block:: python

        grid_bounds, mask = define_latlon_domain(
            latlon_bounds=[[37.91, -91.43], [37.34, -90.63]], grid="conus2"
        )
    """
    _validate_grid(grid)
    _validate_latlon_list(latlon_bounds)
    if len(latlon_bounds) != 2:
        raise ValueError("latlon_bounds must contain exactly two points.")

    grid = grid.lower()
    point0, point1 = [
        hf_hydrodata.to_ij(grid, latlon_pt[0], latlon_pt[1])
        for latlon_pt in latlon_bounds
    ]
    imin, imax = [min(point0[0], point1[0]), max(point0[0], point1[0]) + 1]
    jmin, jmax = [min(point0[1], point1[1]), max(point0[1], point1[1]) + 1]
    grid_bounds = (imin, jmin, imax, jmax)
    try:
        mask = hf_hydrodata.get_gridded_data(
            dataset="huc_mapping",
            grid=grid,
            level="2",
            grid_bounds=grid_bounds,
        )
    except Exception as exc:
        raise ValueError(
            f"Failed to get huc mapping data for the grid {grid}."
        ) from exc
    mask[mask > 0] = 1
    return grid_bounds, mask.astype(int)


def latlon_to_ij(latlon_bounds, grid):
    """This function is deprecated.

    Use define_latlon_domain() instead.
    """
    warnings.warn(
        "This function is deprecated. Use define_latlon_domain() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    _validate_grid(grid)
    _validate_latlon_list(latlon_bounds)
    if len(latlon_bounds) != 2:
        raise ValueError(
            "latlon_bounds must contain exactly two lat-lon points: [[lat1, lon1], [lat2, lon2]]"
        )

    grid = grid.lower()
    point0 = hf_hydrodata.to_ij(grid, latlon_bounds[0][0], latlon_bounds[0][1])
    point1 = hf_hydrodata.to_ij(grid, latlon_bounds[1][0], latlon_bounds[1][1])
    imin, imax = [
        min(point0[0], point1[0]),
        max(point0[0], point1[0]),
    ]
    jmin, jmax = [
        min(point0[1], point1[1]),
        max(point0[1], point1[1]),
    ]
    return (imin, jmin, imax, jmax)


def write_mask_solid(mask, grid, write_dir):
    """Create ParFlow mask and solid files from a mask array.

    Given an integer mask array consisting of 0s and 1s, this function will
    create three files in write_dir.
        - a 2D mask file that indicates which cells inside the box domain are
          part of the selected HUCS.
        - a solid file that defines a 3D domain extending to the depth of
          whichever grid has been selected and tracing the boundaries of the
          selected HUCS.
        - a vtk file, which can be used to visualize the solid file in ParaView.

    Args:
        mask (numpy.ndarray): an integer array such that mask[i, j] == 1 if the
            cell (i, j) is part of the domain, and mask[i, j] == 0 otherwise.
        grid (str): The spatial grid that the ij indices are calculated relative
            to and that the subset data will be returned on. Possible values:
            “conus1” or “conus2”
        write_dir (str): directory path where the mask and solid files will be
            written

    Returns:
        dict: A dictionary mapping the keys ("mask", "mask_vtk", "solid") to the
            corresponding filepaths of the created files.

    Example:

    .. code-block:: python

        filepaths = write_mask_solid(
            mask=np.array([[0, 1], [1, 1]]),
            grid="conus2",
            write_dir="/path/to/your/chosen/directory"
        )
    """
    _validate_mask(mask)
    _validate_grid(grid)
    _validate_dir(write_dir)
    grid = grid.lower()

    if grid == "conus1":
        dz = CONUS1_DZ
        z_top = CONUS1_Z_TOP
    elif grid == "conus2":
        dz = CONUS2_DZ
        z_top = CONUS2_Z_TOP

    nj, ni = mask.shape
    new_mask = mask.reshape((1, nj, ni)).astype(float)
    mask_path = os.path.join(write_dir, "mask.pfb")
    write_pfb(mask_path, new_mask, dx=CONUS_DX, dy=CONUS_DY, dz=dz, dist=False)
    print("Wrote mask.pfb")
    mask_vtk_path = os.path.join(write_dir, "mask_vtk.vtk")
    solid_path = os.path.join(write_dir, "solidfile.pfsol")

    try:
        parflow_dir = os.environ["PARFLOW_DIR"]
    except KeyError:
        raise KeyError(
            "The environment variable PARFLOW_DIR is undefined. Please make "
            'sure you have ParFlow installed and os.environ["PARFLOW_DIR"] '
            "points to that installation."
        )
    script_path = os.path.join(parflow_dir, "bin", "pfmask-to-pfsol")
    if not os.path.exists(script_path):
        raise FileNotFoundError(
            "pfmask-to-pfsol file not found. Please make sure you have ParFlow "
            'installed and os.environ["PARFLOW_DIR"] points to that '
            "installation."
        )
    try:
        subprocess.run(
            [
                script_path,
                "--mask",
                mask_path,
                "--pfsol",
                solid_path,
                "--vtk",
                mask_vtk_path,
                "--z-bottom",
                str(CONUS_Z_BOTTOM),
                "--z-top",
                str(z_top),
            ],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError("pfmask-to-pfsol error:", e.stderr)

    print(f"Wrote solidfile and mask_vtk with total z of {z_top} meters")
    file_paths = {"mask": mask_path, "mask_vtk": mask_vtk_path, "solid": solid_path}
    return file_paths


def create_mask_solid(huc_list, grid, write_dir):
    """This function is deprecated.

    Use write_mask_solid() instead.
    """
    warnings.warn(
        "This function is deprecated. Use write_mask_solid() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    _validate_huc_list(huc_list)
    _validate_grid(grid)
    _validate_dir(write_dir)
    grid = grid.lower()
    _, sel_hucs, indices_j, indices_i = _get_conus_hucs_indices(huc_list, grid)
    if indices_i.size == 0 or indices_j.size == 0:
        raise ValueError(
            f"The area defined by the provided HUCs is not part of the {grid} grid."
        )
    imin, jmin, imax, jmax = _indices_to_ij(indices_j, indices_i)
    nj = jmax - jmin
    ni = imax - imin

    # checks conus1 / 2 grid and assigns appripriate dz and z_total for making the mask and solid file
    if grid == "conus1":
        print("grid is conus1")
        layz = 100
        z_total = str(500)
    else:
        print("grid is conus2")
        layz = 200
        z_total = str(2000)

    # create and write the pfb mask
    mask_clip = np.zeros((1, nj, ni))
    mask_clip[0, :, :] = sel_hucs[jmin:jmax, imin:imax]
    mask_clip = mask_clip.astype(float)
    mask_file_path = os.path.join(write_dir, "mask.pfb")
    write_pfb(mask_file_path, mask_clip, dx=1000, dy=1000, dz=layz, dist=False)
    print("Wrote mask.pfb")
    mask_vtk_path = os.path.join(write_dir, "mask_vtk.vtk")
    solid_file_path = os.path.join(write_dir, "solidfile.pfsol")

    try:
        parflow_dir = os.environ["PARFLOW_DIR"]
    except KeyError:
        raise KeyError(
            "The environment variable PARFLOW_DIR has not been defined. Please make sure you have ParFlow installed "
            'and os.environ["PARFLOW_DIR"] points to that installation.'
        )
    file_path = os.path.join(parflow_dir, "bin", "pfmask-to-pfsol")
    if not os.path.exists(file_path):
        raise FileNotFoundError(
            "pfmask-to-pfsol file not found. Please make sure you have ParFlow installed "
            'and os.environ["PARFLOW_DIR"] points to that installation.'
        )
    try:
        subprocess.run(
            [
                file_path,
                "--mask",
                mask_file_path,
                "--pfsol",
                solid_file_path,
                "--vtk",
                mask_vtk_path,
                "--z-bottom",
                "0.0",
                "--z-top",
                z_total,
            ],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError("pfmask-to-pfsol error:", e.stderr)

    print(f"Wrote solidfile.pfsol and mask_vtk.vtk with total z of {z_total} meters")
    file_paths = {
        "mask": mask_file_path,
        "mask_vtk": mask_vtk_path,
        "solid": solid_file_path,
    }
    return file_paths


def define_upstream_domain(outlets, grid):
    """Define a domain that is the upstream area of the points in outlets.

    The domain is defined by the grid ij bounds of the bounding box that
    encompasses the upstream area of all the points in outlets and a mask for
    that bounding box indicating which cells are part of the selected area.

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

        bounds, mask = define_upstream_domain(
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
