"""This module contains helper functions for the subsetting module."""

import os
from datetime import datetime, timedelta
import pytz
import numpy as np
from parflow.tools.io import read_clm

_VEGM_COLUMNS = 25


def subset_vegm(path, ij_bounds):
    """Read in vegm file and subset it according to ij_bounds.

    Subset a national vegm file to a smaller domain based on ij_bounds that are
    provided relative to the grid that the national file is on.

    Args:
        path (str): path to read vegm file
        ij_bounds (Tuple[int]): bounding box for subset. This should be given as
            i,j index values where 0,0 is the lower left hand corner of a domain.
            ij_bounds are given relative to whatever grid is being used for the
            subset.

    Returns:
        ndarray:
            Subset vegm data.
    """
    vegm_data = read_clm(path, type="vegm")  # returns (j,i,k)
    vegm_data = np.transpose(vegm_data, (2, 0, 1))  # transpose to k,j,i

    imin, jmin, imax, jmax = ij_bounds
    vegm_data = vegm_data[:, jmin:jmax, imin:imax]  # slicing on k,j,i

    return _reshape_ndarray_to_vegm_format(vegm_data)


def _reshape_ndarray_to_vegm_format(data):
    """Reshape ndarray returned by datacatalog to vegm format.

    Args:
        data (ndarray): raw subset vegm data (2d array)

    Returns:
        Ndarray reshaped to vegm format.
    """
    _, nj, ni = data.shape
    indices = np.indices((nj, ni)) + 1
    indices = indices[::-1, :, :]
    data = np.vstack([indices, data])  # stack x,y indices on vegm
    # transpose and reshape back into expected 2D vegm file format for the subset
    return data.transpose(1, 2, 0).reshape(-1, _VEGM_COLUMNS)


def get_utc_time(date_string, time_zone):
    """Convert the given date and time_zone to UTC time.

    Args:
        date_string (str): date in the form 'yyyy-mm-dd'
        time_zone (str): a pytz-supported time zone

    Returns:
        A timezone-unaware datetime object representing the time in UTC.
    """
    date = datetime.strptime(date_string, "%Y-%m-%d")
    if time_zone != "UTC":
        date = (
            date.replace(tzinfo=pytz.timezone(time_zone))
            .astimezone(pytz.UTC)
            .replace(tzinfo=None)
        )
    return date
