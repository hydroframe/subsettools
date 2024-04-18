import pytest
import numpy as np
from subsettools.upstream_area import delin_watershed

# unit test, multiple points, point outside grid, big area test, conus1 tests

@pytest.mark.parametrize(
    "outlets, grid, points, bounds", [([[39.8344, -74.3853]],
                                       "conus2",
                                       np.array([[1914, 4057]]),
                                       (4057, 1914, 4058, 1915)),
                                      ([[39.8522, -74.3786]],
                                       "conus2",
                                       np.array([[1914, 4055],
                                                 [1914, 4056],
                                                 [1914, 4057],
                                                 [1915, 4054],
                                                 [1915, 4055],
                                                 [1915, 4056],
                                                 [1915, 4057],
                                                 [1916, 4054],
                                                 [1916, 4055],
                                                 [1916, 4056],
                                                 [1916, 4057]]),
                                       (4054, 1914, 4058, 1917)),
    ]
)
def test_delin_watershed(outlets, grid, points, bounds):
    ij_bounds, mask = delin_watershed(outlets, grid)
    assert np.array_equal(np.argwhere(mask == 1), points)
    assert ij_bounds == bounds


