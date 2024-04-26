import pytest
import numpy as np
from subsettools.upstream_area import delin_watershed

# unit tests, multiple points, point outside grid, big area test, conus1 tests

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
                                      ([[39.8195, -75.3820]],
                                       "conus2",
                                       np.load("tests/correct_output/upstream_mask_lower_delaware.npy"),
                                       (3874, 1883, 4054, 2180)),
    ]
)
def test_delin_watershed(outlets, grid, points, bounds):
    ij_bounds, mask = delin_watershed(outlets, grid)
    assert np.array_equal(np.argwhere(mask == 1), points)
    assert ij_bounds == bounds


@pytest.mark.parametrize(
    "outlets, grid", [([[57.44, -107.33]], "conus2"),    # point outside the grid
                      ([[22.36, -117.85]], "conus2"),    # point corresponds to [0, 0], where flow_direction is NaN
    ]
)
def test_delin_watershed_errors(outlets, grid):
    with pytest.raises(ValueError) as e:
        delin_watershed(outlets, grid)
