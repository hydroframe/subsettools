import pytest
import numpy as np
from subsettools.upstream_area import upstream_area_to_ij


@pytest.mark.parametrize(
    "outlets, grid, correct_mask, correct_bounds",
    [
        pytest.param(
            [[39.8344, -74.3853]],
            "conus2",
            np.array([[1]]),
            (4057, 1914, 4058, 1915),
            id="single point upstream area",
        ),
        pytest.param(
            [[39.8522, -74.3786]],
            "conus2",
            np.array(
                [[0, 1, 1, 1],
                 [1, 1, 1, 1],
                 [1, 1, 1, 1]
                ]
            ),
            (4054, 1914, 4058, 1917),
            id="small upstream area",
        ),
        pytest.param(
            [[44.1348, -95.5084], [44.1352, -95.4949]],
            "conus2",
            np.array(
                [[1, 0, 0],
                 [1, 1, 0],
                 [1, 1, 1],
                 [1, 1, 1],
                 [1, 1, 1]
                ]
            ),
            (2322, 2111, 2325, 2116),
            id="two outlets with adjacent, non-overlapping upstream areas",
        ),
        pytest.param(
            [[44.1348, -95.5084], [44.1074, -95.5086]],
            "conus2",
            np.array(
                [[1, 0],
                 [1, 1],
                 [1, 1],
                 [1, 1],
                 [1, 1],
                ]
            ),
            (2322, 2111, 2324, 2116),
            id="two outlets, one upstream of the other",
        ),
        pytest.param(
            [[39.8195, -75.3820]],
            "conus2",
            np.load("tests/correct_output/upstream_mask_lower_delaware.npy"),
            (3874, 1883, 4054, 2180),
            id="large upstream area",
        ),
    ],
)
def test_upstream_area_to_ij(outlets, grid, correct_mask, correct_bounds):
    bounds, mask = upstream_area_to_ij(outlets, grid)
    assert np.array_equal(mask, correct_mask)
    assert bounds == correct_bounds


@pytest.mark.parametrize(
    "outlets, grid",
    [
        pytest.param([[57.44, -107.33]], "conus2", id="outlet outside the grid"),
        pytest.param(
            [[22.36, -117.85]],
            "conus2",
            id="outlet [0, 0], where flow_direction is NaN",
        ),
    ],
)
def test_upstream_area_to_ij_errors(outlets, grid):
    with pytest.raises(ValueError) as e:
        upstream_area_to_ij(outlets, grid)
