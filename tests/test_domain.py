import os
import filecmp
import pytest
import numpy as np
import subsettools as st


@pytest.mark.parametrize(
    "hucs, grid, result",
    [
        (["15060202"], "conus1", (375, 239, 487, 329)),
        (["140802"], "conus1", (572, 337, 797, 577)),
        (
            ["14080201", "14080202", "14080203", "14080204", "14080205"],
            "conus1",
            (572, 337, 797, 577),
        ),
    ],
)
def test_define_huc_domain(hucs, grid, result):
    bounds, _ = st.define_huc_domain(hucs, grid)
    assert bounds == result


@pytest.mark.parametrize(
    "hucs, grid",
    [
        (["01010001"], "conus1"),
        (["01010001"], "conus2"),
        (["03130003"], "conus1"),
        (["1710"], "conus1"),
        (["01010002", "01010001"], "conus1"),
    ],
)
def test_define_huc_domain_errors(hucs, grid):
    """Check that a ValueError is raised if a HUC is not part of the grid."""
    with pytest.raises(ValueError):
        st.define_huc_domain(hucs, grid)


@pytest.mark.parametrize(
    "latlon_bounds, grid, correct_bounds, correct_mask",
    [
        (
            [[37.91, -91.43], [37.34, -90.63]],
            "conus1",
            (2285, 436, 2359, 496),
            np.ones((60, 74)),
        ),
        (
            [[37.91, -91.43], [37.34, -90.63]],
            "conus2",
            (2684, 1403, 2758, 1461),
            np.ones((58, 74)),
        ),
        (
            [[40.17, -124.27], [40.19, -124.25]],
            "conus2",
            (4, 2066, 7, 2069),
            np.array([[0, 0, 1], [0, 0, 1], [0, 1, 1]]),
        ),
    ],
)
def test_define_latlon_domain(latlon_bounds, grid, correct_bounds, correct_mask):
    bounds, mask = st.define_latlon_domain(latlon_bounds, grid)
    assert bounds == correct_bounds
    assert np.array_equal(mask, correct_mask)


@pytest.mark.parametrize(
    "latlon_bounds, grid",
    [
        pytest.param(
            [[57.44, -107.33], [57.44, -107.33]],
            "conus2",
            id="lat/lon points outside the grid",
        ),
    ],
)
def test_define_latlon_domain_errors(latlon_bounds, grid):
    with pytest.raises(ValueError):
        st.define_latlon_domain(latlon_bounds, grid)


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
            np.array([[0, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]),
            (4054, 1914, 4058, 1917),
            id="small upstream area",
        ),
        pytest.param(
            [[44.1348, -95.5084], [44.1352, -95.4949]],
            "conus2",
            np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]]),
            (2322, 2111, 2325, 2116),
            id="two outlets with adjacent, non-overlapping upstream areas",
        ),
        pytest.param(
            [[44.1348, -95.5084], [44.1074, -95.5086]],
            "conus2",
            np.array(
                [
                    [1, 0],
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
def test_define_upstream_domain(outlets, grid, correct_mask, correct_bounds):
    bounds, mask = st.define_upstream_domain(outlets, grid)
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
def test_define_upstream_domain_errors(outlets, grid):
    with pytest.raises(ValueError):
        st.define_upstream_domain(outlets, grid)


def test_write_mask_solid(set_parflow_dir, tmp_path):
    hucs = ["15060202"]
    grid = "conus1"
    _, mask = st.define_huc_domain(hucs, grid)
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    st.write_mask_solid(mask, grid, test_dir)
    correct_output_dir = os.path.join(
        os.getcwd(), "tests", "correct_output", "conus1_upper_verde"
    )
    match, mismatch, errors = filecmp.cmpfiles(
        test_dir, correct_output_dir, ["mask.pfb", "solidfile.pfsol"], shallow=False
    )
    assert len(match) == 2 and len(mismatch) == 0 and len(errors) == 0
