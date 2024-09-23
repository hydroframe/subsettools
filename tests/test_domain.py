import os
import filecmp
import pytest
import numpy as np
import subsettools as st
from parflow.tools.io import read_pfb


@pytest.mark.parametrize(
    "hucs, grid, result",
    [
        pytest.param(
            ["15060202"],
            "conus1",
            (375, 239, 487, 329),
            id="simple conus1 level 8 HUC domain",
        ),
        pytest.param(
            ["140802"],
            "conus1",
            (572, 337, 797, 577),
            id="simple conus1 level 6 HUC domain",
        ),
        pytest.param(
            ["14080201", "14080202", "14080203", "14080204", "14080205"],
            "conus1",
            (572, 337, 797, 577),
            id="list of HUCs in conus 1",
        ),
        pytest.param(
            ["02050304"],
            "conus2",
            (3744, 1853, 3854, 1952),
            id="level 8 conus2 HUC domain",
        ),
        pytest.param(
            ["17100306"],
            "conus2",
            (73, 2295, 110, 2368),
            id="level 8 conus2 coastal HUC domain",
        ),
    ],
)
def test_define_huc_domain(hucs, grid, result):
    bounds, _ = st.define_huc_domain(hucs, grid)
    assert bounds == result


@pytest.mark.parametrize(
    "hucs, grid",
    [
        pytest.param(
            ["01010001"], "conus1", id="Level 8 HUC located outside of conus1 grid"
        ),
        pytest.param(
            ["21010005"], "conus2", id="Level 8 HUC located outside of conus2 grid"
        ),
        pytest.param(
            ["03130003"], "conus1", id="Level 8 HUC located outside of conus1 grid"
        ),
        pytest.param(
            ["1710"], "conus1", id="Level 4 HUC located outside of conus1 grid"
        ),
        pytest.param(
            ["01010002", "01010001"],
            "conus1",
            id="list of Level 8 HUCs located outside of conus1 grid",
        ),
        pytest.param(["220102"], "conus2", id="HUC number does not exist"),
        pytest.param(
            ["14080205", "140700"], "conus2", id="HUC IDs of different levels"
        ),
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
            (2284, 435, 2359, 495),
            np.ones((60, 75)),
        ),
        (
            [[37.91, -91.43], [37.34, -90.63]],
            "conus2",
            (2683, 1403, 2757, 1460),
            np.ones((57, 74)),
        ),
        (
            [[40.17, -124.27], [40.19, -124.25]],
            "conus2",
            (3, 2066, 7, 2068),
            np.array([[0, 0, 0, 1], [0, 0, 0, 1]]),
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
        pytest.param([[57.44, -107.33]], "conus2", id="only one lat/lon point"),
        pytest.param(
            [[57.44, -107.33], [57.44, -107.33], [57.55, -108.00]],
            "conus2",
            id="more than two lat/lon points",
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
            [[39.83447238411436, -74.38539661569031]],
            "conus2",
            np.array([[1]]),
            (4057, 1914, 4058, 1915),
            id="single point upstream area",
        ),
        pytest.param(
            [[39.85227324179805, -74.3786676880639]],
            "conus2",
            np.array([[0, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]),
            (4054, 1914, 4058, 1917),
            id="small upstream area",
        ),
        pytest.param(
            [
                [44.13538466649678, -95.50793718973576],
                [44.13521041397351, -95.4949657382491],
            ],
            "conus2",
            np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]]),
            (2322, 2111, 2325, 2116),
            id="two outlets with adjacent, non-overlapping upstream areas",
        ),
        pytest.param(
            [
                [44.13538466649678, -95.50793718973576],
                [44.10762872972563, -95.52162729883625],
            ],
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
            [[39.820145944247464, -75.3815714544093]],
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


@pytest.mark.parametrize(
    "mask",
    [np.array([[0]]), np.array([[1]]), np.array([[0, 1], [0, 1]])],
)
def test_write_mask(mask, set_parflow_dir, tmp_path):
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    st.write_mask_solid(mask=mask, grid="conus2", write_dir=test_dir)
    data = read_pfb(f"{test_dir}/mask.pfb")[0]
    assert np.array_equal(data, mask)


def test_write_solid(set_parflow_dir, tmp_path):
    mask = np.array([[0]])
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    st.write_mask_solid(mask=mask, grid="conus2", write_dir=test_dir)
    expected_solid = ["1\n", "0\n", "1\n", "0\n", "0\n"]
    with open(f"{test_dir}/solidfile.pfsol") as f:
        solid = f.readlines()
    assert solid == expected_solid
