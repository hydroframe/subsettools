import pytest
from subsettools import subsettools, datasets


@pytest.mark.parametrize(
    "huc_id, grid, result",
    [
        ("15060202", "conus1", [375, 239, 487, 329])
    ]
)
def test_get_conus_ij(huc_id, grid, result):
    assert subsettools.get_conus_ij(huc_id, grid) == result
