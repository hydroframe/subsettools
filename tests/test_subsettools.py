import pytest
from subsettools import subsettools, datasets


@pytest.mark.parametrize(
    "huc_list, grid, result", [(["15060202"], "conus1", [375, 239, 487, 329])]
)
def test_huc_to_ij(huc_list, grid, result):
    assert subsettools.huc_to_ij(huc_list, grid) == result
