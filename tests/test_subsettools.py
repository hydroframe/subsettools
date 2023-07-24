from subsettools import subsettools, datasets

def test_get_conus_ij():
    assert subsettools.get_conus_ij("15060202", "conus1") == [375, 239, 487, 329]
