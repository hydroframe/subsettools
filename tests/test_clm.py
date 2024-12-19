"""Tests for clm.py module."""
import pytest
import numpy as np
import subsettools as st


def test_vegm_to_land_cover():
    """Test vegm_to_land_cover conversion function."""
    test_land_cover = st.vegm_to_land_cover("tests/correct_output/test_vegm.dat")
    assert np.array_equal(test_land_cover, [[13., 6.], [18., 1.]])
