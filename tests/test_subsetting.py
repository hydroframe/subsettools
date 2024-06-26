import os
from parflow.tools.io import read_pfb
import numpy as np
import subsettools as st


def test_forcing_timezones(tmp_path):
    "Check if we get the correct forcing (temperature) in EST time."
    utc = tmp_path / "UTC_out"
    utc.mkdir()
    est = tmp_path / "EST_out"
    est.mkdir()
    ij_bounds = (375, 239, 487, 329)
    grid = "conus1"
    start = "2005-10-03"
    dataset = "NLDAS2"
    st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end="2005-10-05",
        dataset=dataset,
        write_dir=utc,
        time_zone="UTC",
    )
    st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end="2005-10-04",
        dataset=dataset,
        write_dir=est,
        time_zone="EST",
    )
    utc_temp1 = read_pfb(os.path.join(utc, "NLDAS.Temp.000001_to_000024.pfb"))
    utc_temp2 = read_pfb(os.path.join(utc, "NLDAS.Temp.000025_to_000048.pfb"))
    est_temp_correct = np.concatenate(
        (utc_temp1[5:, :, :], utc_temp2[:5, :, :]), axis=0
    )
    est_temp = read_pfb(os.path.join(est, "NLDAS.Temp.000001_to_000024.pfb"))
    assert np.array_equal(est_temp_correct, est_temp)


def test_subset_press_init(tmp_path):
    """Check that the call succeeds when it fetches data for the beginning of WY 2003.

    subset_press_init was initially designed to fetch data one hour before midnight on
    the date given. So in this case, it would look at 11pm on 2002-09-30, which does
    not exist.
    """
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    filename = st.subset_press_init(
        ij_bounds=(375, 239, 487, 329),
        dataset="conus1_baseline_mod",
        date="2002-10-01",
        write_dir=test_dir,
    )
    assert filename == os.path.join(
        test_dir, "conus1_baseline_mod_2002.10.01:00.00.00_UTC0_press.pfb"
    )
