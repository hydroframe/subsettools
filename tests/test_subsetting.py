import os
from parflow.tools.io import read_pfb
import numpy as np
import subsettools as st
import hf_hydrodata as hf
import pytest


#  resolution: static or hourly
#  if it hourly
#      hourly press -> 4D -> 3D (t, p_x, p_y, p_z) and t = 0
#      forcing data -> 3D (t, f_x, f_y)
#  if static, 2-d or 3-d, can pretend always 3-d
#
#  have options dict
#  imin, jmin, imax, jmax = options['grid_bounds']
#  ni, nj = imax - imin, jmax - jmin
#  if resolution == 'hourly':
#      nt = end_time - start_time if end_time is not None else 1
#  if variable == 'pressure_head':
#      shape = (nt, 5, nj, ni)    
#  else:
#      shape = (nt, nj, ni)
#  if resolution == 'static':
#      shape = (5, nj, ni)

_SECONDS_PER_HOUR = 3600
_DUMMY_NZ = 5

@pytest.fixture
def mock_hf(monkeypatch):

    def mock_get_data(options):
        grid_bounds = options.get("grid_bounds")
        resolution = options.get("temporal_resolution")
        variable = options.get("variable")
        start_time = options.get("start_time")
        end_time = options.get("end_time")
        imin, jmin, imax, jmax = grid_bounds
        ni, nj = imax - imin, jmax - jmin
        if resolution == "hourly":
            if end_time is not None:
                nt = int((end_time - start_time).total_seconds() / _SECONDS_PER_HOUR)
            else:
                nt = 1
            if variable == "pressure_head":
                shape = (nt, _DUMMY_NZ, ni, nj)
            else:
                shape = (nt, ni, nj)
        elif resolution == "static":
            shape = (_DUMMY_NZ, ni, nj)            
        return np.ones(shape)
    
    monkeypatch.setattr(hf, "get_gridded_data", mock_get_data)


def test_subset_static(tmp_path, mock_hf):
    test_dir = tmp_path / "test_static"
    test_dir.mkdir()
    paths = st.subset_static(
        ij_bounds=(0, 0, 10, 20),
        dataset="my_dataset",
        var_list=("var1", "var2"),
        write_dir=test_dir,
    )
    assert os.path.isfile(paths["var1"])
    data1 = read_pfb(paths["var1"])
    assert np.array_equal(data1, np.ones((_DUMMY_NZ, 10, 20)))
    assert os.path.isfile(paths["var2"])
    data2 = read_pfb(paths["var2"])
    assert np.array_equal(data2, np.ones((_DUMMY_NZ, 10, 20)))

def test_subset_press_init_function(tmp_path, mock_hf):
    test_dir = tmp_path / "test_press_init"
    test_dir.mkdir()
    file_path = st.subset_press_init(
        ij_bounds=(0,0,10,20),
        dataset="my_dataset",
        date="2010-10-02",
        write_dir=test_dir
    )
    assert os.path.isfile(file_path)
    data = read_pfb(file_path)
    print(data.shape)
    assert np.array_equal(data, np.ones((_DUMMY_NZ, 10, 20)))

def test_subset_forcing(tmp_path, mock_hf):
    test_dir = tmp_path / "test_forcing"
    test_dir.mkdir()
    paths = st.subset_forcing(
        ij_bounds=(0,0,10,20),
        grid=""
    )

    

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
