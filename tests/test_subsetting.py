import os
from parflow.tools.io import read_pfb
import numpy as np
import subsettools as st
import hf_hydrodata
import pytest
import glob


_SECONDS_PER_HOUR = 3600
_DUMMY_NZ = 5


@pytest.fixture
def mock_hf_data(monkeypatch):
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
                shape = (nt, _DUMMY_NZ, nj, ni)
            else:
                shape = (nt, nj, ni)
        elif resolution == "static":
            shape = (_DUMMY_NZ, nj, ni)
        return np.ones(shape)

    monkeypatch.setattr(hf_hydrodata, "get_gridded_data", mock_get_data)


@pytest.fixture
def mock_hf_paths(monkeypatch):
    def mock_get_paths(options):
        dataset = options.get("dataset")
        variable = options.get("variable")
        file_path = [f"/folder/{dataset}.{variable}.000000_to_000000.pfb"]
        return file_path

    monkeypatch.setattr(hf_hydrodata, "get_paths", mock_get_paths)


def test_subset_static_filenames(tmp_path, mock_hf_data):
    test_dir = tmp_path / "test_static"
    test_dir.mkdir()
    var_list = ("var1", "var2", "var3")
    static_paths = st.subset_static(
        ij_bounds=(0, 0, 10, 20),
        dataset="dataset",
        var_list=var_list,
        write_dir=test_dir,
    )
    expected_paths = {
        var: glob.glob(os.path.join(test_dir, "*" + var + "*"), root_dir=test_dir)[0]
        for var in var_list
    }
    assert static_paths == expected_paths


def test_subset_static_data(tmp_path, mock_hf_data):
    test_dir = tmp_path / "test_static"
    test_dir.mkdir()
    static_paths = st.subset_static(
        ij_bounds=(0, 0, 10, 20),
        dataset="dataset",
        var_list=("precipitation",),
        write_dir=test_dir,
    )
    data = read_pfb(static_paths["precipitation"])
    assert np.array_equal(data, np.ones((_DUMMY_NZ, 20, 10)))


def test_subset_press_init_filenames(tmp_path, mock_hf_data):
    test_dir = tmp_path / "test_press_filenames"
    test_dir.mkdir()
    file_path = st.subset_press_init(
        ij_bounds=(0, 0, 10, 20),
        dataset="dataset",
        date="2010-10-02",
        write_dir=test_dir,
    )
    assert os.path.join(test_dir, os.listdir(test_dir)[0]) == file_path


def test_subset_press_init_data(tmp_path, mock_hf_data):
    test_dir = tmp_path / "test_press_data"
    test_dir.mkdir()
    file_path = st.subset_press_init(
        ij_bounds=(0, 0, 10, 20),
        dataset="dataset",
        date="2010-10-02",
        write_dir=test_dir,
    )
    data = read_pfb(file_path)
    assert np.array_equal(data, np.ones((_DUMMY_NZ, 20, 10)))


@pytest.mark.parametrize(
    "end_date, num_files",
    [
        ("2005-10-02", 1),
        ("2005-10-05", 4),
        ("2005-11-01", 31),
    ],
)
def test_subset_forcing_num_files(
    end_date, num_files, tmp_path, mock_hf_data, mock_hf_paths
):
    test_dir = tmp_path / "test_forcing"
    test_dir.mkdir()
    paths = st.subset_forcing(
        ij_bounds=(0, 0, 10, 20),
        grid="conus1",
        start="2005-10-01",
        end=end_date,
        dataset="dataset",
        forcing_vars=("var1",),
        write_dir=test_dir,
    )
    assert len(os.listdir(test_dir)) == num_files


def test_subset_forcing_filenames(tmp_path, mock_hf_data, mock_hf_paths):
    test_dir = tmp_path / "test_forcing"
    test_dir.mkdir()
    paths = st.subset_forcing(
        ij_bounds=(0, 0, 10, 20),
        grid="conus1",
        start="2005-09-01",
        end="2005-10-01",
        dataset="my_ds",
        forcing_vars=("var1",),
        write_dir=test_dir,
    )
    expected_files = [
        f"my_ds.var1.{start:06d}_to_{(start + 23):06d}.pfb"
        for start in range(1, 721, 24)
    ]
    assert [os.path.basename(path) for path in paths["var1"]] == expected_files


@pytest.mark.parametrize("time_zone", ["UTC", "EST", "US/Pacific", "US/Central", "US/Mountain", "US/Eastern"])
def test_subset_forcing_timezones(time_zone, tmp_path, mock_hf_data, mock_hf_paths):
    test_dir = tmp_path / "test_forcing"
    test_dir.mkdir()
    paths = st.subset_forcing(
        ij_bounds=(0, 0, 10, 20),
        grid="conus1",
        start="2005-10-02",
        end="2005-10-03",
        dataset="my_ds",
        forcing_vars=("var1",),
        write_dir=test_dir,
        time_zone=time_zone,
    )
    path = paths["var1"][0]
    assert read_pfb(path).shape == (24, 20, 10)


#####################
# Integration tests #
#####################


@pytest.mark.parametrize("time_zone, utc_offset",
                         [
                             ("EST", 5),
                             ("US/Eastern", 5),
                             ("US/Central", 6),
                             ("US/Mountain", 7),                             
                             ("US/Pacific", 8),
                         ]
)
def test_forcing_timezones(tmp_path, time_zone, utc_offset):
    "Check if we get the correct data offset for all US timezones."
    utc = tmp_path / "UTC_out"
    utc.mkdir()
    tz = tmp_path / "tz_out"
    tz.mkdir()
    ij_bounds = (375, 239, 487, 329)
    grid = "conus1"
    start = "2006-01-03"
    end = "2006-01-04"
    dataset = "NLDAS2"
    st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=utc,
        time_zone="UTC",
        forcing_vars=("air_temp",),
    )
    st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=tz,
        time_zone=time_zone,
        forcing_vars=("air_temp",),        
    )
    utc_temp = read_pfb(os.path.join(utc, "NLDAS.Temp.000001_to_000024.pfb"))
    tz_temp = read_pfb(os.path.join(tz, "NLDAS.Temp.000001_to_000024.pfb"))
    assert np.array_equal(tz_temp[0], utc_temp[utc_offset])


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
