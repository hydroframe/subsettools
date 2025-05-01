""""
    Unit tests for subset_forcing function.
"""

# pylint: disable=R0917,C0301,W0621,W0613

import os
import glob
import numpy as np
from parflow.tools.io import read_pfb, read_pfb_sequence
import subsettools as st
import hf_hydrodata
import pytest

_SECONDS_PER_HOUR = 3600
_DUMMY_NZ = 5


@pytest.fixture
def mock_hf_data(monkeypatch):
    """Unit test fixture to use monkeypatch mock calls to get_gridded_data"""

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
        else:
            shape = (1, nj, ni)
        return np.ones(shape)

    monkeypatch.setattr(hf_hydrodata, "get_gridded_data", mock_get_data)


@pytest.fixture
def mock_hf_paths(monkeypatch):
    """Unit test fixture to mock hf_hydrodata.get(paths)"""

    def mock_get_paths(options):
        dataset = options.get("dataset")
        variable = options.get("variable")
        file_path = [f"/folder/{dataset}.{variable}.000000_to_000000.pfb"]
        return file_path

    monkeypatch.setattr(hf_hydrodata, "get_paths", mock_get_paths)


def test_subset_static_filenames(tmp_path, mock_hf_data):
    """Test subset_static with 3 variable names."""

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
    """Test subset_static with precipition."""

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
    """Test subset_press_init"""

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
    """Test subset_press_init and verify numpy array."""

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
    """Test subset_forcing with one forcing variable"""
    test_dir = tmp_path / "test_forcing"
    test_dir.mkdir()
    st.subset_forcing(
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
    """Test subset_forcing and verify returned file names."""

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


@pytest.mark.parametrize(
    "time_zone", ["UTC", "EST", "US/Pacific", "US/Central", "US/Mountain", "US/Eastern"]
)
def test_subset_forcing_timezones(time_zone, tmp_path, mock_hf_data, mock_hf_paths):
    """Test subset forcing with various different timezones."""

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


@pytest.mark.parametrize(
    "time_zone, utc_offset",
    [
        ("EST", 5),
        ("US/Eastern", 5),
        ("US/Central", 6),
        ("US/Mountain", 7),
        ("US/Pacific", 8),
    ],
)
def test_forcing_timezones(tmp_path, time_zone, utc_offset):
    """Check if we get the correct data offset for all US timezones."""

    utc = tmp_path / "UTC_out"
    utc.mkdir()
    tz = tmp_path / "tz_out"
    tz.mkdir()
    ij_bounds = (375, 239, 487, 329)
    grid = "conus1"
    start = "2006-01-03"
    end = "2006-01-04"
    dataset = "NLDAS2"
    paths_utc = st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=utc,
        time_zone="UTC",
        forcing_vars=("air_temp",),
    )
    paths_tz = st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=tz,
        time_zone=time_zone,
        forcing_vars=("air_temp",),
    )
    utc_temp = read_pfb(paths_utc["air_temp"][0])
    tz_temp = read_pfb(paths_tz["air_temp"][0])
    assert np.array_equal(tz_temp[0], utc_temp[utc_offset])


def test_forcing_timezone_time_change(tmp_path):
    """Check for when start and end dates cross the time change (EDT->EST)."""

    write_dir = tmp_path
    ij_bounds = (375, 239, 487, 329)
    grid = "conus1"
    start = "2006-04-02"  # Time change on 02-04-2006 at 2:00am.
    end = "2006-04-03"
    dataset = "NLDAS2"
    paths = st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=write_dir,
        time_zone="US/Eastern",
        forcing_vars=("air_temp",),
    )
    data = read_pfb(paths["air_temp"][0])
    assert data.shape == (24, 90, 112)


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


def test_subset_forcing_dataset_version(tmp_path):
    """Test subset_forcing with dataset version."""

    test_dir_1 = tmp_path / "1.0"
    test_dir_1.mkdir()
    test_dir_latest = tmp_path / "latest"
    test_dir_latest.mkdir()

    paths = {}
    paths["1.0"] = st.subset_forcing(
        ij_bounds=(1225, 1738, 1347, 1811),
        grid="conus2",
        start="2005-10-02",
        end="2005-10-03",
        dataset="CW3E",
        forcing_vars=("air_temp",),
        write_dir=test_dir_1,
        dataset_version="1.0",
    )
    paths["latest"] = st.subset_forcing(
        ij_bounds=(1225, 1738, 1347, 1811),
        grid="conus2",
        start="2005-10-02",
        end="2005-10-03",
        dataset="CW3E",
        forcing_vars=("air_temp",),
        write_dir=test_dir_latest,
        dataset_version=None,
    )

    file_1 = paths["1.0"]["air_temp"][0]
    file_latest = paths["latest"]["air_temp"][0]
    assert np.array_equal(read_pfb(file_1), read_pfb(file_latest))


def test_forcing_conus2_est(tmp_path):
    """Test read_subset_forcing using conus2 for multiple days with a EST timezone."""

    tz = tmp_path / "tz_out"
    tz.mkdir()
    start = "2009-03-01"
    end = "2009-03-04"
    time_zone = "US/Eastern"
    ij_bounds = (813, 1200, 925, 1290)
    grid = "conus2"
    dataset = "CW3E"

    paths_tz = st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=tz,
        time_zone=time_zone,
        forcing_vars=("air_temp",),
    )
    result_np = read_pfb_sequence(paths_tz.get("air_temp"))
    assert len(paths_tz.get("air_temp")) == 3
    assert result_np.shape == (3, 24, 90, 112)
    assert np.nansum(result_np) == pytest.approx(206304706.1352234)


def test_forcing_conus2_utc(tmp_path):
    """Test read_subset_forcing using conus2 for multiple days with a UTC timezone."""

    tz = tmp_path / "tz_out"
    tz.mkdir()
    start = "2009-03-01"
    end = "2009-03-04"
    ij_bounds = (813, 1200, 925, 1290)
    grid = "conus2"
    dataset = "CW3E"

    paths_tz = st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=tz,
        forcing_vars=("air_temp",),
    )
    result_np = read_pfb_sequence(paths_tz.get("air_temp"))
    assert len(paths_tz.get("air_temp")) == 3
    assert result_np.shape == (3, 24, 90, 112)
    assert np.nansum(result_np) == pytest.approx(206204619.5949707)


def test_forcing_all_variables_timezone(tmp_path):
    """Test read_subset_forcing using conus2 for all variables with a timezone."""

    tz = tmp_path / "tz_out"
    tz.mkdir()
    start = "2009-03-01"
    end = "2009-03-04"
    time_zone = "US/Eastern"
    ij_bounds = (813, 1200, 925, 1290)
    grid = "conus2"
    dataset = "CW3E"

    paths_tz = st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=tz,
        time_zone=time_zone,
    )

    assert len(paths_tz) == 8
    result_np = read_pfb_sequence(paths_tz.get("air_temp"))
    assert len(paths_tz.get("air_temp")) == 3
    assert len(paths_tz.get("downward_shortwave")) == 3
    assert result_np.shape == (3, 24, 90, 112)
    assert np.nansum(result_np) == pytest.approx(206304706.1352234)


def test_forcing_full_conus2(tmp_path):
    """
    Test ability to subset tools to download a forcing file file
    with a full conus2 grid.
    This test can only be executed on a machine with access to /hydrodata files
    to be able to compare the result with the raw data files. This would not
    work remotely for older version of subset tools because files were too big.
    """

    raw_file = "/hydrodata/forcing/processed_data/CONUS2/CW3E_v1.0/hourly/WY2009/CW3E.Temp.003625_to_003648.pfb"
    if os.path.exists(raw_file):

        tz = tmp_path / "tz_out"
        tz.mkdir()
        start = "2009-03-01"
        end = "2009-03-02"
        ij_bounds = (0, 0, 4442, 3256)
        grid = "conus2"
        dataset = "CW3E"

        paths_tz = st.subset_forcing(
            ij_bounds=ij_bounds,
            grid=grid,
            start=start,
            end=end,
            dataset=dataset,
            write_dir=tz,
            forcing_vars=("air_temp",),
        )
        result_np = read_pfb_sequence(paths_tz.get("air_temp"))
        raw_np = read_pfb_sequence([raw_file])
        assert len(paths_tz.get("air_temp")) == 1
        assert result_np.shape == (1, 24, 3256, 4442)
        assert np.array_equal(result_np, raw_np)


def test_forcing_1pt(tmp_path):
    """
    Test ability to subset tools to download a forcing file file
    with a 1x1 subgrid size. This caused a bug before.
    """

    tz = tmp_path / "tz_out"
    tz.mkdir()
    start = "2009-10-01"
    end = "2009-10-06"
    ij_bounds = (1000, 1000, 1001, 1001)
    grid = "conus2"
    dataset = "CW3E"

    paths_tz = st.subset_forcing(
        ij_bounds=ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=dataset,
        write_dir=tz,
        forcing_vars=("air_temp",),
    )
    result_np = read_pfb_sequence([paths_tz.get("air_temp")[0]])
    assert len(paths_tz.get("air_temp")) == 5
    assert result_np.shape == (1, 24, 1, 1)
