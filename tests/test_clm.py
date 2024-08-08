import os
import pytest
import numpy as np
import subsettools as st
import hf_hydrodata
import shutil
from test_subsetting import mock_hf_data

@pytest.fixture
def mock_hf_get_file(monkeypatch):
    def mock_get_raw_file(filepath, **kwargs):
        file_type = kwargs.get("file_type")
        write_dir = filepath
        if file_type == "vegp":
            file_path = os.path.join(os.getcwd(), "tests/clm_mocking_files/drv_vegp.dat")
            shutil.copy(file_path, write_dir)
        elif file_type == "drv_clm":
            file_path = os.path.join(os.getcwd(), "tests/clm_mocking_files/drv_clmin.dat")
            shutil.copy(file_path, write_dir)
    monkeypatch.setattr(hf_hydrodata, "get_raw_file", mock_get_raw_file)



def test_config_clm_filepaths(mock_hf_get_file, mock_hf_data, tmp_path):
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    file_paths = st.config_clm(
    ij_bounds=(0,0,10,10),
    start="2005-10-01",
    end="2006-10-01",
    dataset="conus1_baseline_mod",
    write_dir=write_dir
    )
    assert file_paths["vegp"] == f"{write_dir}/drv_vegp.dat"
    assert file_paths["pfb"] == f"{write_dir}/drv_vegm.dat"
    assert file_paths["drv_clm"] == f"{write_dir}/drv_clmin.dat"

def test_config_clm_data_vegm(mock_hf_get_file, mock_hf_data, tmp_path):
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    file_paths = st.config_clm(
    ij_bounds=(0,0,10,10),
    start="2005-10-01",
    end="2006-10-01",
    dataset="conus1_baseline_mod",
    write_dir=write_dir
    )
    with open(f"{write_dir}/drv_vegm.dat") as my_file:
        file_lines = my_file.readlines()
    assert file_lines[2] == '1 1 1.000000 1.000000 1.00 1.00 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n'
    assert len(file_lines) == 3


@pytest.mark.parametrize(
    "start, end, expected_start, expected_end",
    [
        pytest.param("2005-10-01","2005-10-06", "2005101", "2005105", id="time window of one year"),
        pytest.param("2005-10-01", "2006-10-01", "2005101", "2006930", id="time window of 5 days")
    ],
)
def test_config_clm_data_drv_clm(start, end, expected_start, expected_end, mock_hf_get_file, mock_hf_data, tmp_path):
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    file_paths = st.config_clm(
    ij_bounds=(0,0,10,10),
    start=start,
    end=end,
    dataset="conus1_baseline_mod",
    write_dir=write_dir
    )
    # open file and save lines
    with open(f"{write_dir}/drv_clmin.dat") as my_file:
        file_lines = my_file.readlines()
    # start and end dates are written as three different lines
    # so I use this function to pull out the numbers from each line and put them together in an array
    start_numbers = []
    end_numbers = []
    def pull_numbers(line_no, my_list):
        line=file_lines[line_no]
        stripped = line.replace(" ", "")
        for char in stripped:
            if char.isdigit():
                my_list.append(char)
    # pull out the start date
    pull_numbers(60, start_numbers)
    pull_numbers(59, start_numbers)
    pull_numbers(58, start_numbers)
    # pull out the end date
    pull_numbers(67, end_numbers)
    pull_numbers(66, end_numbers)
    pull_numbers(65, end_numbers)

    # check that it matches up with provided expected values
    assert expected_start == "".join(start_numbers)
    assert expected_end == "".join(end_numbers)