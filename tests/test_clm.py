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
    

    
