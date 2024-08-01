import os
import pytest
import numpy as np
import subsettools as st
import hf_hydrodata
import shutil

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


@pytest.fixture
def mock_hf_get_data(monkeypatch):
    def mock_get_gridded_data(options):
        return np.ones((10,10,10))
    monkeypatch.setattr(hf_hydrodata, "get_gridded_data", mock_get_gridded_data)

def test_config_clm_filepaths(mock_hf_get_file, mock_hf_get_data, tmp_path):
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

def test_config_clm_data
    

    
