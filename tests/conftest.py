import pytest
import os
from parflow.tools.fs import get_absolute_path, mkdir, rm

@pytest.fixture
def setup_dir_structure():
    home = os.path.expanduser("~")
    base_dir = f"{home}/subsettools_test_output/"
    input_dir = base_dir + "inputs/"
    output_dir = base_dir + "outputs/"
    temp_pf_dir = os.environ.get("PARFLOW_DIR")
    os.environ["PARFLOW_DIR"] = "/home/SHARED/software/parflow/3.10.0"
    correct_output_dir = get_absolute_path('tests/correct_output')

    def setup(run_name):
        static_write_dir = input_dir + run_name + "/static/"
        mkdir(static_write_dir)
        forcing_dir = input_dir + run_name + "/forcing/"
        mkdir(forcing_dir)
        pf_out_dir = output_dir + run_name + "/"
        mkdir(pf_out_dir)
        target_runscript = pf_out_dir + run_name + '.yaml'
        return  static_write_dir, forcing_dir, pf_out_dir, correct_output_dir, target_runscript

    yield setup
    
    if temp_pf_dir is not None:
        os.environ["PARFLOW_DIR"]
    else:
        del os.environ["PARFLOW_DIR"]
        
    rm(base_dir)
