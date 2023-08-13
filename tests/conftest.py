import os
from pathlib import Path
import pytest
from parflow.tools.fs import get_absolute_path, mkdir, rm


@pytest.fixture(scope="module")
def setup_dir_structure():
    home = os.path.expanduser("~")
    base_dir = os.path.join(home, "subsettools_test_output")
    input_dir = os.path.join(base_dir, "inputs")
    output_dir = os.path.join(base_dir, "outputs")
    temp_pf_dir = os.environ.get("PARFLOW_DIR")
    os.environ["PARFLOW_DIR"] = "/home/SHARED/software/parflow/3.10.0"
    correct_output_dir = get_absolute_path(os.path.join("tests", "correct_output"))

    def setup(run_name):
        static_write_dir = os.path.join(input_dir, run_name, "static")
        mkdir(static_write_dir)
        forcing_dir = os.path.join(input_dir, run_name, "forcing")
        mkdir(forcing_dir)
        pf_out_dir = os.path.join(output_dir, run_name)
        mkdir(pf_out_dir)
        target_runscript = os.path.join(pf_out_dir, run_name + ".yaml")
        return (
            static_write_dir,
            forcing_dir,
            pf_out_dir,
            correct_output_dir,
            target_runscript,
        )

    yield setup

    if temp_pf_dir is not None:
        os.environ["PARFLOW_DIR"] = temp_pf_dir
    else:
        del os.environ["PARFLOW_DIR"]

    rm(base_dir)


@pytest.fixture
def remove_output_files():
    def ret_fun(dir_path):
        for p in Path(dir_path).glob("*.out.*"):
            p.unlink()

    return ret_fun
