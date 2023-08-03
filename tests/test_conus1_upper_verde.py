from datetime import datetime, timedelta
import numpy as np
import os
from parflow import Run
from parflow.tools.io import read_pfb, write_pfb, read_clm, read_pfb_sequence
from parflow.tools.fs import cp, mkdir, get_absolute_path
from parflow.tools.settings import set_working_directory
from subsettools.subsettools import *
from subsettools.datasets import *
import pathlib
import shutil
import xarray as xr
import hydrodata.data_catalog.data_access
from pathlib import Path

from testutils import pf_test_file

home = os.path.expanduser("~")
run_name = "conus1_upper_verde"
subset_target = "15060202"
start = "2005-10-01" 
end = "2005-10-03"
grid = "conus1"  
run_ds = "conus1_baseline_mod"
var_ds = "conus1_domain"
forcing_ds = "NLDAS2"
base_dir = f"{home}/subsettools_test_output/"
input_dir = base_dir + "inputs/"
output_dir = base_dir + "outputs/"
static_write_dir = input_dir + run_name + "/static/"
mkdir(static_write_dir) 
forcing_dir = input_dir + run_name + "/forcing/"
mkdir(forcing_dir)
pf_out_dir = output_dir + run_name + "/"
mkdir(pf_out_dir) 
os.environ["PARFLOW_DIR"] = "/home/SHARED/software/parflow/3.10.0"
reference_run = get_ref_yaml_path(grid, "transient", "solid")
target_runscript = pf_out_dir + run_name + '.yaml'
correct_output_dir = get_absolute_path('tests/correct_output')

ij_bounds = get_conus_ij(domain = subset_target, grid = grid)

create_mask_solid(huc_id = subset_target, grid = grid, write_dir = static_write_dir)

subset_static(ij_bounds, dataset = var_ds, write_dir = static_write_dir)

subset_press_init(ij_bounds, dataset = run_ds, date = start, write_dir = static_write_dir, time_zone = 'UTC')

config_clm(ij_bounds, start = start, end = end, dataset = run_ds, write_dir = static_write_dir)

subset_forcing(ij_bounds, grid = grid, start = start, end = end, dataset = forcing_ds, write_dir = forcing_dir)

edit_runscript_for_subset(ij_bounds, runscript_path = reference_run, write_dir = pf_out_dir, runname = run_name, forcing_dir = forcing_dir)

copy_static_files(static_input_dir = static_write_dir, pf_dir = pf_out_dir)

change_filename_values(runscript_path = target_runscript, ip = "conus1_baseline_mod_2005.09.30:23.00.00_UTC0_press.pfb")

dist_run(P=4, Q=4, runscript_path = target_runscript, write_dir = pf_out_dir, dist_clim_forcing = True)


set_working_directory(f'{pf_out_dir}')
run = Run.from_definition(target_runscript)
run.TimingInfo.StopTime = 10
run.run(working_directory=pf_out_dir)

def test_conus1_upper_verde():
    filename = f"/{run_name}.out.perm_x.pfb"
    assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_x")

    filename = f"/{run_name}.out.perm_y.pfb"
    assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_y")

    filename = f"/{run_name}.out.perm_z.pfb"
    assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_z")


    timesteps = ["00000", "00001", "00002", "00003", "00004", "00005"]
    for i in timesteps:
        filename = f"/{run_name}.out.press.{i}.pfb"
        assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, f"Max difference in Pressure for timestep {i}")

        filename = f"/{run_name}.out.satur.{i}.pfb"
        assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, f"Max difference in Saturation for timestep {i}")

