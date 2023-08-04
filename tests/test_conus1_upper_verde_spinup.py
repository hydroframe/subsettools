from parflow import Run
from parflow.tools.settings import set_working_directory
from subsettools.subsettools import *
from subsettools.datasets import get_ref_yaml_path

from testutils import pf_test_file

def test_conus1_upper_verde_spinup(setup):
    run_name = "conus1_upper_verde_spinup"
    static_write_dir, forcing_dir, pf_out_dir, correct_output_dir, target_runscript = setup(run_name)

    subset_target = "15060202"
    start = "2005-10-01" 
    wy = 2006
    grid = "conus1"  
    run_ds = "conus1_baseline_mod"
    var_ds = "conus1_domain"
    P = 1 
    Q = 1
    reference_run = get_ref_yaml_path(grid, "spinup", "solid")

    ij_bounds = get_conus_ij(domain = subset_target, grid = grid)
    create_mask_solid(huc_id = subset_target, grid = grid, write_dir = static_write_dir)
    subset_static(ij_bounds, dataset = var_ds, write_dir = static_write_dir)
    subset_press_init(ij_bounds, dataset = run_ds, date = start, write_dir = static_write_dir, time_zone = 'UTC')
    edit_runscript_for_subset(ij_bounds, runscript_path = reference_run, write_dir = pf_out_dir, runname = run_name)
    copy_static_files(static_input_dir = static_write_dir, pf_dir = pf_out_dir)
    change_filename_values(runscript_path = target_runscript, ip = "conus1_baseline_mod_2005.09.30:23.00.00_UTC0_press.pfb")
    dist_run(P=P, Q=Q, runscript_path = target_runscript, write_dir = pf_out_dir, dist_clim_forcing = False)

    set_working_directory(f'{pf_out_dir}')
    run = Run.from_definition(target_runscript)
    run.TimingInfo.StopTime = 48.0 
    run.TimingInfo.DumpInterval = 12.0
    run.run(working_directory=pf_out_dir)

    filename = f"/{run_name}.out.perm_x.pfb"
    assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_x")
    filename = f"/{run_name}.out.perm_y.pfb"
    assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_y")
    filename = f"/{run_name}.out.perm_z.pfb"
    assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_z")


    for i in range(5):
        filename = f"/{run_name}.out.press.0000{i}.pfb"
        assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, f"Max difference in Pressure for timestep {i}")
        filename = f"/{run_name}.out.satur.0000{i}.pfb"
        assert pf_test_file(pf_out_dir + filename, correct_output_dir + filename, f"Max difference in Saturation for timestep {i}")
