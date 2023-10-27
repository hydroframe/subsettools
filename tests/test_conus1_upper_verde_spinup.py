from parflow import Run
from parflow.tools.settings import set_working_directory
from subsettools.subsettools import *
from subsettools.datasets import get_ref_yaml_path

from testutils import pf_test_file


def test_conus1_upper_verde_spinup(setup_dir_structure, remove_output_files):
    run_name = "conus1_upper_verde_spinup"
    (
        static_write_dir,
        _,
        pf_out_dir,
        correct_output_dir,
        target_runscript,
    ) = setup_dir_structure(run_name)

    huc_list = ["15060202"]
    start = "2005-10-01"
    grid = "conus1"
    run_ds = "conus1_baseline_mod"
    var_ds = "conus1_domain"
    P = 1
    Q = 1
    reference_run = get_ref_yaml_path(grid, "spinup", "solid", static_write_dir)

    ij_bounds = huc_to_ij(huc_list=huc_list, grid=grid)
    create_mask_solid(huc_list=huc_list, grid=grid, write_dir=static_write_dir)
    subset_static(ij_bounds, dataset=var_ds, write_dir=static_write_dir)
    init_press_filename = subset_press_init(
        ij_bounds,
        dataset=run_ds,
        date=start,
        write_dir=static_write_dir,
        time_zone="UTC",
    )
    target_runscript = edit_runscript_for_subset(
        ij_bounds, runscript_path=reference_run, write_dir=pf_out_dir, runname=run_name
    )
    copy_files(read_dir=static_write_dir, write_dir=pf_out_dir)
    target_runscript = change_filename_values(
        runscript_path=target_runscript,
        write_dir=pf_out_dir,
        init_press=init_press_filename,
    )
    dist_run(
        P=P,
        Q=Q,
        runscript_path=target_runscript,
        dist_clim_forcing=False,
    )

    set_working_directory(f"{pf_out_dir}")
    run = Run.from_definition(target_runscript)
    run.TimingInfo.StopTime = 48.0
    run.TimingInfo.DumpInterval = 12.0
    run.run(working_directory=pf_out_dir)

    test_vars = ["perm_x", "perm_y", "perm_z"]
    for var in test_vars:
        filename = f"{run_name}.out.{var}.pfb"
        assert pf_test_file(
            os.path.join(pf_out_dir, filename),
            os.path.join(correct_output_dir, filename),
            "Max difference in perm_x",
        )
    for i in range(5):
        filename = f"{run_name}.out.press.0000{i}.pfb"
        assert pf_test_file(
            os.path.join(pf_out_dir, filename),
            os.path.join(correct_output_dir, filename),
            f"Max difference in Pressure for timestep {i}",
        )
        filename = f"{run_name}.out.satur.0000{i}.pfb"
        assert pf_test_file(
            os.path.join(pf_out_dir, filename),
            os.path.join(correct_output_dir, filename),
            f"Max difference in Saturation for timestep {i}",
        )

    remove_output_files(pf_out_dir)
