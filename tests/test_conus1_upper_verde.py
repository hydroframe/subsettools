from parflow import Run
from parflow.tools.settings import set_working_directory
from subsettools.subsettools import *
from subsettools.datasets import get_ref_yaml_path
import pytest

# remove this for newer parflow versions
from testutils import pf_test_file


@pytest.fixture(scope="module")
def setup_run(setup_dir_structure):
    run_name = "conus1_upper_verde"
    (
        static_write_dir,
        forcing_dir,
        pf_out_dir,
        correct_output_dir,
        target_runscript,
    ) = setup_dir_structure(run_name)

    huc_list = ["15060202"]
    start = "2005-10-01"
    end = "2005-10-03"
    grid = "conus1"
    run_ds = "conus1_baseline_mod"
    var_ds = "conus1_domain"
    forcing_ds = "NLDAS2"
    reference_run = get_ref_yaml_path(grid, "transient", "solid")

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

    config_clm(
        ij_bounds, start=start, end=end, dataset=run_ds, write_dir=static_write_dir
    )

    subset_forcing(
        ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=forcing_ds,
        write_dir=forcing_dir,
    )

    target_runscript = edit_runscript_for_subset(
        ij_bounds,
        runscript_path=reference_run,
        write_dir=pf_out_dir,
        runname=run_name,
        forcing_dir=forcing_dir,
    )

    copy_static_files(read_dir=static_write_dir, write_dir=pf_out_dir)

    target_runscript = change_filename_values(
        runscript_path=target_runscript, write_dir=pf_out_dir, init_press=init_press_filename
    )

    return run_name, target_runscript, pf_out_dir, correct_output_dir


@pytest.mark.parametrize("P, Q", [(1, 1), (2, 2), (4, 4)])
def test_conus1_upper_verde(setup_run, remove_output_files, P, Q):
    run_name, target_runscript, pf_out_dir, correct_output_dir = setup_run

    # overwrite run_name
    run_name = run_name + f"_{P}_{Q}"
    target_runscript = change_filename_values(
        runscript_path=target_runscript, write_dir=pf_out_dir, runname=run_name
    )

    dist_run(
        P,
        Q,
        runscript_path=target_runscript,
        write_dir=pf_out_dir,
        dist_clim_forcing=True,
    )

    set_working_directory(pf_out_dir)
    run = Run.from_definition(target_runscript)
    run.TimingInfo.StopTime = 10
    run.run(working_directory=pf_out_dir)

    filename = f"/{run_name}.out.perm_x.pfb"
    assert pf_test_file(
        pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_x"
    )

    filename = f"/{run_name}.out.perm_y.pfb"
    assert pf_test_file(
        pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_y"
    )

    filename = f"/{run_name}.out.perm_z.pfb"
    assert pf_test_file(
        pf_out_dir + filename, correct_output_dir + filename, "Max difference in perm_z"
    )

    for i in range(5):
        filename = f"/{run_name}.out.press.0000{i}.pfb"
        assert pf_test_file(
            pf_out_dir + filename,
            correct_output_dir + filename,
            f"Max difference in Pressure for timestep {i}",
        )

        filename = f"/{run_name}.out.satur.0000{i}.pfb"
        assert pf_test_file(
            pf_out_dir + filename,
            correct_output_dir + filename,
            f"Max difference in Saturation for timestep {i}",
        )

    remove_output_files(pf_out_dir)
