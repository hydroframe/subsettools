import pytest
from parflow import Run
from parflow.tools.settings import set_working_directory
import subsettools as st
import os
import pathlib

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

    hucs = ["15060202"]
    start = "2005-10-01"
    end = "2005-10-03"
    grid = "conus1"
    run_ds = "conus1_baseline_mod"
    var_ds = "conus1_domain"
    forcing_ds = "NLDAS2"
    reference_run = st.get_template_runscript(
        grid, "transient", "solid", static_write_dir
    )

    ij_bounds, mask = st.define_huc_domain(hucs=hucs, grid=grid)
    assert ij_bounds == (375, 239, 487, 329)

    st.write_mask_solid(mask=mask, grid=grid, write_dir=static_write_dir)

    st.subset_static(
        ij_bounds,
        dataset=var_ds,
        write_dir=static_write_dir,
        var_list=(
            "slope_x",
            "slope_y",
            "pf_indicator",
            "pme",
            "ss_pressure_head",
        ),
    )

    init_press_filepath = st.subset_press_init(
        ij_bounds,
        dataset=run_ds,
        date=start,
        write_dir=static_write_dir,
        time_zone="UTC",
    )
    for path in pathlib.Path(static_write_dir).glob("*.pfb"):
        filename = os.path.basename(path)
        print(">>>>> checking: " + filename)
        assert pf_test_file(
            os.path.join(static_write_dir, filename),
            os.path.join(correct_output_dir, filename),
            "Max difference in " + filename,
        )

    st.config_clm(
        ij_bounds, start=start, end=end, dataset=run_ds, write_dir=static_write_dir
    )

    st.subset_forcing(
        ij_bounds,
        grid=grid,
        start=start,
        end=end,
        dataset=forcing_ds,
        write_dir=forcing_dir,
    )
    for filename in os.listdir(forcing_dir):
        print(">>>>> checking: " + filename)
        assert pf_test_file(
            os.path.join(forcing_dir, filename),
            os.path.join(correct_output_dir, filename),
            "Max difference in " + filename,
        )

    target_runscript = st.edit_runscript_for_subset(
        ij_bounds,
        runscript_path=reference_run,
        runname=run_name,
        forcing_dir=forcing_dir,
    )

    st.copy_files(read_dir=static_write_dir, write_dir=pf_out_dir)

    target_runscript = st.change_filename_values(
        runscript_path=target_runscript,
        init_press=os.path.basename(init_press_filepath),
    )

    return run_name, target_runscript, pf_out_dir, correct_output_dir


@pytest.mark.parametrize("P, Q", [(1, 1), (2, 2), (4, 4)])
def test_conus1_upper_verde(setup_run, remove_output_files, P, Q):
    run_name, target_runscript, pf_out_dir, correct_output_dir = setup_run

    # overwrite run_name
    run_name = run_name + f"_{P}_{Q}"
    target_runscript = st.change_filename_values(
        runscript_path=target_runscript, write_dir=pf_out_dir, runname=run_name
    )

    st.dist_run(
        P,
        Q,
        runscript_path=target_runscript,
        dist_clim_forcing=True,
    )

    set_working_directory(pf_out_dir)
    run = Run.from_definition(target_runscript)
    run.TimingInfo.StopTime = 10
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
