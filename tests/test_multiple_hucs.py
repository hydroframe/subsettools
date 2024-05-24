from parflow import Run
from parflow.tools.settings import set_working_directory
import subsettools as st
import pathlib
import os

from testutils import pf_test_file


def test_multiple_hucs(setup_dir_structure, remove_output_files):
    run_name = "conus1_multiple_hucs"
    (
        static_write_dir,
        forcing_dir,
        pf_out_dir,
        correct_output_dir,
        target_runscript,
    ) = setup_dir_structure(run_name)

    hucs = [
        "1302020701",
        "1302020702",
        "1408010604",
        "1408010608",
        "1502000601",
        "1502000602",
    ]
    start = "2005-11-14"
    end = "2005-11-16"
    grid = "conus1"
    run_ds = "conus1_baseline_mod"
    var_ds = "conus1_domain"
    forcing_ds = "NLDAS2"
    P = 1
    Q = 1
    reference_run = st.get_template_runscript(
        grid, "transient", "solid", static_write_dir
    )

    ij_bounds, mask = st.define_huc_domain(hucs=hucs, grid=grid)
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
    target_runscript = st.edit_runscript_for_subset(
        ij_bounds,
        runscript_path=reference_run,
        write_dir=pf_out_dir,
        runname=run_name,
        forcing_dir=forcing_dir,
    )
    st.copy_files(read_dir=static_write_dir, write_dir=pf_out_dir)
    target_runscript = st.change_filename_values(
        runscript_path=target_runscript,
        write_dir=pf_out_dir,
        init_press=os.path.basename(init_press_filepath),
    )
    st.dist_run(
        topo_p=P,
        topo_q=Q,
        runscript_path=target_runscript,
        working_dir=pf_out_dir,
        dist_clim_forcing=True,
    )

    set_working_directory(f"{pf_out_dir}")
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
    for i in range(11):
        timestep = str(i).rjust(5, "0")
        filename = f"{run_name}.out.press.{timestep}.pfb"
        assert pf_test_file(
            os.path.join(pf_out_dir, filename),
            os.path.join(correct_output_dir, filename),
            f"Max difference in Pressure for timestep {timestep}",
        )
        filename = f"{run_name}.out.satur.{timestep}.pfb"
        assert pf_test_file(
            os.path.join(pf_out_dir, filename),
            os.path.join(correct_output_dir, filename),
            f"Max difference in Saturation for timestep {timestep}",
        )

    remove_output_files(pf_out_dir)
