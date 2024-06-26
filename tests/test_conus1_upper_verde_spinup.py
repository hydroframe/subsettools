from parflow import Run
from parflow.tools.settings import set_working_directory
import subsettools as st
import os
import pathlib

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

    hucs = ["15060202"]
    start = "2005-10-01"
    grid = "conus1"
    run_ds = "conus1_baseline_mod"
    var_ds = "conus1_domain"
    P = 1
    Q = 1
    reference_run = st.get_template_runscript(grid, "spinup", "solid", static_write_dir)

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
    target_runscript = st.edit_runscript_for_subset(
        ij_bounds, runscript_path=reference_run, write_dir=pf_out_dir, runname=run_name
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
        dist_clim_forcing=False,
    )

    set_working_directory(f"{pf_out_dir}")
    run = Run.from_definition(target_runscript)
    run.TimingInfo.StopTime = 48.0
    run.TimingInfo.DumpInterval = 12.0
    run.run(working_directory=pf_out_dir)
    remove_output_files(pf_out_dir)
