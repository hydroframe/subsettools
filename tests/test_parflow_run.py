import os
import pytest
import numpy as np
import shutil
import re
from netCDF4 import Dataset
from parflow import Run
from parflow.tools.io import read_pfb, write_pfb
import subsettools as st


_NUM_TIMESTEPS = 28


def test_get_ref_yaml_path(tmp_path):
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    assert runscript == os.path.join(test_dir, "conus1_transient_box.yaml")
    run = Run.from_definition(runscript)
    assert run.ComputationalGrid.NX == 3342
    assert run.GeomInput.domaininput.InputType == "Box"


def test_edit_runscript_for_subset_1(tmp_path):
    """Check the edited fiedls of the runscript file."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    runscript = st.edit_runscript_for_subset(
        (10, 10, 25, 25), runscript, runname="my_new_run", forcing_dir=str(forcing_dir)
    )
    run = Run.from_definition(runscript)
    assert run.get_name() == "my_new_run"
    assert run.Solver.CLM.MetFilePath == str(forcing_dir)
    assert run.ComputationalGrid.NX == 15
    assert run.ComputationalGrid.NY == 15
    assert run.Geom.domain.Upper.X == 15000
    assert run.Geom.domain.Upper.Y == 15000


def test_edit_runscript_for_subset_2(tmp_path):
    """Check the runscript file in write_dir"""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    filename = os.path.basename(runscript)
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    runscript = st.edit_runscript_for_subset(
        (10, 10, 25, 25), runscript, write_dir=write_dir, forcing_dir=str(forcing_dir)
    )
    assert runscript == os.path.join(write_dir, filename)
    assert os.path.exists(runscript)
    run = Run.from_definition(runscript)
    assert run.Solver.CLM.MetFilePath == str(forcing_dir)
    assert run.ComputationalGrid.NX == 15
    assert run.ComputationalGrid.NY == 15
    assert run.Geom.domain.Upper.X == 15000
    assert run.Geom.domain.Upper.Y == 15000


def test_edit_runscript_for_subset_3(tmp_path):
    """Check that exception is raised if forcing_dir is invalid."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    forcing_dir = os.path.join(tmp_path, "forcing")
    with pytest.raises(Exception):
        runscript = st.edit_runscript_for_subset(
            (10, 10, 25, 25), runscript, runname="my_new_run", forcing_dir=forcing_dir
        )


def test_change_filename_values_1(tmp_path):
    """Check the edited fiedls of the runscript file."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    old_runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    test_file = test_dir / "slope_x.pfb"
    test_file.write_bytes(b"Binary file contents")
    new_runscript = st.change_filename_values(
        old_runscript, runname="my_new_run", slopex=str(test_file)
    )
    # check that both files exist in the folder:
    assert os.path.exists(old_runscript)
    assert os.path.exists(new_runscript)
    assert old_runscript != new_runscript
    # check the edited fields of the new runscript:
    run = Run.from_definition(new_runscript)
    assert run.get_name() == "my_new_run"
    assert run.TopoSlopesX.FileName == str(test_file)


def test_change_filename_values_2(tmp_path):
    """Check that the file is replaced if write_dir==None and runname==None"""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    old_runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    test_file = test_dir / "slope_x.pfb"
    test_file.write_bytes(b"Binary file contents")
    new_runscript = st.change_filename_values(old_runscript, slopex=str(test_file))
    # check that old file has been replaced:
    assert os.path.exists(new_runscript)
    assert old_runscript == new_runscript
    # check the edited fields of the new runscript:
    run = Run.from_definition(new_runscript)
    assert run.TopoSlopesX.FileName == str(test_file)


@pytest.fixture
def setup_dummy_run(tmp_path):
    def setup(clm_on=False):
        old_dir = tmp_path / "old"
        old_dir.mkdir()
        slope_x = old_dir / "slope_x.pfb"
        slope_x.write_text("Dummy slope_x input file")
        slope_y = old_dir / "slope_y.pfb"
        slope_y.write_text("Dummy slope_y input file")
        ic_pressure = old_dir / "old_ic_pressure.pfb"
        ic_pressure.write_text("Dummy ic_pressure input file")
        run = Run("dummy_run")
        run.TopoSlopesX.FileName = os.path.basename(slope_x)
        run.TopoSlopesY.FileName = os.path.basename(slope_y)
        run.GeomInput.Names = "domain_input"
        run.GeomInput.domain_input.InputType = "Box"
        run.GeomInput.domain_input.GeomName = "domain"
        run.Geom.domain.ICPressure.FileName = os.path.basename(ic_pressure)
        _create_nc_output(old_dir)
        _create_pfb_output(old_dir)
        if clm_on:
            run.Solver.LSM = "CLM"
            _write_clm_restart_files(old_dir)
            _write_clm_driver_files(old_dir)
        runscript_path, _ = run.write(file_format="yaml", working_directory=old_dir)
        return runscript_path

    return setup


def _write_clm_driver_files(working_directory):
    file_path = os.path.join(working_directory, "drv_vegm.dat")
    with open(file_path, "w") as file:
        file.write("Dummy vegm file.")
    file_path = os.path.join(working_directory, "drv_vegp.dat")
    with open(file_path, "w") as file:
        file.write("Dummy vegp file.")
    # Copy real drv_clmin.dat for testing.
    file_path = os.path.join("tests", "correct_output", "drv_clmin.dat")
    shutil.copy(file_path, working_directory)


def _write_clm_restart_files(working_directory):
    file_path = os.path.join(working_directory, "clm_restart.tcl")
    restart_time = _NUM_TIMESTEPS - _NUM_TIMESTEPS % 24
    with open(file_path, "w") as file:
        file.write(f"set istep           {restart_time}")
    file_path = os.path.join(working_directory, "clm.rst.00000.0")
    with open(file_path, "w") as file:
        file.write("This is a dummy CLM restart file.")


def _create_nc_output(working_directory):
    value = 0
    num_days = (_NUM_TIMESTEPS + 23) // 24
    for idx in range(num_days):
        ncfile = Dataset(
            os.path.join(working_directory, f"dummy_run.out.{idx:05d}.nc"),
            "w",
            format="NETCDF4",
        )
        ncfile.createDimension("time", None)
        ncfile.createDimension("x", 5)
        ncfile.createDimension("y", 3)
        ncfile.createDimension("z", 1)
        pressure = ncfile.createVariable("pressure", "f8", ("time", "z", "y", "x"))
        for i in range(24):
            pressure[i, :, :, :] = value
            value = value + 1
            if value >= _NUM_TIMESTEPS:
                break
        ncfile.close()


def _create_pfb_output(working_directory):
    for idx in range(_NUM_TIMESTEPS):
        data = np.ones((1, 3, 5)) * idx
        file_path = os.path.join(
            working_directory, f"dummy_run.out.press.{idx:05d}.pfb"
        )
        write_pfb(file_path, data)


def test_restart_run_new_directory(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run()
    new_dir = tmp_path / "new"
    assert not os.path.isdir(new_dir)
    st.restart_run(runscript_path=old_runscript, new_dir=new_dir)
    assert os.path.isdir(new_dir)


def test_restart_run_valid_runscript(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run()
    new_dir = tmp_path / "new"
    new_runscript = st.restart_run(
        runscript_path=old_runscript,
        new_dir=new_dir,
    )
    try:
        run = Run.from_definition(new_runscript)
    except Exception as e:
        pytest.fail("Failed to load new run from returned runscript.")


@pytest.mark.parametrize("dir_name", [None, "new"])
def test_restart_run_runscript_path(dir_name, setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run()
    new_dir = None if dir_name is None else tmp_path / dir_name
    new_runscript = st.restart_run(
        runscript_path=old_runscript,
        new_dir=new_dir,
    )
    working_directory = (
        os.path.dirname(old_runscript) if new_dir is None else str(new_dir)
    )
    assert os.path.dirname(new_runscript) == working_directory


@pytest.mark.parametrize("runname", [None, "new_run"])
def test_restart_run_filename(runname, setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run()
    new_runscript = st.restart_run(
        runscript_path=old_runscript,
        runname=runname,
    )
    filename = "dummy_run.yaml" if runname is None else "new_run.yaml"
    assert os.path.basename(new_runscript) == filename


def test_restart_run_copy_inputs(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run()
    new_dir = tmp_path / "new"
    st.restart_run(runscript_path=old_runscript, new_dir=new_dir)
    assert "slope_x.pfb" in os.listdir(new_dir)
    assert "slope_y.pfb" in os.listdir(new_dir)
    assert "old_ic_pressure.pfb" not in os.listdir(new_dir)


def test_restart_run_forcing_directory(setup_dummy_run):
    old_runscript = setup_dummy_run()
    runscript_path = st.restart_run(
        runscript_path=old_runscript,
        forcing_dir="forcing",
    )
    run = Run.from_definition(runscript_path)
    assert run.Solver.CLM.MetFilePath == "forcing"


@pytest.mark.parametrize(
    "output_type, clm_on",
    [("netcdf", False), ("netcdf", True), ("pfb", False), ("pfb", True)],
)
def test_restart_run_initial_pressure(output_type, clm_on, setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run(clm_on)
    new_dir = tmp_path / "new"
    _ = st.restart_run(
        runscript_path=old_runscript,
        new_dir=new_dir,
        output_type=output_type,
    )
    ic_data = read_pfb(os.path.join(new_dir, "initial_pressure.pfb"))
    assert ic_data.shape == (1, 3, 5)
    assert ic_data[0, 0, 0] == 24 if clm_on else _NUM_TIMESTEPS - 1


@pytest.mark.parametrize(
    "clm_on, create_new_dir, result",
    [
        (True, True, (0, 0, 1)),
        (True, False, (24, 24, 25)),
        (False, True, (0, 0, 1)),
        (False, False, (27, 27, 1)),
    ],
)
def test_restart_run_restart_timestep(
    clm_on, create_new_dir, result, setup_dummy_run, tmp_path
):
    old_runscript = setup_dummy_run(clm_on)
    new_dir = (tmp_path / "new") if create_new_dir else None
    new_runscript = st.restart_run(
        runscript_path=old_runscript,
        new_dir=new_dir,
    )
    run = Run.from_definition(new_runscript)
    assert (
        run.TimingInfo.StartCount,
        run.TimingInfo.StartTime,
        run.Solver.CLM.IstepStart,
    ) == result


def test_restart_run_stoptime(setup_dummy_run):
    old_runscript = setup_dummy_run()
    stop_time = _NUM_TIMESTEPS + 24
    new_runscript = st.restart_run(runscript_path=old_runscript, stop_time=stop_time)
    run = Run.from_definition(new_runscript)
    assert run.TimingInfo.StopTime == stop_time


def test_restart_run_stoptime_errors(setup_dummy_run):
    old_runscript = setup_dummy_run()
    stop_time = _NUM_TIMESTEPS - 2
    with pytest.raises(ValueError):
        st.restart_run(runscript_path=old_runscript, stop_time=stop_time)


def test_restart_run_copy_restart_files(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run(clm_on=True)
    new_dir = tmp_path / "new"
    new_runscript = st.restart_run(runscript_path=old_runscript, new_dir=new_dir)
    assert os.path.exists(os.path.join(new_dir, "clm.rst.00000.0"))


def test_restart_run_rename_clm_restart_files(setup_dummy_run):
    old_runscript = setup_dummy_run(clm_on=True)
    working_directory = os.path.dirname(old_runscript)
    new_runscript = st.restart_run(runscript_path=old_runscript)
    clm_restart = _NUM_TIMESTEPS - _NUM_TIMESTEPS % 24
    clm_restart = str(clm_restart).rjust(5, '0')
    assert os.path.exists(os.path.join(working_directory, f"clm.rst.{clm_restart}.0"))


@pytest.mark.parametrize(
    "stop_time, end_day, end_hour",
    [
        (24, 2, 0),
        (28, 2, 4),
        (48, 3, 0),
        (50, 3, 2)
    ],
)
def test_restart_run_drvclm(setup_dummy_run, stop_time, end_day, end_hour):
    stop_time = _NUM_TIMESTEPS + 1
    old_runscript = setup_dummy_run()
    working_directory = os.path.dirname(old_runscript)
    file_path = os.path.join("tests", "correct_output", "drv_clmin.dat")
    shutil.copy(file_path, working_directory)
    _ = st.restart_run(runscript_path=old_runscript, stop_time=stop_time)
    with open(file_path, 'r') as file:
        content = file.read()
    startcode = int(re.search(r'(startcode\s+)(\d{1})', content).group(2))
    start_hour = int(re.search(r'(shr\s+)(\d{1,2})', content).group(2))
    eda = int(re.search(r'(eda\s+)(\d{1,2})', content).group(2))
    ehr = int(re.search(r'(ehr\s+)(\d{1,2})', content).group(2))
    assert startcode == 1
    assert start_hour == 0
    assert (eda, ehr) == (end_day, end_hour)


def test_restart_run_copy_clm_files(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run()
    run = Run.from_definition(old_runscript)
    run.Solver.LSM = "CLM"
    runscript_path, _ = run.write(file_format="yaml", working_directory=os.path.dirname(old_runscript))
    new_dir = tmp_path / "new"
    _ = st.restart_run(runscript_path=runscript_path, new_dir=new_dir)
    assert os.path.exists(os.path.join(new_dir, "drv_clmin.dat"))
    assert os.path.exists(os.path.join(new_dir, "drv_vegm.dat"))
    assert os.path.exists(os.path.join(new_dir, "drv_vegp.dat"))
    
