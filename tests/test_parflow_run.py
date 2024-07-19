import os
import pytest
from parflow import Run
import subsettools as st


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
    old_dir = tmp_path / "old"
    old_dir.mkdir()
    slope_x = old_dir / "slope_x.pfb"
    slope_x.write_text("Dummy slope_x input file")
    slope_y = old_dir / "slope_y.pfb"
    slope_y.write_text("Dummy slope_y input file")
    run = Run("dummy_run")
    run.TopoSlopesX.FileName = str(slope_x)
    run.TopoSlopesY.FileName = str(slope_y)
    runscript_path, _ = run.write(
        file_format="yaml",
        working_directory=old_dir
    )
    return runscript_path
    

def test_restart_run_new_directory(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run
    new_dir = tmp_path / "new"
    assert not os.path.isdir(new_dir)
    st.restart_run(runscript_path=old_runscript, new_dir=new_dir)
    assert os.path.isdir(new_dir)

    
def test_restart_run_valid_runscript(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run
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
    old_runscript = setup_dummy_run
    new_dir = None if dir_name is None else tmp_path / dir_name
    new_runscript = st.restart_run(
        runscript_path=old_runscript,
        new_dir=new_dir,
    )
    working_directory = os.path.dirname(old_runscript) if new_dir is None else str(new_dir)
    assert os.path.dirname(new_runscript) == working_directory


@pytest.mark.parametrize("new_runname", [None, "new_run"])
def test_restart_run_same_name(new_runname, setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run
    new_runscript = st.restart_run(
        runscript_path=old_runscript,
        new_runname=new_runname,
    )
    filename = "dummy_run.yaml" if new_runname is None else "new_run.yaml"
    assert os.path.basename(new_runscript) == filename


def test_restart_run_copy_inputs(setup_dummy_run, tmp_path):
    old_runscript = setup_dummy_run
    new_dir = tmp_path / "new"
    st.restart_run(
        runscript_path=old_runscript,
        new_dir=new_dir
    )
    assert "slope_x.pfb" in os.listdir(new_dir)
    assert "slope_y.pfb" in os.listdir(new_dir)
