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


@pytest.mark.parametrize(
    "grid, mode, input_file_type",
    [
        pytest.param("conus1", "spinup", "box"),
        pytest.param("conus2", "spinup", "box"),
        pytest.param("conus1", "transient", "box"),
        pytest.param("conus2", "transient", "box"),
        pytest.param("conus1", "spinup", "solid"),
        pytest.param("conus2", "spinup", "solid"),
        pytest.param("conus1", "transient", "solid"),
        pytest.param("conus2", "transient", "solid"),

    ]
)
def test_get_template_runscript_filepaths(grid, mode, input_file_type, tmp_path):
    test_dir = tmp_path / "test_runscript"
    test_dir.mkdir()
    runscript_path = st.get_template_runscript(
        grid=grid,
        mode=mode,
        input_file_type=input_file_type,
        write_dir=test_dir,
)
    output=f"{test_dir}/{grid}_{mode}_{input_file_type}.yaml"
    print(os.path.exists(runscript_path))
    assert runscript_path == output
    assert os.path.exists(runscript_path)


@pytest.mark.parametrize(
    "grid, mode, input_file_type",
    [
        pytest.param("conus1", "spinup", "box"),
        pytest.param("conus2", "spinup", "box"),
        pytest.param("conus1", "transient", "box"),
        pytest.param("conus2", "transient", "box"),
        pytest.param("conus1", "spinup", "solid"),
        pytest.param("conus2", "spinup", "solid"),
        pytest.param("conus1", "transient", "solid"),
        pytest.param("conus2", "transient", "solid"),

    ]
)
def test_get_template_runscript_data(grid, mode, input_file_type, tmp_path):
    test_dir = tmp_path / "test_runscript_data"
    test_dir.mkdir()
    runscript_path = st.get_template_runscript(
    grid=grid,
    mode=mode,
    input_file_type=input_file_type,
    write_dir=test_dir
)
    run = Run.from_definition(runscript_path)

    #check mode
    if mode == "transient":
        assert run.Solver.LSM == "CLM"
    else:
        assert run.Solver.LSM == "none"
    #check grid
    if grid == "conus1":
        assert run.ComputationalGrid.NX == 3342
        assert run.ComputationalGrid.NY == 1888
    else:
        assert run.ComputationalGrid.NX == 4442
        assert run.ComputationalGrid.NY == 3256
    #check input_file_type
    if input_file_type == "box":
        assert run.Geom.domain.Upper.X == run.ComputationalGrid.NX*1000.0
        assert run.Geom.domain.Upper.Y == run.ComputationalGrid.NY*1000.0
    else:
        assert run.Geom.domain.Upper.X == None
        assert run.Geom.domain.Upper.Y == None


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
