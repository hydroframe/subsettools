import os
import pytest
from parflow import Run
import subsettools as st
import filecmp


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
    assert os.path.exists(runscript_path)
    assert runscript_path == os.path.join(test_dir, os.listdir(test_dir)[0])


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


def test_edit_runscript_filepaths(tmp_path):
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    runscript_path = st.get_template_runscript("conus1", "transient", "box",test_dir)
    ij_bounds = (10, 10, 25, 25)
    runname = "test-edit-runscript"

    new_runscript_path = st.edit_runscript_for_subset(
        ij_bounds=ij_bounds, 
        runscript_path=runscript_path, 
        write_dir=write_dir,
        runname=runname, 
        forcing_dir=forcing_dir,
    )
    assert os.path.exists(new_runscript_path)
    assert os.path.join(write_dir, os.listdir(write_dir)[0]) == new_runscript_path

@pytest.mark.parametrize(
    "runname, ij_bounds",
    [
        pytest.param("new_runscript", (10,10,25,25)),
        pytest.param(None, (10,10,25,23))
    ]
)
def test_edit_runscript_data(runname, ij_bounds, tmp_path):
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    old_runscript_path = st.get_template_runscript("conus1", "transient", "box", test_dir)
    new_runscript_path = st.edit_runscript_for_subset(
        ij_bounds=ij_bounds,
        runscript_path=old_runscript_path,
        write_dir=write_dir,
        runname=runname,
        forcing_dir=str(forcing_dir),
    )
    i_min, j_min, i_max, j_max = ij_bounds
    nx, ny = (i_max-i_min, j_max-j_min)
    new_run = Run.from_definition(new_runscript_path)
    old_run = Run.from_definition(old_runscript_path)

    if runname:
        assert new_run.get_name() == runname
    else:
        assert new_run.get_name() == old_run.get_name()
    assert new_run.ComputationalGrid.NX == nx
    assert new_run.ComputationalGrid.NY == ny
    assert new_run.Geom.domain.Upper.X == nx*1000
    assert new_run.Geom.domain.Upper.Y == ny*1000


def test_copy_single_file(tmp_path):
    read_dir = tmp_path / "read"
    read_dir.mkdir()
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    runscript_path = st.get_template_runscript(
        grid="conus1",
        mode="transient",
        input_file_type="box",
        write_dir=read_dir,
    )
    st.copy_files(read_dir, write_dir)
    f1 = runscript_path
    f2 = f"{write_dir}/{os.path.basename(runscript_path)}"
    assert os.path.exists(f2)
    assert filecmp.cmp(f1, f2, shallow=False)


def test_copy_multiple_files(tmp_path):
    write_dir = tmp_path / "write"
    write_dir.mkdir()
    correct_output_dir = os.path.join(
        os.getcwd(), "tests", "correct_output", "conus1_upper_verde"
    )
    st.copy_files(read_dir=correct_output_dir, write_dir=write_dir)
    match, mismatch, errors = filecmp.cmpfiles(
        write_dir, correct_output_dir, ["mask.pfb", "solidfile.pfsol", "drv_clmin.dat"], shallow=True
    )
    assert len(match) == 3
    assert len(mismatch) == 0
    assert len(errors) == 0
    


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
