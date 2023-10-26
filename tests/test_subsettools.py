import os
import pytest
from subsettools import subsettools
from subsettools.datasets import get_ref_yaml_path
from hf_hydrodata.grid import to_latlon
from parflow import Run
from parflow.tools.fs import mkdir, rm

@pytest.mark.parametrize(
    "huc_list, grid, result", [(["15060202"], "conus1", (375, 239, 487, 329))]
)
def test_huc_to_ij(huc_list, grid, result):
    assert subsettools.huc_to_ij(huc_list, grid) == result


def test_latlon_to_ij():
    latlon_points = to_latlon("conus1", *[375, 239, 487, 329])
    print(latlon_points)
    latlon_points = [latlon_points[:2], latlon_points[2:]]
    print(latlon_points)
    assert subsettools.latlon_to_ij(latlon_points, "conus1") == (375, 239, 487, 329)

    
def test_edit_runscript_for_subset_1(tmp_path):
    """Check the edited fiedls of the runscript file."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = get_ref_yaml_path("conus1", "transient", "box", test_dir)
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    runscript = subsettools.edit_runscript_for_subset((10, 10, 25, 25),
                                                      runscript,
                                                      runname="my_new_run",
                                                      forcing_dir=str(forcing_dir)
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
    runscript = get_ref_yaml_path("conus1", "transient", "box", test_dir)
    filename = os.path.basename(runscript)
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    runscript = subsettools.edit_runscript_for_subset((10, 10, 25, 25),
                                                      runscript,
                                                      write_dir=write_dir,
                                                      forcing_dir=str(forcing_dir)
    )
    assert runscript == os.path.join(write_dir, filename)
    assert os.path.exists(runscript)
    run = Run.from_definition(runscript)
    assert run.Solver.CLM.MetFilePath == str(forcing_dir)
    assert run.ComputationalGrid.NX == 15
    assert run.ComputationalGrid.NY == 15
    assert run.Geom.domain.Upper.X == 15000
    assert run.Geom.domain.Upper.Y == 15000


def test_change_filename_values_1(tmp_path):
    """Check the edited fiedls of the runscript file."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    old_runscript = get_ref_yaml_path("conus1", "transient", "box", test_dir)
    test_file = test_dir / "slope_x.pfb"
    test_file.write_bytes(b'Binary file contents')    
    new_runscript = subsettools.change_filename_values(old_runscript,
                                                       runname="my_new_run",
                                                       slopex=str(test_file)
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
    old_runscript = get_ref_yaml_path("conus1", "transient", "box", test_dir)
    test_file = test_dir / "slope_x.pfb"
    test_file.write_bytes(b'Binary file contents')
    new_runscript = subsettools.change_filename_values(old_runscript,
                                                       slopex=str(test_file)
    )
    # check that old file has been replaced:
    assert os.path.exists(new_runscript)
    assert old_runscript == new_runscript
    # check the edited fields of the new runscript:
    run = Run.from_definition(new_runscript)
    assert run.TopoSlopesX.FileName == str(test_file)
