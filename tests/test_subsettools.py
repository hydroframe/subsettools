import os
import pytest
from subsettools import subsettools
from parflow.tools.io import read_pfb
import numpy as np
import os
from subsettools.datasets import get_ref_yaml_path
from hf_hydrodata.grid import to_latlon
from parflow import Run


@pytest.mark.parametrize(
    "huc_list, grid, result", [(["15060202"], "conus1", (375, 239, 487, 329)),
                               (["140802"], "conus1", (572, 337, 797, 577)),
                               (["14080201", "14080202", "14080203", "14080204", "14080205"], "conus1", (572, 337, 797, 577))]
)
def test_huc_to_ij(huc_list, grid, result):
    assert subsettools.huc_to_ij(huc_list, grid) == result

    
def test_forcing_timezones(tmp_path):
    "Check if we get the correct forcing (temperature) in EST time."
    utc = tmp_path / 'UTC_out'
    utc.mkdir()
    est = tmp_path / 'EST_out'
    est.mkdir()
    ij_bounds = (375, 239, 487, 329)
    grid = 'conus1'
    start = '2005-10-03' 
    dataset = 'NLDAS2'
    subsettools.subset_forcing(ij_bounds=ij_bounds, grid=grid, start=start, end='2005-10-05', dataset=dataset, write_dir=utc, time_zone='UTC')
    subsettools.subset_forcing(ij_bounds=ij_bounds, grid=grid, start=start, end='2005-10-04', dataset=dataset, write_dir=est, time_zone='EST')
    utc_temp1 = read_pfb(os.path.join(utc, "NLDAS.Temp.000001_to_000024.pfb"))
    utc_temp2 = read_pfb(os.path.join(utc, "NLDAS.Temp.000025_to_000048.pfb"))
    est_temp_correct = np.concatenate((utc_temp1[5:, :, :], utc_temp2[:5, :, :]), axis=0)
    est_temp = read_pfb(os.path.join(est, "NLDAS.Temp.000001_to_000024.pfb"))
    assert np.array_equal(est_temp_correct, est_temp)

    
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


def test_edit_runscript_for_subset_3(tmp_path):
    """Check that exception is raised if forcing_dir is invalid."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = get_ref_yaml_path("conus1", "transient", "box", test_dir)
    forcing_dir = os.path.join(tmp_path, "forcing")
    with pytest.raises(Exception) as e_info:
        runscript = subsettools.edit_runscript_for_subset((10, 10, 25, 25),
                                                          runscript,
                                                          runname="my_new_run",
                                                          forcing_dir=forcing_dir
        )
    

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
