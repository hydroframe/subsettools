import os
import pytest
import subsettools as st
from parflow.tools.io import read_pfb
import numpy as np
import os
import hf_hydrodata
from parflow import Run


@pytest.mark.parametrize(
    "huc_list, grid, result", [(["15060202"], "conus1", (375, 239, 487, 329)),
                               (["140802"], "conus1", (572, 337, 797, 577)),
                               (["14080201", "14080202", "14080203", "14080204", "14080205"], "conus1", (572, 337, 797, 577))]
)
def test_huc_to_ij(huc_list, grid, result):
    bounds, _ = st.huc_to_ij(huc_list, grid)
    assert bounds == result

    
@pytest.mark.parametrize(
    "huc_list, grid", [(["01010001"], "conus1"),
                       (["01010001"], "conus2"),
                       (["03130003"], "conus1"),
                       (["1710"], "conus1"),
                       (["01010002", "01010001"], "conus1"),
    ]
)
def test_huc_to_ij_errors(huc_list, grid):
    """Check that a ValueError is raised if a HUC is not part of the grid."""
    with pytest.raises(ValueError) as e:
        st.huc_to_ij(huc_list, grid)
    

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
    st.subset_forcing(ij_bounds=ij_bounds, grid=grid, start=start, end='2005-10-05', dataset=dataset, write_dir=utc, time_zone='UTC')
    st.subset_forcing(ij_bounds=ij_bounds, grid=grid, start=start, end='2005-10-04', dataset=dataset, write_dir=est, time_zone='EST')
    utc_temp1 = read_pfb(os.path.join(utc, "NLDAS.Temp.000001_to_000024.pfb"))
    utc_temp2 = read_pfb(os.path.join(utc, "NLDAS.Temp.000025_to_000048.pfb"))
    est_temp_correct = np.concatenate((utc_temp1[5:, :, :], utc_temp2[:5, :, :]), axis=0)
    est_temp = read_pfb(os.path.join(est, "NLDAS.Temp.000001_to_000024.pfb"))
    assert np.array_equal(est_temp_correct, est_temp)

    
def test_latlon_to_ij():
    """Check that hf_hydrodata.grid.to_latlon and subsettools.latlon_to_ij are inverse to each other."""
    latlon_points = hf_hydrodata.grid.to_latlon("conus1", *[375, 239, 487, 329])
    latlon_points = [latlon_points[:2], latlon_points[2:]]
    assert st.latlon_to_ij(latlon_points, "conus1") == (375, 239, 487, 329)

    
def test_edit_runscript_for_subset_1(tmp_path):
    """Check the edited fiedls of the runscript file."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    runscript = st.edit_runscript_for_subset((10, 10, 25, 25),
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
    runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    filename = os.path.basename(runscript)
    forcing_dir = tmp_path / "forcing"
    forcing_dir.mkdir()
    runscript = st.edit_runscript_for_subset((10, 10, 25, 25),
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
    runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    forcing_dir = os.path.join(tmp_path, "forcing")
    with pytest.raises(Exception) as e_info:
        runscript = st.edit_runscript_for_subset((10, 10, 25, 25),
                                                          runscript,
                                                          runname="my_new_run",
                                                          forcing_dir=forcing_dir
        )
    

def test_change_filename_values_1(tmp_path):
    """Check the edited fiedls of the runscript file."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    old_runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    test_file = test_dir / "slope_x.pfb"
    test_file.write_bytes(b'Binary file contents')    
    new_runscript = st.change_filename_values(old_runscript,
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
    old_runscript = st.get_template_runscript("conus1", "transient", "box", test_dir)
    test_file = test_dir / "slope_x.pfb"
    test_file.write_bytes(b'Binary file contents')
    new_runscript = st.change_filename_values(old_runscript,
                                                       slopex=str(test_file)
    )
    # check that old file has been replaced:
    assert os.path.exists(new_runscript)
    assert old_runscript == new_runscript
    # check the edited fields of the new runscript:
    run = Run.from_definition(new_runscript)
    assert run.TopoSlopesX.FileName == str(test_file)


def test_subset_press_init(tmp_path):
    """Check that the call succeeds when it fetches data for the beginning of WY 2003.

    subset_press_init was initially designed to fetch data one hour before midnight on
    the date given. So in this case, it would look at 11pm on 2002-09-30, which does 
    not exist.
    """
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    filename = st.subset_press_init(
        ij_bounds=(375, 239, 487, 329),
        dataset="conus1_baseline_mod",
        date="2002-10-01",
        write_dir=test_dir,
    )
    assert filename == os.path.join(test_dir, "conus1_baseline_mod_2002.10.01:00.00.00_UTC0_press.pfb")
