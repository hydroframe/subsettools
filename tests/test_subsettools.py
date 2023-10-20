import os
import pytest
from subsettools import subsettools
from subsettools.datasets import get_ref_yaml_path
from parflow import Run
from parflow.tools.fs import mkdir, get_absolute_path, rm

@pytest.mark.parametrize(
    "huc_list, grid, result", [(["15060202"], "conus1", (375, 239, 487, 329))]
)
def test_huc_to_ij(huc_list, grid, result):
    assert subsettools.huc_to_ij(huc_list, grid) == result


def test_edit_runscript_for_subset_1(tmp_path):
    """Check the edited fiedls of the runscript file."""
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = get_ref_yaml_path("conus1", "transient", "box", test_dir)
    runscript = subsettools.edit_runscript_for_subset((10, 10, 25, 25),
                                                      runscript,
                                                      runname="my_new_run",
                                                      forcing_dir="my_forcing_dir"
    )
    run = Run.from_definition(runscript)
    assert run.get_name() == "my_new_run"
    assert run.Solver.CLM.MetFilePath == "my_forcing_dir"
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
    runscript = subsettools.edit_runscript_for_subset((10, 10, 25, 25),
                                                      runscript,
                                                      write_dir=write_dir,
                                                      forcing_dir="my_forcing_dir"
    )
    assert runscript == os.path.join(write_dir, filename)
    run = Run.from_definition(runscript)
    assert run.Solver.CLM.MetFilePath == "my_forcing_dir"
    assert run.ComputationalGrid.NX == 15
    assert run.ComputationalGrid.NY == 15
    assert run.Geom.domain.Upper.X == 15000
    assert run.Geom.domain.Upper.Y == 15000
