import os
from parflow import Run
from subsettools.datasets import get_ref_yaml_path

def test_get_ref_yaml_path(tmp_path):
    test_dir = tmp_path / "test"
    test_dir.mkdir()
    runscript = get_ref_yaml_path("conus1", "transient", "box", test_dir)
    assert runscript == os.path.join(test_dir, "conus1_pfclm_transient_box.yaml")
    run = Run.from_definition(runscript)
    assert run.ComputationalGrid.NX == 3342
    assert run.GeomInput.domaininput.InputType == "Box"
