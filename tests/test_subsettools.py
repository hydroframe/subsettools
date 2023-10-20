import pytest
from subsettools import subsettools
from parflow.tools.fs import mkdir, rm
from parflow.tools.io import read_pfb
import numpy as np
import os

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
