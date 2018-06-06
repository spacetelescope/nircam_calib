import numpy as np
import pytest

def pytest_generate_tests(metafunc):
    if 'cases' in metafunc.funcargnames:
        metafunc.parametrize("cases", [
            ("nrca1_47Tuc_subpix_dither1_newpos_rate.fits"),
            ("NRCNRCA1-DARK-60012216201_1_481_SE_2016-01-02T02h34m28_rate.fits")
        ])

    if 'reffile' in metafunc.funcargnames:
        metafunc.parametrize("reffile", [
            ("None"),
            ("/grp/crds/cache/references/jwst/jwst_nircam_flat_0054.fits")
        ])
