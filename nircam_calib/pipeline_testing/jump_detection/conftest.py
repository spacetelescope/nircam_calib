import numpy as np
import pytest

def pytest_generate_tests(metafunc):
    if 'cases' in metafunc.funcargnames:
        metafunc.parametrize("cases", [
            ("nrca1_47Tuc_rapid_dark8_sat_superbias_refpix_linearity.fits"),
            ("nrca1_47Tuc_rapid_dark1_sat_superbias_refpix_linearity.fits")
        ])

    if 'xpix' in metafunc.funcargnames:
        metafunc.parametrize("xpix", [
            (100),
            (1284)
        ])

    if 'ypix' in metafunc.funcargnames:
        metafunc.parametrize("ypix", [
            (1320),
            (567)
        ])

    if 'thresholds' in metafunc.funcargnames:
        metafunc.parametrize("thresholds", [
            (4),
            (8),
            (15)
        ])
