import numpy as np
import pytest

def pytest_generate_tests(metafunc):
    if 'darkcases' in metafunc.funcargnames:
        metafunc.parametrize("darkcases", [
            ("NRCNRCA1-DARK-60012216201_1_481_SE_2016-01-02T02h34m28_uncal_sliced.fits")
        ])

    if 'simcases' in metafunc.funcargnames:
        metafunc.parametrize("simcases", [
            ("nrca1_47Tuc_rapid_dark8_sat_superbias_refpix_linearity.fits")
        ])

    if 'rates' in metafunc.funcargnames:
        metafunc.parametrize("rates", [
            (800)
        ])

    if 'pedestals' in metafunc.funcargnames:
        metafunc.parametrize("pedestals", [
            (50)
        ])

    #
    # if 'rates' in metafunc.funcargnames:
    #     metafunc.parametrize("rates", [
    #         (20),
    #         (800),
    #         (3500)
    #     ])
    #
    #
    # if 'pedestals' in metafunc.funcargnames:
    #     metafunc.parametrize("pedestals", [
    #         (0),
    #         (12),
    #         (50)
    #     ])
