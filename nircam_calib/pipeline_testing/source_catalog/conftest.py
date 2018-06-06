import numpy as np
import pytest

def pytest_generate_tests(metafunc):
    if 'cases' in metafunc.funcargnames:
        metafunc.parametrize("cases", [
            ("image3_dithered_modA_allSW_asn.json")
        ])

    if 'RAtol' in metafunc.funcargnames:
        metafunc.parametrize("RAtol", [
            (0.00002)
        ])

    if 'Dectol' in metafunc.funcargnames:
        metafunc.parametrize("Dectol", [
            (0.0002)
        ])

    if 'kernel_fwhm' in metafunc.funcargnames:
        metafunc.parametrize("kernel_fwhm", [
            (3.)
        ])

    if 'kernel_xsize' in metafunc.funcargnames:
        metafunc.parametrize("kernel_xsize", [
            (5.)
        ])

    if 'kernel_ysize' in metafunc.funcargnames:
        metafunc.parametrize("kernel_ysize", [
            (5.)
        ])

    if 'npixels' in metafunc.funcargnames:
        metafunc.parametrize("npixels", [
            (50)
        ])

    if 'snr_threshold' in metafunc.funcargnames:
        metafunc.parametrize("snr_threshold", [
            (5.),
            (10.),
            (20.)
        ])
