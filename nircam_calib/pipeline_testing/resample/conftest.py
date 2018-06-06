import numpy as np
import pytest

def pytest_generate_tests(metafunc):

    if 'asnfile' in metafunc.funcargnames:
        metafunc.parametrize("asnfile", [
            ("image3_directImage_asn.json")
        ])

    if 'RAtol' in metafunc.funcargnames:
        metafunc.parametrize("RAtol", [
            (0.00002)
        ])

    if 'Dectol' in metafunc.funcargnames:
        metafunc.parametrize("Dectol", [
            (0.0002)
        ])

    if 'kernel' in metafunc.funcargnames:
        metafunc.parametrize("kernel", [
            ('square')
        ])

    if 'good_bits' in metafunc.funcargnames:
        metafunc.parametrize("good_bits", [
            (4)
        ])
