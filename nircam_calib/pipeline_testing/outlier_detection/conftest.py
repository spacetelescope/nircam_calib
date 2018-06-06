import numpy as np
import pytest

def pytest_generate_tests(metafunc):
    if 'cases' in metafunc.funcargnames:
        metafunc.parametrize("cases", [
            ("image3_forOutlierImageTest_asn.json")
        ])

    if 'median' in metafunc.funcargnames:
        metafunc.parametrize("median", [
            ("image3_forMedianTest_asn.json")
        ])

    if 'outlier' in metafunc.funcargnames:
        metafunc.parametrize("outlier", [
            ("image3_forOutlierTest_asn.json")
        ])
