import numpy as np
import pytest

def pytest_generate_tests(metafunc):
    if 'cases' in metafunc.funcargnames:
        metafunc.parametrize("cases", [
            ("spec2_redo_dispersion.json")
        ])

    if 'direct_image_file' in metafunc.funcargnames:
        metafunc.parametrize("direct_image_file", [
            ("V54321001002P000000000110d_A5_F444W_directImage_i2d.fits")
        ])

    # if 'object_list' in metafunc.funcargnames:
    #     metafunc.parametrize("object_list", [
    #         ("test_grism_obj_list.list")
    #     ])
