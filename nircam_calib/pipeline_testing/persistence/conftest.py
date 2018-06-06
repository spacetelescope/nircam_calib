#! /usr/bin/env python

'''
This file contains pytest fixtures that can be called by multiple
test*py files. Good for tests that will be common to more than one
pipeline step
'''

import pytest


def pytest_generate_tests(metafunc):
    if 'input_file' in metafunc.funcargnames:
        metafunc.parametrize("input_file", [
            ("dither1_B5_F250M_linear.fits")
        ])

    if 'trapsfilled' in metafunc.funcargnames:
        metafunc.parametrize("trapsfilled", [
            (""),
            ("dither1_B5_F250M_linear_trapsfilled.fits")
        ])

    if 'flagthresh' in metafunc.funcargnames:
        metafunc.parametrize("flagthresh", [
            (None),
            (40)
        ])


