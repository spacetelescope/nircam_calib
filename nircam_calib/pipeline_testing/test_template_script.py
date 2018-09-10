"""Template for formatting pytest code for the pipeline steps.

This module demonstrates the format you should use to build scripts to test the
calibration steps of the pipeline. For more information on Python's pytest, see
the link below.


Helpful links:
    Pipeline code: https://github.com/spacetelescope/jwst/tree/master/jwst
    Pytest: https://docs.pytest.org/en/documentation-restructure/how-to/index.html


Dependencies (for reporting and parallelization):
    pytest-cov
    pytest-xdist
    pytest-html


Example command to run pytests:
    This command depends on your specific directory structure, but this gives
    you the main arguments:

        $ pytest tests/ --html=report.html --self-contained-html --cov=tests --cov-report=html

    - "html" keyword tells pytest to display the test results in an html
      file (it will be saved in the current directory).
    - "self-contained-html" keyword creates a self-contained html report instead
      of assets such as CSS and images getting stored separately by default to
      respect the Content Security Policy.
    - "cov" keyword tells pytest to show coverage for the routines in the listed
      directory.
    - "cov-report" keyword tells pytest to generate an html report that goes
      into a lower directory. The html for a given module is annotated to show
      what is and is not covered. To get the default formatting, you have to
      display the html report from a directory that has the needed javascript
      files (they end in .js) that are placed in the directory.


Notes:
    - unit tests should not require or depend on input files (use fake
      reference data for unit tests)
    - all test function names must have a 'test_' prefix
    - functions without the prefix won't be run as pytests (but other functions
      can use them, e.g., setup_cube below)
    - test names should be descriptive for easy identification in the report

Some example scripts were taken from Mike Regan's jump test directory in the
JWST pipeline GitHub repository: https://github.com/spacetelescope/jwst

If you have any questions about the formatting, please talk to Alicia Canipe or
Bryan Hilbert.

"""


# Imports
import sys
import pytest
import numpy as np

# Pipeline-specific imports
from jwst.jump.twopoint_difference import find_CRs
from jwst.datamodels import dqflags




# To add a test that is not relevant at the moment, you can
#   mark it so that it's skipped
@pytest.mark.skip(reason="no way of currently testing this")
def test_the_unknown():
    '''This function tests a future pipeline step.'''

    print('No way to test.')


# You can also skip a test based on a condition, for example:
@pytest.mark.skipif(sys.version_info < (3, 3),
                    reason="requires python3.3")
def test_python3_function():
    '''This function is skipped if the wrong Python version is used.'''

    print('Only print if Python 3.3 is used.')


# You can use the xfail marker to indicate that you expect a test to fail:
@pytest.mark.xfail
def test_failing_function():
    '''This test is designed to fail.'''

    print('This fails, as expected.')
    # This test will run but no traceback will be reported when it fails.
    # Terminal reporting will list it in the “expected to fail” (XFAIL)
    # or “unexpectedly passing” (XPASS) sections.


# Below are two of Mike Regan's tests for jump detection
def test_5grps_cr3_noFlux():
    '''Test 5 groups: CR in Group 2.'''

    # create fake test data
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    data[0, 0:2, 100, 100] = 10.0
    data[0, 2:5, 100, 100] = 1000

    # calculate CR-cleaned median slope for this pixel
    median_diff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print(median_diff.shape)
    assert np.max(gdq) == 4 # a CR was found
    assert np.argmax(gdq[0, :, 100, 100]) == 2 # find the CR in the expected group


def test_6grps_Satat6_CRat1():
    '''Test 6 groups: CR in Group 2, Saturation in Group 6.'''

    # create fake test data
    ngroups = 6
    #crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 35000  # CR
    data[0, 2, 100, 100] = 40005
    data[0, 3, 100, 100] = 45029
    data[0, 4, 100, 100] = 50014
    data[0, 5, 100, 101] = 61000  # Saturation
    data[0, 0, 100, 101] = 10000
    data[0, 1, 100, 101] = 15001
    data[0, 2, 100, 101] = 20003
    data[0, 3, 100, 101] = 25006
    data[0, 4, 100, 101] = 30010
    data[0, 5, 100, 101] = 35015
    gdq[0, 5, 100, 100] = dqflags.group['SATURATED']

    # calculate CR-cleaned median slope for this pixel
    median_diff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print(median_diff.shape)
    # assert(4 == np.max(gdq))  # no CR was found
    assert np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0,
                           dqflags.group['SATURATED']], gdq[0, :, 100, 100])


# This function will be used by the pytest functions, but it isn't considered
# a test by the pytest tool because there is no 'test_' prefix.
def setup_cube(ngroups, readnoise=10):
    '''Convenience function to set up fake data cube.'''

    nints = 1
    nrows = 2048
    ncols = 2048
    rej_threshold = 3
    nframes = 1
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    read_noise = np.zeros(shape=(nrows, ncols), dtype=np.float32)
    read_noise[:, :] = readnoise
    #primary_gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    return data, gdq, nframes, read_noise, rej_threshold
