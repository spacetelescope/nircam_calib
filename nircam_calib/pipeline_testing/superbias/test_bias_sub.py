"""

Unit tests for superbias subtraction

(A. Canipe)

"""

import pytest
import numpy as np

from jwst.superbias.bias_sub import do_correction
from jwst.datamodels import RampModel, SuperBiasModel


def test_superbias_subtraction():
    '''Test subtraction of the bias.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    blevel = 2000.

    data, bias, ngroups = setup_cube(ngroups, nrows, ncols)
    data.data[:] = blevel
    bias.data[:] = blevel
    output = do_correction(data, bias)
    manual = data.data - bias.data
    assert np.array_equal(output.data, manual)


def test_subarray_correction():
    '''Test subarray extraction and correction.'''

    ngroups = 5
    nrows = 500
    ncols = 500
    blevel = 2000.

    data, bias, ngroups = setup_cube(ngroups, nrows, ncols)
    data.data[:] = blevel
    bias.data[:] = blevel
    output = do_correction(data, bias)
    manual = data.data - bias.data[:ncols, :nrows]
    assert np.array_equal(output.data, manual)


def test_dq_propagation():
    '''Test pixel DQ propagation.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    dqval1 = 5
    dqval2 = 10

    data, bias, ngroups = setup_cube(ngroups, nrows, ncols)
    data.pixeldq[5, 5] = dqval1
    bias.dq[5, 5] = dqval2
    output = do_correction(data, bias)
    assert output.pixeldq[5, 5] == dqval1 + dqval2


@pytest.fixture(scope='function')
def setup_cube(ngroups, nrows, ncols):
    ''' Set up fake data to test.'''

    nints = 1

    data_model = RampModel()
    data_model.data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    data_model.pixeldq = np.zeros(shape=(nrows, ncols), dtype=np.int32)
    data_model.meta.subarray.xstart = 1
    data_model.meta.subarray.ystart = 1
    data_model.meta.subarray.xsize = ncols
    data_model.meta.subarray.ysize = nrows
    data_model.meta.instrument.name = 'NIRCAM'

    bias_model = SuperBiasModel()
    bias_model.data = np.zeros(shape=(2048, 2048), dtype=np.float32)
    bias_model.dq = np.zeros(shape=(2048, 2048), dtype=np.int32)
    bias_model.meta.subarray.xstart = 1
    bias_model.meta.subarray.ystart = 1
    bias_model.meta.subarray.xsize = 2048
    bias_model.meta.subarray.ysize = 2048
    bias_model.meta.instrument.name = 'NIRCAM'
    bias_model.meta.description = 'Fake data.'
    bias_model.meta.telescope = 'JWST'
    bias_model.meta.reftype = 'SuperBiasModel'
    bias_model.meta.author = 'Alicia'
    bias_model.meta.pedigree = 'Dummy'
    bias_model.meta.useafter = '2015-10-01T00:00:00'

    return data_model, bias_model, ngroups
