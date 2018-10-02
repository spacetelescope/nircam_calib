"""

Unit tests for saturation flagging

"""

import pytest
import numpy as np

from jwst.saturation.saturation import do_correction, correct_for_NaN
from jwst.datamodels import RampModel, SaturationModel, dqflags


def test_basic_saturation_flagging():
    '''Test basic saturation flagging.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    satvalue = 60000

    data, satmap = setup_cube(ngroups, nrows, ncols)
    data.data[0,0,500,500] = 0
    data.data[0,1,500,500] = 20000
    data.data[0,2,500,500] = 40000
    data.data[0,3,500,500] = 60000   # signal reaches saturation limit
    data.data[0,4,500,500] = 80000
    satmap.data[500,500] = satvalue
    output = do_correction(data, satmap)
    satindex = np.argmax(output.data[0,:,500,500] == satvalue)
    assert np.all(output.groupdq[0,satindex:,500,500] == dqflags.group['SATURATED'])


def test_signal_fluctuation_flagging():
    '''Test case where signal drops after reaching saturation.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    satvalue = 60000

    data, satmap = setup_cube(ngroups, nrows, ncols)
    data.data[0,0,500,500] = 0
    data.data[0,1,500,500] = 20000
    data.data[0,2,500,500] = 40000
    data.data[0,3,500,500] = 60000
    data.data[0,4,500,500] = 40000  # signal dips back below saturation limit
    satmap.data[500,500] = satvalue
    output = do_correction(data, satmap)
    satindex = np.argmax(output.data[0,:,500,500] == satvalue)
    assert np.all(output.groupdq[0,satindex:,500,500] == dqflags.group['SATURATED'])


def test_first_group_saturation():
    '''Test case where all groups are saturated.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    satvalue = 60000

    data, satmap = setup_cube(ngroups, nrows, ncols)
    data.data[0,0,500,500] = 60000
    data.data[0,1,500,500] = 70000
    data.data[0,2,500,500] = 70000
    data.data[0,3,500,500] = 60000
    data.data[0,4,500,500] = 70000
    satmap.data[500,500] = satvalue
    output = do_correction(data, satmap)
    satindex = np.argmax(output.data[0,:,500,500] == satvalue)
    assert np.all(output.groupdq[0,satindex:,500,500] == dqflags.group['SATURATED'])


def test_subarray_flagging():
    '''Test subarray extraction and saturation flagging.'''

    ngroups = 5
    nrows = 200
    ncols = 200
    satvalue = 60000

    data, satmap = setup_cube(ngroups, nrows, ncols)
    print(data.data.shape,satmap.data.shape)
    data.data[0,0,150,150] = 0
    data.data[0,1,150,150] = 20000
    data.data[0,2,150,150] = 40000
    data.data[0,3,150,150] = 60000
    data.data[0,4,150,150] = 80000
    satmap.data[150,150] = satvalue
    output = do_correction(data, satmap)
    satindex = np.argmax(output.data[0,:,150,150] == satvalue)
    assert np.all(output.groupdq[0,satindex:,150,150] == dqflags.group['SATURATED'])


def test_dq_propagation():
    '''Test pixel DQ propagation.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    dqval1 = 5
    dqval2 = 10

    data, satmap = setup_cube(ngroups, nrows, ncols)
    data.pixeldq[5, 5] = dqval1
    satmap.dq[5, 5] = dqval2
    output = do_correction(data, satmap)
    assert output.pixeldq[5, 5] == dqval1 + dqval2


def test_no_sat_check():
    '''Test pixels with NO_SAT_CHECK are not flagged.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    satvalue = 60000

    data, satmap = setup_cube(ngroups, nrows, ncols)
    satmap.dq[500,500] = dqflags.pixel['NO_SAT_CHECK']
    data.data[0,0,500,500] = 0
    data.data[0,1,500,500] = 20000
    data.data[0,2,500,500] = 40000
    data.data[0,3,500,500] = 60000
    data.data[0,4,500,500] = 80000
    output = do_correction(data, satmap)
    assert np.all(output.groupdq[0,:,500,500] != dqflags.group['SATURATED'])
    assert output.pixeldq[500,500] == dqflags.pixel['NO_SAT_CHECK']


def test_nans_in_mask():
    '''Test pixels with with NaN saturation values are not flagged.'''

    ngroups = 5
    nrows = 2048
    ncols = 2048
    HUGE_NUM = 100000.

    data, satmap = setup_cube(ngroups, nrows, ncols)
    satmap.data[500,500] = np.nan
    correct_for_NaN(satmap.data, satmap.dq)
    assert satmap.data[500,500] == HUGE_NUM
    assert satmap.dq[500,500] == dqflags.pixel['NO_SAT_CHECK']


@pytest.fixture(scope='function')
def setup_cube(ngroups, nrows, ncols):
    ''' Set up fake data to test.'''

    nints = 1

    data_model = RampModel()
    data_model.data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    data_model.pixeldq = np.zeros(shape=(nrows, ncols), dtype=np.int32)
    data_model.groupdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    data_model.meta.subarray.xstart = 1
    data_model.meta.subarray.ystart = 1
    data_model.meta.subarray.xsize = ncols
    data_model.meta.subarray.ysize = nrows
    data_model.meta.exposure.ngroups = ngroups
    data_model.meta.instrument.name = 'NIRCAM'

    saturation_model = SaturationModel()
    saturation_model.data = np.zeros(shape=(2048, 2048), dtype=np.float32)
    saturation_model.dq = np.zeros(shape=(2048, 2048), dtype=np.int32)
    saturation_model.meta.subarray.xstart = 1
    saturation_model.meta.subarray.ystart = 1
    saturation_model.meta.subarray.xsize = 2048
    saturation_model.meta.subarray.ysize = 2048
    saturation_model.meta.instrument.name = 'NIRCAM'
    saturation_model.meta.description = 'Fake data.'
    saturation_model.meta.telescope = 'JWST'
    saturation_model.meta.reftype = 'SaturationModel'
    saturation_model.meta.author = 'Alicia'
    saturation_model.meta.pedigree = 'Dummy'
    saturation_model.meta.useafter = '2015-10-01T00:00:00'

    return data_model, saturation_model
