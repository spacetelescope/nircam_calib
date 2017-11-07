import pytest
import numpy as np
from jwst.datamodels import RampModel
from jwst.dq_init import DQInitStep
from jwst.jump import JumpStep
from jwst.superbias import SuperBiasStep
from astropy.io import fits



def test_fake_data(thresholds):
    '''Mike Regan's basic jump test.'''
    with fits.open('NRCNRCA1-DARK-60012216201_1_481_SE_2016-01-02T02h34m28_uncal_sliced.fits') as hduraw:
        nints=1
        ngroups=10
        nrows=2048
        ncols=2048
        new_data =np.zeros((nints, ngroups, nrows, ncols), dtype=np.float32)
        new_data[0, 0:2, 100, 100] = 10.0
        new_data[0, 2:10, 100, 100] = 1000
        hduraw['SCI'].data = new_data
        hduraw.writeto("fake_test.fits",overwrite=True)
    JumpStep.call('fake_test.fits',output_file='fake_test_jump.fits',rejection_threshold=thresholds)
    with fits.open('fake_test_jump.fits') as hduout:
        dqdata = hduout['GROUPDQ'].data
        assert 4 == np.max(dqdata)  # a CR was found
        assert 2 == np.argmax(dqdata[0, :, 100, 100])  # find the CR in the expected group



def test_cr_threshold(cases,thresholds):
    '''Test to check the default threshold and that the user can change it.'''

    default = JumpStep.call(cases)
    grpdq_arr_default = default.groupdq
    outname = cases[:-5]+'_jump_threshold'+str(thresholds)+'.fits'

    if thresholds == 3:
        # make sure that the default is actually 3 (groupdq should be same)
        calibdata = JumpStep.call(cases,rejection_threshold=thresholds)
        calibdata.save(outname)
        grpdq_arr_out = calibdata.groupdq
        assert np.array_equal(grpdq_arr_default,grpdq_arr_out) == True
    else:
        # if threshold isn't 3, groupdq should be different than the default
        calibdata = JumpStep.call(cases,rejection_threshold=thresholds)
        calibdata.save(outname)
        grpdq_arr_out = calibdata.groupdq
        assert np.array_equal(grpdq_arr_default,grpdq_arr_out) == False



def test_dq_array(cases,xpix,ypix,thresholds):
    '''Test to check that pixel DQ values are not changed after jump step.'''

    dq_arr_in = RampModel(cases).pixeldq
    dq_arr_out = RampModel(cases[:-5]+'_jump_threshold'+str(thresholds)+'.fits').pixeldq
    assert dq_arr_in[xpix,ypix] == dq_arr_out[xpix,ypix]



def test_sat_pix(cases,thresholds):
    '''Test to make sure saturated pixels don't get flagged.'''

    pixdq_arr = RampModel(cases[:-5]+'_jump_threshold'+str(thresholds)+'.fits').pixeldq
    groupdq_arr = RampModel(cases[:-5]+'_jump_threshold'+str(thresholds)+'.fits').groupdq

    for inds in np.arange(0,len(groupdq_arr[0,:,:,:])):
        for x in np.arange(0,2048):
            for y in np.arange(0,2048):
                if pixdq_arr[y,x] == 2:
                    assert groupdq_arr[inds,y,x] == 0
