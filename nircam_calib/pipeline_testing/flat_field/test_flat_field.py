import pytest
import numpy as np
from jwst.datamodels import RampModel
from jwst.flatfield import FlatFieldStep
from astropy.io import fits



def test_fake_data():
    '''Test flat field step with fake data.'''

    # open ramp to get data shape and headers, fill it with fake data
    with fits.open('nrca1_47Tuc_subpix_dither1_newpos_rate.fits') as hduraw:
        nrows=2048
        ncols=2048
        new_data =np.zeros((nrows, ncols), dtype=np.float32)
        new_data[:, :] = 50
        hduraw['SCI'].data = new_data
        hduraw['ERR'].data = new_data/10.
        hduraw['DQ'].data[500,500] = 50
        hduraw.writeto("fake_test.fits",overwrite=True)

    # call flat field step
    FlatFieldStep.call('fake_test.fits',output_file='fake_test_flat.fits')

    # grab reference file
    with fits.open('fake_test_flat.fits') as hduout:
        hdr = hduout['PRIMARY'].header['R_FLAT']
        print('header',hdr)
        jwst = np.int(hdr.find('jwst'))
        print('jwst',jwst)
        if hdr[:4] == 'crds':
            reffile = '/grp/crds/cache/references/jwst/'+hdr[jwst:]
        else:
            reffile = hdr
        print(reffile)

        # get output data -- check SCI, ERR, and DQ extensions
        refdata = fits.getdata(reffile,1)
        scidata = hduout['SCI'].data
        errdata = hduout['ERR'].data
        dqdata = hduout['DQ'].data

        # SCI data should be SCI/REF, ERR data should be ERR/ref, DQout == DQin
        assert np.array_equal(scidata,hduraw['SCI'].data/refdata) == True
        assert np.array_equal(errdata,hduraw['ERR'].data/refdata) == True
        assert dqdata[500,500] == hduraw['DQ'].data[500,500]



def test_flat_cases(cases,reffile):
    '''Test to check flat field step for different files.'''

    # get input data
    with fits.open(cases) as hduraw:
        in_data = hduraw['SCI'].data
        in_err = hduraw['ERR'].data
        in_dq = hduraw['DQ'].data[500,500]

    # call flat field step
    outname = cases[:-5]+'_flat_field.fits'
    if reffile is not "None":
        print(reffile)
        default = FlatFieldStep.call(cases,override_flat=reffile,output_file=outname)
    else:
        print(reffile)
        default = FlatFieldStep.call(cases,output_file=outname)

    # grab reference file
    with fits.open(outname) as hduout:
        hdr = hduout['PRIMARY'].header['R_FLAT']
        print('header',hdr)
        jwst = np.int(hdr.find('jwst'))
        print('jwst',jwst)
        if hdr[:4] == 'crds':
            reffile = '/grp/crds/cache/references/jwst/'+hdr[jwst:]
        else:
            reffile = hdr
        print(reffile)

        # get output data -- check SCI, ERR, and DQ extensions
        ref_data = fits.getdata(reffile,1)
        out_data = hduout['SCI'].data
        out_err = hduout['ERR'].data
        out_dq = hduout['DQ'].data[500,500]

        # SCI data should be SCI/REF, ERR data should be ERR/ref, DQout == DQin
        assert np.array_equal(out_data,in_data/ref_data) == True
        assert np.array_equal(out_err,in_err/ref_data) == True
        assert in_dq == out_dq
