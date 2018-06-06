import pytest
import numpy as np
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.refpix import RefPixStep
from jwst.linearity import LinearityStep
from jwst.dark_current import DarkCurrentStep
from jwst.jump import JumpStep
from jwst.ramp_fitting import RampFitStep
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import RampModel
from astropy.io import fits
from astropy.stats import sigma_clip


def test_extensions(darkcases):
    '''Test ramp-fit step with fake data.'''

    # open ramp to get data shape and headers
    m = RampModel(darkcases)

    fake_data_outname = darkcases[:-5]+"_test_extensions_uncal.fits"
    output,outint = RampFitStep.call(m,save_opt=True,opt_name=fake_data_outname[:-5]+"_rate_opt.fits")

    with fits.open(fake_data_outname[:-5]+"_rate_opt.fits") as h:

        # check all optional outputs are there
        assert h[1] == h['SLOPE']
        assert h[2] == h['SIGSLOPE']
        assert h[3] == h['YINT']
        assert h[4] == h['SIGYINT']
        assert h[5] == h['PEDESTAL']
        assert h[6] == h['WEIGHTS']
        assert h[7] == h['CRMAG']


def test_dq(darkcases):
    '''Test ramp-fit step with fake data.'''

    # open ramp to get data shape and headers
    m = RampModel(darkcases)

    # create one saturated pixel
    m.data[0,:,250,250] = m.data[0,:,250,250]+57000
    m.pixeldq[250,250] = 2.0

    # add in jump to one of the pixels
    m.data[0,2:,500,500] = m.data[0,2:,500,500]+5000
    m.groupdq[0,2,500,500] = 4.0

    output,outint = RampFitStep.call(m)

    # check DQ array
    before_satdq = m.pixeldq[250,250]
    after_satdq = output.dq[250,250]
    before_CRdq = m.groupdq[0,2,500,500]
    after_CRdq = output.dq[500,500]

    assert before_satdq == after_satdq
    assert before_CRdq == after_CRdq


def test_fake_rates_singleInt(darkcases,rates,pedestals):
    '''Test ramp-fit step with fake data.'''

    # open ramp to get data shape and headers
    m = RampModel(darkcases)
    tgroup = m.meta.exposure.group_time
    rates = np.float(rates)
    pedestals = np.float(pedestals)

    nrows = int(m.meta.subarray.xsize)
    ncols = int(m.meta.subarray.ysize)
    ngroups = int(m.meta.exposure.ngroups)
    nints = int(m.meta.exposure.nints)

    # create fake ramps with known slope and pedestal
    new_data =np.zeros((nints,ngroups,nrows, ncols), dtype=np.float32)
    for i in np.arange(0,ngroups):
        for j in np.arange(4,2044):
            for k in np.arange(4,2044):
                new_data[0,i,j, k] = pedestals+rates*((i+1)*tgroup)

    m.data = new_data
    m.err = np.zeros((nints,ngroups,nrows, ncols), dtype=np.float32)
    fake_data_outname = darkcases[:-5]+"_rate"+str(rates)+"_pedestal"+str(pedestals)+"_test_fake_rates_singleInt_uncal.fits"
    m.save(fake_data_outname,overwrite=True)

    output,outint = RampFitStep.call(m,output_file=fake_data_outname[:-5]+"_rate.fits",save_opt=True,opt_name=fake_data_outname[:-5]+"_rate_opt.fits")

    # check output rates in rate.fits file
    clip = sigma_clip(output.data)
    clip.data[clip.mask] = np.nan
    clip.data[output.dq != 0] = np.nan
    meanrate = np.nanmean(clip.data)
    assert np.allclose(meanrate,rates,rtol=8,atol=8) == True


def test_fake_pedestals(darkcases,rates,pedestals):
    '''Test ramp-fit step with fake data.'''

    # open ramp to get data shape and headers
    m = RampModel(darkcases)
    tgroup = m.meta.exposure.group_time
    rates = np.float(rates)
    pedestals = np.float(pedestals)

    nrows = int(m.meta.subarray.xsize)
    ncols = int(m.meta.subarray.ysize)
    ngroups = int(m.meta.exposure.ngroups)
    nints = int(m.meta.exposure.nints)

    # create fake ramps with known slope and pedestal
    new_data =np.zeros((nints,ngroups,nrows, ncols), dtype=np.float32)
    for i in np.arange(0,ngroups):
        for j in np.arange(4,2044):
            for k in np.arange(4,2044):
                new_data[0,i,j, k] = pedestals+rates*((i+1)*tgroup)

    # save it
    m.data = new_data
    m.err = np.zeros((nints,ngroups,nrows, ncols), dtype=np.float32)
    fake_data_outname = darkcases[:-5]+"_rate"+str(rates)+"_pedestal"+str(pedestals)+"_test_fake_pedestals_uncal.fits"
    m.save(fake_data_outname,overwrite=True)

    output,outint = RampFitStep.call(m,output_file=fake_data_outname[:-5]+"_rate.fits",save_opt=True,opt_name=fake_data_outname[:-5]+"_rate_opt.fits")
    optoutput = fits.open(fake_data_outname[:-5]+"rate_opt.fits")

    # check pedestal
    clip = sigma_clip(optoutput['PEDESTAL'].data)
    clip.data[clip.mask] = np.nan
    meanped = np.nanmean(clip.data)
    assert np.allclose(pedestals,meanped,rtol=2,atol=2) == True

    optoutput.close()

def test_CR_handling(darkcases,rates,pedestals):
    '''Test ramp-fit step with fake data.'''

    # open ramp to get data shape and headers
    m = RampModel(darkcases)
    tgroup = m.meta.exposure.group_time
    rates = np.float(rates)
    pedestals = np.float(pedestals)

    ngroups = 10
    nints = 1
    nrows = int(m.meta.subarray.xsize)
    ncols = int(m.meta.subarray.ysize)
    m.meta.exposure.ngroups = ngroups
    m.meta.exposure.ngroup = ngroups
    m.meta.exposure.nints = nints
    m.err = np.zeros((nints,ngroups,nrows, ncols), dtype=np.float32)
    m.groupdq = np.zeros((nints,ngroups,nrows, ncols), dtype=np.float32)

    # create fake ramps with known slope and pedestal
    new_data =np.zeros((nints,ngroups,nrows, ncols), dtype=np.float32)
    for ints in np.arange(0,nints):
        for i in np.arange(0,ngroups):
            for j in np.arange(4,2044):
                for k in np.arange(4,2044):
                    new_data[ints,i,j, k] = pedestals+rates*((i+1)*tgroup)

    # add in jump to one of the pixels
    new_data[0,2:,500,500] = new_data[0,2:,500,500]+(rates*5)
    m.groupdq[0,2,500,500] = 4.0

    # add in two jumps to one of the pixels
    new_data[0,2:,740,740] = new_data[0,2:,740,740]+(rates*5)
    new_data[0,6:,740,740] = new_data[0,6:,740,740]+(rates*6)
    m.groupdq[0,2,740,740] = 4.0
    m.groupdq[0,6,740,740] = 4.0

    # save it
    m.data = new_data
    fake_data_outname = darkcases[:-5]+"_rate"+str(rates)+"_pedestal"+str(pedestals)+"_test2_uncal.fits"
    # m.save(fake_data_outname,overwrite=True)

    output,outint = RampFitStep.call(m,output_file=fake_data_outname[:-5]+"rate.fits",save_opt=True,opt_name=fake_data_outname[:-5]+"rate_opt.fits")
    optoutput = fits.open(fake_data_outname[:-5]+"rate_opt.fits")

    # check output rates in regular output
    clip = sigma_clip(output.data)
    clip.data[clip.mask] = np.nan
    clip.data[output.dq != 0] = np.nan
    meanrate = np.nanmean(clip.data)
    assert np.allclose(meanrate,rates,rtol=8,atol=8) == True

    # check output rates in INTS output
    if nints > 1:
        for i in np.arange(0,nints):
            clip = sigma_clip(outint.data[nints,:,:])
            clip.data[clip.mask] = np.nan
            clip.data[output.dq != 0] = np.nan
            meanrate = np.nanmean(clip.data)
            assert np.allclose(meanrate,rates,rtol=8,atol=8) == True

    # CR rates from rate_opt.fits file
    ratebeforeCR1 = optoutput['SLOPE'].data[0,0,740,740]
    rateafterCR1 = optoutput['SLOPE'].data[0,1,740,740]
    ratebeforeCR2 = optoutput['SLOPE'].data[0,1,740,740]
    rateafterCR2 = optoutput['SLOPE'].data[0,2,740,740]
    assert np.allclose(ratebeforeCR1,rateafterCR1,rtol=1e-2,atol=1e-2) == True
    assert np.allclose(ratebeforeCR2,rateafterCR2,rtol=1e-2,atol=1e-2) == True

    # # check to make sure slope is weighted average of intervals
    # weights = optoutput['WEIGHTS'].data
    # interval1 = optoutput['SLOPE'].data[0,0,740,740]*weights[0,0,740,740]
    # interval2 = optoutput['SLOPE'].data[0,1,740,740]*weights[0,1,740,740]
    # interval3 = optoutput['SLOPE'].data[0,2,740,740]*weights[0,2,740,740]
    # calc = (interval1 + interval2 + interval3)/(weights[0,0,740,740] + weights[0,1,740,740] +weights[0,2,740,740])
    # print(calc)

    # other integrations shouldn't have CR hit
    if nints > 1:
        int2_noCRbefore = optoutput['SLOPE'].data[1,0,740,740]
        int2_noCRafter1 = optoutput['SLOPE'].data[1,1,740,740]
        int2_noCRafter2 = optoutput['SLOPE'].data[1,2,740,740]
        assert int2_noCRbefore == rates
        assert int2_noCRafter1 == 0.0
        assert int2_noCRafter2 == 0.0

    # CR rates for pix with no CR hit
    ratebefore = optoutput['SLOPE'].data[0,0,800,800]
    rateafter = optoutput['SLOPE'].data[0,1,800,800]
    assert ratebefore == output.data[800,800]
    assert rateafter == 0.0

    # Check CR magnitude
    # right now this is just calculated as the difference
    # between the two group values for the pixel. Is that right?
    manualCRmag = new_data[0,2,500,500]-new_data[0,1,500,500]
    pipeCRmag = optoutput['CRMAG'].data[0,0,500,500]
    assert np.allclose(manualCRmag,pipeCRmag,rtol=1,atol=1) == True

    manualCRmag = new_data[0,2,740,740]-new_data[0,1,740,740]
    pipeCRmag = optoutput['CRMAG'].data[0,0,740,740]
    assert np.allclose(manualCRmag,pipeCRmag,rtol=1,atol=1) == True

    manualCRmag = new_data[0,6,740,740]-new_data[0,5,740,740]
    pipeCRmag = optoutput['CRMAG'].data[0,1,740,740]
    assert np.allclose(manualCRmag,pipeCRmag,rtol=1,atol=1) == True


    optoutput.close()
