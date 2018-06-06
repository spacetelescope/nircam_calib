import copy,os,sys
import pytest
import numpy as np
import pandas as pd
from tools import *
from astropy.io import fits
from astropy.io import ascii
from astropy.convolution import convolve
from astropy.convolution import Box1DKernel
from astropy.convolution import Box2DKernel
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from jwst.refpix import RefPixStep
from jwst.superbias import SuperBiasStep
from jwst.saturation import SaturationStep
from jwst.dq_init import DQInitStep
from jwst.datamodels import RampModel,dqflags


######
# note
######
# To run: pytest test_refpix.py --html=report.html --self-contained-html

def display_multi_image(diff,number,rtol,test_name):
        '''Function to generate subplot for all groups.'''

        plt.figure(figsize=(15,30))
        columns = 2
        titles = ['g'+str(i) for i in np.arange(0,number)]
        images = [diff[0,i,:,:] for i in np.arange(0,number)]
        for i, image in enumerate(images):
            plt.subplot(len(images) / columns + 1, columns, i + 1)
            plt.imshow(image,vmin=-rtol,vmax=rtol, cmap=plt.cm.gray)
            plt.colorbar(orientation='horizontal',pad=0.1)
            plt.title(titles[i])
        plt.savefig(test_name+'_differences.png')


def oneoverf_medians(refpix,minn,maxx,pixeldq,smoothing_length):
    '''Function to calculate the 1/f medians.'''

    # calculate the array of medians of a set of side reference pixels
    nrows,ncols = refpix.shape
    dq = pixeldq[:,minn:maxx]

    # mirror the refpix on top and bottom so that we can use a moving box median
    addlen = np.int((smoothing_length-1)/2)

    #tmp = refpix[0:addlen,:]
    tmp = refpix[1:addlen+1,:] #to mirror without the 0th row
    mirrorbottom = tmp[::-1,:]

    tmpdq = dq[1:addlen+1,:] #to mirror without the 0th row
    mirrorbottomdq = tmpdq[::-1,:]

    #tmp = refpix[nrows-addlen:nrows,:]
    tmp = refpix[nrows-addlen-1:nrows-1,:] #to mirror without top row
    mirrortop = tmp[::-1,:]

    #tmp = refpix[nrows-addlen:nrows,:]
    tmpdq = dq[nrows-addlen-1:nrows-1,:] #to mirror without top row
    mirrortopdq = tmpdq[::-1,:]

    extra_rpix = np.zeros((nrows+smoothing_length-1,ncols))
    extra_rpix[0:addlen,:] = mirrorbottom
    extra_rpix[nrows+addlen:nrows+2*addlen,:] = mirrortop
    extra_rpix[addlen:nrows+addlen,:] = refpix

    extra_rpixdq = np.zeros((nrows+smoothing_length-1,ncols),dtype='int32')
    extra_rpixdq[0:addlen,:] = mirrorbottomdq
    extra_rpixdq[nrows+addlen:nrows+2*addlen,:] = mirrortopdq
    extra_rpixdq[addlen:nrows+addlen,:] = dq

    exshape = extra_rpix.shape

    numrows,numcol = extra_rpix.shape
    medians = np.zeros(nrows)
    for row in range(nrows):
        goodpixels = np.where(np.bitwise_and(extra_rpixdq[(row+addlen)-addlen:(row+addlen)+(addlen+1),:], dqflags.pixel['DO_NOT_USE']) == 0)
        medians[row] = np.median(extra_rpix[(row+addlen)-addlen:(row+addlen)+(addlen+1),:][goodpixels])

    return medians.astype(float)


def amp_only_refpix_corr(data,sigma,left,right,bottom,top,goodpix):
    '''Function to manually calculate the amp-only refpix correction.'''

    # check number of amps used
    if left > 0:
        if right > 0:
            amps = 4
        else:
            amps = 1
    else:
        amps = 1

    # get data shape and amp boundaries
    nint,ngroup,ys,xs = data.shape

    if amps == 4:
        bounds = [0,left+508,left+1020,left+1532,xs]
    else:
        bounds = [0,xs]

    vbounds = [0,ys]

    for integ in range(nint):
        hold = []

        for group in range(ngroup):

            # get sigma-clipped top and bottom averages, then average those
            for i in range(len(bounds)-1):
                if bottom > 0:
                    tmp = data[integ,group,0:bottom,bounds[i]:bounds[i+1]]
                    good = goodpix[0:bottom,bounds[i]:bounds[i+1]]
                    elembottom,clower,cupper = sigmaclip(tmp[good],low=sigma,high=sigma)
                    meanbottom = np.mean(elembottom)
                if top > 0:
                    tmp = data[integ,group,ys-top:ys,bounds[i]:bounds[i+1]]
                    good = goodpix[ys-top:ys,bounds[i]:bounds[i+1]]
                    elemtop,clower,cupper = sigmaclip(tmp[good],low=sigma,high=sigma)
                    meantop = np.mean(elemtop)
                if ((bottom>0) & (top>0)):
                    amean = (meanbottom + meantop) * 0.5
                elif ((bottom > 0) & (top == 0)):
                    amean = meanbottom
                elif ((bottom == 0) & (top > 0)):
                    amean = meantop

                # subtract top and bottom average
                data[integ,group,:,bounds[i]:bounds[i+1]] -= amean

                # put values into a table to compare groups and amps
                cdat,datlow,dathigh = sigmaclip(data[integ,group,:,bounds[i]:bounds[i+1]],low=sigma,high=sigma)
                meandat = np.mean(cdat)
                hold.append([group,i,meandat,amean])

    return data.astype(float), hold


def even_odd_refpix_corr(data,sigma,left,right,bottom,top,goodpix):
    '''Function to manually calculate the even/odd refpix correction.'''

    # check number of amps used
    if left > 0:
        if right > 0:
            amps = 4
        else:
            amps = 1
    else:
        amps = 1

    # get data shape and amp boundaries
    nint,ngroup,ys,xs = data.shape

    if amps == 4:
        bounds = [0,left+508,left+1020,left+1532,xs]
    else:
        bounds = [0,xs]

    vbounds = [0,ys]

    for integ in range(nint):
        holde = []
        holdo = []

        for group in range(ngroup):

            # get sigma-clipped even/odd column averages for the top and bottom
            # then average those
            for i in range(len(bounds)-1):
                if bottom > 0:
                    tmpeven = data[integ,group,0:bottom,bounds[i]:bounds[i+1]:2]
                    tmpodd = data[integ,group,0:bottom,bounds[i]+1:bounds[i+1]:2]
                    evenelembottom,clower,cupper = sigmaclip(tmpeven,low=sigma,high=sigma)
                    oddelembottom,clower,cupper = sigmaclip(tmpodd,low=sigma,high=sigma)
                    evenbottom = np.mean(evenelembottom)
                    oddbottom = np.mean(oddelembottom)
                if top > 0:
                    tmpeven = data[integ,group,ys-top:ys,bounds[i]:bounds[i+1]:2]
                    tmpodd = data[integ,group,ys-top:ys,bounds[i]+1:bounds[i+1]:2]
                    evenelemtop,clower,cupper = sigmaclip(tmpeven,low=sigma,high=sigma)
                    oddelemtop,clower,cupper = sigmaclip(tmpodd,low=sigma,high=sigma)
                    eventop = np.mean(evenelemtop)
                    oddtop = np.mean(oddelemtop)
                if ((bottom>0) & (top>0)):
                    emean = (evenbottom + eventop) * 0.5
                    omean = (oddbottom + oddtop) * 0.5
                elif ((bottom > 0) & (top == 0)):
                    emean = evenbottom
                    omean = oddbottom
                elif ((bottom == 0) & (top > 0)):
                    emean = eventop
                    omean = oddtop

                # subtract top and bottom average for even and odd columns
                data[integ,group,:,bounds[i]:bounds[i+1]:2] -= emean
                data[integ,group,:,bounds[i]+1:bounds[i+1]:2] -= omean

                # put even values into a table to compare groups and amps
                cedat,datelow,datehigh = sigmaclip(data[integ,group,:,bounds[i]:bounds[i+1]:2],low=sigma,high=sigma)
                meanedat = np.mean(cedat)
                holde.append([group,i,meanedat,emean])

                # put odd values into a table to compare groups and amps
                codat,datolow,datohigh = sigmaclip(data[integ,group,:,bounds[i]+1:bounds[i+1]:2],low=sigma,high=sigma)
                meanodat = np.mean(codat)
                holdo.append([group,i,meanodat,omean])


    return data.astype(float),holde,holdo


def include_1overf(data,sigma,left,right,bottom,top,pixeldq,smoothing_length,side_gain):
    '''Function to remove 1/f noise using side reference pixels.'''

    # check number of amps used
    if left > 0:
        if right > 0:
            amps = 4
        else:
            amps = 1
    else:
        amps = 1

    # get data shape and amp boundaries
    nint,ngroup,ys,xs = data.shape

    if amps == 4:
        bounds = [0,left+508,left+1020,left+1532,xs]
    else:
        bounds = [0,xs]

    vbounds = [0,ys]

    # get left and right medians, then get the average of those
    for integ in range(nint):
        hold = []

        for group in range(ngroup):
            lmeds = None
            rmeds = None
            if left > 0:
                lpix = data[integ,group,:,0:left]
                lmeds = oneoverf_medians(lpix,0,left,pixeldq,smoothing_length)
            if right > 0:
                rpix = data[integ,group,:,xs-right:xs]
                rmeds = oneoverf_medians(rpix,xs-right,xs,pixeldq,smoothing_length)
            if (lmeds is not None):
                if (rmeds is not None):
                    meds = np.mean([lmeds,rmeds],axis=0)
                else:
                    meds = lmeds
            else:
                if (rmeds is not None):
                    meds = rmeds

            # remove the 1/f noise
            for row in range(ys):
                data[integ,group,row,:] -= side_gain * meds[row]

                # put values into a table to compare groups and amps
                cdat,datlow,dathigh = sigmaclip(data[integ,group,row,:],low=sigma,high=sigma)
                meandat = np.mean(cdat)
                hold.append([group,row,meandat,meds[row]])

    return data.astype(float), hold


def test_amps_only(cases,sigmas,tolerances):
    '''Test amp-only average subtraction.'''

    test_name = "amp_only"

    # pipeline refpix correction results
    # ----------------
    refq = RefPixStep.call(cases,
                            odd_even_columns=False,
                            use_side_ref_pixels=False,
                            odd_even_rows = False)
    outname = cases[:-5]+'_refpix_ampsonly_pipeline.fits'
    refq.save(outname)


    # manual refpix correction results
    # --------------

    # read in input file
    ramp = RampModel(cases)
    rampdata = np.copy(ramp.data)
    pixeldq = ramp.pixeldq
    goodpix = pixeldq == 0

    #extract subarray if necessary
    xs,xe,ys,ye = get_coords_rampmodel(ramp)

    # get bounds
    num_left = np.max([4-xs,0])
    num_right = np.max([xe-2044,0])
    num_bottom = np.max([4-ys,0])
    num_top = np.max([ye-2044,0])

    # do the manual subtraction
    refq_manual, table = amp_only_refpix_corr(rampdata,sigmas,num_left,num_right,num_bottom,num_top,goodpix)

    # save out table of results to compare groups and amps
    save_df_table(table,cases[:-5]+'amp_only_refpix_corr.dat')

    # compare manual to pipeline
    diff = refq_manual - refq.data

    # save a fits file of the differences between manual and pipeline subtraction
    images = RampModel()
    images.data = diff
    images.save(cases[:-5]+'_refpix_ampsonly_differences.fits',overwrite=True)

    # check some values
    print("Group, Diffs: Amp1 through Amp4")
    for i in np.arange(0,diff.shape[1]):
        print('')
        print('Pipeline: '+str(i)+','+str(refq.data[0,i,12,12])+','+str(refq.data[0,i,12,600])+','+str(refq.data[0,i,12,1030])+','+str(refq.data[0,i,12,1600]))
        print('Manual: '+str(i)+','+str(refq_manual[0,i,12,12])+','+str(refq_manual[0,i,12,600])+','+str(refq_manual[0,i,12,1030])+','+str(refq_manual[0,i,12,1600]))

    # save out data
    outramp = RampModel()
    outramp.data = refq_manual
    outramp.save(cases[:-5]+'_refpix_ampsonly_manual.fits',overwrite=True)

    # pytest to make sure pipeline ~= manual
    if np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True) == False:

        # if test fails, get image of (manual - pipeline) for each group
        display_multi_image(diff,np.shape(diff)[1],tolerances[1],test_name)

        print('')
        print("Group, Max Difference")
        for i in np.arange(0,diff.shape[1]):
            print('Max difference between pipeline and manual: '+str(i)+','+str(np.max(np.absolute(diff[0,i,:,:]))))
        print('')

    assert np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True)== True


def test_amps_even_odd(cases,sigmas,tolerances):
    '''Test even/odd column average subtraction.'''

    test_name = 'amps_even_odd'

    # pipeline refpix correction results
    # ----------------
    refq = RefPixStep.call(cases,
                            odd_even_columns=True,
                            use_side_ref_pixels=False,
                            odd_even_rows = False)
    outname = cases[:-5]+'_refpix_ampsevenodd_pipeline.fits'
    refq.save(outname)


    # manual refpix correction results
    # --------------

    # read in input file
    ramp = RampModel(cases)
    rampdata = np.copy(ramp.data)
    pixeldq = ramp.pixeldq
    goodpix = pixeldq == 0

    #extract subarray if necessary
    xs,xe,ys,ye = get_coords_rampmodel(ramp)

    # get bounds
    num_left = np.max([4-xs,0])
    num_right = np.max([xe-2044,0])
    num_bottom = np.max([4-ys,0])
    num_top = np.max([ye-2044,0])

    # do the manual subtraction
    refq_manual,table_even,table_odd = even_odd_refpix_corr(rampdata,sigmas,num_left,num_right,num_bottom,num_top,goodpix)

    # save tables to compare groups and amps
    save_df_table(table_even,cases[:-5]+'_even_column_refpix_corr.dat')
    save_df_table(table_odd,cases[:-5]+'_odd_column_refpix_corr.dat')

    # compare manual to pipeline
    diff = refq_manual - refq.data

    # save a fits file of the differences between manual and pipeline subtraction
    images = RampModel()
    images.data = diff
    images.save(cases[:-5]+'_refpix_ampsevenodd_differences.fits',overwrite=True)

    # check some values
    print("Group, Diffs: Amp1 through Amp4")
    for i in np.arange(0,diff.shape[1]):
        print('')
        print('Pipeline: '+str(i)+','+str(refq.data[0,i,12,12])+','+str(refq.data[0,i,12,600])+','+str(refq.data[0,i,12,1030])+','+str(refq.data[0,i,12,1600]))
        print('Manual: '+str(i)+','+str(refq_manual[0,i,12,12])+','+str(refq_manual[0,i,12,600])+','+str(refq_manual[0,i,12,1030])+','+str(refq_manual[0,i,12,1600]))

    # save out data
    outramp = RampModel()
    outramp.data = refq_manual
    outramp.save(cases[:-5]+'_refpix_ampsevenodd_manual.fits',overwrite=True)

    # pytest to make sure pipeline ~= manual
    if np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True) == False:

        # if test fails, get image of (manual - pipeline) for each group
        display_multi_image(diff,np.shape(diff)[1],tolerances[1],test_name)

        print('')
        print("Group, Max Difference")
        for i in np.arange(0,diff.shape[1]):
            print('Max difference between pipeline and manual: '+str(i)+','+str(np.max(np.absolute(diff[0,i,:,:]))))
        print('')

    assert np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True)== True


def test_pipeline_even_odd_1overf(cases,sigmas,smoothing_lengths, tolerances,side_gains):
    '''Test 1/f noise subtraction with pipeline even/odd column averages.'''

    test_name = 'pipeline_even_odd_1overf'

    # pipeline refpix correction results
    # ----------------
    refq = RefPixStep.call(cases,
                            odd_even_columns=True,
                            use_side_ref_pixels=True,
                            side_smoothing_length=smoothing_lengths,
                            side_gain=side_gains,
                            odd_even_rows = False)
    outname = cases[:-5]+'_refpix_include1overf_pipeline_pipeevenodd.fits'
    refq.save(outname)


    # manual refpix correction results
    # --------------

    # read in input file
    ramp = RampModel(cases)
    rampdata = np.copy(ramp.data)
    pixeldq = ramp.pixeldq
    goodpix = pixeldq == 0

    #extract subarray if necessary
    xs,xe,ys,ye = get_coords_rampmodel(ramp)

    # get bounds
    num_left = np.max([4-xs,0])
    num_right = np.max([xe-2044,0])
    num_bottom = np.max([4-ys,0])
    num_top = np.max([ye-2044,0])

    # do the manual 1/f subtraction using pipeline even/odd averages
    refq_manual = RefPixStep.call(cases,
                                odd_even_columns=True,
                                use_side_ref_pixels=False,
                                odd_even_rows = False)
    refq_manual,table = include_1overf(refq_manual.data,sigmas,num_left,num_right,num_bottom,num_top,pixeldq,smoothing_lengths,side_gains)

    # save table to compare groups and amps
    save_df_table(table,cases[:-5]+'_include_1overf_pipeline_evenodd.dat')

    # compare manual to pipeline
    diff = refq_manual - refq.data

    # save a fits file of the differences between manual and pipeline subtraction
    images = RampModel()
    images.data = diff
    images.save(cases[:-5]+'_refpix_include1overf_differences_pipeevenodd.fits',overwrite=True)

    # check some values
    print("Group, Diffs: Amp1 through Amp4")
    for i in np.arange(0,diff.shape[1]):
        print('')
        print('Pipeline: '+str(i)+','+str(refq.data[0,i,12,12])+','+str(refq.data[0,i,12,600])+','+str(refq.data[0,i,12,1030])+','+str(refq.data[0,i,12,1600]))
        print('Manual: '+str(i)+','+str(refq_manual[0,i,12,12])+','+str(refq_manual[0,i,12,600])+','+str(refq_manual[0,i,12,1030])+','+str(refq_manual[0,i,12,1600]))

    # save out data
    outramp = RampModel()
    outramp.data = refq_manual
    outramp.save(cases[:-5]+'_refpix_include1overf_manual_pipeevenodd.fits',overwrite=True)

    # pytest to make sure pipeline ~= manual
    if np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True) == False:

        # if test fails, get image of (manual - pipeline) for each group
        display_multi_image(diff,np.shape(diff)[1],tolerances[1],test_name)

        print('')
        print("Group, Max Difference")
        for i in np.arange(0,diff.shape[1]):
            print('Max difference between pipeline and manual: '+str(i)+','+str(np.max(np.absolute(diff[0,i,:,:]))))
        print('')

    assert np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True)== True


def test_manual_even_odd_1overf(cases,sigmas, smoothing_lengths, tolerances,side_gains):
    '''Test 1/f noise subtraction with manual even/odd column averages.'''

    test_name = 'manual_even_odd_1overf'

    # pipeline refpix correction results
    # ----------------
    refq = RefPixStep.call(cases,
                            odd_even_columns=True,
                            use_side_ref_pixels=True,
                            side_smoothing_length=smoothing_lengths,
                            side_gain=side_gains,
                            odd_even_rows = False)
    outname = cases[:-5]+'_refpix_include1overf_pipeline_manevenodd.fits'
    refq.save(outname)


    # manual refpix correction results
    # --------------

    # read in input file
    ramp = RampModel(cases)
    rampdata = np.copy(ramp.data)
    pixeldq = ramp.pixeldq
    goodpix = pixeldq == 0

    #extract subarray if necessary
    xs,xe,ys,ye = get_coords_rampmodel(ramp)

    # get bounds
    num_left = np.max([4-xs,0])
    num_right = np.max([xe-2044,0])
    num_bottom = np.max([4-ys,0])
    num_top = np.max([ye-2044,0])

    # do the manual subtraction
    refq_manual,table_even,table_odd = even_odd_refpix_corr(rampdata,sigmas,num_left,num_right,num_bottom,num_top,goodpix)
    refq_manual,table = include_1overf(refq_manual,sigmas,num_left,num_right,num_bottom,num_top,pixeldq,smoothing_lengths,side_gains)

    # save table to compare groups and amps
    save_df_table(table, cases[:-5]+'_include_1overf_manual_evenodd.dat')

    # compare manual to pipeline
    diff = refq_manual - refq.data

    # save an image of the differences between manual and pipeline subtraction
    images = RampModel()
    images.data = diff
    images.save(cases[:-5]+'_refpix_include1overf_differences_manevenodd.fits',overwrite=True)

    # check some values
    print("Group, Diffs: Amp1 through Amp4")
    for i in np.arange(0,diff.shape[1]):
        print('')
        print('Pipeline: '+str(i)+','+str(refq.data[0,i,12,12])+','+str(refq.data[0,i,12,600])+','+str(refq.data[0,i,12,1030])+','+str(refq.data[0,i,12,1600]))
        print('Manual: '+str(i)+','+str(refq_manual[0,i,12,12])+','+str(refq_manual[0,i,12,600])+','+str(refq_manual[0,i,12,1030])+','+str(refq_manual[0,i,12,1600]))

    # save out data
    outramp = RampModel()
    outramp.data = refq_manual
    outramp.save(cases[:-5]+'_refpix_include1overf_manual_manevenodd.fits',overwrite=True)

    # pytest to make sure pipeline ~= manual
    if np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True) == False:

        # if test fails, get image of (manual - pipeline) for each group
        display_multi_image(diff,np.shape(diff)[1],tolerances[1],test_name)

        print('')
        print("Group, Max Difference")
        for i in np.arange(0,diff.shape[1]):
            print('Max difference between pipeline and manual: '+str(i)+','+str(np.max(np.absolute(diff[0,i,:,:]))))
        print('')

    assert np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True)== True


def test_1overf(cases,sigmas, smoothing_lengths, tolerances,side_gains):
    '''Test amp average and 1/f noise subtraction.'''

    test_name = 'only_1overf'

    # pipeline refpix correction results
    # ----------------
    refq = RefPixStep.call(cases,
                            odd_even_columns=False,
                            use_side_ref_pixels=True,
                            side_smoothing_length=smoothing_lengths,
                            side_gain=side_gains,
                            odd_even_rows = False)
    outname = cases[:-5]+'_refpix_only1overf_pipeline.fits'
    refq.save(outname)


    # manual refpix correction results
    # --------------

    # read in input file
    ramp = RampModel(cases)
    rampdata = np.copy(ramp.data)
    pixeldq = ramp.pixeldq
    goodpix = pixeldq == 0

    #extract subarray if necessary
    xs,xe,ys,ye = get_coords_rampmodel(ramp)

    # get bounds
    num_left = np.max([4-xs,0])
    num_right = np.max([xe-2044,0])
    num_bottom = np.max([4-ys,0])
    num_top = np.max([ye-2044,0])

    # do the manual subtraction
    refq_manual,table = amp_only_refpix_corr(rampdata,sigmas,num_left,num_right,num_bottom,num_top,goodpix)
    refq_manual,outtable = include_1overf(refq_manual,sigmas,num_left,num_right,num_bottom,num_top,pixeldq,smoothing_lengths,side_gains)

    # save table to compare groups and amps
    save_df_table(outtable, cases[:-5]+'_include_1overf_amponly.dat')

    # compare manual to pipeline
    diff = refq_manual - refq.data

    # save an image of the differences between manual and pipeline subtraction
    images = RampModel()
    images.data = diff
    images.save(cases[:-5]+'_refpix_only1overf_differences.fits',overwrite=True)

    # check some values
    print("Group, Diffs: Amp1 through Amp4")
    for i in np.arange(0,diff.shape[1]):
        print('')
        print('Pipeline: '+str(i)+','+str(refq.data[0,i,12,12])+','+str(refq.data[0,i,12,600])+','+str(refq.data[0,i,12,1030])+','+str(refq.data[0,i,12,1600]))
        print('Manual: '+str(i)+','+str(refq_manual[0,i,12,12])+','+str(refq_manual[0,i,12,600])+','+str(refq_manual[0,i,12,1030])+','+str(refq_manual[0,i,12,1600]))

    # save out data
    outramp = RampModel()
    outramp.data = refq_manual
    outramp.save(cases[:-5]+'_refpix_only1overf_manual.fits',overwrite=True)

    # pytest to make sure pipeline ~= manual
    if np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True) == False:

        # if test fails, get image of (manual - pipeline) for each group
        display_multi_image(diff,np.shape(diff)[1],tolerances[1],test_name)

        print('')
        print("Group, Max Difference")
        for i in np.arange(0,diff.shape[1]):
            print('Max difference between pipeline and manual: '+str(i)+','+str(np.max(np.absolute(diff[0,i,:,:]))))
        print('')

    assert np.allclose(refq_manual,refq.data,rtol=tolerances[1],atol=tolerances[0],equal_nan=True)== True
