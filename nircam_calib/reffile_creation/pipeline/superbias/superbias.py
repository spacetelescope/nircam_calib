#!/usr/bin/env python

'''try a new method for creating a NIRCam superbias image

This is the updated version that works with data processed using Build 5
of the SSB pipeline. The difference is that in Build 4, the bias_drift step
included a 0th read subtraction, refpix correction, then 0th read was added back
in. In Build 5 the software assumes that a superbias has already been subtracted
and so the bias_drift step (which has been renamed to refpix) performs the 
refpix correction without doing anything special with the 0th read.




UPDATING FOR MKREFS.PY. ADDING ABILITY TO CALL FROM WITHIN PYTHON
MAKING INTO A CLASS
BNH - 7 Nov 2017

Conversion to class complete. Now need to streamline...

'''


import os
import sys
import copy
import glob
import time
import shutil
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits,ascii
from jwst.datamodels import SuperBiasModel
from jwst.datamodels import dqflags
from jwst.ipc import IPCStep

#quad boundaries, in pixel units not including refpix
qstart = [0,508,1020,1532,2040]

checkx = 1015
checky = 1037

class Superbias:
    def __init__(self):
        self.something = something

    def calculate_means(self,data,boxwidth,edge_truncate=True):
        #given a 2D array of many differences for a single row,
        #calulate and return the means and uncertainties for every
        #boxwidth pixel. 
        #NOTE: if boxwidth=100, then we will calculate the mean in a
        #box that is +/-50 pixels around the central pixel. Later, the
        #solution will be applied within a box +/-25 pixels around the
        #central pixel. This implies that we need overlap between 
        #boxes. So for boxwidth of 100, we perform the calculation
        #centered around every 50th pixel.
        targets = np.arange(0,data.shape[1],boxwidth/2)
        means = np.zeros((data.shape[0],len(targets)))
        uncs = np.zeros_like(means)
        for group in xrange(data.shape[0]):
            for i,center in enumerate(targets):
                if center < boxwidth/2:
                    pix = data[group,0:boxwidth/2+1]
                elif data.shape[1]-boxwidth/2 < (boxwidth/2):
                    pix = data[group,data.shape[1]-boxwidth/2:]
                else:
                    pix = data[group,center-(boxwidth/2.):center+(boxwidth/2.)+1]
                    
                #sigma-clip
                goodpix, toolow, toohigh = sigmaclip(pix,low=3.,high=3.)
                means[group,i] = np.mean(goodpix)
                uncs[group,i] = np.std(goodpix)
        return(means,uncs,targets)
        

    def correct_1overf(self,data,filename,boxwidth,appwidth):
        print('Beginning 1/f correction, using NN2 method.')
        # Check to see if the output file already exists
        fdir, fname = os.path.split(filename)
        voutname = os.path.join(fdir,'NN2_V_for_',fname)
        coutname = os.path.join(fdir,'Corrected_data_for_',fname)

        # Output variables
        nn2corrdata = np.zeros_like(data)
        nn2corrdataerr = np.zeros_like(data)
        nn2corrmean = np.zeros((data.shape[1],data.shape[2]))
        nn2corrmeanunc = np.zeros_like(nn2corrmean)

        # Save some of the means, for debugging
        meancheck = np.zeros((2040,44))

        # We don't want to mix pixels from different amplifiers,
        # so work on each quad separately
        for amp in xrange(4):
            qdata = data[:,:,qstart[amp]:qstart[amp+1]]

            # Calculate mean values here. expect one mean value
            # per box of boxidth in each row

            # Work one row at a time to save memory
            zd,yd,xd = qdata.shape
            for rownum in range(yd):
                if rownum % 100 == 0:
                    print("Amp: {}, Row: {}".format(amp+1,rownum))
                rowdata = qdata[:,rownum,:]

                # Calculate the matrix of differences
                diffrowdata = np.zeros((zd*(zd-1)/2,xd))
                outidx = 0
                for idx in xrange(zd):
                    for idx2 in xrange(idx+1,zd):
                        diffrowdata[outidx,:] = rowdata[idx,:] - rowdata[idx2,:]
                        outidx = outidx + 1

                # Calculate the means in overlapping boxes of size boxwidth
                # output means,uncs are 2d arrays, groups x number of boxes in the row
                means,uncs,targx = self.calculate_means(diffrowdata,boxwidth,edge_truncate=True)

                # Save the means as a check for debugging
                meancheck[rownum,11*amp:11*amp+11] = means[0,:]

                # Now we need to work on each meanbox of the row.
                # For each box, construct the NxN matrix of all
                # possible difference values. This matrix is then fed
                #into the NN2 machinery.
                for box in xrange(means.shape[1]):
                    allmeans = np.zeros((zd,zd))
                    alluncs = np.zeros_like(allmeans)
                    idx = 0
                    for y in range(0,zd):
                        for x in range(y+1,zd):
                            allmeans[y,x] = means[idx,box]
                            allmeans[x,y] = 0. - means[idx,box]
                            alluncs[y,x] = uncs[idx,box]
                            alluncs[x,y] = uncs[idx,box]
                            idx = idx + 1

                    # Pass the means and uncertainties to the NN2 engine
                    nn2inst = nn2.nn2()
                    nn2inst.A = allmeans
                    nn2inst.E = alluncs
                    nn2inst.solveMatrix()

                    # Save the NN2 output vector and uncertainty vector
                    # for each pixel. We also need the mean vector, for
                    # later subtraction
                    for xpt in range(qstart[amp]+targx[box]-(appwidth/2),qstart[amp]+targx[box]+(appwidth/2)):
                        if xpt >= (qstart[amp]) and xpt < (qstart[amp+1]):
                            nn2corrdata[:,rownum,xpt] = nn2inst.V
                            nn2corrdataerr[:,rownum,xpt] = nn2inst.Verr
                            nn2corrmean[rownum,xpt] = np.mean(nn2inst.V)
                            nn2corrmeanunc[rownum,xpt] = np.std(nn2inst.V) / np.sqrt(len(nn2inst.V))

        print("Pixel tracking, nn2 corrected data. Pixel ({},{}) = {}".format(checkx,checky,nn2corrdata[0,checkx-4,checky-4]))

        # Save nn2 outputs for examination later
        hdu0 = fits.PrimaryHDU()
        hdu1 = fits.ImageHDU(nn2corrdata)
        hdu2 = fits.ImageHDU(nn2corrdataerr)
        hdu3 = fits.ImageHDU(nn2corrmean)
        hdu4 = fits.ImageHDU(nn2corrmeanunc)
        hdu = fits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4])
        hdu.writeto(voutname,overwrite=True)
        return nn2corrdata,nn2corrmean,nn2corrdataerr,nn2corrmeanunc

    
    def mean_refpixframe(self,data):
        # Compute the simple mean of each reference pixel
        # for use in the superbias with only 30(?) reads,
        # we probably don't have enough samples to make a
        # sigma-clipped mean, right?
        rpixframe = np.zeros(2048,2048)
        top = data[:,2044:,:]
        top_mean = np.mean(top,axis=0)
        bottom = data[:,0:4,:]
        bottom_mean = np.mean(bottom,axis=0)
        left = data[:,4:2044,0:4]
        left_mean = np.mean(left,axis=0)
        right = data[:,4:2044,2044:]
        right_mean = np.mean(right,axis=0)

        rpixframe[:,2044:,:] = top_mean
        rpixframe[:,0:4,:] = bottom_mean
        rpixframe[:,4:2044,0:4] = left_mean
        rpixframe[:,4:2044,2044:] = right_mean
        return rpixframe
        

    def bias_from_file(self,filename,boxwidth,appwidth,groupstart,groupend,integration,skipnn2,ipc,xtalk,showplots=False,saveplots=False):

        # Read in fits file
        hdu = fits.open(filename)
        data = hdu[1].data
        dq = hdu[2].data
        detector = hdu[0].header['DETECTOR'].strip()

        # Exptime per group. To be used later in line-fitting of
        # pixels with signal
        tgroup = hdu[0].header['TGROUP']
        hdu.close()

        # Just grab the first integration for the moment
        if len(data.shape) == 4:
            if data.shape[0] > 1:
                print("Extracting integration {} from the input file.".format(integration))
                data = data[integration,:,:,:]
            else:
                data = data[0,:,:,:]
    
        # Chop to save time in testing
        numgroups = data.shape[0]
        if ((groupstart >= 0) & (groupstart < numgroups)):
            if ((groupend > groupstart) & (groupend < numgroups)):
                print("========================")
                print(("Using groups {} through {} to construct the superbias."
                       .format(groupstart,groupend)))
                print("========================")
                data = data[groupstart:groupend+1,:,:]
            else:
                print(("WARNING: groupend is either less than groupstart, "
                       "or too large for the input file. Quitting."))
                sys.exit
        else:
            print(("WARNING: groupstart is either negative or too large "
                   "for the input file. Quitting."))
            print(groupstart)
            sys.exit
                
        # Remove reference pixels
        if data.shape[1] == 2048:
            # Make a copy of the original data (so that we have the refpix)
            # for testing
            sb_data = copy.deepcopy(data)
            withrefpixy = data.shape[1]
            withrefpixx = data.shape[2]
            data = data[:,4:withrefpixy-4,4:withrefpixx-4]
            dq = dq[4:withrefpixy-4,4:withrefpixx-4]
        else:
            print("WARNING: INPUT DATA ARE NOT 2048x2048! Quitting!")
            return 0

        # For DQ, all we care about is the DO_NOT_USE bit.
        # Extract that bit and throw out the rest
        mask = (dq & dqflags.pixel['DO_NOT_USE'] > 0)
        print(("{} pixels have the DO NOT USE flag set"
               .format(np.sum(mask))))
        
        # Set all pixels that have the DO_NOT_USE flag set
        # to NaN so they are not used in future calculations
        print(("Masking of bad pixels turned off because the number"
               "of NaNs was multiplying within the IPC and XTALK steps."))

        # Uncertainty in the data. At the moment, pipline-output
        # files have 0s in the error arrays, so we need to create
        # our own errors. Since the dark current signal
        #is negligible, let's just use single frame readnoise. 
        ron = 6.

        # Poor man's bias correction
        data = self.poor_mans_bias(data)

        # 1/f correction using NN2
        if skipnn2 == False:
            nn2corrdata, nn2corrmean, nn2corrdataerr, nn2corrmeanunc = self.correct_1overf(data,filename,boxwidth,appwidth)
            data = data + (nn2corrdata-nn2corrmean)
            sb_data[:,4:2044,4:2044] = sb_data[:,4:2044,4:2044] + (nn2corrdata-nn2corrmean)
            sb_unc = np.sqrt(np.absolute(sb_data[:,4:2044,4:2044]) + nn2corrdataerr*nn2corrdataerr + nn2corrmeanunc*nn2corrmeanunc + ron*ron)
            uncnan = np.isnan(sb_unc[0,:,:])
            print('Subtracting 1/f correction from data')
        else:
            print(("1/f correction using NN2 skipped. Need to "
                   "figure out true uncertainties in this case!!!!!!!!"))
            sb_unc = np.sqrt(np.absolute(sb_data[:,4:2044,4:2044]) + ron*ron)
        
        # Perform IPC correction if asked
        if ipc:
            ipc_image = ipc_correction_image(data,detector)
            print('Subtracting IPC image from data')
            sb_data[:,4:2044,4:2044] = sb_data[:,4:2044,4:2044] + ipc_image
            # Save ipc and image here, so that we can include the filename
            ipchdu = fits.PrimaryHDU(ipc_image)
            ipchdu.writeto('ipc_correction_image_for_integration'+str(integration)+'_of_'+filename,overwrite=True)

        if xtalk:
            # If IPC correction is done, remove IPC effects before attempting crosstalk correction
            data = data + ipc_image
            xtalk_corr = self.xtalk_correction_image(data,detector)
            print('Subtracting Xtalk image from data')
            sb_data[:,4:2044,4:2044] = sb_data[:,4:2044,4:2044]  - xtalk_corr
            # Save xtalk image
            xhdu = fits.PrimaryHDU(xtalk_corr)
            xhdu.writeto('xtalk_correction_image_for_integration'+str(integration)+'_of_'+filename,clobber=True)
    
        # Add reference pixels back in to sb_unc
        sb_unc_full = np.zeros((sb_unc.shape[0],2048,2048))
        sb_unc_full[:,2044:,:] = ron 
        sb_unc_full[:,0:4,:] = ron
        sb_unc_full[:,4:2044,0:4] = ron
        sb_unc_full[:,4:2044,2044:] = ron
        sb_unc_full[:,4:2044,4:2044] = sb_unc
        return (sb_data,sb_unc_full)
    
    def poor_mans_bias(self,orig):
        # Quick bias correction. Subtract 0th read, and
        # make an adjustment for the signal between reset
        # and read 0, approximating it as the sigal between
        # reads 0 and 1.
        print("Performing poor man's bias correction.")
        temp = 0. - orig[0,:,:] + (orig[1,:,:]-orig[0,:,:])
        orig = orig - orig[0,:,:] + (orig[1,:,:]-orig[0,:,:])
        return orig


    def ipc_correction_image(self,orig,detector,save=False):
        print('Generating IPC correction image')
        ipc_dir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/cv3_reffile_conversion/ipc/'
        ipc_prefix = {'NRCA1':'NRCA1_17004','NRCA2':'NRCA2_17006','NRCA3':'NRCA3_17012','NRCA4':'NRCA4_17048','NRCALONG':'NRCA5_17158','NRCB1':'NRCB1_16991','NRCB2':'NRCB2_17005','NRCB3':'NRCB3_17011','NRCB4':'NRCB4_17047','NRCBLONG':'NRCB5_17161'}
        ipc_base = '_IPCDeconvolutionKernel_2016-03-18_ssbipc.fits'
        
        ipc_kernel_file = ipc_dir + ipc_prefix[detector] + ipc_base
        ipchdu = fits.open(ipc_kernel_file)
        kernel = ipchdu[1].data
        ipchdu.close()

        shift_dict = {'0':1,'1':0,'2':-1}
        kernel[1,1] = kernel[1,1] - 1.
        print("IPC kernel being used is: {}".format(kernel))
        corr_image = np.zeros_like(orig)
    
        for amp in xrange(4):
            subframe = orig[:,:,qstart[amp]:qstart[amp+1]]

            nansubframe = np.isnan(subframe)
            if len(np.where(nansubframe == True)[0]) > 0:
                print("nan subframe!!")
                print(len(np.where(nansubframe == True)[0]))
            
            app_kernel = kernel
            if ((amp == 1) | (amp == 3)):
                app_kernel = np.fliplr(kernel)

            for group in range(orig.shape[0]):
                summed_image = np.zeros((subframe.shape[1],subframe.shape[2]))

                for j in range(kernel.shape[0]):
                    for i in range(kernel.shape[1]):
                        prod = subframe[group,:,:] * kernel[j,i]

                        nanlen = self.nancheck(prod)
                        if nanlen > 0:
                            print("NAN 1st Prod {},{},{},{}".format(amp,group,j,i))

                        prod = np.roll(prod,shift_dict[str(j)],axis=0)
                        nanlen = self.nancheck(prod)
                        if nanlen > 0:
                            print("NAN 2nd Prod {},{},{},{}".format(amp,group,j,i))

                        prod = np.roll(prod,shift_dict[str(i)],axis=1)

                        nanlen = self.nancheck(prod)
                        if nanlen > 0:
                            print("NAN 3rd Prod {},{},{},{}".format(amp,group,j,i))

                        summed_image = summed_image + prod

                corr_image[group,:,qstart[amp]:qstart[amp+1]] = summed_image
        nancorr = np.isnan(corr_image)
        if len(np.where(nancorr==True)[0]) > 0:
            print("nan corr!!")
            print(len(np.where(nancorr==True)[0]))

        # Save the correction image
        if save == True:
            hdu = fits.PrimaryHDU(corr_image)
            ipcname = 'ipc_correction_image_for_integration'+str(integration)+'_of_'+listfile+'.fits'
            hdu.writeto(ipcname,clobber=True)
        return corr_image

    
    def nancheck(self,data):
        # Get number of nans in data
        nan = np.isnan(np.ravel(data))
        num = len(np.where(nan == True)[0])
        return num

    def read_xtalk_file(self,file,detector):
        # Read in appropriate line from the xtalk
        # file and return the coeffs
        xtcoeffs = ascii.read(file,header_start=0,data_start=1)

        indet = detector[3:5]
        if indet == 'BL':
            indet = 'B5'
        if indet == 'AL':
            indet = 'A5'

        coeffs = []
        for i,det in enumerate(xtcoeffs['Det']):
            if det == indet:
                coeffs = xtcoeffs[i]

        if len(coeffs) == 0:
            print('Detector not found in xtalk file: {}'.format(indet))
            sys.exit()    
        return coeffs
        

    def xtalk_correction_image(self,orig,detector,save=False):
        print('Generating crosstalk correction image')
        xtfile = 'xtalk20150303g0.errorcut.txt'
        coeffs = read_xtalk_file(xtfile,detector)
        
        # Need to tack refpix back on in order to have
        # all four quadrants be the same width
        tmp = np.zeros((orig.shape[0],orig.shape[1],orig.shape[2]+8))
        tmp[:,:,4:2044] = orig
        orig = tmp
        xtqstart = [0,512,1024,1536,2048]

        # Work amp by amp, multiply by appropriate coefficient,
        # and subtract from appropriate amp
        xtalk_corr_im = np.zeros_like(orig)
        subamp_shift = {"0":1,"1":-1,"2":1,"3":-1}
        for group in range(orig.shape[0]):
            for amp in xrange(4):
                to_mult = orig[group,:,xtqstart[amp]:xtqstart[amp+1]]
                receivers = []
                for i in range(4):
                    if i != amp:
                        receivers.append(i)
                # Reverse the values to multply if the amps
                # being used are adjacent or 3 amps apart
                for subamp in receivers:
                    index = 'xt'+str(amp+1)+str(subamp+1)
                    if ((np.absolute(amp-subamp) == 1) | (np.absolute(amp-subamp) == 3)):
                        corr_amp = np.fliplr(to_mult) * coeffs[index]
                    if (np.absolute(amp-subamp) == 2):
                        corr_amp = to_mult * coeffs[index]
            
                    xtalk_corr_im[group,:,xtqstart[subamp]:xtqstart[subamp+1]] = self.xtalk_corr_im[group,:,xtqstart[subamp]:xtqstart[subamp+1]] + corr_amp

                #per Armin's instructions, now repeat the process
                # using his xt??post coefficients, but shift the arrays
                #by one pixel according to readout direction.
                for subamp in receivers:
                    index = 'xt'+str(amp+1)+str(subamp+1)+'post'
                    if ((np.absolute(amp-subamp) == 1) | (np.absolute(amp-subamp) == 3)):
                        corr_amp = np.fliplr(to_mult) * coeffs[index] 
                        corr_amp = np.roll(corr_amp,subamp_shift[str(subamp)],axis=1)
                    if (np.absolute(amp-subamp) == 2):
                        corr_amp = to_mult * coeffs[index]
                        corr_amp = np.roll(corr_amp,subamp_shift[str(subamp)])

                    xtalk_corr_im[group,:,xtqstart[subamp]:xtqstart[subamp+1]] = self.xtalk_corr_im[group,:,xtqstart[subamp]:xtqstart[subamp+1]] + corr_amp
                    
        # Strip the refpix back off
        xtalk_corr_im = xtalk_corr_im[:,:,4:2044]

        # Save the crosstalk correction image
        if save == True:
            phdu = fits.PrimaryHDU(xtalk_corr_im)
            phdu.writeto('xtalk_correction_image.fits',clobber=True)

        print("Xtalk correction image created")
        return xtalk_corr_im
                        
                    
    
    def collapse(self,ramp,err):
        # Take 3D array of original data - NN2
        # output (calculated in bias_from_file)
        # and take the mean in each pixel along
        # all reads. This will create a mean bias
        # frame.
        print("Averaging together the bias frames from all reads into a single superbias")
        mns = np.average(ramp,axis=0,weights=1./err)
        
        # Error propagation
        uncs = np.std(ramp,axis=0)/np.sqrt(ramp.shape[0])
        return (mns,uncs)


    def collapse_median(self,ramp,err):
        #same as collapse above, but uses a median
        # rather than a weighted average
        mns = np.median(ramp,axis=0)
        
        # Error propagation
        uncs = np.std(ramp,axis=0)/np.sqrt(ramp.shape[0])
        return (mns,uncs)
        

    def unreliable_bias(self,biaserr,sigmacut=5):
        # Find pixels with large variances in the
        # bias, and flag these for the DQ extension

        # Calculate the mean variance across the detector
        goodpix, toolow, toohigh = sigmaclip(biaserr)
        largevar = biaserr > toohigh
        dq = np.zeros(2048,2048)
        dq[largevar] = 1
        return dq


    def linfit_yerror(self,x,y,yerr):
        # From Bevington page 104
        # Any points with yerr == 0 will have yerr
        # set to a large number. Similarly, any points
        # with x or y == np.nan will be ignored
        baderr = np.where(yerr == 0)
        if len(baderr[0]) > 0:
            yerr[baderr] = 1e6

        xnan = np.isnan(x)
        y[xnan] = np.nan
        yerr[xnan] = np.nan

        ynan = np.isnan(y)
        x[ynan] = np.nan
        yerr[ynan] = np.nan
    
        # Sum over only the non-nan points
        good = ~np.isnan(x)
        sumis2 = np.sum(1.0/(yerr[good] * yerr[good]))
        sumx2 = np.sum(x[good]*x[good]/(yerr[good] * yerr[good]))
        sumx = np.sum(x[good]/(yerr[good] * yerr[good]))
        sumy = np.sum(y[good]/(yerr[good] * yerr[good]))
        sumxy = np.sum(x[good]*y[good]/(yerr[good] * yerr[good]))

        delta = sumis2*sumx2 - sumx*sumx
        offset_a     = 1.0/delta * (sumx2*sumy - sumx*sumxy)
        slope_b      = 1.0/delta * (sumis2*sumxy - sumx*sumy)
        offset_a_err = np.sqrt(1.0/delta * sumx2)
        slope_b_err  = np.sqrt(1.0/delta * sumis2)
        return offset_a, slope_b, offset_a_err, slope_b_err

    
    def pixels_with_signal(self,biasramp,biasramperr,biasframe,biasframeerr,detector,tgroup,nsig=9):
        # Identify pixels with significant signal in them
        # line-fit to these pixels and use the intercept
        # as the value in the superbias
        print(('Dealing with pixels that conatin significant signal,'
               'through CR checking and line-fitting, etc.'))

        float_formatter = lambda x: "%.2f" %x

        # Determine the appropriate linearity correction
        # coefficient file. We're going to need to
        # apply the linearity correction to the bright
        # pixels before line-fitting.
        linfiles = glob.glob('/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/cv3_reffile_conversion/linearity/*ADU0*DMSorient.fits')
        lin_reffile = [s for s in linfiles if detector in s][0]
        print('Using '+lin_reffile+' for linearity correction of pixels with signal.')

        hdu = fits.open(lin_reffile)
        lin_coeff = hdu[1].data
        hdu.close()

        # Remove reference pixels from non-lin coeffs
        lin_coeff = lin_coeff[:,4:2044,4:2044]

        # What is 'significant' signal? work amp by amp,
        # calculate the amp mean and stdev. 
        # significant signal is > N-sigma above the mean?
        # or just > N DN above the mean?
        times = tgroup * np.arange(biasramp.shape[0])

        for amp in range(4):
            biasframe_amp = biasramp[-1,:,qstart[amp]:qstart[amp+1]]
            nanpts = np.isnan(biasframe_amp)
            goodsigpix, toolow, toohigh = sigmaclip(biasframe,low=3.,high=3.)
            sigmn = np.mean(goodsigpix)
            sigdev = np.std(goodsigpix)
            sigpix = np.where((biasframe_amp >= (sigmn+nsig*sigdev)) | (biasframe_amp <= (sigmn-nsig*siddev)))
            perc = float_formatter(len(sigpix[0])/(510.*2040.)*100.)
            print(("Amp {}: {} pixels ({}% of total pixels)"
                   .format(amp+1,len(sigpix[0]),perc)))
            print("have appreciable signal and will be line-fit.")

            # Loop over pixels with signal so that we can line-fit
            for i in range(len(sigpix[0])):
                signal = biasramp[:,sigpix[0][i],sigpix[1][i]]
                uncertainty = biasramperr[:,sigpix[0][i],sigpix[1][i]]

                # Check the range of signal levels here. If
                # the range is high enough, we need to
                # linearity correct
                sigrange = np.max(signal) - np.min(signal)

                if sigrange > 500:
                    # Linearity correct here, but only if the
                    # linearity coeffs are not nan
                    lin = lin_coeff[:,sigpix[0][i],sigpix[1][i]]
                    linnan = ~np.isnan(lin)
                    if all(linnan):
                        signal = lin[0] + lin[1]*signal + lin[2]*signal**2 + lin[3]*signal**3 + lin[4]*signal**4 + lin[5]*signal**5
                        # CR check
                        cr = self.cr_check(times,signal,5)
                        if cr > 2:
                            print('linefit 1')
                            offset, slope, offset_err, slope_err = self.linfit_yerror(times[0:cr],signal[0:cr],uncertainty[0:cr])
                            if np.isnan(offset):
                                print("NaN returned from line-fit")
                                print("coordinates:",sigpix[0][i],sigpix[1][i])
                                print("signal:",signal)

                        if cr == 2:
                            slope = (signal[cr-1] - signal[0]) / (times[cr-1] - times[0])
                            offset = signal[0] - slope*times[0]
                            if np.isnan(offset):
                                print("NaN returned from 2-group fit")
                                print("coordinates:",sigpix[0][i],sigpix[1][i])
                                print("signal:",signal)

                        if cr == 1:
                            slope = 0
                            offset = signal[0]
                            if np.isnan(offset):
                                print("NaN returned from one-group fit")
                                print("coordinates:",sigpix[0][i],sigpix[1][i])
                                print("signal:",signal)

                        cflag = 'red'
                    else:
                        # If the linearity coefficients are NaN, check
                        # for CR hits, and then do a 3rd order fit to
                        # the uncorrected data
                        print('cr_check2')
                        cr = cr_check(times,signal,5)
                        if cr > 4:
                            coeffs = np.polyfit(times[0:cr],signal[0:cr],3)
                            offset = coeffs[3]
                            slope = coeffs[2]
                            if np.isnan(offset):
                                print(("NaN returned from 3rd order "
                                       "polynomial fit, pix with nan for "
                                       "lin coeffs"))
                                print("coordinates:",sigpix[0][i],sigpix[1][i])
                                print("signal:",signal)

                        if ((cr<=4) & (cr > 1)):
                            slope = (signal[cr-1] - signal[0]) / (times[cr-1] - times[0])
                            offset = signal[0] - slope*times[0]
                            if np.isnan(offset):
                                print(("NaN returned from 2-group fit in "
                                       "case where lin coeffs are nan"))
                                print("coordinates:",sigpix[0][i],sigpix[1][i])
                                print("signal:",signal)

                        if cr == 1:
                            slope = 0
                            offset = signal[0]
                            if np.isnan(offset):
                                print(("NaN returned from 1-group fit in case"
                                       " where nonlin coeffs are nan"))
                                print("coordinates:",sigpix[0][i],sigpix[1][i])
                                print("signal:",signal)


                        cflag = 'blue'
                else:
                    # Quick check for cosmic ray hits
                    cr = cr_check(times,signal,5)

                    if cr > 2:
                        # Now line-fit
                        offset, slope, offset_err, slope_err = self.linfit_yerror(times[0:cr],signal[0:cr],uncertainty[0:cr])
                        if np.isnan(offset):
                            print(("NaN returned from line-fit in case of "
                                   "pixel with <500 counts"))
                            print("coordinates:",sigpix[0][i],sigpix[1][i])
                            print("signal:",signal[0:cr])
                            print("error:",uncertainty[0:cr])
                            print("time:",times[0:cr])
                    if cr == 2:
                        slope = (signal[cr-1] - signal[0]) / (times[cr-1] - times[0])
                        offset = signal[0] - slope*times[0]
                        if np.isnan(offset):
                            print(("NaN returned from 2-group fit in pixel "
                                   "with signal under 500 counts"))
                            print("coordinates:",sigpix[0][i],sigpix[1][i])
                            print("signal:",signal)

                    if cr == 1:
                        slope = 0
                        offset = signal[0]
                        if np.isnan(offset):
                            print(("NaN returned from 1-group fit for pixel "
                                   "with signal under 500 counts"))
                            print("coordinates:",sigpix[0][i],sigpix[1][i])
                            print("signal:",signal)
                    cflag = 'green'

                f,a = plt.subplots()
                a.plot(times,signal,'ro')
                a.plot(times[0:cr],signal[0:cr],'bo')
                a.plot(times,offset+times*slope,linestyle='--',color=cflag)
                a.set_xlabel('Time (sec)')
                a.set_ylabel('Signal (DN)')
                a.set_title('offset: '+str(offset)+'  slope: '+str(slope))

                # Stick these plots in a subdirectory, since
                # there can be a lot of them
                if os.path.isdir('pixels_with_signal') == False:
                    os.mkdir('pixels_with_signal')

                f.savefig('pixels_with_signal/linfit_pixel_amp_'+str(amp)+'_'+str(sigpix[0][i])+'_'+str(sigpix[1][i])+'.png')
                plt.close(f)
                
            biasframe[sigpix[0][i],sigpix[1][i]] = offset
            biasframeerr[sigpix[0][i],sigpix[1][i]] = offset_err
            print("Pixel tracking, pixel with signal.")
            print(("Pixel ({},{}) = {}"
                   .format(checkx,checky,biasframe[checkx-4,checky-4])))

            nancheck = np.ravel(np.isnan(biasframe))
            print("after pixels_with_signal, number of nans: ",len(np.where(nancheck == True)[0]))

        return biasframe,biasframeerr

      
    def cr_check(self,time,signal,nsig):
        # Given a list of signals and times, do a crude CR check
        # Return the read number in which the CR appears, or
        # the total number of reads if none.
        rates = (np.roll(signal,-1) - signal) / (np.roll(time,-1) - time)
        rates = rates[0:-1]
        goodpix, toolow, toohigh = sigmaclip(rates)
        mnrate = np.mean(goodpix)
        devrate = np.std(goodpix)
    
        cr = np.where((rates < (mnrate-nsig*devrate)) | (rates > (mnrate+nsig*devrate)))[0]
        if len(cr) > 0:
            # Add 1 to go from rates indices to signal indices
            cr = cr[0] + 1
        else:
            cr = len(signal)
        return cr
        

    def create_superbias(self,tup):
        # 1.read in list file
        # 2.for each file in the listfile, read in the (fits) file
        # 3.go row by row through the image
        # 4.create all group_x-group_y differences
        # 5.calculate the mean and uncertainty within an array of boxsize pixels,for each difference.
        # 6.collect means and uncertainties as a matrix and pass to NN2 engine.
        # 7.apply NN2 results to pixels within a width box in all groups
        # 8.average together all groups
        # 9.average together all files
        listfile,boxsize,width,showplots,saveplots,outdir,outfile,fulledge,skipnn2,groupstart,groupend,integration,ipc,xtalk = tup

        # Read in list of files to use
        files = []
        with open(listfile,'r') as fin:
            for line in fin:
                files.append(line.strip())
        
        # Check dimensions of files we'll be working with so we can set up final variables
        ext = files[0][-4:]
        if ext == 'fits':
            h = fits.open(files[0])
            xin = h[1].header['NAXIS1']
            yin = h[1].header['NAXIS2']
            detector = h[0].header['DETECTOR']
            tgroup = h[0].header['TGROUP']
            h.close()
            allbias = np.zeros((1,yin,xin))
            allbiaserr = np.zeros_like(allbias)

        # Create a superbias from each file
        for file in files:
            ext = file[-4:]
            if ext == 'fits':
                biasramp,biasramperr = self.bias_from_file(file,boxsize,width,groupstart,groupend,integration,skipnn2,ipc,xtalk,showplots=showplots,saveplots=saveplots)
                biasframe,biasframeerr = self.collapse(biasramp,biasramperr)

                # Only send the science pixels through pixels_with_signal
                tmpframe,tmpframeerr = self.pixels_with_signal(biasramp[:,4:2044,4:2044],biasramperr[:,4:2044,4:2044],biasframe[4:2044,4:2044],biasframeerr[4:2044,4:2044],detector,tgroup)
                # Add the updated science pixels back to the refpix
                biasframe[4:2044,4:2044] = tmpframe
                biasframeerr[4:2044,4:2044] = tmpframeerr
                del biasramp,biasramperr
                biasframe = np.expand_dims(biasframe,axis=0)
                biasframeerr = np.expand_dims(biasframeerr,axis=0)
                allbias = np.append(allbias,biasframe,axis=0)
                allbiaserr = np.append(allbiaserr,biasframeerr,axis=0)
                del biasframe,biasframeerr
            else:
                print("WARNING. Expecting fits files. Quitting.")
                sys.exit()

        # Remove the initial read from allbias and
        # allbiaserr, which is all zeros
        allbias = allbias[1:,:,:]
        allbiaserr = allbiaserr[1:,:,:]

        # If more than one fits file was used, calculate the mean across all of the bias files
        if len(files) > 1:
            finalbias,finalbiaserr = self.collapse(allbias,allbiaserr)
            finaldq = self.unreliable_bias(finalbiaserr,sigmacut=sigmacut)
        else:
            finalbias = allbias[0,:,:]
            finalbiaserr = allbiaserr[0,:,:]
            finaldq = np.zeros((2048,2048))

        # Create dq_def 
        finaldqdef = []
        dqssb={'DO_NOT_USE':np.uint8(1),'UNRELIABLE_BIAS':np.uint8(2)}
        for bitname in dqssb:
            bitvalue = dqssb[bitname]
            bitnumber = int(np.log(bitvalue)/np.log(2))
            newrow = (bitnumber,bitvalue,bitname,'')
            finaldqdef.append(newrow)
        del allbias,allbiaserr
        
        # Save results
        slashpos = outdir.rfind('/')
        if slashpos != len(outdir)-1:
            outdir = outdir + '/'
        
        if outfile == None:
            outfile = 'superbias_from_integration'+str(integration)+'_of_'+listfile+'_fromGroups'+str(groupstart)+'to'+str(groupend)+'.fits'
            if ((ipc == True) & (xtalk == True)):
                outfile = 'superbias_withIPCXTALK_from_integration'+str(integration)+'_of_'+listfile+'_fromGroups'+str(groupstart)+'to'+str(groupend)+'.fits' 
            if ((ipc == True) & (xtalk == False)):
                outfile = 'superbias_withIPC_from_integration'+str(integration)+'_of_'+listfile+'_fromGroups'+str(groupstart)+'to'+str(groupend)+'.fits' 
            if ((ipc == False) & (xtalk == True)):
                outfile = 'superbias_withXTALK_from_integration'+str(integration)+'_of_'+listfile+'_fromGroups'+str(groupstart)+'to'+str(groupend)+'.fits' 

        # Create an instance of a SuperbiasModel and place
        #the data into it.
        self.save_superbias(finalbias,finalbiaserr,finaldq,finaldqdef,files,outdir,outfile)

        # Fits file format checks
        check_ssb = fits.open(outdir+outfile)
        print(check_ssb.info())
        print(check_ssb['DQ_DEF'].data)

        # Redcat team checks
        subprocess.call(['fitsverify',outdir+outfile])
            

    def save_superbias(bias,err,dq,dqdef,files,outdir,outfile):
        # From jwst_lib.models import SuperBiasModel
        # Save the new superbias in an SSB-formatted reference file

        # Get header of one of the input files for reference
        h = fits.open(files[0])
        header0 = h[0].header
        header1 = h[1].header
        h.close()

        # Create the superbias model instance to use
        sbmodel = SuperBiasModel()
        sbmodel.data = bias
        sbmodel.err = err
        sbmodel.dq = dq
        sbmodel.dq_def = dqdef

        # Other required keywords
        sbmodel.meta.reffile.type = 'SUPERBIAS'
        sbmodel.meta.subarray.name = 'FULL'
        sbmodel.meta.subarray.xstart = 1
        sbmodel.meta.subarray.xsize = header1['NAXIS1']
        sbmodel.meta.subarray.ystart = 1
        sbmodel.meta.subarray.ysize = header1['NAXIS2']
        sbmodel.meta.instrument.name = 'NIRCAM'
        detector = str(header0['DETECTOR'])
        if detector == 'NRCA5':
            detector = 'NRCALONG'
        if detector == 'NRCB5':
            detector = 'NRCBLONG'
        sbmodel.meta.instrument.detector = str(detector)

        # Look for the fastaxis and slowaxis keywords in the input data.
        # if they are present propogate those values into the bad pixel
        # mask. If they are not present, then you must be working with 
        # native orientation data, so use the appropriate values
        inhdu = fits.open(files[0])
        try:
            sbmodel.meta.subarray.fastaxis = inhdu[0].header['FASTAXIS']
            sbmodel.meta.subarray.slowaxis = inhdu[0].header['SLOWAXIS']
        except KeyError:
            print('===============================================')
            print(("FASTAXIS and SLOWAXIS header keywords not found "
                   "in the input data."))
            print(("Assuming they are in native (fitswriter) "
                   "orientation, and adding the"))
            print(("native orientation values for those keywords "
                   "to the static pixel mask."))
            print('===============================================')
            sbmodel.meta.subarray.fastaxis = 1
            sbmodel.meta.subarray.slowaxis = 2

        sbmodel.meta.reffile.author = 'Hilbert'
        sbmodel.meta.reffile.description = 'Superbias reffile from CV3 data'
        sbmodel.meta.reffile.pedigree = 'GROUND'
        sbmodel.meta.reffile.useafter = '2015-10-01'
        sbmodel.meta.exposure.readpatt = fits.getval(files[0],'READPATT')

        # HISTORY keyword
        sbmodel.history.append('Description of Reference File Creation')
        sbmodel.history.append('DOCUMENT:')
        sbmodel.history.append('JWST-STScI-TR-XXXX')
        sbmodel.history.append('SOFTWARE:')
        sbmodel.history.append('/ifs/jwst/wit/witserv/data4/nrc/')
        sbmodel.history.append('hilbert/superbias/cv3/B1/superbias_create_ngroups.py')
        # Put the list of input files into the HISTORY keyword
        sbmodel.history.append('DATA USED:')
        for file in files:
            totlen = len(file)
            div = np.arange(0,totlen,60)
            for val in div:
                if totlen > (val+60):
                    sbmodel.history.append(file[val:val+60])
                else:
                    sbmodel.history.append(file[val:])
        
        sbmodel.history.append('DIFFERENCES:')
        sbmodel.history.append('N/A. No previous version.')

        # Save
        sbmodel.save(outdir+outfile)


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description="Using the NN2 method, create a superbias from a list of input ramps")
        parser.add_argument("listfile",help="file containing a list of fits ramps to process")
        parser.add_argument("-v","--verbose",help="detailed outputs displayed",action="store_true",default="False")
        parser.add_argument("-b","--boxwidth",help="width of the box (in pixels) over which to calculate average values",default=100,type=int)
        parser.add_argument("-w","--applybox",help="width of the box (in pixels) over which the mean value is subtracted",type=int,default=50)
        parser.add_argument("-e","--fulledge",help="if True, points closer than boxwidth to the edge of the detector will have the mean taken over a full boxwidth box. If false, only pixels within boxwidth of the pixel will be used",action="store_true",default=False)
        parser.add_argument("-o","--outdir",help="output directory",default='.')
        parser.add_argument("-s","--saveplots",help="save output plots",action="store_true",default=False)
        parser.add_argument("-p","--showplots",help="display plots to the screen",action="store_true",default=False)
        parser.add_argument("-f","--outfile",help="file to save bias frame and uncertainty to.",default=None)
        parser.add_argument("-n","--skipnn2",help="if True, the use of NN2 to remove the 1/f noise will be skipped.",action="store_true")
        parser.add_argument("-g","--groupstart",help="Number of the first group to extract to make the superbias",default=0,type=int)
        parser.add_argument("-i","--groupend",help="Number of the final group to extract to make the superbias",default=30,type=int)
        parser.add_argument("-t","--integration",help="For files containing multiple integrations, the integration number to extract and work on.",default=0)
        parser.add_argument("--ipc",help="Perform IPC correction on the input data. Default = False.",action = 'store_true', default=False)
        parser.add_argument("--xtalk",help="Perform cross-talk correction on the input data. Default=False",action='store_true',default=False)
        return parser
        

if __name__ == '__main__':
        
    usagestring = 'USAGE: superbias.py files_to_use.list'
    
    superbias = Superbias()
    parser = add_options(usage=usagestring)
    args = parser.parse_args(namespace=superbias)

    #set up variables
    verbose=args.verbose
    boxsize = args.boxwidth
    appboxsize = args.applybox
    outdir = args.outdir
    saveplots = args.saveplots
    showplots = args.showplots
    outfile = args.outfile
    fulledge = args.fulledge
    listfile = args.listfile
    skipnn2 = args.skipnn2
    groupstart = args.groupstart
    groupend = args.groupend
    integration = args.integration
    ipc = args.ipc
    xtalk = args.xtalk

    if verbose == True:
        verbosefile = 'superbias_verbose_output_for'+listfile+'.txt'
        verb = open(verbosefile,'w')
        verb.write('Calculating a superbias from files: \n')
        verb.write(listfile)
        verb.write('\n')
        verb.write('Parameters: \n')
        verb.write('Boxsize: {}\n'.format(boxsize))
        verb.write('Boxsize means are applied to: {}\n'.format(appboxsize))
        verb.write('Output directory : {}\n'.format(outdir))
        verb.write('Output file: {}\n'.format(outfile))
        verb.write('Skip NN2: {}\n'.format(skipnn2))
        

    output = superbias.create_superbias((args.listfile,args.boxwidth,args.applybox,args.showplots,args.saveplots,args.outdir,args.outfile,args.fulledge,args.skipnn2,args.groupstart,args.groupend,args.integration,args.ipc,args.xtalk))

