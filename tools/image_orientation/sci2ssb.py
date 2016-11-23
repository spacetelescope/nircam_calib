#! /usr/bin/env python

# This script converts the fits files from the NIRCam CRYO runs at Lockheed
# into ssb-conform fits files.

import sys, os,re,math
import optparse,scipy
from astropy.io import fits as pyfits
from jwst_lib import models
import numpy as np 

sys.path.append(os.environ['UAZCONVDIR'])
from nircam2ssb import nircam2ssbclass

class sci2ssbclass(nircam2ssbclass):
    def __init__(self):
        nircam2ssbclass.__init__(self)
        
    def image2ssb(self,inputfilename, outfilebasename='auto',outdir=None,outsuffix=None,outsubdir=None):
        outfilebasename = self.mkoutfilebasename(inputfilename, outfilebasename=outfilebasename,outdir=outdir,outsuffix=outsuffix,outsubdir=outsubdir)
        (self.data,self.hdr)=pyfits.getdata(inputfilename, 0, header=True)
        print 'Input shape:',self.data.shape

        self.runID = self.getRunID(filename=inputfilename)
        
        # How many groups and integrations?
        if self.runID=='TUCSONNEW':
            Nint=1
        else:
            Nint = int(self.hdr['NINT'])
        Ngroup =int(self.hdr['NGROUP'])
        Nframe = int(self.hdr['NFRAME'])
        print 'NGROUP:',Ngroup
        print 'NINT:',Nint
        if (Ngroup*Nint)!=self.data.shape[0]:
            raise RuntimeError,'something is wrong! NGROUP=%d, NINT=%d, sum is not shape[0]=%d' % (Ngroup,Nint,self.data.shape[0])

        # rearrange the data: 4th coordinate is integration
        scinew=scipy.zeros((Nint,Ngroup,self.data.shape[1],self.data.shape[2]), dtype=float)
        for i in xrange(Nint):
            scinew[i,:,:,:]=self.data[i*Ngroup:(i+1)*Ngroup,:,:]
        print 'Output shape:',scinew.shape

        self.outputmodel = models.RampModel(data=scinew)

        #Mask dead pixels 
        #mask = scipy.where(np.isnan(scinew))
        # mask = np.where(np.isnan(a))
        #mask = np.isnan(scinew)
        # a=scinew[mask]
        #print mask
        #sys.exit(0)

        #updates the date string
        self.updatemetadata(inputfilename,reffileflag=False)

        #update detector
        self.cryo_update_meta_detector(reffileflag=False)
        
        #flip the data around to place it in the science orientation expected by SSB
        if self.runID in ['CV2','CV3']:
            self.native_to_science_image_flip()

        if self.runID in ['CV3']:
            self.outputmodel.meta.exposure.readpatt = self.hdr['READOUT']



        self.outputmodel.meta.exposure.nints = Nint
        self.outputmodel.meta.exposure.nframes = Nframe
        #self.outputmodel.meta.exposure.ngroups = self.hdr['NAXIS3']
        self.outputmodel.meta.exposure.ngroups = Ngroup
        self.outputmodel.meta.subarray.name = 'FULL'
        self.outputmodel.meta.subarray.xstart = 1
        self.outputmodel.meta.subarray.xsize = self.hdr['NAXIS1']
        self.outputmodel.meta.subarray.ystart = 1
        self.outputmodel.meta.subarray.ysize = self.hdr['NAXIS2']

        # put the version of nircam2ssb into header
        # self.outputmodel.meta.channel = self.version

        outfilename = outfilebasename
        if not re.search('fits$',outfilename): outfilename += '_uncal.fits'

        print 'Saving %s' % outfilename
        self.outputmodel.save(outfilename)

        return(outfilename)


    def native_to_science_image_flip(self):
        #flip the data to match the orientation expected by SSB
        data = self.outputmodel.data
        
        if self.hdr['DETECTOR'] in ['NRCA2','NRCA4','NRCB1','NRCB3','NRCBLONG']:
            #vertical flip
            flip = data[:,:,::-1]
            self.outputmodel.data = flip
            self.outputmodel.meta.subarray.fastaxis = 1
            self.outputmodel.meta.subarray.slowaxis = -2

        elif self.hdr['DETECTOR'] in ['NRCA1','NRCA3','NRCB2','NRCB4','NRCALONG']:
            #horizontal flip
            flip = data[:,:,:,::-1]
            self.outputmodel.data = flip
            self.outputmodel.meta.subarray.fastaxis = -1
            self.outputmodel.meta.subarray.slowaxis = 2
        else:
            print("WARNING! I don't recognize {} as a valid detector!".format(self.detector))
            sys.exit(0)



if __name__=='__main__':

    usagestring='USAGE: sci3ssb.py infile1 infile2 ...'

    sci2ssb=sci2ssbclass()
    parser = sci2ssb.add_options(usage=usagestring)
    options,  args = parser.parse_args()

    if len(args)<1:
        parser.parse_args(['--help'])
        sys.exit(0)                        

    sci2ssb.verbose=options.verbose
    
    for infile in args:
        sci2ssb.image2ssb(infile,outfilebasename=options.outfilebasename,outdir=options.outdir,outsuffix=options.outsuffix,outsubdir=options.outsubdir)
        if options.verbose>=3:
            print 'output image attributes:'
            for i in  sci2ssb.outputmodel.iteritems():
                print i
           
