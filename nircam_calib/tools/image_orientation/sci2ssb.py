#! /usr/bin/env python

# This script converts the fits files from the NIRCam CRYO runs
# into ssb-conform fits files.

import sys, os,re,math
import optparse,scipy
from astropy.io import fits as pyfits
from jwst import datamodels as models
from astropy.io import ascii
import numpy as np
from nircam2ssb import nircam2ssbclass



# Before running it, make sure to set environment variables:
cmd1 = "export UAZCONVDIR='/grp/jwst/wit/nircam/nircam-tools/pythonmodules/'"
cmd2 = "export JWSTTOOLS_PYTHONMODULES='$JWSTTOOLS_ROOTDIR/pythonmodules'"
cmd3 = "export JWSTTOOLS_ROOTDIR='/grp/jwst/wit/nircam/nircam-tools'"
cmd4 = "export JWSTTOOLS_INSTRUMENT='NIRCAM'"

try:
    sys.path.append(os.environ['UAZCONVDIR'])
except KeyError:
    for cmd in [cmd1,cmd2,cmd3,cmd4]:
        os.system(cmd)


# Subarray list from latest OSS summary
# has format: AperName Detector ColCorner RowCorner Ncols Nrows
subarrays = ascii.read("NIRCam_raw_coord_subarray_definitions.dat")


class sci2ssbclass(nircam2ssbclass):
    def __init__(self):
        nircam2ssbclass.__init__(self)

    def image2ssb(self,inputfilename, outfilebasename='auto',outdir=None,outsuffix=None,outsubdir=None):
        outfilebasename = self.mkoutfilebasename(inputfilename, outfilebasename=outfilebasename,outdir=outdir,outsuffix=outsuffix,outsubdir=outsubdir)
        (self.data,self.hdr)=pyfits.getdata(inputfilename, 0, header=True)
        print('input shape:',self.data.shape)

        self.runID = self.getRunID(filename=inputfilename)
        self.hdr['SUBARRAY'] = 'FULL'

        # How many groups and integrations?
        if self.runID=='TUCSONNEW':
            Nint=1
        else:
            Nint = int(self.hdr['NINT'])
        Ngroup =int(self.hdr['NGROUP'])
        Nframe = int(self.hdr['NFRAME'])
        print('NGROUP:',Ngroup)
        print('NINT:',Nint)
        if (Ngroup*Nint)!=self.data.shape[0]:
            print('something is wrong! NGROUP=%d, NINT=%d, sum is not shape[0]=%d' % (Ngroup,Nint,self.data.shape[0]))
            raise RuntimeError
        # rearrange the data: 4th coordinate is integration
        scinew=scipy.zeros((Nint,Ngroup,self.data.shape[1],self.data.shape[2]), dtype=float)
        for i in range(Nint):
            scinew[i,:,:,:]=self.data[i*Ngroup:(i+1)*Ngroup,:,:]
        print('output shape:',scinew.shape)

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
        print('meta data updated')
        #update detector
        self.cryo_update_meta_detector(reffileflag=False)
        print('cryo meta data updated')
        #flip the data around to place it in the science orientation expected by SSB
        if self.runID in ['CV2','CV3', 'OTIS']:
            self.native_to_science_image_flip()

        if self.runID in ['CV3','OTIS']:
            self.outputmodel.meta.exposure.readpatt = self.hdr['READOUT']


        self.outputmodel.meta.exposure.nints = Nint
        self.outputmodel.meta.exposure.nframes = Nframe
        #self.outputmodel.meta.exposure.ngroups = self.hdr['NAXIS3']
        self.outputmodel.meta.exposure.ngroups = Ngroup

        # put the version of nircam2ssb into header
        # self.outputmodel.meta.channel = self.version

        outfilename = outfilebasename
        if not re.search('fits$',outfilename): outfilename += '_uncal.fits'

        print('Saving %s' % outfilename)
        self.outputmodel.save(outfilename)

        return(outfilename)


    def native_to_science_image_flip(self):
        #flip the data to match the orientation expected by SSB
        data = self.outputmodel.data
        print('trying to flip the image')

        ncols = self.hdr['NAXIS1']
        nrows = self.hdr['NAXIS2']
        endcol = int(self.hdr['COLCORNR']) + int(ncols) - 1
        endrow = int(self.hdr['ROWCORNR']) + int(nrows) - 1

        print('detector is: ',self.hdr['DETECTOR'])

        # NIRCAM A2, A4, B1, B3, BLONG
        # Flip vertically
        if self.hdr['DETECTOR'] in ['NRCA2','NRCA4','NRCB1','NRCB3','NRCBLONG']:
            flip = data[:,:,::-1]
            self.outputmodel.data = flip
            self.outputmodel.meta.subarray.fastaxis = 1
            self.outputmodel.meta.subarray.slowaxis = -2
            if ncols == 2048 and nrows == 2048:
                xstrt = 1
                ystrt = 2048
                xend = 1
                yend = 2048
            else:
                xstrt = int(self.hdr['COLCORNR'])
                ystrt = int(2048 - self.hdr['ROWCORNR'])
                xend = int(endcol)
                yend = int(2048 - endrow)


        # NIRCAM A1, A3, ALONG, B2, B4
        # Flip horizontally
        elif self.hdr['DETECTOR'] in ['NRCA1','NRCA3','NRCB2','NRCB4','NRCALONG']:
            print('original rowcorner,colcorner: ',self.hdr['ROWCORNR'],self.hdr['COLCORNR'])
            flip = data[:,:,:,::-1]
            self.outputmodel.data = flip
            self.outputmodel.meta.subarray.fastaxis = -1
            self.outputmodel.meta.subarray.slowaxis = 2
            if ncols == 2048 and nrows == 2048:
                xstrt = 1
                ystrt = 2048
                xend = 1
                yend = 2048
            else:
                xstrt = int(2048 - self.hdr['COLCORNR'])
                ystrt = int(self.hdr['ROWCORNR'])
                xend = int(2048 - endcol)
                yend = int(endrow)

        else:
            print("WARNING! I don't recognize {} as a valid detector!".format(self.detector))
            sys.exit(0)

        if xstrt > xend:
            subxstrt = xend
            subxend = xstrt

        elif xstrt < xend:
            subxstrt = xstrt
            subxend = xend

        if ystrt > yend:
            subystrt = yend
            subyend = ystrt

        elif ystrt < yend:
            subystrt = ystrt
            subyend = yend

        print('final rowcorner,colcorner',subxstrt,subystrt)
        self.outputmodel.meta.subarray.xstart = subxstrt
        self.outputmodel.meta.subarray.xsize = ncols
        self.outputmodel.meta.subarray.ystart = subystrt
        self.outputmodel.meta.subarray.ysize = nrows


        #  Update the subarray parameters
        print('trying to get the subarray name')        
        subarray_name = self.get_subarray_name(subarrays,self.hdr['DETECTOR'],self.hdr['COLCORNR'], self.hdr['ROWCORNR'])
        self.outputmodel.meta.subarray.name = subarray_name


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
