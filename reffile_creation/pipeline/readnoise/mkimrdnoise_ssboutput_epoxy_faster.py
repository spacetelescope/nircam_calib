#!/usr/bin/env python

import sys, os,re
from astropy.io import fits
#import pyfits
import numpy as np
import subprocess
import argparse
import scipy
import math
import datetime, glob
import copy
from itertools import izip

# put the tools directory into the path
# add pythonmodules to PATH
if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
else:
    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
    sys.exit(0)

from sigmacut import calcaverageclass
from texttable import txttableclass
from tools import rmfile
from nircam2ssb import nircam2ssbclass
from jwst_lib import models
#nircam_read import NRC_RampModel

lowedir = '/grp/jwst/wit/nircam/hilbert/epoxy_thinned_area_definition/'
lowepoxy = glob.glob(lowedir+'*CV3_low_epoxy_region_DMSorientation.fits')

class mkimrdnoiseclass(nircam2ssbclass):
    def __init__(self):
        self.filename = None
        nircam2ssbclass.__init__(self)

    #def add_options_orig(self, parser=None, usage=None):
    #    if parser == None:
    #        parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

    #    parser.add_option('-v','--verbose', default=False, action="store_true",
    #                      help='Verbose output')
    #    parser.add_option('-d','--debug', default=False, action="store_true",
    #                      help='Debugging output')
    #    parser.add_option('-o','--outfile'  , default='rdnoise.fits' , type="string",
    #                      help='file name of output file (default=%default)')
    #    parser.add_option('--bpm'  , default=None , type="string",
    #                      help='file name of bad pixel mask (default=%default)')
    #    parser.add_option('--bpmval'  , default=0xffff , type="int",
    #                      help='pixels with bpm&bpmval>0 are masked (default=%default)')
    #    parser.add_option('--xmin'  , default=None , type="int",
    #                      help='xmin for stats (default=%default)')
    #    parser.add_option('--xmax'  , default=None , type="int",
    #                      help='xmax for stats (default=%default)')
    #    parser.add_option('--ymin'  , default=None , type="int",
    #                      help='ymin for stats (default=%default)')
    #    parser.add_option('--ymax'  , default=None , type="int",
    #                      help='ymax for stats (default=%default)')
    #    parser.add_option('--forcexylimits', default=False, action="store_true",
    #                      help='Do not use any pixels outside the xylimits, even if possible')
    #    parser.add_option('--boxsizex'  , default=128 , type="int",
    #                      help='boxsize for stats in x direction (default=%default)')
    #    parser.add_option('--boxsizey'  , default=128 , type="int",
    #                      help='boxsize for stats in y direction (default=%default)')
    #    parser.add_option('--stepsizex'  , default=None , type="int",
    #                      help='stepsize between positions in x direction (default=%default)')
    #    parser.add_option('--stepsizey'  , default=None , type="int",
    #                      help='stepsize between positions in y direction (default=%default)')
    #    parser.add_option('--gmin'  , default=None , type="int",
    #                      help='minimum group used. If None, then group 0 (default=%default)')
    #    parser.add_option('--gmax'  , default=None , type="int",
    #                      help='maximum group used. If None, then last group (default=%default)')
    #    parser.add_option('-f','--fill', default=False, action="store_true",
    #                      help='if stepsize>1: fill out the output matrix')
    #    parser.add_option('--Pclipmax'  , default=10.0 , type="float",
    #                      help='maximum % (0-100) of pixels clipped for statistic (default=%default)')
    #    parser.add_option('--Npixmin'  , default=32 , type="int",
    #                      help='minimum #  of valid pixels in box required. If None, then  (default=%default)')
    #    parser.add_option('--invert_epoxy', default = False, action='store_true', 
    #                      help="Invert the epoxy void mask, in order to calculate readnoise values inside the void region.")
    #    parser.add_option('--fill_iterations', default = 4, type="int", help="Number of times to run the filling function that calculates readnoise values for boxes with too few pixels by calculating the mean of surrounding good boxes.")
    #
    #    return(parser)

    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='calculate readnoise')
        #parser.add_argument("listfile",help="Name of file containing the list of files to calculate readnoise for.")
        parser.add_argument("infile",help="Name of file to process")
        parser.add_argument("outfilebasename",help="Name of output file", default = 'NRC_readnoise.fits')
        parser.add_argument('-v','--verbose', default=False, action="store_true",help='Verbose output')
        parser.add_argument('-d','--debug', default=False, action="store_true",help='Debugging output')
        parser.add_argument('-o','--outfile', default='rdnoise.fits',help='file name of output file')
        parser.add_argument('--bpm'  , default=None ,help='file name of bad pixel mask')
        parser.add_argument('--bpmval'  , default=0xffff ,type=int,help='pixels with bpm&bpmval>0 are masked.')
        parser.add_argument('--xmin'  , default=None , type=int,help='xmin for stats.')
        parser.add_argument('--xmax'  , default=None , type=int,help='xmax for stats.')
        parser.add_argument('--ymin'  , default=None , type=int,help='ymin for stats.')
        parser.add_argument('--ymax'  , default=None , type=int,help='ymax for stats.')
        parser.add_argument('--forcexylimits', default=False, action="store_true",help='Do not use any pixels outside the xylimits, even if possible')
        parser.add_argument('--boxsizex'  , default=128 , type=int,help='boxsize for stats in x direction.')
        parser.add_argument('--boxsizey'  , default=128 , type=int,help='boxsize for stats in y direction.')
        parser.add_argument('--stepsizex'  , default=None , type=int,help='stepsize between positions in x direction.')
        parser.add_argument('--stepsizey'  , default=None , type=int,help='stepsize between positions in y direction.')
        parser.add_argument('--gmin'  , default=None , type=int,help='minimum group used. If None, then group 0.')
        parser.add_argument('--gmax'  , default=None , type=int,help='maximum group used. If None, then last group.')
        parser.add_argument('-f','--fill', default=False, action="store_true",help='if stepsize>1: fill out the output matrix')
        parser.add_argument('--Pclipmax'  , default=10.0 , type=float,help='maximum % (0-100) of pixels clipped for statistic.')
        parser.add_argument('--Npixmin'  , default=32 , type=int,help='minimum # of valid pixels in box required. If None, then .')
        parser.add_argument('--invert_epoxy', default = False, action='store_true',help="Invert the epoxy void mask, in order to calculate readnoise values inside the void region.")
        parser.add_argument('--fill_iterations', default = 4, type=int, help="Number of times to run the filling function that calculates readnoise values for boxes with too few pixels by calculating the mean of surrounding good boxes.")

        return parser


    def mkimrdnoise(self,infile,outfilebasename,pattern='.fits',format='g%03d'):
        if self.verbose:
            print 'Loading ',infile
        #self.data, self.hdr = pyfits.getdata(infile, 0, header=True)

        #astropy.fits
        with fits.open(infile) as h:
            self.data = h[1].data
            self.hdr = h[1].header
            self.hdr0 = h[0].header

        #read in via RampModel
        #inramp = NRCRampModel()
        #indata = inramp.loadimage(infile)
        #self.data = indata.data
        #self.hdr = ??indata.header
        #self.hdr0 = ??indata.header

        if len(self.data.shape)==4:
            print('WARNING: Input data cube has 4 dimensions.')
            print('Extracting the first integration only and continuing.')
            self.data = self.data[0,:,:,:]
        elif len(self.data.shape) < 3:
            print 'ERROR: data cube has less than 3 dimensions!!! Quitting!'
            sys.exit(0)

        Ngroups = self.data.shape[0]
        if self.verbose:
             print 'Ngroups ',Ngroups
        #outfilebasename=re.sub('\.fits$','',outfilebasename)


        stepsizex = self.stepsizex
        stepsizey = self.stepsizey
        if stepsizex == None: stepsizex = self.boxsizex
        if stepsizey == None: stepsizey = self.boxsizey

        halfstepsizex = int(0.5*stepsizex)
        halfstepsizey = int(0.5*stepsizey)

        # get the xy limits.
        if (self.xmin==None): xmin=0
        else: xmin=self.xmin
        if (self.xmax==None): xmax=self.data.shape[2]
        else: xmax=self.xmax
        if (self.ymin==None): ymin=0
        else: ymin=self.ymin
        if (self.ymax==None): ymax=self.data.shape[1]
        else: ymax=self.ymax
        if self.verbose:
            print 'xmin, xmax, ymin, ymax:',xmin,xmax,ymin,ymax

        sigmacut = calcaverageclass()
        diffim = self.data[1:Ngroups,:,:]- self.data[0:Ngroups-1,:,:]

        # define output matrices
        im_rdnoise = scipy.zeros((self.data.shape[1],self.data.shape[2]), dtype=float)
        im_rdnoise_err = scipy.zeros(im_rdnoise.shape, dtype=float)
        im_Nused = scipy.zeros(im_rdnoise.shape, dtype=int)
        im_Pskipped = scipy.zeros(im_rdnoise.shape, dtype=float)

        # load mask
        if self.bpm != None:
            print 'loading ',self.bpm
            print('BPM read in via astropy for the moment, to avoid MaskModel bug!!!!!')
            #self.bpmdata = pyfits.getdata(self.bpm, 0).astype('int_')
            with fits.open(self.bpm) as bpmhdu:
                self.bpmdata = bpmhdu[1].data

            #bpminstance = MaskModel(self.bpm)
            #self.bpmdata = bpminstance.dq
            if self.bpmdata.shape != (2048,2048):
                print("WARNING! BPM shape is incorrect or was improperly read!")
                sys.exit()

            #print self.bpmdata
            mask = scipy.where(scipy.logical_and(self.bpmdata,self.bpmval)>0 ,True,False)

            # mask the overscan
            mask[0:4,:]=1
            mask[2044:2048,:]=1
            mask[:,0:4]=1
            mask[:,2044:2048]=1

            #print mask
        else:
            mask =  scipy.zeros((self.data.shape[1],self.data.shape[2]),dtype=np.uint16)


        #load epoxy mask if needed
        detector = self.hdr0['DETECTOR']
        if detector in ['NRCA1','NRCA3','NRCA4','NRCALONG','NRCB1','NRCB4']:
            epoxymask = self.load_epoxy_mask(detector)

            #invert the epoxy mask if requested
            if self.invert_epoxy:
                print("Inverting the epoxy mask to work within the epoxy void region.")
                epoxymask = epoxymask - 1
                epoxymask[epoxymask == -1] = 1

            #add the epoxy mask to the bpm
            mask = scipy.logical_or(mask,epoxymask)

            #create an inverted epoxy mask for later when filling in readnoise values
            inv_epoxymask = epoxymask.astype(bool)
            inv_epoxymask = np.invert(inv_epoxymask)


        x=xmin + halfstepsizex
        OneoverSQR2 = 1.0/math.sqrt(2)
        xtransitions = [512,2*512,3*512]

        if self.gmin==None:
            gmin=0
        else:
            gmin=self.gmin
        if self.gmax==None:
            gmax=diffim.shape[0]
        else:
            gmax=self.gmax
            if gmax>diffim.shape[0]:
                gmax=diffim.shape[0]



        #lists of starting and stopping indexes for the filling later
        xfills = []
        xfille = []
        yfills = []
        yfille = []

        grange = xrange(gmin,gmax)

        while x<xmax:

            # define the x1,x2 of the box for stats
            x1 = x-int(0.5*self.boxsizex)
            if x1<0: x1=0
            if self.forcexylimits and (x1<xmin): x1=xmin


            x2 = int(x+max(1,0.5*self.boxsizex))
            if x2>self.data.shape[2]: x2=self.data.shape[2]
            if self.forcexylimits and (x2>xmax): x2=xmax

            # make sure that the box contains only data from the same amp!!!
            for xtransition in xtransitions:
                if x>=xtransition and x1<xtransition:x1=xtransition
                if x< xtransition and x2>=xtransition:x2=xtransition
            print '!!!x',x,0.5*self.boxsizex,x1,x2

            y=ymin + halfstepsizey
            while y<ymax:

                # define the y1,y2 of the box for stats
                y1 = y-int(0.5*self.boxsizey)
                if y1<0: y1=0
                if self.forcexylimits and (y1<ymin): y1=ymin

                y2 = int(y+max(1,0.5*self.boxsizey))
                if y2>self.data.shape[2]: y2=self.data.shape[1]
                if self.forcexylimits and (y2>ymax): y2=ymax

                if (x % 64) == 0 and (y == 0):
                    print '(x,y)=(%d,%d) box:%d:%d,%d:%d' % (x,y,x1,x2,y1,y2)

                stdevs = []
                Nused = []
                Pskipped = []
                for g in grange:
                    sigmacut.calcaverage_sigmacutloop(diffim[g,y1:y2,x1:x2],mask=mask[y1:y2,x1:x2],Nsigma=3.0,verbose=0)
                    if self.verbose:
                        print 'x:%d y:%d g:%d' % (x,y,g),sigmacut.__str__()
                    if sigmacut.converged and sigmacut.Nused>self.Npixmin and 100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)<self.Pclipmax:
                        stdevs.append(sigmacut.stdev)
                        Nused.append(sigmacut.Nused)
                        Pskipped.append(100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped))
                #if len(stdevs)>0:
                if len(stdevs)>1:
                    sigmacut.calcaverage_sigmacutloop(np.array(stdevs),Nsigma=3.0,verbose=0)
                    if sigmacut.converged:
                        im_rdnoise[y,x]=sigmacut.mean*OneoverSQR2
                        im_rdnoise_err[y,x]=sigmacut.mean_err*OneoverSQR2
                    if self.verbose:
                        print 'x:%d y:%d average' % (x,y),sigmacut.__str__()
                    im_Nused[y,x]=scipy.median(Nused)
                    im_Pskipped[y,x]=scipy.median(Pskipped)
                elif len(stdevs) == 1:
                    im_rdnoise[y,x]=sigmacut.stdev*OneoverSQR2
                    im_rdnoise_err[y,x]=sigmacut.stdev_err*OneoverSQR2
                    im_Nused[y,x]=1
                    im_Pskipped[y,x]=0.0

                if self.fill:

                    xx1=x-halfstepsizex
                    for xtransition in xtransitions:
                        if x-stepsizex<xtransition and x>=xtransition:xx1=xtransition
                    if x-stepsizex<0:xx1=0
                    if xx1<x1:xx1=x1
                    if xx1<0:xx1=0
                    xx2=x+max(1,halfstepsizex)
                    for xtransition in xtransitions:
                        if x+stepsizex>=xtransition and x<xtransition:
                            xx2=xtransition
                    if x+stepsizex>self.data.shape[2]: xx2=self.data.shape[2]
                    if xx2>x2:xx2=x2
                    if xx2>self.data.shape[2]: xx2=self.data.shape[2]

                    yy1=y-halfstepsizey
                    if y-stepsizey<0:yy1=0
                    if yy1<y1:yy1=y1
                    if yy1<0:yy1=0
                    yy2=y+max(1,halfstepsizey)
                    if y+stepsizey>self.data.shape[1]: yy2=self.data.shape[1]
                    if yy2>y2:yy2=y2
                    if yy2>self.data.shape[1]: yy2=self.data.shape[1]

                    #save the x and y coordinates for filling in missing data later
                    xfills.append(xx1)
                    xfille.append(xx2)
                    yfills.append(yy1)
                    yfille.append(yy2)

                    if len(stdevs)>0:
                        #submask = inv_epoxymask[yy1:yy2,xx1:xx2]
                        im_rdnoise[yy1:yy2,xx1:xx2]=im_rdnoise[y,x]#*submask
                        im_rdnoise_err[yy1:yy2,xx1:xx2]=im_rdnoise_err[y,x]#*submask
                        im_Nused[yy1:yy2,xx1:xx2]=im_Nused[y,x]
                        im_Pskipped[yy1:yy2,xx1:xx2]=im_Pskipped[y,x]

                #####
                y+=stepsizey
            x+=stepsizex


        #fill in the gaps in the map (i.e. the partial boxes that contain less than the minimum number
        #of pixels needed for calculating the readnoise
        
        #save pre-filled image
        #h0=fits.PrimaryHDU(im_rdnoise)
        #hhlist = fits.HDUList([h0])
        #tmpfile = 'before_filling_'+infile
        #if os.path.isfile(tmpfile):
        #    os.rename(tmpfile,'outside_'+tmpfile)
        #hhlist.writeto(tmpfile)
        #goo = im_rdnoise > 0
        #print(np.min(im_rdnoise[goo]),np.max(im_rdnoise[goo]))
        #print("BEFORE FILLING EDGES, FILE SAVED TO before_filling*fits")

        for iter in xrange(self.fill_iterations):
            print('BEGINNING FILL ITERATION {}'.format(iter))
            im_rdnoise = self.fill_empty_regions(im_rdnoise,xfills,xfille,yfills,yfille)
            im_rdnoise_err = self.fill_empty_regions(im_rdnoise_err,xfills,xfille,yfills,yfille)

        #print("AFTER FILLING EDGES, FILE SAVED TO after_filling*fits")
        #h0=fits.PrimaryHDU(im_rdnoise)
        #hhlist = fits.HDUList([h0])
        #tmpfile = 'after_filling_'+infile
        #if os.path.isfile(tmpfile):
        #    os.rename(tmpfile,'outside_'+tmpfile)
        #hhlist.writeto(tmpfile)
        #goo = im_rdnoise > 0
        #print(np.min(im_rdnoise[goo]),np.max(im_rdnoise[goo]))


        #now mask out the pixels in the epoxymask
        #print("why is this not being done????")
        #l0=fits.PrimaryHDU(epoxymask)
        #ll=fits.HDUList([l0])
        #ll.writeto('test.fits',clobber=True)
        
        if detector in ['NRCA1','NRCA3','NRCA4','NRCALONG','NRCB1','NRCB4']:
            im_rdnoise[epoxymask.astype(bool)] = 0.
            im_rdnoise_err[epoxymask.astype(bool)] = 0.
        
        #l0=fits.PrimaryHDU(im_rdnoise)
        #ll=fits.HDUList([l0])
        #ll.writeto('test2.fits',clobber=True)
        #stop
        

        #save output file
        outfilename = self.save_readnoise_file(im_rdnoise,infile,outfilebasename)

        #fits format checks
        check_ssb = fits.open(outfilename)
        print(check_ssb.info())
        print(check_ssb['SCI'].data)
        print(check_ssb['SCI'].data[500,500])
        print(im_rdnoise[500,500])

        #redcat team checks
        subprocess.call(['fitsverify',outfilename])


        #hdu_rdnoise = pyfits.PrimaryHDU(im_rdnoise)
        #hdu_rdnoise.header.set('object','readnoise')
        #hdu_rdnoise_err = pyfits.ImageHDU(im_rdnoise_err)
        #hdu_rdnoise_err.header.set('object','readnoise_err')
        #hdu_Nused = pyfits.ImageHDU(im_Nused)
        #hdu_Nused.header.set('object','Nused')
        #hdu_Pskipped = pyfits.ImageHDU(im_Pskipped)
        #hdu_Pskipped.header.set('object','Pskipped')

        #outfilename = re.sub('fits','bx%d.by%d.fits' % (self.boxsizex,self.boxsizey),outfilebasename)

        #hdulist = pyfits.HDUList([hdu_rdnoise,hdu_rdnoise_err,hdu_Nused,hdu_Pskipped])
        #print 'Saving ',outfilename
        #rmfile(outfilename)
        #hdulist.writeto(outfilename)


    def fill_empty_regions(self,image,xstart,xend,ystart,yend):
        '''fill in the boxes that did not have readnoise values calculated because they 
        have less than the minimum number of pixels'''
        fixedimage = copy.deepcopy(image)


        print('INSIDE FILL_EMPTY_REGIONS: ',np.nanmin(image),np.nanmax(image))
        #low = np.where((image>0) & (image < 2))
        #print("LOW: ",len(low[0]))


        #for each empty box, just fill the whole thing with the mean of the nearby boxes
        #repeat several times (set as a tunable parameter)
        #print("Filling empty boxes")

        xstart = np.unique(xstart)
        ystart = np.unique(ystart)
        xend = np.unique(xend)
        yend = np.unique(yend)

        #print(xstart)
        #print(xend)
        #print(ystart)
        #print(yend)

        #sanity check
        if len(xstart) > 66:
            stop

        #loop over boxes
        numx = len(xend)
        numy = len(yend)
        for i in xrange(numx):
            for j in xrange(numy):
                xs = xstart[i]
                xe = xend[i]
                ys = ystart[j]
                ye = yend[j]

                if ((xs > 1120) & (xs < 1130) & (ys > 540) & (ys < 550)):
                    nzero = np.where(image[ys:ye,xs:xe] == 0)[0]
                    print("box check, number of zeros:".format(len(nzero)))
                    


                #print("looking at box {},{}".format(i,j))
                #print("with xstart,xend,ystart,yend of {},{},{},{}".format(xs,xe,ys,ye))

                #if the box has no valid readnoise values....
                if np.nanmax(image[ys:ye,xs:xe]) <= 0.:
        
                    #print("Empty box found. x,y start is {},{}".format(xs,ys))
                    #print("x,y end is {},{}.".format(xe,ye))

                    above = -999
                    below = -999
                    right = -999
                    left = -999

                    if j+1 < numy:
                        if image[ystart[j+1]+1,xs+1] != 0:
                            above = image[ystart[j+1]+1,xs+1]
                            print('ABOVE MINMAX: ',np.nanmin(above),np.nanmax(above))
                    if j-1 > 0:
                        if image[ystart[j-1]+1,xs+1] != 0:
                            below = image[ystart[j-1]+1,xs+1]
                            print('BELOW MINMAX:',np.nanmin(below),np.nanmax(below))
                    if i+1 < numx:
                        if image[ys+1,xstart[i+1]+1] != 0:
                            right = image[ys+1,xstart[i+1]+1]
                            print('right MINMAX:',np.nanmin(right),np.nanmax(right))
                    if i-1 > 0:
                        if image[ys+1,xstart[i-1]+1] != 0:
                            left = image[ys+1,xstart[i-1]+1]
                            print('left MINMAX:',np.nanmin(left),np.nanmax(left))
                        
                    neighbors = np.array([above,below,right,left])
                    print('neighbors:',neighbors)
                    
                    goodpts = np.where(neighbors != -999)
                    if len(goodpts[0]) > 0:
                        print("GOODPoINTS: ",goodpts)
                        mnval = np.nanmean(neighbors[goodpts[0]])
                    else:
                        mnval = 0
                    print('mean neighbors: ',mnval)
                    fixedimage[ys:ye,xs:xe] = mnval
                    
        return fixedimage


    def load_epoxy_mask(self,detector):
        '''Add the epoxy void region to the grid of boxes'''
        
        #read in epoxy void mask
        #efile = lowepoxy[detector]
        efile = [s for s in lowepoxy if detector[3:] in s]

        if len(efile) > 1:
            print("WARNING: Found multiple low-epoxy masks for detector {}. Quitting.".format(detector))
        efile = efile[0]
        print("Using low epoxy mask: {}".format(efile))

        ehdu = fits.open(efile)
        epoxymask = ehdu[1].data
        ehdu.close()
        return epoxymask


    def save_readnoise_file(self,im,infile,outfilebasename):

        #create readnoise model instance and add data
        self.outputmodel = models.ReadnoiseModel()
        self.outputmodel.data = im
        
        #update info on the detector
        try:
            self.outputmodel.meta.instrument.detector = self.hdr0['DETECTOR']
            self.outputmodel.meta.instrument.channel = self.hdr0['CHANNEL']
            self.outputmodel.meta.instrument.module = self.hdr0['MODULE']
        except KeyError:
            print('Keyword detector not found in header of {}. Assuming this file is NOT in SSB format, and attempting to proceed.'.format(infile))
            #get the runID so we know how to populate the detector-related header keywords
            runID = self.getRunID(filename=infile)

            #update the info on the detector
            self.cryo_update_meta_detector(runID=runID,reffileflag=False)
            

        #other required keywords
        self.outputmodel.meta.reffile.type = 'READNOISE'
        self.outputmodel.meta.subarray.name = 'FULL'
        self.outputmodel.meta.subarray.xstart = 1
        self.outputmodel.meta.subarray.xsize = self.hdr['NAXIS1']
        self.outputmodel.meta.subarray.ystart = 1
        self.outputmodel.meta.subarray.ysize = self.hdr['NAXIS2']
        self.outputmodel.meta.reffile.pedigree = 'GROUND'
        self.outputmodel.meta.instrument.name = 'NIRCAM'
        self.outputmodel.meta.telescope = 'JWST'
        self.outputmodel.meta.reffile.author = 'NIRCam Instrument Team'
        self.outputmodel.meta.reffile.description = 'Detector Readnoise maps.'
        try:
            #for fitswriter format files
            self.outputmodel.meta.exposure.readpatt = self.hdr['READOUT']
        except KeyError:
            #for SSB-format files
            self.outputmodel.meta.exposure.readpatt = self.hdr0['READPATT']

        self.outputmodel.meta.reffile.author = 'Hilbert'
        self.outputmodel.meta.reffile.description = 'Readnoise reffile from CV3 data'
        self.outputmodel.meta.reffile.pedigree = 'GROUND'
        self.outputmodel.meta.reffile.useafter = '2015-10-01'

        #look for the fastaxis and slowaxis keywords in the input data.
        #if they are present propogate those values into the bad pixel
        #mask. If they are not present, then you must be working with 
        #native orientation data, so use the appropriate values
        #inhdu = fits.open(infile)
        try:
            self.outputmodel.meta.subarray.fastaxis = self.hdr0['FASTAXIS']
            self.outputmodel.meta.subarray.slowaxis = self.hdr0['SLOWAXIS']
        except KeyError:
            print('===============================================')
            print("FASTAXIS and SLOWAXIS header keywords not found in the input data.")
            print("Assuming they are in native (fitswriter) orientation, and adding the")
            print("native orientation values for those keywords to the static pixel mask.")
            print('===============================================')
            model.meta.subarray.fastaxis = 1
            model.meta.subarray.slowaxis = 2

        #HISTORY keyword
        self.outputmodel.history.append('Description of Reference File Creation')

        self.outputmodel.history.append('DOCUMENT:')
        self.outputmodel.history.append('JWST-STScI-TR-XXXX')

        self.outputmodel.history.append('SOFTWARE:')
        self.outputmodel.history.append('/grp/jwst/wit/nircam/nircam-tools/pythonmodules/')
        self.outputmodel.history.append('mkimrdnoise_ssboutput.py')
        
        #put the list of input files into the HISTORY keyword
        self.outputmodel.history.append('DATA USED:')
        self.outputmodel.history.append(infile)
        
        self.outputmodel.history.append('DIFFERENCES:')
        self.outputmodel.history.append('N/A. No previous version.')

        #save output file
        outfilename = re.sub('fits','bx%d.by%d_ssbreadnoise.fits' % (self.boxsizex,self.boxsizey),outfilebasename)
        print(outfilename+' inside save_readnoise_file')
        if not re.search('fits$',outfilename): outfilename += '_ssbreadnoise.fits'
        print 'Saving %s' % outfilename
        self.outputmodel.save(outfilename)
        return outfilename


if __name__=='__main__':

    print("Starting: "),datetime.datetime.now()

    usagestring='USAGE: mkimrdnoise.py infile outfilebasename'

    mkimrdnoise=mkimrdnoiseclass()
    parser = mkimrdnoise.add_options(usage=usagestring,parser=None)
    args = parser.parse_args(namespace=mkimrdnoise)
    #mkimrdnoise.options,  args = parser.parse_args()

    #if len(args)!=2:
    #    sys.argv.append('--help')
    #    options,  args = parser.parse_args()
    #    sys.exit(0)

    infile = args.infile
    outfilebasename = args.outfilebasename
    mkimrdnoise.mkimrdnoise(infile,outfilebasename)
    print "SUCCESS calcmkimrdnoise.py",datetime.datetime.now()

