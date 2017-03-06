#! /usr/bin/env python

'''
Create a readnoise reference file. 

mkrefs.py version

This version is for creating a single readnoise reference file from
a single input ramp.

Ready for testing. Would probably be cleaner to replace sigmacut references
with the scipy or numpy sigma clipped mean calcultor, although I'm not sure
if those have all the functionality of sigmacut that is used here...

'''

import sys, os
from astropy.io import fits
import numpy as np
import subprocess
import argparse
import scipy
import math
import copy
import glob
from jwst.datamodels import ReadnoiseModel

# put the tools directory into the path
# add pythonmodules to PATH
if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
else:
    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
    sys.exit(0)

from sigmacut import calcaverageclass


class Readnoise():
    def __init__(self):
        self.infile = ''
        self.verbose = False
        self.outfile = None
        self.bpm = None
        self.bpmval = 0xffff
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.forcexylimits = False
        self.boxsizex = 128
        self.boxsizey = 128
        self.stepsizex = None
        self.stepsizey = None
        self.gmin = None
        self.gmax = None
        self.fill = False
        self.Pclipmax = 10.0
        self.Npixmin = 32
        self.invert_epoxy = False
        self.fill_iterations = 4
        self.useafter = ''
        self.description = 'NIRCam Readnoise file'
        self.author = 'NIRCam Instrument Team'
        self.pedigree = 'GROUND'


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='calculate readnoise')
        parser.add_argument("infile",help="Name of file to process")
        parser.add_argument('-v','--verbose', default=False, action="store_true",help='Verbose output')
        parser.add_argument('-o','--outfile', default=None,help='file name of output file')
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
        parser.add_argument('--Npixmin'  , default=32 , type=int,help='minimum # of valid pixels in box required.')
        parser.add_argument('--invert_epoxy', default = False, action='store_true',help="Invert the epoxy void mask, in order to calculate readnoise values inside the void region.")
        parser.add_argument('--fill_iterations', default = 4, type=int, help="Number of times to run the filling function that calculates readnoise values for boxes with too few pixels by calculating the mean of surrounding good boxes.")
        parser.add_argument('--useafter',help='Use after date to put into readnoise reference file',default='')
        parser.add_argument('--description',help='String to place into the DESCRIP header keyword of output file',default='NIRCam Readnoise file')
        parser.add_argument('--author',help='Name of person creating the file',default='NIRCam Instrument Team')
        parser.add_argument('--pedigree',help='Source of the data used to make reffile',default='GROUND')
        return parser



    def make_rdnoise(self):
        if self.verbose:
            print('Loading {}'.format(infile))

        #Read in input fits file using astropy
        with fits.open(self.infile) as h:
            self.data = h[1].data
            self.hdr = h[1].header
            self.hdr0 = h[0].header

        #Keep only the first integration
        if len(self.data.shape)==4:
            print('WARNING: Input data cube has 4 dimensions.')
            print('Extracting the first integration only and continuing.')
            self.data = self.data[0,:,:,:]
        elif len(self.data.shape) < 3:
            print 'ERROR: data cube has less than 3 dimensions!!! Quitting!'
            sys.exit(0)

        #Make sure the data are full frame
        Ngroups,ny,nx = self.data.shape
        if (ny,nx) != (2048,2048):
            print("WARNING: input data from {} appear not to be full frame!".format(self.infile))
            print("x,y dimensions are {}x{}".format(nx,ny))
            sys.exit()
            
        if self.verbose:
             print('Number of groups: {}'.format(Ngroups))

        #Set step sizes
        stepsizex = self.stepsizex
        stepsizey = self.stepsizey
        if stepsizex is None:
            stepsizex = self.boxsizex
        if stepsizey is None:
            stepsizey = self.boxsizey

        halfstepsizex = int(0.5*stepsizex)
        halfstepsizey = int(0.5*stepsizey)

        # get the xy limits.
        if (self.xmin is None):
            xmin=0
        else:
            xmin=self.xmin
        if (self.xmax is None):
            xmax=self.data.shape[2]
        else:
            xmax=self.xmax
        if (self.ymin is None):
            ymin=0
        else:
            ymin=self.ymin
        if (self.ymax is None):
            ymax=self.data.shape[1]
        else:
            ymax=self.ymax
        if self.verbose:
            print('xmin, xmax, ymin, ymax: {}, {}, {}, {}'.format(xmin,xmax,ymin,ymax))

        sigmacut = calcaverageclass()
        #Create a series of CDS frames
        diffim = self.data[1:Ngroups,:,:]- self.data[0:Ngroups-1,:,:]

        # define output matrices
        self.im_rdnoise = scipy.zeros((self.data.shape[1],self.data.shape[2]), dtype=float)
        self.im_rdnoise_err = scipy.zeros(self.im_rdnoise.shape, dtype=float)
        im_Nused = scipy.zeros(self.im_rdnoise.shape, dtype=int)
        im_Pskipped = scipy.zeros(self.im_rdnoise.shape, dtype=float)

        # load mask
        if self.bpm is not None:
            print('Loading mask file {}'.format(self.bpm))

            with fits.open(self.bpm) as bpmhdu:
                self.bpmdata = bpmhdu[1].data

            if self.bpmdata.shape != (2048,2048):
                print("WARNING! Bad pixel mask shape is incorrect or was improperly read!")
                sys.exit()

            mask = scipy.where(scipy.logical_and(self.bpmdata,self.bpmval)>0 ,True,False)

            # mask the overscan (THIS ASSUMES WE ARE WORKING ON FULL FRAME DATA)
            mask[0:4,:]=1
            mask[2044:2048,:]=1
            mask[:,0:4]=1
            mask[:,2044:2048]=1

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
            #inv_epoxymask = epoxymask.astype(bool)
            #inv_epoxymask = np.invert(inv_epoxymask)


        x = xmin + halfstepsizex
        OneoverSQR2 = 1.0/math.sqrt(2)
        xtransitions = [512,2*512,3*512]

        if self.gmin is None:
            gmin=0
        else:
            gmin=self.gmin
        if self.gmax is None:
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

        grange = range(gmin,gmax)

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
            #print('!!!x',x,0.5*self.boxsizex,x1,x2)

            y=ymin + halfstepsizey
            while y<ymax:

                # define the y1,y2 of the box for stats
                y1 = y-int(0.5*self.boxsizey)
                if y1<0: y1=0
                if self.forcexylimits and (y1<ymin): y1=ymin

                y2 = int(y+max(1,0.5*self.boxsizey))
                if y2>self.data.shape[2]: y2=self.data.shape[1]
                if self.forcexylimits and (y2>ymax): y2=ymax

                if self.verbose:
                    if (x % 64) == 0 and (y == 0):
                        print('(x,y)=(%d,%d) box:%d:%d,%d:%d' % (x,y,x1,x2,y1,y2))

                stdevs = []
                Nused = []
                Pskipped = []
                for g in grange:
                    sigmacut.calcaverage_sigmacutloop(diffim[g,y1:y2,x1:x2],mask=mask[y1:y2,x1:x2],Nsigma=3.0,verbose=0)
                    if self.verbose:
                        print('x:%d y:%d g:%d' % (x,y,g),sigmacut.__str__())
                    if sigmacut.converged and sigmacut.Nused>self.Npixmin and 100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)<self.Pclipmax:
                        stdevs.append(sigmacut.stdev)
                        Nused.append(sigmacut.Nused)
                        Pskipped.append(100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped))

                if len(stdevs)>1:
                    sigmacut.calcaverage_sigmacutloop(np.array(stdevs),Nsigma=3.0,verbose=0)
                    if sigmacut.converged:
                        self.im_rdnoise[y,x]=sigmacut.mean*OneoverSQR2
                        self.im_rdnoise_err[y,x]=sigmacut.mean_err*OneoverSQR2
                    if self.verbose:
                        print('x:%d y:%d average' % (x,y),sigmacut.__str__())
                    im_Nused[y,x]=scipy.median(Nused)
                    im_Pskipped[y,x]=scipy.median(Pskipped)
                elif len(stdevs) == 1:
                    self.im_rdnoise[y,x]=sigmacut.stdev*OneoverSQR2
                    self.im_rdnoise_err[y,x]=sigmacut.stdev_err*OneoverSQR2
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
                        self.im_rdnoise[yy1:yy2,xx1:xx2] = self.im_rdnoise[y,x]
                        self.im_rdnoise_err[yy1:yy2,xx1:xx2] = self.im_rdnoise_err[y,x]
                        im_Nused[yy1:yy2,xx1:xx2]=im_Nused[y,x]
                        im_Pskipped[yy1:yy2,xx1:xx2]=im_Pskipped[y,x]


                y+=stepsizey
            x+=stepsizex


        #fill in the gaps in the map (i.e. the partial boxes that contain less than the minimum number
        #of pixels needed for calculating the readnoise
        
        for iter in range(self.fill_iterations):
            #print('BEGINNING FILL ITERATION {}'.format(iter))
            self.im_rdnoise = self.fill_empty_regions(self.im_rdnoise,xfills,xfille,yfills,yfille)
            self.im_rdnoise_err = self.fill_empty_regions(self.im_rdnoise_err,xfills,xfille,yfills,yfille)

        
        #mask out the pixels in the low epoxy area
        if detector in ['NRCA1','NRCA3','NRCA4','NRCALONG','NRCB1','NRCB4']:
            self.im_rdnoise[epoxymask.astype(bool)] = 0.
            self.im_rdnoise_err[epoxymask.astype(bool)] = 0.
        
        #save output file
        outfilename = self.save_readnoise_file(self.im_rdnoise)

        #redcat team checks
        #subprocess.call(['fitsverify',outfilename])
        
        return outfilename,self.im_rdnoise
        

    def fill_empty_regions(self,image,xstart,xend,ystart,yend):
        '''fill in the boxes that did not have readnoise values calculated because they 
        have less than the minimum number of pixels'''
        fixedimage = copy.deepcopy(image)

        #for each empty box, just fill the whole thing with the mean of the nearby boxes
        #repeat several times (set as a tunable parameter)

        xstart = np.unique(xstart)
        ystart = np.unique(ystart)
        xend = np.unique(xend)
        yend = np.unique(yend)

        #loop over boxes
        numx = len(xend)
        numy = len(yend)
        for i in xrange(numx):
            for j in xrange(numy):
                xs = xstart[i]
                xe = xend[i]
                ys = ystart[j]
                ye = yend[j]

                #if ((xs > 1120) & (xs < 1130) & (ys > 540) & (ys < 550)):
                #    nzero = np.where(image[ys:ye,xs:xe] == 0)[0]
                #    print("box check, number of zeros:".format(len(nzero)))

                #if the box has no valid readnoise values....
                if np.nanmax(image[ys:ye,xs:xe]) <= 0.:
        
                    above = -999
                    below = -999
                    right = -999
                    left = -999

                    if j+1 < numy:
                        if image[ystart[j+1]+1,xs+1] != 0:
                            above = image[ystart[j+1]+1,xs+1]
                            #print('ABOVE MINMAX: ',np.nanmin(above),np.nanmax(above))
                    if j-1 > 0:
                        if image[ystart[j-1]+1,xs+1] != 0:
                            below = image[ystart[j-1]+1,xs+1]
                            #print('BELOW MINMAX:',np.nanmin(below),np.nanmax(below))
                    if i+1 < numx:
                        if image[ys+1,xstart[i+1]+1] != 0:
                            right = image[ys+1,xstart[i+1]+1]
                            #print('right MINMAX:',np.nanmin(right),np.nanmax(right))
                    if i-1 > 0:
                        if image[ys+1,xstart[i-1]+1] != 0:
                            left = image[ys+1,xstart[i-1]+1]
                            #print('left MINMAX:',np.nanmin(left),np.nanmax(left))
                        
                    neighbors = np.array([above,below,right,left])
                    #print('neighbors:',neighbors)
                    
                    goodpts = np.where(neighbors != -999)
                    if len(goodpts[0]) > 0:
                        #print("GOODPoINTS: ",goodpts)
                        mnval = np.nanmean(neighbors[goodpts[0]])
                    else:
                        mnval = 0
                    #print('mean neighbors: ',mnval)
                    fixedimage[ys:ye,xs:xe] = mnval
                    
        return fixedimage


    def load_epoxy_mask(self,detector):
        '''Add the epoxy void region to the grid of boxes'''

        lowedir = '/grp/jwst/wit/nircam/hilbert/epoxy_thinned_area_definition/'
        lowepoxy = glob.glob(lowedir+'*CV3_low_epoxy_region_DMSorientation.fits')
        
        #read in epoxy void mask
        #efile = lowepoxy[detector]
        efile = [s for s in lowepoxy if detector[3:] in s]

        if len(efile) > 1:
            print("WARNING: Found multiple low-epoxy masks for detector {}. Quitting.".format(detector))
        efile = efile[0]
        print("Using low epoxy mask: {}".format(efile))

        try:
            with fits.open(efile) as ehdu:
                epoxymask = ehdu[1].data
        except:
            print("Unable to read in exposy mask {}.".format(efile))

        return epoxymask


    def save_readnoise_file(self,im):

        #create readnoise model instance and add data
        self.outputmodel = ReadnoiseModel()
        self.outputmodel.data = im
        
        #update info on the detector
        try:
            self.outputmodel.meta.instrument.detector = self.hdr0['DETECTOR']
            self.outputmodel.meta.instrument.channel = self.hdr0['CHANNEL']
            self.outputmodel.meta.instrument.module = self.hdr0['MODULE']
        except KeyError:
            print('WARNING: Unable to find needed header keywords in {}.'.format(infile))
            sys.exit()
            
        #other required keywords
        self.outputmodel.meta.reffile.type = 'READNOISE'
        self.outputmodel.meta.subarray.name = 'GENERIC'
        self.outputmodel.meta.subarray.xstart = 1
        self.outputmodel.meta.subarray.xsize = self.hdr['NAXIS1']
        self.outputmodel.meta.subarray.ystart = 1
        self.outputmodel.meta.subarray.ysize = self.hdr['NAXIS2']
        self.outputmodel.meta.reffile.pedigree = self.pedigree
        self.outputmodel.meta.instrument.name = 'NIRCAM'
        self.outputmodel.meta.telescope = 'JWST'
        self.outputmodel.meta.reffile.author = self.author
        self.outputmodel.meta.reffile.description = self.description
        self.outputmodel.meta.reffile.useafter = self.useafter

        #look for the fastaxis and slowaxis keywords in the input data.
        #if they are present propogate those values into the bad pixel
        #mask. If they are not present, then you must be working with 
        #native orientation data, so use the appropriate values
        try:
            self.outputmodel.meta.subarray.fastaxis = self.hdr0['FASTAXIS']
            self.outputmodel.meta.subarray.slowaxis = self.hdr0['SLOWAXIS']
        except KeyError:
            print("WARNING: FASTAXIS and SLOWAXIS header keywords not found in the input data.")
            sys.exit()

        #HISTORY keyword
        self.outputmodel.history.append('Description of Reference File Creation')
        self.outputmodel.history.append('DOCUMENT:')
        self.outputmodel.history.append('JWST-STScI-TR-XXXX')
        self.outputmodel.history.append('SOFTWARE:')
        self.outputmodel.history.append('readnoise_reffile.py, a module within')
        self.outputmodel.history.append('makeref.py')
        
        #put the list of input files into the HISTORY keyword
        self.outputmodel.history.append('DATA USED:')
        self.outputmodel.history.append(self.infile)
        self.outputmodel.history.append('DIFFERENCES:')
        self.outputmodel.history.append('N/A. No previous version.')

        #save output file
        if self.outfile is None:
            self.outfile = self.infile[0:-5] + '_readnoise_reffile.fits'

        print('Saving readnoise reference file as {}'.format(self.outfile))
        self.outputmodel.save(self.outfile)
        return self.outfile

        


if __name__=='__main__':

    usagestring='USAGE: readnoise.py infile --outfile readnoise.fits'

    rdnoise = Readnoise()
    parser = rdnoise.add_options(usage=usagestring,parser=None)
    args = parser.parse_args(namespace=rdnoise)

    rdnoise.make_rdnoise()

