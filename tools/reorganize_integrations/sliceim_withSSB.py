#!/usr/bin/env python

# Purpose:  This Python script slices up an image based on the optional group, x, and y parameters

import sys, os,re, copy
import scipy
from astropy.io import fits as pyfits
import numpy as np
import optparse
from tools import rmfile

class sliceimclass:
    def __init__(self):
        self.filename = None
        self.hdulist = None

        self.data=None
        self.hdr=None
        self.outdata=None
        self.verbose=0

        self.ssbformat = None
        self.xindex = self.yindex = self.gindex = self.intindex = None

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
        parser.add_option('-v', '--verbose', action='count', dest='verbose',default=0)
        parser.add_option('--intmin'  , default=None, type="int",
                          help='minimum integration')
        parser.add_option('--intmax'  , default=None, type="int",
                          help='maximum integration')
        parser.add_option('--gmin'  , default=None, type="int",
                          help='minimum group')
        parser.add_option('--gmax'  , default=None, type="int",
                          help='maximum group')
        parser.add_option('--xmin'  , default=None, type="int",
                          help='minimum x')
        parser.add_option('--xmax'  , default=None, type="int",
                          help='maximum x')
        parser.add_option('--ymin'  , default=None, type="int",
                          help='minimum y')
        parser.add_option('--ymax'  , default=None, type="int",
                          help='maximum y')
        #parser.add_option('--ssbformat', default=True, help='type of data')
        parser.add_option('-o','--outfilename'  , default='auto' , type="string",
                          help='file name of output file (default=%default)')
        return(parser)

    def loadimage(self, filename):
        if self.verbose: print 'Loading ',filename
        self.filename = filename

        self.hdulist=pyfits.open(filename)
        if self.verbose>1: 
            print self.hdulist.info()
        if self.hdulist[0].header['bitpix'] in [8,-32]:
            print("SSB format!!")
            # ssb format: 0th extension is just header, so first extension is data...
            shape = self.hdulist[1].data.shape
            if len(shape)!=4:
                raise RuntimeError,"SSB format expected, len(shape) should be 4, but is %d" % len(shape)
            if self.verbose: print 'SSB format!'
            self.ssbformat = True
        else:
            # this must be fitswriter data
            print("FITSWRITER DATA!")
            shape = self.hdulist[0].data.shape
            if len(shape)!=3:
                raise RuntimeError,"Fitswriter format expected, len(shape) should be 3, but is %d" % len(shape)
            if self.verbose: print 'Fitswriter format!'
            self.ssbformat = False

        if self.verbose: print 'image shape:',shape

    def sliceim(self, intmin=None, intmax=None, gmin=None, gmax=None, xmin=None, xmax=None, ymin=None, ymax=None):
        if self.hdulist is None:
            raise RuntimeError,'No input data ...'

        if self.ssbformat:
            data = self.hdulist[1].data
            hdr0 = self.hdulist[0].header
            hdr  = self.hdulist[1].header
            #print hdr
            #print 'bbbbbb'
            #print  self.hdulist[0].header
            Ngroups = hdr0['NGROUPS']
            #self.outdata=None # not used in ssb format
            #self.outhdr=None # not used in ssb format
        else:
            data = self.hdulist[0].data
            hdr0 = self.hdulist[0].header
            hdr  = self.hdulist[0].header
            Ngroups = hdr0['NGROUP']
            self.outhdr=copy.deepcopy(self.hdulist[0].header)

        Nint = hdr0['NINT']
        if self.verbose: 
            print 'Nint: %d' % Nint
            print 'Ngroups: %d' % Ngroups
            
        if gmin==None: gmin=0
        if xmin==None: xmin=0
        if ymin==None: ymin=0
        if intmin==None: intmin=0

        #Ngroups_per_int = int(hdr0['NGROUPS']/hdr0['NINT'])
        #print 'Ngroups per integration:',Ngroups_per_int

        if xmax==None: xmax = hdr['NAXIS1']
        if ymax==None: ymax = hdr['NAXIS2']
        if gmax==None: gmax = Ngroups
        if intmax==None: intmax = Nint

        if gmax>Ngroups:
            raise RuntimeError,'gmax=%d exceeds limit of %d groups per integration!' % (gmax,Ngroups)
 
        if self.verbose:
            print 'i: %d-%d' % (intmin,intmax)
            print 'g: %d-%d' % (gmin,gmax)
            print 'y: %d-%d' % (ymin,ymax)
            print 'x: %d-%d' % (xmin,xmax)

        #Ngroups_per_int_new = gmax-gmin
        Ngroups_new = gmax-gmin
        Nint_new = intmax-intmin

        #print 'Ngroups_per_int_new = %s' % Ngroups_per_int_new
        print 'Nint_new = %s' % Nint_new
        print 'Ngroups_new = %s' % Ngroups_new

        if self.ssbformat:
            for ext in xrange(1,len(self.hdulist)):
                extname = self.hdulist[ext].header['EXTNAME']
                if extname in ['SCI','GROUPDQ','ERR']:
                    print 'Extension %d:%s' % (ext,self.hdulist[ext].header['EXTNAME'])
                    self.hdulist[ext].data = self.hdulist[ext].data[intmin:intmax,gmin:gmax,ymin:ymax,xmin:xmax]
                elif extname=='PIXELDQ':
                    print 'pixel DQ Extension %d:%s' % (ext,self.hdulist[ext].header['EXTNAME'])
                    self.hdulist[ext].data = self.hdulist[ext].data[ymin:ymax,xmin:xmax]
                else:
                    print 'skipping extension %d:%s' % (ext,self.hdulist[ext].header['EXTNAME'])
            self.outdata=None

            self.hdulist[0].header['NGROUPS']=Ngroups_new
            self.hdulist[0].header['NINTS']=Nint_new
           
            if xmax-xmin<2048 or ymax-ymin<2048:
                self.hdulist[0].header['SUBARRAY']= 'SUB128' #Is this correct???
                self.hdulist[0].header['SUBSTRT1']= xmin+1 
                self.hdulist[0].header['SUBSTRT2']= ymin+1
                self.hdulist[0].header['SUBSIZE1']= xmax-xmin
                self.hdulist[0].header['SUBSIZE2']= ymax-ymin
     
        else:
            self.outdata = np.zeros(((intmax-intmin)*(gmax-gmin),ymax-ymin,xmax-xmin),dtype=np.uint16)
            print 'outdata shape:',self.outdata.shape
            iout=0
            for i in xrange(intmin,intmax):
                print 'integration %d' %i
                print 'gs:',i*Ngroups_new,(i+1)*Ngroups_new,i*Ngroups+gmin,i*Ngroups+gmax
                self.outdata[iout*Ngroups_new:(iout+1)*Ngroups_new,ymin:ymax,xmin:xmax] = data[i*Ngroups+gmin:i*Ngroups+gmax,ymin:ymax,xmin:xmax]
                iout+=1

            self.outhdr['NGROUP']=Ngroups_new
            self.outhdr['NINT']=Nint_new

            #SUBARRAY= 'FULL    '           / Subarray used
            #SUBSTRT1=                    1 / Starting pixel in axis 1 direction             
            #SUBSTRT2=                    1 / Starting pixel in axis 2 direction             
            #SUBSIZE1=                 2048 / Number of pixels in axis 1 direction           
            #SUBSIZE2=                 2048 / Number of pixels in axis 2 direction           
            if xmax-xmin<2048 or ymax-ymin<2048:
                self.outhdr['SUBARRAY']= 'SUB128' #Is this correct???
                self.outhdr['SUBSTRT1']= xmin+1 
                self.outhdr['SUBSTRT2']= ymin+1
                self.outhdr['SUBSIZE1']= xmax-xmin
                self.outhdr['SUBSIZE2']= ymax-ymin

        return(0)


    def old(self):    

        if gmin==None: gmin=0
        if xmin==None: xmin=0
        if ymin==None: ymin=0
        if intmin==None: intmin=0

        if xmax==None: xmax = self.data.shape[self.xindex]
        if ymax==None: ymax = self.data.shape[self.yindex]
        if gmax==None: gmax = self.hdr['NGROUPS']
        if intmax==None: intmax = self.hdr['NINT']

        if len(self.data.shape)==3:
            print '3-dim: fitswrite data'
            self.ssbformat = False
            self.xindex = 2
            self.yindex = 1
            self.gindex = 0
            self.intindex = None
        elif len(self.data.shape)==4:
            print '4-dim: SSB format'
            self.ssbformat = True
            self.xindex = 3
            self.yindex = 2
            self.gindex = 1
            self.intindex = 0
        else:
            print 'SHAPE:',self.data.shape
            raise RuntimeError,'wrong shape!!!'
            

        if self.verbose:
            if self.ssbformat:
                print '(i,g,x,y)=',self.data.shape
            else:
                print '(g,x,y)=',self.data.shape


        if gmin==None: gmin=0
        if xmin==None: xmin=0
        if ymin==None: ymin=0
        if intmin==None: intmin=0

        if xmax==None: xmax = self.data.shape[self.xindex]
        if ymax==None: ymax = self.data.shape[self.yindex]
        if gmax==None: gmax = self.hdr['NGROUPS']
        if intmax==None: intmax = self.hdr['NINT']
            

        if self.verbose:
            print 'i: %d-%d' % (intmin,intmax)
            print 'g: %d-%d' % (gmin,gmax)
            print 'y: %d-%d' % (ymin,ymax)
            print 'x: %d-%d' % (xmin,xmax)

        Ngroups_new = gmax-gmin
        Nint_new = intmax-intmin

        if self.ssbformat:
            self.outdata = self.data[intmin:intmax,gmin:gmax,ymin:ymax,xmin:xmax]
        else:
            self.outdata = np.zeros(((intmax-intmin)*(gmax-gmin),ymax-ymin,xmax-xmin),dtype=np.uint16)
            Ngroups = self.hdr['NGROUPS']
            for i in xrange(intmin,intmax):
                self.outdata[i*Ngroups_new:(i+1)*Ngroups_new,ymin:ymax,xmin:xmax] = self.outdata[i*Ngroups+gmin:i*Ngroups+gmax,ymin:ymax,xmin:xmax]

        self.hdr['NGROUP']=Ngroups_new
        self.hdr['NINT']=Nint_new

        #SUBARRAY= 'FULL    '           / Subarray used
        #SUBSTRT1=                    1 / Starting pixel in axis 1 direction             
        #SUBSTRT2=                    1 / Starting pixel in axis 2 direction             
        #SUBSIZE1=                 2048 / Number of pixels in axis 1 direction           
        #SUBSIZE2=                 2048 / Number of pixels in axis 2 direction           
        if xmax-xmin<2048 or ymax-ymin<2048:
            self.hdr['SUBARRAY']= 'SUBARRAY' #Is this correct???
            self.hdr['SUBSTRT1']= xmin+1 
            self.hdr['SUBSTRT2']= ymin+1
            self.hdr['SUBSIZE1']= xmax-xmin
            self.hdr['SUBSIZE2']= ymax-ymin


    def saveimage(self, outfilename):
        if not self.ssbformat:

            if self.outdata is None:
                raise RuntimeError,'No output data to save...'

            if outfilename.lower() == 'auto':
                outfilename = re.sub('fits$','sliced.fits',self.filename)
                if outfilename == self.filename:
                    raise RuntimeError,'Could not determine automatic filename for %s' % (self.filename)
            if self.verbose:
                print 'Saving ',outfilename
            rmfile(outfilename)
            pyfits.writeto(outfilename,self.outdata,self.outhdr)
        else:
            rmfile(outfilename)
            self.hdulist.writeto(outfilename)
            print('Sliced file saved to {}'.format(outfilename))
            #raise RuntimeError,'Not yet implemented!!'

if __name__=='__main__':

    usagestring='USAGE: sliceim.py infile'

    sliceim=sliceimclass()
    parser = sliceim.add_options(usage=usagestring)
    options,  args = parser.parse_args()
    (infile,)=args

    sliceim.verbose=options.verbose

    sliceim.loadimage(infile)

    sliceim.sliceim(intmin=options.intmin,
                    intmax=options.intmax,
                    gmin=options.gmin,
                    gmax=options.gmax,
                    xmin=options.xmin,
                    xmax=options.xmax,
                    ymin=options.ymin,
                    ymax=options.ymax)

    if options.outfilename == 'auto':
        options.outfilename = infile[0:-5] + '_sliced.fits'
    sliceim.saveimage(options.outfilename)

