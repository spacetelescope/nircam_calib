#! /usr/bin/env python

'''
Given a series of individual gain maps, create a mean gain
map and error array

If the indiviual gain files were created using gain.py, then
this script will check the headers of those files
and create a list of flats and darks used to create the gain.
These lists will be appended to the history string.

'''

import sys
import copy
import numpy as np
from astropy.stats import sigma_clip
from astropy.io import fits
from . import gain_reffile as gainref

class Map():
    def __init__(self):
        self.gainlist = []
        self.author = ''
        self.descrip = ''
        self.pedigree = ''
        self.useafter = ''
        self.history = ''
        self.save_output = True
        self.outfile = ''
        self.outdir = './'
        
        
    def check_headers(self,hdrlist):
        '''Make sure inputs have consistent header info'''
        dets = []
        for i,h in enumerate(hdrlist):
            det = h['DETECTOR']
            if i == 0:
                compdet = copy.deepcopy(det)
            if det != compdet:
                print(("WARNING: inconsistent detector values"
                       "in input file headers. Quitting."))
                sys.exit()
        

    def create_map(self):
        '''Main function'''

        # Read in individual input files
        if isinstance(self.gainlist[0],str):
            hdrs = []
            flatfiles = []
            darkfiles = []
            for i,file in enumerate(self.gainlist):
                im,err,hdr0 = self.read_data(file)
                if i == 0:
                    gainims = copy.deepcopy(im)
                    gainerrs = copy.deepcopy(err)
                else:
                    if im.shape[1:] != gainims.shape[1:]:
                        print(("WARNING: inconsistent array sizes "
                               "in input files. Quitting."))
                        sys.exit
                    gainims = np.vstack([gainims,im])
                    gainerrs = np.vstack([gainerrs,err])
                hdrs.append(hdr0)

                # Get a list of files used to create the
                # mean gain map
                flats = self.get_ffiles(hdr0)
                flatfiles += flats
                darks = self.get_dfiles(hdr0)
                if darks[0] not in [None,'']:
                    darkfiles += darks

        # If the lists of flat and dark files have entries
        # in them, add them to the HISTORY
        if len(flatfiles) > 0:
            newtext = ("Flatfield files used to create this "
                       "gain map:  ")
            for f in flatfiles:
                newtext += "{}, ".format(f)
            self.history += newtext
        if len(darkfiles) > 0:
            newtext = ("Dark current files used to create this "
                       "gain map: ")
            for f in darkfiles:
                newtext += "{}, ".format(f)
            self.history += newtext
        
        # Check headers for consistency
        self.check_headers(hdrs)

        # Calculate the sigma-clipped mean through the stack
        # astropy.stats.SigmaClip gives confusing results
        # so let's fall back to scipy's version
        yd,xd = gainims.shape[1:]
        self.gain = np.zeros((yd,xd))
        self.gainerr = np.zeros((yd,xd))

        p = sigma_clip(gainims,sigma=3,axis=0)
        for x in range(xd):
            for y in range(yd):
                ppix = p[:,y,x]
                mn = np.mean(ppix.data[~ppix.mask])
                self.gain[y,x] = mn
                gs = gainims[:,y,x]
                errs = gainerrs[:,y,x]
                perr = np.sqrt(np.sum((gs[~ppix.mask]*errs[~ppix.mask])**2)
                               / np.sum(~ppix.mask)**2)
                self.gainerr[y,x] = perr
                
        if self.save_output:
            # Save in CRDS format
            det = hdrs[0]['DETECTOR']
            if '5' in det:
                det = det.replace('5','LONG')

            self.get_axes_values(det)
            gfile = gainref.GainFile()
            gfile.detector = det
            gfile.author = self.author
            gfile.descrip = self.descrip
            gfile.useafter = self.useafter
            gfile.pedigree = self.pedigree
            gfile.history = self.history
            gfile.outdir = self.outdir
            gfile.outfile = self.outfile
            gfile.fastaxis = self.fastaxis
            gfile.slowaxis = self.slowaxis
            gfile.save(self.gain,self.gainerr)
            

    def get_axes_values(self,detector):
        '''Retrieve the appropriate fastaxis and slowaxis
        keywords for the given detector'''
        if detector in ['NRCA2','NRCA4','NRCB1','NRCB3','NRCBLONG']:
            self.fastaxis = 1
            self.slowaxis = -2
        elif detector in ['NRCA1','NRCA3','NRCB2','NRCB4','NRCALONG']:
            self.fastaxis = -1
            self.slowaxis = 2

            
    def get_dfiles(self,hdu):
        '''Return the dark current files used to
        create an individual gain map'''
        d1 = hdu['DARKFIL1']
        d2 = hdu['DARKFIL2']
        return [d1,d2]


    def get_ffiles(self,hdu):
        '''Return the flat field files used to
        create an individual gain map'''
        f1 = hdu['FLATFIL1']
        f2 = hdu['FLATFIL2']
        return [f1,f2]

    
    def read_data(self,infile):
        '''Read in an individual gain file, return
        gain map, error map, and header'''
        with fits.open(infile) as h:
            im = h['gain'].data
            err = h['gain_err'].data
            hdr0 = h[0].header
            im = np.expand_dims(im,0)
            err = np.expand_dims(err,0)
        return im, err, hdr0
