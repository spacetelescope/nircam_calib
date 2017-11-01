#! /usr/bin/env python

'''
Create gain files from already existing arrays
'''

import os
from glob import glob
from astropy.io import fits
import gain_reffile


files = glob('/ifs/jwst/wit/witserv/data4/nrc/hilbert/gain/cv3/armins_method/finalgain/NRC*dsubij.fits')
files = ['/ifs/jwst/wit/witserv/data4/nrc/hilbert/gain/cv3/armins_method/finalgain/NRCN815A-LIN-6015180318_5_481_SE_2016-01-15T18h44m39_uncal_dq_init_saturation_refpixg0.doublediff.gain.dsubij.fits']


detectors = {'481':'NRCA1','482':'NRCA2','483':'NRCA3',
             '484':'NRCA4','485':'NRCA5','486':'NRCB1',
             '487':'NRCB2','488':'NRCB3','489':'NRCB4',
             '490':'NRCB5'}

author = 'Bryan Hilbert'
useafter = '2015-10-01T00:00:00'
pedigree = 'GROUND'
descrip = 'Gain file created with CV3 data and PTC method'  
history = "This file was created using Armin Rest's new method for gain calculations. (REFERENCE DOCUMENT) Inputs were a pair of CV3 flat field integrations. Gain arrays were saved as reference files using gain_reffile.py."


for file in files:

    filename = os.path.split(file)[1]
    detnum = filename[26:29]
    detector = detectors[detnum]
    outfile = 'NIRCam_CV3_gain_{}.fits'.format(detector)
    
    with fits.open(file) as h:
        data = h[1].data

    gain = gain_reffile.GainFile()
    gain.detector = detector
    gain.outfile = outfile
    gain.author = author
    gain.descrip = descrip
    gain.useafter = useafter
    gain.pedigree = pedigree
    gain.history = history
    gain.save(data)
    
