#! /usr/bin/env python

'''
Create persistence trap density reference file
'''

from jwst.datamodels import TrapDenisityModel
from astropy.io import fits

class TrapDen():

    def __init__(self):
        self.input = None
        self.input_file = False
        self.detector = ''
        self.outfile = ''
        self.author = ''
        self.descrip = ''
        self.useafter = ''
        self.pedigree = ''
        self.history = ''
        self.outdir = '.'
        self.fastaxis = 0
        self.slowaxis = 0
        
    def save(self):
        # If the input is a file, read in
        if self.input_file:
            with fits.open(self.input) as h:
                traps = h[1].data
        else:
            traps = self.input

        mod = TrapDensityModel(traps)

        # Metadata
        mod.meta.telescope = 'JWST'
        mod.meta.instrument.name = 'NIRCAM'
        mod.meta.instrument.detector = self.detector
        mod.meta.author = self.author
        mod.meta.description = self.descrip
        mod.meta.reftype = 'TRAPDENSITY'
        mod.meta.pedigree = self.pedigree
        mod.meta.useafter = self.useafter
        mod.meta.subarray.name = 'GENERIC'
        mod.meta.subarray.slowaxis = self.slowaxis
        mod.meta.subarray.fastaxis = self.fastaxis
        mod.meta.xsize = 2048
        mod.meta.ysize = 2048
        mod.meta.xstart = 1
        mod.meta.ystart = 1
        mod.history.append()
        mod.save(os.path.join(self.outdir,self.outfile))
