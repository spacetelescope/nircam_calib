#! /usr/bin/env python

'''
Save a readnoise reference file given the
data arrays and metadata
'''

import os
from jwst.datamodels import ReadnoiseModel

class GainFile:
    def __init__(self):
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
        
    def save(self,array,error):
        # Save using datamodel
        mod = ReadnoiseModel()
        mod.data = array
        mod.meta.telescope = 'JWST'
        mod.meta.author = self.author
        mod.meta.description = self.descrip
        mod.meta.useafter = self.useafter
        mod.meta.instrument.name = 'NIRCAM'
        mod.meta.instrument.detector = self.detector
        mod.meta.pedigree = self.pedigree
        mod.meta.reftype = 'GAIN'
        mod.meta.subarray.fastaxis = self.fastaxis
        mod.meta.subarray.slowaxis = self.slowaxis
        mod.meta.subarray.name = 'GENERIC'
        mod.meta.subarray.xsize = 2048
        mod.meta.subarray.ysize = 2048
        mod.meta.subarray.xstart = 1
        mod.meta.subarray.ystart = 1
        mod.history.append(self.history)
        mod.save(os.path.join(self.outdir,self.outfile))

