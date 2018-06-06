#! /usr/bin/env python

'''
Modify the input bad pixel table such that every possible type of bad
pixel flag is present, so we can then check that all are propagated
correctly into an observation
'''

from jwst import datamodels
import argparse
from astropy.io import fits

class DQ():
    def __init__(self):
        self.infile = None
        self.outfile = None

        
    def modify(self):
        #read in existing bad pixel mask
        #mask = datamodels.open(self.infile)
        h = fits.open(self.infile)
        
        #bad pixel flags can have values of any
        #power of 2 up through 1073741824 (2**30)
        for i in range(31):
            h[1].data[5,i+400] = 2**(i)

        #now add in some combinations of bad pix
        for i in range(0,30,2):
            h[1].data[5,i+431] = 2**(i) + 2**(i+1)

        #save the result
        if self.outfile is None:
            self.outfile = self.infile[0:-5] + '_modified.fits'
        h.writeto(self.outfile,overwrite=True)


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Simulate JWST ramp')
        parser.add_argument("infile",help='Input bad pixel mask filename')
        parser.add_argument("--outfile",help='Output name for modified bad pixel mask.',default=None)
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: modify_badpix_table.py bpm.fits'

    dq = DQ()
    parser = dq.add_options(usage = usagestring)
    args = parser.parse_args(namespace=dq)
    dq.modify()


