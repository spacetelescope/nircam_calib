#!/usr/bin/env python

'''reference pixel subtraction using the G0 method.
Subtract 0th read from all reads, calculate and subtract
bias drifts, add 0th read back in.

To speed up the writing of this script, let's just call SSB's 
refpix step after subtracting the 0th read.
'''

#import numpy as np
#from astropy.io import fits,ascii
import argparse,sys
import os,copy
from jwst.datamodels import RampModel
from jwst.refpix import RefPixStep
#import sigmacut

class refpix_g0:
    def __init__(self):
        self.verbose = False

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument("infile",help="Name of the fits file to re-group.")
        parser.add_argument("--outfile",help="Name of file to output results into.",default=None)
        return parser

    def run(self):
        #check the proposed output name. If it exists, remove it.
        if self.outfile == None:
            dot = self.infile.rfind('.')
            self.outfile = self.infile[0:dot] + '_refpixg0.fits'

        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

        #read in data
        ramp = RampModel(self.infile)
        data = ramp.data

        #get data shape
        nint,ngroup,ny,nx = data.shape

        #make a copy of the 0th read
        zero_read = copy.deepcopy(ramp.data[:,0,:,:])

        #subtract the zeroth read from all subsequent reads
        for integration in range(nint):
            data[integration,:,:,:] -= zero_read[integration,:,:]
        ramp.data = data

        #run the SSB pipeline's refpix step
        ramp = RefPixStep.call(ramp,use_side_ref_pixels=True,odd_even_columns=True,odd_even_rows=False,config_file='refpix.cfg')

        #now add the original 0th read back in
        data = ramp.data
        for integration in range(nint):
            data[integration,:,:,:] += zero_read[integration,:,:]
            #dd = data[0,0,:,:] - zero_read[0,:,:]
        ramp.data = data

        #save the result
        ramp.save(self.outfile)


if __name__ == '__main__':
    usagestring = 'USAGE: refpixcorr_g0.py filename.fits'

    refsub = refpix_g0()
    parser = refsub.add_options(usage=usagestring)
    args = parser.parse_args(namespace=refsub)

    refsub.run()
