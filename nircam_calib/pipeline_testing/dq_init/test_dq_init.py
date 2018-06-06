#! /usr/bin/env python

'''
Test the dq_init step of the pipeline
'''

from jwst.datamodels import RampModel
from jwst.dq_init import DQInitStep
import argparse
import numpy as np
from astropy.io import fits

class DQTest():
    def __init__(self):
        self.infile = None
        self.maskfile = None
        self.outfile = None
        
    def run_dq_step(self,file,maskfile,outfile):
        #run the dq_init pipeline step
        m = DQInitStep.call(file,override_mask=maskfile,output_file=outfile)
        return m

    def compare(self):
        #set up name of output file if the user didn't
        if self.outfile is None:
            uncal = self.infile.rfind('uncal')
            suffix = 'dq_init.fits'
            if uncal == -1:
                uncal = -5
                suffix = '_' + suffix
            self.outfile = self.infile[0:uncal] + suffix

        
        #Run the pipeline on the input file
        pipeout = self.run_dq_step(self.infile,self.maskfile,self.outfile)

        #Read in the mask file
        with fits.open(self.maskfile) as msk:
            mask = msk[1].data
        
        #Now compare the relevant pixels' values with what they should be
        self.test_dq(pipeout.pixeldq,mask)

        print("DQ test complete.")
        
        
    def test_dq(self,pipeline,truth):
        np.testing.assert_allclose(pipeline,truth,atol=0)


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Test DQ mask propagation')
        parser.add_argument("infile",help='Input ramp to test')
        parser.add_argument("maskfile",help='Mask file to use for testing (fits)')
        parser.add_argument("--outfile",help='Output name for ramp after applying DQ mask.',default=None)
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: test_dq_init.py myramp.fits mymask.fits'

    dq = DQTest()
    parser = dq.add_options(usage = usagestring)
    args = parser.parse_args(namespace=dq)
    dq.compare()

