#! /usr/bin/env python

'''
Test the saturation flagging step of the pipeline
'''

from jwst.datamodels import RampModel,dqflags
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
import argparse
import numpy as np
from astropy.io import fits
from copy import copy

class SatTest():
    def __init__(self):
        self.infile = None
        self.satfile = None
        self.outfile = None
        self.run_dq = False
        self.maskfile = None
        
    def run_sat_step(self,file,maskfile,satfile,outfile,run_dq=True):
        #run the saturation pipeline step

        #if the dq_init step needs to be run, do that first
        if run_dq:
            if maskfile is not None:
                m = DQInitStep.call(file,override_mask=maskfile)
            else:
                m = DQInitStep.call(file)
            m = SaturationStep.call(m,override_saturation=satfile,output_file=outfile)    

        else:
            #run the saturation step
            m = SaturationStep.call(file,override_saturation=satfile,output_file=outfile)
        return m

    
    def compare(self):
        #set up name of output file if the user didn't
        if self.outfile is None:
            dot = self.infile.rfind('.')
            suffix = '_saturation.fits'
            self.outfile = self.infile[0:dot] + suffix

        #Run the pipeline on the input file
        pipeout = self.run_sat_step(self.infile,self.maskfile,self.satfile,self.outfile,run_dq=self.run_dq)

        #Read in the saturation file
        with fits.open(self.satfile) as sat:
            satdata = sat[1].data
        
        #Now compare the relevant pixels' values with what they should be
        self.test_sat(pipeout,satdata)

        print("Saturation test complete.")


    def test_sat(self,pipeline,truth):
        #get a list of points for which no sat check should be performed
        nocheck = (pipeline.pixeldq & dqflags.pixel['NO_SAT_CHECK'] > 0)
        nocheck = np.expand_dims(nocheck,axis=0)
        
        #create manually an array of saturated and unsaturated data points
        nint,ngroup,ny,nx = pipeline.data.shape
        for integ in range(nint):
            my_saturated = pipeline.data[integ,:,:,:] >= truth
        
            nochecks = copy(nocheck)
            for i in range(ngroup-1):
                nochecks = np.vstack([nochecks,nocheck])

            #set no_sat_check pix to false
            my_saturated[nochecks] = False

            #now we need to deal with pix that saturate but
            #then have signal decrease to below the threshold
            #Need to make sure these are still flagged as saturated
            for group in range(ngroup):
                satframe = my_saturated[group,:,:]
                saty,satx = np.where(satframe == True)
                my_saturated[group:,saty,satx] = True

            #create a boolean array of the saturated data points from the pipeline
            pipe_saturated = (pipeline.groupdq[integ,:,:,:] & dqflags.pixel['SATURATED'] > 0)
        
            np.testing.assert_allclose(my_saturated,pipe_saturated)
            


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Test saturation flagging')
        parser.add_argument("infile",help='Input ramp to test')
        parser.add_argument("satfile",help='Satruation reference file to use')
        parser.add_argument("--maskfile",help='Mask file to use for testing (fits)',default=None)
        parser.add_argument("--outfile",help='Output name for ramp after applying DQ mask.',default=None)
        parser.add_argument("--run_dq",help='If True, run DQ_init before saturation step.',action='store_true',default=False)
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: test_saturation.py myramp.fits satreffile.fits'

    sat = SatTest()
    parser = sat.add_options(usage = usagestring)
    args = parser.parse_args(namespace=sat)
    sat.compare()
