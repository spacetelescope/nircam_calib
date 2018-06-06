#! /usr/bin/env python


## ---------------
## Modification of Bryan's test_superbias.py to include DQ check
## DQ check is done with Matt's MESA tool scripts (modified for this script)
## ---------------


'''
Test the superbias subtraction step of the pipeline
'''

from jwst.datamodels import RampModel,dqflags,SuperBiasModel
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
import argparse
import numpy as np
import extract_subarray as ex

class SuperbiasTest():
    def __init__(self):
        self.infile = None
        self.sbfile = None
        self.outfile = None

    def run_superbias_step(self):
        #run the superbias step (whether dq_init and saturation
        #flagging have been run is irrelelvant)
        m = SuperBiasStep.call(self.infile,override_superbias=self.sbfile)
        if self.outfile is not None:
            m.save(self.outfile)
        return m


    def pixeldq_propagation(self,pipeout, check):
        #check that all DQ flags are propogated from ref file to output PIXDQ array
        input_dq = np.zeros_like(pipeout.pixeldq)
        result = np.all(self.bitwise_propagate(check, input_dq) == pipeout.pixeldq)
        return result

    def bitwise_propagate(self,refhdu, pixeldq):
        # find pixels with bit sets and propogate them into PIXDQ extension
        dq_dict = {
        'DO_NOT_USE' : 0,
        'SATURATED' : 1,
        'JUMP_DET' : 2,
        'DROPOUT' : 3,
        'RESERVED' : 4,
        'RESERVED' : 5,
        'RESERVED' : 6,
        'RESERVED' : 7,
        'UNRELIABLE_ERROR' : 8,
        'NON_SCIENCE' : 9,
        'DEAD' : 10,
        'HOT' : 11,
        'WARM' : 12,
        'LOW_QE' : 13,
        'RC' : 14,
        'TELEGRAPH' : 15,
        'NONLINEAR' : 16,
        'BAD_REF_PIXEL' : 17,
        'NO_FLAT_FIELD' : 18,
        'NO_GAIN_VALUE' : 19,
        'NO_LIN_CORR' : 20,
        'NO_SAT_CHECK' : 21,
        'UNRELIABLE_BIAS' : 22,
        'UNRELIABLE_DARK' : 23,
        'UNRELIABLE_SLOPE' : 24,
        'UNRELIABLE_FLAT' : 25,
        'OPEN' : 26,
        'ADJ_OPEN' : 27,
        'UNRELIABLE_RESET' : 28,
        'MSA_FAILED_OPEN' : 29,
        'OTHER_BAD_PIXEL' : 30,
        }

        for row in refhdu.dq_def:
            try:
                # find which pixels have the bit set
                flagged = (np.bitwise_and(1, np.right_shift(refhdu.dq.astype(np.uint32), row['BIT'])))
                # shift them to the correct bit for PIXELDQ
                flagged = np.left_shift(flagged, dq_dict[row['NAME']])
                # propagate into the PIXELDQ extension
                pixeldq = np.bitwise_or(pixeldq, flagged)
            except KeyError:
                print("No DQ mnemonic "+row['NAME'])
            return pixeldq


    def compare(self):
        #run the pipeline on the input ramp
        pipeout = self.run_superbias_step()

        #Manual superbias subtraction
        sb = SuperBiasModel(self.sbfile)
        ramp = RampModel(self.infile)
        superbias = sb.data

        #extract the appropriate subarray from the reference
        #file if necessary
        if ramp.data.shape[-2] != 2048:
            xs,xe,ys,ye = ex.get_coords_rampmodel(ramp)
            superbias = sb.data[ys:ye,xs:xe]

        #subtract superbias
        ramp.data -= superbias

        #Compare pipeline output with manual output
        self.test_superbias_sub(pipeout.data,ramp.data)
        self.test_pixeldq_propagation(pipeout,sb)
        print("Superbias testing complete")


    def test_superbias_sub(self,pipe,check):
        # function to compare pipeline output with manual output
        np.testing.assert_allclose(pipe,check,atol=1e-7)

    def test_pixeldq_propagation(self,pipe,check):
        # function to make sure DQ flags get propogated
        assert self.pixeldq_propagation(pipe, check)


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Test saturation flagging')
        parser.add_argument("infile",help='Input ramp to test')
        parser.add_argument("sbfile",help='Superbias reference file to use')
        parser.add_argument("--outfile",help='Output file from pipeline',default=None)
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: test_superbias.py myramp.fits sbreffile.fits'

    sb = SuperbiasTest()
    parser = sb.add_options(usage = usagestring)
    args = parser.parse_args(namespace=sb)
    sb.compare()
