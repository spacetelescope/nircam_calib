#! /usr/bin/env python

'''
Test the jump detection step of the pipeline
'''

from jwst.datamodels import RampModel
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.refpix import RefPixStep
from jwst.linearity import LinearityStep
from jwst.dark_current import DarkCurrentStep
from jwst.jump import JumpStep
from astropy.io import ascii
from astropy.table import Table
from astropy.io import fits
from astrotable import astrotableclass
from scipy import stats
import argparse
import numpy as np
import yaml
import sys

class JumpTest():
    """This class runs and tests the output of the jump detection step.

    It does several things:
        - runs jump detection and (optionally) the steps prior
        - if simulated data, loads the input CR list for comparison
        - looks for detections and notes the CR energy and pixel spread
        - counts how many input CRs were flagged
        - counts total CR flags in each group after jump step
        - returns masks showing which pixels were flagged
        - returns tables showing counts of flags
        - runs pytest to make sure > 80 percent of input CRs were detected

    Attributes:
        run_jump_step:           Runs the pipeline.
        count_jumps:             Counts # flags for different cases.
        test_totals:             Pytest to determine how well pipeline did.
        (more tests to be added)

    Usage (in JWST pipeline environment):
        $ python test_jump.py uncal.fits --crlist="cosmicrays.list" --paramfile="file.yaml" --threshold=5 --run_steps'

    """



    def __init__(self):
        self.infile = None
        self.run_steps = False
        self.crlist = None
        self.threshold = None
        self.paramfile = None



    def run_jump_step(self,infile,threshold,run_steps):
        '''Function to run the jump detection step.'''

        # output file name
        out = infile[:-5]+"_jump_CRthresh"+str(threshold)+".fits"

        # if run_steps, run all steps prior to jump
        if run_steps:

            m = DQInitStep.call(infile)
            m = SaturationStep.call(m)
            m = SuperBiasStep.call(m)
            m_ref = RefPixStep.call(m,config_file='refpix.cfg')
            m_lin = LinearityStep.call(m)
            m_dark = DarkCurrentStep.call(m)

            # if threshold is given, use that rejection threshold
            if threshold is not None:
                m = JumpStep.call(m,output_file=out,rejection_threshold=threshold)
            else:
                m = JumpStep.call(m,output_file=out)

        # else, run only jump_step
        else:
            if threshold is not None:
                m = JumpStep.call(infile,output_file=out,rejection_threshold=threshold)
            else:
                m = JumpStep.call(infile,output_file=out)

        return m



    def count_jumps(self):
        '''Main function to count number of flags after pipeline.'''

        # run the jump step
        pipeline = self.run_jump_step(self.infile,self.threshold,self.run_steps)

        # get the detector to find the simulated CR input file
        detector = pipeline.meta.instrument.detector[3:]

        if detector == "NRCALONG":
            detector = "A5"
        if detector == "NRCBLONG":
            detector = "B5"

        # get the dq values to check what was flagged
        pixdq_arr = pipeline.pixeldq
        groupdq_arr = pipeline.groupdq
        groupdq_jump = groupdq_arr[0,:,:,:]

        # now compare input CR locations to jump step output CR flags
        if self.crlist is not None:

            # if yaml file used to generate simulated data is given, load it
            # get CR library input file
            if self.paramfile is not None:
                with open(self.paramfile,'r') as y:
                    self.params = yaml.load(y)
                crlibrary = self.params['cosmicRay']['library']
                crpath = self.params['cosmicRay']['path']
                crsuffix = self.params['cosmicRay']['suffix']
            else:
                print('ERROR, need parameter file!')
                sys.exit(0)

            # load simulated cosmicrays.list file
            crx, cry, groups, frame, cr_index, cr_frame, max_cr_sig = np.loadtxt(
                                              self.crlist,unpack=True,skiprows=3)

            # make sure data types are integers for slicing later
            crx = crx.astype(int)
            cry = cry.astype(int)
            cr_index = cr_index.astype(int)
            groups = groups.astype(int)
            cr_frame = cr_frame.astype(int)

            # set up arrays to hold results
            res = np.empty(len(crx))
            multipix = np.empty(len(crx))

            # loop over data
            notfound = np.zeros((2048,2048),dtype='bool')
            notfound_energy = np.empty(len(crx))
            found = np.zeros((2048,2048),dtype='bool')
            for i in np.arange(0,len(crx)):

                # if CR hit was flagged in groupdq, mark True
                if groupdq_arr[0,groups[i],cry[i],crx[i]] == 4:
                    found[cry[i],crx[i]] = True
                    res[i] = True

                # else, mark false
                elif groupdq_arr[0,groups[i],cry[i],crx[i]] is not 4:
                    notfound[cry[i],crx[i]] = True
                    notfound_energy[i] = max_cr_sig[i].astype(float)
                    res[i] = False

                # get the input CR file from the library
                crinfile = fits.getdata(crpath+"CRs_MCD1.7_"+crlibrary+"_0"+str(cr_index[i])+"_"+crsuffix+".fits",1)

                # look for CRs that impact more than 2 pixels
                if np.shape(np.where(crinfile[cr_frame[i],:,:] > 1))[1] >= 2:
                    multipix[i] = np.shape(np.where(crinfile[cr_frame[i],:,:] > 1))[1]
                if np.shape(np.where(crinfile[cr_frame[i],:,:] > 1))[1] < 2:
                    multipix[i] = 1


            # get total number of CRs input, number found, number missed
            totalcrs = np.float(len(crx))
            totalfound = np.float(len(crx[res == 1]))
            totalnotfound = np.float(len(crx[res == 0]))
            print('Total input CRS:',totalcrs)
            print('Total input found:',totalfound)
            print('Total input not found:',totalnotfound)
            print('Percent input found:',(totalfound/totalcrs)*100)

            # save information on CR flags based on group, number of pixels affected,
            # and their energies
            buildtable=astrotableclass()
            buildtable.t['x_Pixel'] =crx.astype(int)
            buildtable.t['y_Pixel']=cry.astype(int)
            buildtable.t['Group']=groups.astype(int)
            buildtable.t['Energy']=max_cr_sig.astype(float)
            buildtable.t['Spread']=multipix.astype(int)
            buildtable.t['Found']=res.astype(bool)
            outfilename = self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_energies.dat"
            ascii.write(buildtable.t,outfilename,format='fixed_width_two_line',overwrite=True)

            # save counts in an ASCII table
            buildtable=astrotableclass()
            buildtable.t['n_input'] =np.asfarray([totalcrs])
            buildtable.t['n_found']=np.asfarray([totalfound])
            buildtable.t['n_missed']=np.asfarray([totalnotfound])
            buildtable.t['perc_found']=np.asfarray([np.round(np.divide(totalfound,totalcrs)*100,decimals=2)])
            outfilename = self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_stats.dat"
            ascii.write(buildtable.t,outfilename,format='fixed_width_two_line',overwrite=True)

            # write out energies for CRs that were missed
            buildtable=astrotableclass()
            buildtable.t['missed_CR_energies'] = notfound_energy
            outfilename = self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_notfound_energies.dat"
            ascii.write(buildtable.t,outfilename,format='fixed_width_two_line',overwrite=True)

            # write out mask file showing which CRs were found
            hduf = fits.PrimaryHDU()
            hduf.data = np.asfarray(found)
            hduf.header['DATA'] = 'Input CRs found'
            hduf.writeto(self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_found.fits",overwrite=True)

            # write out mask file showing which CRs were not found
            hdunf = fits.PrimaryHDU()
            hdunf.data = np.asfarray(notfound)
            hdunf.header['DATA'] = 'Input CRs not found'
            hdunf.writeto(self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_notfound.fits",overwrite=True)


        # set up counter and arrays to check number of pixels flagged based
        # on groups, even if no input CR list is supplied
        counter = 0.
        grpcount = []
        slopes = []
        flagged = np.zeros((len(groupdq_jump),2048,2048),dtype='bool')

        # loop over groups and pixels
        for inds in np.arange(0,len(groupdq_jump)):
            print('Total flagged CRs in group ',inds,':',np.count_nonzero(groupdq_jump[inds,:,:] == 4))
            for x in np.arange(0,2048):
                for y in np.arange(0,2048):

                    # if GROUPDQ is flagged for a CR hit, flag pixel in mask
                    if groupdq_jump[inds,y,x] == 4:
                        counter += 1
                        flagged[inds,y,x] = True
                        slope, intercept, r_val, p_val, std_err = stats.linregress(np.arange(0,2), pipeline.data[0,inds-1:inds+1,y,x])
                        slopes.append(slope)

            grpcount.append(np.count_nonzero(groupdq_jump[inds,:,:] == 4))

        # save it as an ASCII table
        buildtable=astrotableclass()
        buildtable.t['group'] = np.arange(0,len(groupdq_jump)).astype(int)
        buildtable.t['n_CRflags'] = grpcount
        outfilename = self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_groupstats.dat"
        ascii.write(buildtable.t,outfilename,format='fixed_width_two_line',overwrite=True)

        # save out slope information
        buildtable=astrotableclass()
        buildtable.t['jump_slopes'] = slopes
        outfilename = self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_slopes.dat"
        ascii.write(buildtable.t,outfilename,format='fixed_width_two_line',overwrite=True)

        # write out mask showing CR flags in each group
        hdufl = fits.PrimaryHDU()
        hdufl.data = np.asfarray(flagged)
        hdufl.header['DATA'] = 'Overall CRs flagged'
        hdufl.writeto(self.infile[:-5]+"_jump_CRthresh"+str(self.threshold)+"_group_CRflags.fits",overwrite=True)

        # begin actual pytests
        if self.crlist is not None:
            self.test_jump_totals(totalcrs,totalfound)

        print("\nJump detection validation test complete.")



    # beginning of tests to validation jump step
    def test_jump_totals(self,totalcrs,totalfound):
        '''Test to see how well the pipeline did.'''

        # compare difference between pipeline output and input CRs to
        # atol + rtol * abs(input).
        print('\nTesting to see if pipeline flags more than 80% of input CRs...')
        np.testing.assert_allclose(totalfound,totalcrs,rtol=0.2,atol=0.0)





    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Test saturation flagging')
        parser.add_argument("infile",help='Input ramp to test')
        parser.add_argument("--crlist",help='Mask file to use for testing (fits)',default=None)
        parser.add_argument("--threshold",help='Mask file to use for testing (fits)',default=None)
        parser.add_argument("--paramfile",help='Mask file to use for testing (fits)',default=None)
        parser.add_argument("--run_steps",help='If True, run DQ_init before saturation step.',action='store_true',default=False)
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: python jump_validation_tests.py uncal.fits --crlist="cosmicrays.list" --paramfile="sim_params.yaml" --threshold=5 --run_steps'

    jump = JumpTest()
    parser = jump.add_options(usage = usagestring)
    args = parser.parse_args(namespace=jump)
    jump.count_jumps()
