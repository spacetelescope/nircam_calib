#!/usr/bin/env python

'''
Run the dq test from within a python script
'''

import test_dq_init as test


file = '../exposures/nrca1_47Tuc_subpix_dither1_newpos_uncal.fits'
maskfile = 'NRCA1_17004_BPM_ISIMCV3_2016-01-21_ssbspmask_DMSorient_modified.fits'
outfile = 'nrca1_47Tuc_subpix_dither1_newpos_dq_init.fits'

dq = test.DQTest()
dq.infile = file
dq.maskfile = maskfile
dq.outfile = outfile
dq.compare()
