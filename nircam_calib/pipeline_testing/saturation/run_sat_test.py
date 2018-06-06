#! /usr/bin/env python

'''
Test run of the saturation test script
'''

import test_saturation

#Include dq_init run
sat = test_saturation.SatTest()
sat.infile = '../exposures/nrca1_47Tuc_subpix_dither1_newpos_uncal.fits'
sat.satfile = 'nrca1_modified_saturation_reffile.fits'
sat.run_dq = True
sat.maskfile = '../reffiles/NRCA1_17004_BPM_ISIMCV3_2016-01-21_ssbspmask_DMSorient.fits'
sat.outfile = 'nrca1_47Tuc_subpix_dither1_newpos_dq_init_saturation.fits'
sat.compare()

#Without dq_init run
sat = test_saturation.SatTest()
sat.infile = '../exposures/nrca1_47Tuc_subpix_dither1_newpos_uncal.fits'
sat.satfile = 'nrca1_modified_saturation_reffile.fits'
sat.run_dq = False
sat.outfile = 'nrca1_47Tuc_subpix_dither1_newpos_dq_init_saturation.fits'
sat.compare()
