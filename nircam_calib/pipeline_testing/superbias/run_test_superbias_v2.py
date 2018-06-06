#!/usr/bin/env python

'''
Test run the superbias testing step from within a python script
'''

# import test_superbias as test
import test_superbias_v2 as test

sb = test.SuperbiasTest()
sb.infile = 'NRCNRCA1-DARK-60012216201_1_481_SE_2016-01-02T02h34m28_uncal_sliced.fits'
sb.sbfile = 'jwst_nircam_superbias_0026.fits'
sb.outfile = 'NRCNRCA1-DARK-60012216201_1_481_SE_2016-01-02T02h34m28_uncal_sliced_dq_init_saturation_superbias.fits'
sb.compare()
