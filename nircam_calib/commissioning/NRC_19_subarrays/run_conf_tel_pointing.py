#! /usrbin/env python

from nircam_calib.commissioning.NRC_19_subarrays import confirm_telescope_pointing as c

filenames = ['Pipeline_Level2/jw01068005001_01101_00001_nrcb1_cal.fits',
             'Pipeline_Level2/jw01068005001_01101_00002_nrcb5_cal.fits',
             'Pipeline_Level2/jw01068006001_01101_00001_nrcb1_cal.fits',
             'Pipeline_Level2/jw01068006001_01101_00002_nrcb5_cal.fits']

for filename in filenames:
    c.check_pointing_target_star(filename)

for filename in filenames:
    c.check_pointing_using_2mass_catalog(filename)

