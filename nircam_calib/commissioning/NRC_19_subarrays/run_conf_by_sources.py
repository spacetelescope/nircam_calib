#! /usrbin/env python

from nircam_calib.commissioning.NRC_19_subarrays import confirm_subarray_location_via_sources as c

# B5 - seems to be workig well
#full = 'Pipeline_Level1/jw01068004001_01101_00001_nrcb5_rate.fits'
#subs = ['Pipeline_Level1/jw01068001001_01101_00001_nrcb5_rate.fits',
#        'Pipeline_Level1/jw01068002001_01101_00001_nrcb5_rate.fits',
#        'Pipeline_Level1/jw01068003001_01101_00001_nrcb5_rate.fits']
#c.run_using_fits(full, subs)

# B1
#full = 'Pipeline_Level2/jw01068004001_01101_00002_nrcb1_cal.fits'
#subs = ['Pipeline_Level2/jw01068001001_01101_00001_nrcb1_cal.fits',
#        'Pipeline_Level2/jw01068002001_01101_00001_nrcb1_cal.fits',
#        'Pipeline_Level2/jw01068003001_01101_00001_nrcb1_cal.fits']
#c.run_using_fits(full, subs)

# B1 - use datamodels
full = 'Pipeline_Level2/jw01068004001_01101_00002_nrcb1_cal.fits'
subs = ['Pipeline_Level2/jw01068001001_01101_00001_nrcb1_cal.fits',
        'Pipeline_Level2/jw01068002001_01101_00001_nrcb1_cal.fits',
        'Pipeline_Level2/jw01068003001_01101_00001_nrcb1_cal.fits']
c.run(full, subs)
