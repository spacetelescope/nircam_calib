#! /usr/bin/env python

"""Run nircam_calib's confirm_count_rates.py
"""
from glob import glob
from nircam_calib.commissioning.NRC_19_subarrays import confirm_count_rates as c
import os

# Using rate images
fullframe_file = 'Pipeline_Level2/jw01068004001_01101_00001_nrcb5_cal.fits'
#subarray_file = 'Pipeline_Level2/jw01068001001_01101_00001_nrcb1_cal.fits'


# Result for B5 Mirage data - subarray sources are 6-7% brighter than in full frame
#subfiles = sorted(glob('Pipeline_Level2/jw01068001*nrcb5_cal.fits'))
#print('subfiles:', subfiles)
#for subarray_file in subfiles:
#    print('Trying: {}'.format(subarray_file))
#    c.compare_rate_images(subarray_file, fullframe_file)

# Using level3 source catalogs
ff_cat_file = ''
sub_cat_file = ''
outfile = 'Source_comparison_{}_{}.tab'.format(os.path.basename(ff_cat_file), os.path.basename(sub_cat_file))
c.compare_level3_catalogs(sub_cat_file, ff_cat_file, outfile)
