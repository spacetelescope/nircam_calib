#! /usr/bin/env python

"""Call darks_with_zeros.py to make darks for SW detectors.
"""
from glob import glob
import os

from darks_with_zeroes import create_dark

sigma_threshold = 5

descrip = 'CV3 based dark with zero for nominal pixels'
use_after = '2015-10-31'
history = ('This dark reference file was created from CV3 dark current observations. '
           'We first created a slope image for each input dark exposure, and then created '
           'a mean slope image. Pixels with signal rates more than {}-sigma above the mean '
           'are treated as hot and have dark signal values equal to their mean dark rate times the '
           'exposure time of each group. Pixels with signal rates below this threshold have '
           'their dark signals set to zero. This strategy was agreed upon by members of the '
           'NIRCam IDT. '.format(sigma_threshold))


# Development
#a1_output_dir = '/grp/jwst/wit/nircam/reference_files/darks/full_frame/for_commissioning/A1/'
#a1_output = os.path.join(a1_output_dir, 'NRCA1_dark_zeros_reffile.fits')
#a1_linearized_darks = sorted(glob(os.path.join(a1_output_dir, '*linearity.fits')))
#a1_rate_files = sorted(glob(os.path.join(a1_output_dir, '*rate.fits')))

#create_dark(a1_linearized_darks, input_dark_slope_files=a1_rate_files, sigma_threshold_for_hot=sigma_threshold,
#            descrip=descrip, use_after=use_after, history=history, output_file=a1_output,
#            output_slope_file_dir=a1_output_dir)

short_list_dets = ['A1', 'A2', 'A4']
long_list_dets = ['A3', 'B1', 'B2', 'B3', 'B4']


detectors = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4']
#detectors = ['A3', 'A4', 'B1', 'B2', 'B3', 'B4']
for det in detectors:
    output_dir = '/grp/jwst/wit/nircam/reference_files/darks/full_frame/for_commissioning/{}/'.format(det)
    output = os.path.join(output_dir, 'NRC{}_dark_zeros_reffile.fits'.format(det))

    #if det in short_list_dets:
    #    raw_darks = sorted(glob('/ifs/jwst/wit/nircam/isim_cv3_files_for_calibrations/darks/{}/*2016-01-0*level1b_uncal.fits'.format(det)))
    #elif det in long_list_dets:
    #    raw_darks1 = sorted(glob('/ifs/jwst/wit/nircam/isim_cv3_files_for_calibrations/darks/{}/*2015*level1b_uncal.fits'.format(det)))
    #    raw_darks2 = sorted(glob('/ifs/jwst/wit/nircam/isim_cv3_files_for_calibrations/darks/{}/*2016-01-0*level1b_uncal.fits'.format(det)))
    #    raw_darks = raw_darks1 + raw_darks2

    linearized_darks = sorted(glob(os.path.join(output_dir, '*linearity.fits')))
    rate_files = sorted(glob(os.path.join(output_dir, '*rate.fits')))

    #create_dark(raw_darks, sigma_threshold_for_hot=sigma_threshold,
    #        descrip=descrip, use_after=use_after, history=history, output_file=output,
    #        output_slope_file_dir=output_dir)

    create_dark(linearized_darks, input_dark_slope_files=rate_files, sigma_threshold_for_hot=sigma_threshold,
            descrip=descrip, use_after=use_after, history=history, output_file=output,
            output_slope_file_dir=output_dir)




