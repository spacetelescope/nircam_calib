#! /usr/bin/env python

"""Verify that sources are properly registered by the Stage 3 pipeline. Do this by comparing the FWHM of sources in individual
(distortion-corrected?) images to that in the level 3 combined images

Example call of this script:

python verify_source_registration.py individual_files_test.list jw09996001001_01101_00001_nrcb5_rate.fits
"""
import argparse
from astropy.convolution import Gaussian2DKernel
from astropy.io import ascii, fits
from astropy.stats import gaussian_fwhm_to_sigma
import numpy as np
import os
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources, detect_sources, detect_threshold
from photutils.background import Background2D, MedianBackground



def process_image(data, threshold=15, numpix=5):
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                       bkg_estimator=bkg_estimator)
    data -= bkg.background  # subtract the background
    threshold = threshold * bkg.background_rms  # above the background
    #sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    #kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    #kernel.normalize()
    segm = detect_sources(data, threshold, npixels=numpix) # kernel=kernel)
    segm_deblend = deblend_sources(data, segm, npixels=numpix, # kernel=kernel,
                               nlevels=32, contrast=0.001)

    cat = SourceCatalog(data, segm_deblend)
    #>>> labels = [1, 5, 20, 50, 75, 80]
    #>>> cat_subset = cat.get_labels(labels)
    columns = ['label', 'xcentroid', 'ycentroid', 'area', 'segment_flux', 'fwhm', 'elongation']
    tbl = cat.to_table(columns=columns)

    # Flux cut to keep only the bright sources
    flux_min = 1000
    flux_max = 500000
    flux_range = np.where((tbl['segment_flux'] > flux_min) & (tbl['segment_flux'] < flux_max))[0]
    tbl = tbl[flux_range]

    tbl['xcentroid'].info.format = '.4f'  # optional format
    tbl['ycentroid'].info.format = '.4f'
    tbl['segment_flux'].info.format = '.4f'
    print(tbl)
    return tbl


def run(individual_files, combined_image_file):
    for file in individual_files:
        img = fits.getdata(file)
        srcs = process_image(img)
        outfile = f'fwhm_measurements_{os.path.basename(file)}.tab'
        ascii.write(srcs, outfile, overwrite=True)
        median_fwhm = np.median(srcs['fwhm'])
        median_elongation = np.median(srcs['elongation'])
        print(f'{file}: {median_fwhm}, {median_elongation}')

    mosaic = fits.getdata(combined_image_file)
    mosaic_src = process_image(mosaic)
    outfile = f'fwhm_measurements_{os.path.basename(combined_image_file)}.tab'
    ascii.write(mosaic_src, outfile, overwrite=True)
    median_fwhm = np.median(mosaic_src['fwhm'])
    median_elongation = np.median(mosaic_src['elongation'])
    print(f'MOSAIC {file}: {median_fwhm}, {median_elongation}')


def define_options(parser=None, usage=None, conflict_handler='resolve'):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('listfile', default=None,
                        help='File continaing a list of individual images to measure PSF properties')
    parser.add_argument('mosaic_file', default=None,
                        help='Fits file containing mosaic image to measure PSF properties')
    return parser

if __name__ == '__main__':
    parser = define_options()
    args = parser.parse_args()

    individual_files = []
    with open(args.listfile) as fobj:
        individual_files.append(fobj.read().strip('\n'))

    print(individual_files)



    run(individual_files, args.mosaic_file)

