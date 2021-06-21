#! /usr/bin/env python

"""Use source locations to confirm that the correct pixels are being extracted from the full detector when
subarray observations are made.

Input:
subarray slope images
full frame slope images

Method:
Detect sources and create source catalog for each exposure
Get expected offset between full frame and subarray coordinates from pysiaf
Compare subarray source locations + coordinate offsets (+ any dithers) to source locations on full frame data

Note that the sub320, sub640, and full data all have the same set of 4 dithers, so we should be able to
compare exposure 1 of each to each other, and exposure 2, etc.
sub160 has more dithers, but the dither positions of the full frame observations are a subset of those
in the sub160 observations. So we should be able to pull our the appropriate sub160 exposures and compare
those directly as well. The observations with sub64p and sub400p have no dithers, and should be comparable
to the initial exposure in the full data.

For the GRISM256 data, we would need to back out the location of the source in the direct data in order to compare,
since the grism shifts the source in y by 1.449", and we'd need to calculate location in x using the dispersion solution

========================================
* sub160 RAPID (Obs 1)
** Visit 1:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -87.081 -491.730   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.160   +0.144  -87.241 -491.586   +0.160   +0.144 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1    3 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.080   +0.224  -87.161 -491.506   +0.080   +0.224 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1    4 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.240   +0.048  -87.321 -491.682   +0.240   +0.048 SUBDITHER  SCIENCE        0      0  0.238
  1    1   1    5 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -95.274 -483.543   +8.190   +8.190 DITHER     SCIENCE        0      0 11.380
  1    1   1    6 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.350   +8.334  -95.434 -483.399   +8.350   +8.334 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1    7 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.270   +8.414  -95.354 -483.319   +8.270   +8.414 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1    8 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.430   +8.238  -95.514 -483.495   +8.430   +8.238 SUBDITHER  SCIENCE        0      0  0.238
  1    1   1    9 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -95.269 -499.923   +8.190   -8.190 DITHER     SCIENCE        0      0 16.430
  1    1   1   10 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.350   -8.046  -95.429 -499.779   +8.350   -8.046 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1   11 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.270   -7.966  -95.349 -499.699   +8.270   -7.966 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1   12 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.430   -8.142  -95.509 -499.875   +8.430   -8.142 SUBDITHER  SCIENCE        0      0  0.238
  1    1   1   13 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -78.889 -499.918   -8.190   -8.190 DITHER     SCIENCE        0      0 16.620
  1    1   1   14 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.030   -8.046  -79.049 -499.774   -8.030   -8.046 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1   15 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.110   -7.966  -78.969 -499.694   -8.110   -7.966 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1   16 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -7.950   -8.142  -79.129 -499.870   -7.950   -8.142 SUBDITHER  SCIENCE        0      0  0.238

========================================
* sub320 RAPID (Obs 2)
** Visit 2:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -87.081 -491.730   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -95.274 -483.543   +8.190   +8.190 DITHER     SCIENCE        0      0 11.582
  1    1   1    3 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -95.269 -499.923   +8.190   -8.190 DITHER     SCIENCE        0      0 16.380
  1    1   1    4 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -78.889 -499.918   -8.190   -8.190 DITHER     SCIENCE        0      0 16.380

========================================
* sub640 RAPID (Obs 3)
** Visit 3:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -87.081 -491.730   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -95.274 -483.543   +8.190   +8.190 DITHER     SCIENCE        0      0 11.582
  1    1   1    3 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -95.269 -499.923   +8.190   -8.190 DITHER     SCIENCE        0      0 16.380
  1    1   1    4 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -78.889 -499.918   -8.190   -8.190 DITHER     SCIENCE        0      0 16.380

========================================
* full (Obs 4)
** Visit 4:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -82.287 -496.206   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -90.486 -488.025   +8.190   +8.190 DITHER     SCIENCE        0      0 11.582
  1    1   1    3 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -90.468 -504.405   +8.190   -8.190 DITHER     SCIENCE        0      0 16.380
  1    1   1    4 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -74.088 -504.387   -8.190   -8.190 DITHER     SCIENCE        0      0 16.380

========================================
* sub400p (Obs 5)
** Visit 5:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  0    0   0    0 NRCB5_TAPSIMG32           2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -152.003 -440.340   +0.000   +0.000 TARGET     T_ACQ          0      0  0.000 (base)
  1    1   1    1 NRCB1_SUB400P             2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -146.051 -432.206   +0.000   +0.000 TARGET     SCIENCE        0      0 10.079

========================================
* sub64p (Obs 6)
** Visit 6:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  0    0   0    0 NRCB5_TAPSIMG32           2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -152.003 -440.340   +0.000   +0.000 TARGET     T_ACQ          0      0  0.000 (base)
  1    1   1    1 NRCB1_SUB64P              2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -151.093 -427.403   +0.000   +0.000 TARGET     SCIENCE        0      0 12.969

========================================
* subgrism256 (Obs 7)
** Visit 7:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  0    0   0    0 NRCA5_TAGRISMTS32         2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000  +73.262 -551.572   +0.000   +0.000 TARGET     T_ACQ          0      0  0.000 (base)
  1    1   1    1 NRCA5_GRISM256_F322W2     2 LMC-ASTR  +80.48258  -69.49396   +0.000   +1.449   +0.000   +0.000  +50.773 -555.028   +0.000   +1.449 TARGET     SCIENCE        0      0 22.752


"""
from astropy.io import ascii, fits
import numpy as np
import os
#import pysiaf
#import yaml

#from mirage.utils import siaf_interface
from jwst import datamodels

from nircam_calib.commissioning.utils.astrometry import RADec_To_XY, XY_To_RADec
from nircam_calib.commissioning.utils.photometry import find_sources, fwhm, get_fwhm


def run_using_fits(full_frame_file, subarray_files, output_dir='./'):
    """MAIN FUNCTION

    THIS VERSION DOES NOT USE DATAMODELS, AND RELIES ON SIAF FOR COORDINATE
    TRANSFORMS. IT IS MUCH EASIER TO USE THE run() FUNCTION BELOW.

    Paramters
    ---------
    full_frame_file : str

    subarray_files : list
    """
    # Determine which full frame image to compare with which subarray images
    # Detect sources in full-frame image
    # Detect sources in subarray image
    # Get subarray location from pysiaf
    # Add pysiaf result to subarray catalog
    # Compare full frame catalog to shifted subarray catalog to see how well they match.
    # For subarrays not centered on the detector (sub64p/sub400p), this will fail because the reference locations are different!!!
        # In that case, you'll have to use knowledge of the refrence locations to predict source locations
    # Compare
    siaf_instance = pysiaf.Siaf('nircam')

    # Read in full frame data and locate sources
    with fits.open(full_frame_file) as ffhdu:
        full_frame = ffhdu[1].data
        ff_header = ffhdu[0].header
        ff_header1 = ffhdu[1].header

    ff_det = ff_header['DETECTOR']
    ff_ap_name = ff_header['SUBARRAY']
    if 'LONG' in ff_det:
        ff_det = ff_det.replace('LONG', '5')
    ff_aperture = '{}_{}'.format(ff_det, ff_ap_name)
    ff_filter = ff_header['FILTER']

    # Calculate the FWHM in pixels to input to the source finder
    ff_fwhm = get_fwhm(ff_filter)

    ffbasename = os.path.basename(full_frame_file)
    full_frame_sources = find_sources(full_frame, threshold=500, fwhm=ff_fwhm, plot_name='{}_full_frame_source_map.png'.format(ffbasename))
    ascii.write(full_frame_sources, '{}_ff_sources.txt'.format(ffbasename), overwrite=True)

    # Read in full frame file's WCS
    ff_filebase = os.path.join('yaml_files', os.path.basename(full_frame_file))
    ff_local_roll, ff_attitude_matrix, ff_ffsize, ff_subarray_bounds = siaf_interface.get_siaf_information(siaf_instance, ff_aperture,
                                                                                                           ff_header1['RA_REF'],
                                                                                                           ff_header1['DEC_REF'],
                                                                                                           ff_header1['PA_V3'], v2_arcsec=None,
                                                                                                           v3_arcsec=None, verbose=False)

    # Read in subarray data and locate sources
    for subfile in subarray_files:
        subarray, start_coords, end_coords, det, ap_name, filt, ra_ref, dec_ref, pa_v3 = get_data(subfile)
        if 'LONG' in det:
            det = det.replace('LONG', '5')
        aperture = '{}_{}'.format(det, ap_name)

        # Calculate the FWHM in pixels to input to the source finder
        sub_fwhm = get_fwhm(filt)
        subarray_sources = find_sources(subarray, threshold=500, fwhm=sub_fwhm, plot_name='{}_source_map.png'.format(os.path.basename(subfile)))

        filebase = os.path.join('yaml_files', os.path.basename(subfile))
        local_roll, attitude_matrix, ffsize, subarray_bounds = siaf_interface.get_siaf_information(siaf_instance, aperture,
                                                                                                   ra_ref, dec_ref, pa_v3, v2_arcsec=None,
                                                                                                   v3_arcsec=None, verbose=False)
        # Add RA, Dec for each source in subarray image
        ra_vals = []
        dec_vals = []
        ff_equiv_x = []
        ff_equiv_y = []

        # Add an empty column to the table listing the distance to the nearest
        # full frame source
        subarray_sources['delta_ff'] = np.zeros(len(subarray_sources['xcentroid'])) + 99.
        for i, source in enumerate(subarray_sources):
            ra, dec = XY_To_RADec(source['xcentroid'], source['ycentroid'], aperture, attitude_matrix)
            ra_vals.append(ra)
            dec_vals.append(dec)

            # Now calculate the pixel value for that RA, Dec on the full frame image
            ffx, ffy = RADec_To_XY(ra, dec, ff_aperture, ff_attitude_matrix)
            ff_equiv_x.append(ffx)
            ff_equiv_y.append(ffy)

            # Check to see if there is a source in the full frame image at this location
            dx = ffx - full_frame_sources['xcentroid']
            dy = ffy - full_frame_sources['ycentroid']
            deltas = np.sqrt(dx**2 + dy**2)
            subarray_sources['delta_ff'][i] = np.nanmin(deltas)

        subarray_sources['ra'] = ra_vals
        subarray_sources['dec'] = dec_vals
        subarray_sources['fullframe_x_equivalent'] = ff_equiv_x
        subarray_sources['fullframe_y_equivalent'] = ff_equiv_y

        # Save the updated table
        outfile = os.path.basename(subfile)
        outname = os.path.join(output_dir, outfile.replace('.fits', '_fullframe_source_comparison.txt'))
        print('Writing out results to {}'.format(outname))
        ascii.write(subarray_sources, outname, overwrite=True)

        # Now check the distribution of differences between full frame source locations,
        # and calculated full frame equivalent locations from the subarray
        # First, throw out the sources that obviously have no match
        diffs = subarray_sources['delta_ff'].data
        good = diffs < 20.
        med = np.median(diffs[good])
        dev = np.std(diffs[good])
        print('\n\nFile: {}'.format(subfile))
        print('Median difference between subarray source locations and full frame locations: {} pixels'.format(med))
        print('Standard deviation of differences: {} pixels\n\n\n'.format(dev))


def run(full_frame_file, subarray_files, output_dir='./'):
    """MAIN FUNCTION

    Paramters
    ---------
    full_frame_file : str

    subarray_files : list
    """
    # Determine which full frame image to compare with which subarray images
    # Detect sources in full-frame image
    # Detect sources in subarray image
    # Get subarray location from pysiaf
    # Add pysiaf result to subarray catalog
    # Compare full frame catalog to shifted subarray catalog to see how well they match.
    # For subarrays not centered on the detector (sub64p/sub400p), this will fail because the reference locations are different!!!
        # In that case, you'll have to use knowledge of the refrence locations to predict source locations

    # Read in full frame data
    ff_model = datamodels.open(full_frame_file)
    ff_radec_to_xy = ff_model.meta.wcs.get_transform('world', 'detector')

    # Calculate the FWHM in pixels to input to the source finder
    ff_fwhm = get_fwhm(ff_model.meta.instrument.filter)

    ffbasename = os.path.basename(full_frame_file)
    full_frame_sources = find_sources(ff_model.data, threshold=500, fwhm=ff_fwhm, plot_name='{}_full_frame_source_map_datamodels.png'.format(ffbasename))
    ascii.write(full_frame_sources, os.path.join(output_dir, '{}_ff_sources_datamodels.txt'.format(ffbasename)), overwrite=True)

    # Read in subarray data and locate sources
    print("Locating sources in subarray files...")
    for subfile in subarray_files:
        sub_model = datamodels.open(subfile)
        sub_xy_to_radec = sub_model.meta.wcs.get_transform('detector', 'world')

        # Calculate the FWHM in pixels to input to the source finder
        sub_fwhm = get_fwhm(sub_model.meta.instrument.filter)
        plotname = os.path.join(output_dir, '{}_source_map_datamodels.png'.format(os.path.basename(subfile)))
        subarray_sources = find_sources(sub_model.data, threshold=500, fwhm=sub_fwhm, plot_name=plotname)

        # If no sources are found in the subarray, then move on to the next
        if subarray_sources is None:
            print('No sources found in {}. Skipping.\n'.format(subfile))
            continue

        # Add RA, Dec for each source in subarray image
        ra_vals = []
        dec_vals = []
        ff_equiv_x = []
        ff_equiv_y = []

        # Add an empty column to the table listing the distance to the nearest
        # full frame source
        subarray_sources['delta_ff'] = np.zeros(len(subarray_sources['xcentroid'])) + 99.
        for i, source in enumerate(subarray_sources):
            ra, dec = sub_xy_to_radec(source['xcentroid'], source['ycentroid'])
            ra_vals.append(ra)
            dec_vals.append(dec)

            # Now calculate the pixel value for that RA, Dec on the full frame image
            ffx, ffy = ff_radec_to_xy(ra, dec)
            ff_equiv_x.append(ffx)
            ff_equiv_y.append(ffy)

            # Check to see if there is a source in the full frame image at this location
            dx = ffx - full_frame_sources['xcentroid']
            dy = ffy - full_frame_sources['ycentroid']
            deltas = np.sqrt(dx**2 + dy**2)
            subarray_sources['delta_ff'][i] = np.nanmin(deltas)

        subarray_sources['ra'] = ra_vals
        subarray_sources['dec'] = dec_vals
        subarray_sources['fullframe_x_equivalent'] = ff_equiv_x
        subarray_sources['fullframe_y_equivalent'] = ff_equiv_y

        # Save the updated table
        outfile = os.path.basename(subfile)
        outname = os.path.join(output_dir, outfile.replace('.fits', '_fullframe_source_comparison_datamodels.txt'))
        print('Writing out results to {}'.format(outname))
        ascii.write(subarray_sources, outname, overwrite=True)

        # Now check the distribution of differences between full frame source locations,
        # and calculated full frame equivalent locations from the subarray
        # First, throw out the sources that obviously have no match
        diffs = subarray_sources['delta_ff'].data
        good = diffs < 20.
        if len(good) == 0:
            print("No matching sources found with an error of less than 20 pixels. Differences "
                  "(in pixels) are: {}. Skipping.\n".format(diffs))
        else:
            print('{} matching sources found between this file and the full frame file.'.format(len(good)))
            med = np.median(diffs[good])
            dev = np.std(diffs[good])
            print('File: {}'.format(subfile))
            print('Median difference between subarray source locations and full frame locations: {} pixels'.format(med))
            print('Standard deviation of differences: {} pixels\n\n\n'.format(dev))


def read_yaml_file(filename):
    try:
        with open(filename, 'r') as infile:
            data = yaml.safe_load(infile)
    except (ScannerError, FileNotFoundError, IOError) as e:
        print(e)
    return data


def get_data(filename):
    """
    Retrieve the data from the given file, along with basic aperture information

    Parameters
    ----------
    filename : str
        Name of fits file

    Returns
    -------
    data : numpy.ndarray
        Multidimensional array of data values

    substrt : tuple
        Tuple of the (x,y) coordinates of the lower left corner of the aperture

    subend : tuple
        Tuple of the (x,y) coordinates of the upper right corner of the aperture

    detector : str
        Detector name

    aperture : str
        Aperture name (e.g. 'SUB160')

    filter_name : str
        Filter name

    header1['RA_REF'] : float
        RA at the reference location

    header1['DEC_REF'] : float
        Dec at the reference location

    header1['PA_V3'] : float
        Telescope roll angle at (V1, V2) = (0,0)
    """
    with fits.open(filename) as hdulist:
        data = hdulist[1].data
        header = hdulist[0].header
        header1 = hdulist[1].header

    data_shape = data.shape
    if len(data_shape) == 3:
        ylen, xlen = data_shape[1], data_shape[2]
    elif len(data_shape) == 2:
        ylen, xlen = data_shape
    else:
        raise ValueError("Expecting 2D or 3D data, but {} dimensions are {}".format(filename, data_shape))

    detector = header['DETECTOR']
    aperture = header['SUBARRAY']
    filter_name = header['FILTER']

    substrt = (header['SUBSTRT1'] - 1, header['SUBSTRT2'] - 1)
    subend = (substrt[0] + xlen, substrt[1] + ylen)
    return data, substrt, subend, detector, aperture, filter_name, header1['RA_REF'], header1['DEC_REF'], header1['PA_V3']
