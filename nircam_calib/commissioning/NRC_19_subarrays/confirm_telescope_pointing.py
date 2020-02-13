#! /usr/bin/env python

"""This module is intended to be used to confirm that telescope pointing is
correct when taking subarray data. In order to do this, it checks that the
target is placed at the expected position on the detector.

FROM THE CAP notes:
When observing with module B and the extended source subarrays, the telescope
should point (before any dithering) such that the target RA and Dec are centered
in the LW aperture. This corresponds to a location between the 4 SW detectors.
(V2, V3 = (-87.081103, -491.730129). This is the ref location for B5.) When using
the point source subarrays, the target RA and Dec should be centered in the SW
aperture, which corresponds to the upper right quadrant of the LW aperture (in DMS
orientation) in SUB400P, and the upper center of SUB64P.

For extended source subarrays (SUB160, SUB320, SUB640, FULL), the target RA, Dec are:
05:21:57.00, -69:29:51.00, which is: 80.4875, -69.4975

For point source subarrays (SUB400P:, S?Ub64P), the target RA, Dec are:
05:21:55.8187, -69:29:38.25, which is: 80.482577917, -69.493958333


"""

import numpy as np

extended_subarr_ref_ra = 80.4875
extended_subarr_ref_dec = -69.4975
ptsrc_subarr_ref_ra = 80.482577917
ptsrc_subarr_ref_dec = -69.493958333


"""
for a given file:
    1) identify the reference location (x,y)
    2) calculate RA, Dec at the reference location
    3) check that RA, Dec matches the target RA, Dec

for point source files:
    1) check to see that the intended target is located at the reference location,
    using source finding.

"""

def check_via_wcs(filename):
    """This checks that the WCS in the file agrees with what is expected.
    It doesn't check that the actual pointing is correct.

    LW: have target RA, Dec.
        have RA, Dec at reference location
        should be the same after accounting for dithers

    SW: have target RA, Dec (whcih can be outside detector fov)
        have reference location RA, Dec
        calc v2,v3 of target, confirm it's the same as in LW?
        calculate delta bet. target and ref loc, and compare with v2,v3 of ref loc of pointing (between sw detectors)

    """
    with fits.open(filename) as hdulist:
        h0 = hdulist[0].header
        h1 = hdulist[1].header

    # Get the reported RA, Dec at the reference location
    ra_at_ref = h1['RA_REF']
    dec_at_ref = h1['DEC_REF']

    # Get the target RA, Dec
    targ_ra = h0['TARG_RA']
    targ_dec = h0['TARG_DEC']

    # Account for dither
    ditherx = header1['X_OFFSET']
    dithery = header1['Y_OFFSET']

    translate dithers to ra, dec

    # Calculate the pixel location of the target
    if 'LONG' in h0['DETECTOR']:
        delta = ra_at_ref - (targ_ra + ditherra)




def check_via_source_location(filename, location_threshold=2):
    """This function identifies sources in the input image,
    and checks that there is a source located at the reference
    location (after accounting for dithers), and that this corresponds
    to the target RA, Dec. This will only work for the point source
    subarrays, where the target in APT is a particular star
    """
    with fits.open(file) as hdulist:
        data = hdulist[1].data
        header0 = hdulist[0].header
        header1 = hdulist[1].header

    # FIND_SOURCES from confirm_subarray_location_via_sources -- move to utils
    sources = find_sources(data, threshold=30, show_sources=False)

    # Account for dither
    ditherx = header1['X_OFFSET']
    dithery = header1['Y_OFFSET']

    # Get the reference location
    refx = header0['CRPIX1']
    refy = header0['CRPIX2']

    # Check to see if the source is at the reference location
    dx = sources['xcentroid'].data - (refx + ditherx)
    dy = sources['ycentroid'].data - (refy + dithery)
    delta = np.sqrt(dx**2 + dy**2)
    closest = np.nanmin(delta)

    if closest > location_threshold:
        print(("WARNING: In {}: Source closest to reference location is {} pixels away. This is "
               "more than the threshold of {} pixels.".format(os.path.basename(filename), closest, location_threshold)))
    else:
        print("In {}: found a source {} pixels from the reference location, as expected.".format(os.path.basename(filename), closest))




def check_targ_ra_dec(hdu, expected_ra, expected_dec):
    """Confirm that the values of the TARG_RA and TARG_DEC are as expected

    Parameters
    ----------
    hdu : astropy.io.fits.header
        Header to be checked. Should be extension 0 in JWST files

    expected_ra : float
        Expected target RA in decimal degrees

    expected_dec : float
        Expected target Dec in decimal degrees
    """
    assert np.isclose(float(hdu['TARG_RA']), expected_ra, rtol=0., atol=2.8e-7)
    assert np.isclose(float(hdu['TARG_DEC']), expected_dec, rtol=0, atol=2.8e-7)

def ra_dec_at_ref_loc()