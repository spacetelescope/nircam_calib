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

So really maybe all we need to do is confirm for a given input file
that there is a star located at the RA, Dec of the target star, and that
that RA, Dec correspond to the reference location (+dither)

"""
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from jwst.datamodels import ImageModel
import numpy as np
import os

from nircam_calib.commissioning.utils.photometry import find_sources, fwhm, get_fwhm

extended_subarr_ref_ra = 80.4875
extended_subarr_ref_dec = -69.4975

# This is the RA, Dec of the J=14 star as input in APT
# (NOTE THAT THIS IS ACTUALLY NOT CORRECT. IT'S OFF BY 0.13".
# HOPEFULLY WE WILL CORRECT THE COORDINATES IN THE APT FILE)
STAR_RA = 80.482577917
STAR_DEC = -69.493958333

# Real RA and Dec, from the 2MASS catalog
real_star_ra = 80.482541
real_star_dec = -69.4939583
ra_err = 0.08  # arcsec, from 2MASS
dec_err = 0.06  # arcsec, from 2MASS


"""
for a given file:
    1) identify the reference location (x,y)
    2) calculate RA, Dec at the reference location
    3) Calculate pixel that corresponds to RA, Dec of star
    4) Confirm that the star is at that pixel



for point source files:
    1) check to see that the intended target is located at the reference location,
    using source finding.

"""

def check_pointing_target_star(filename, out_dir='./'):
    """Check that the target star is present at the expected location in the image

    Parameters
    ----------
    filename : str
        Name of fits file containing the observation to check. This file must
        contain a valid GWCS (i.e. it must have gone through the assign_wcs
        pipeline step.)

    out_dir : str
        Name of directory into which source catalogs are written
    """
    model = ImageModel(filename)
    pix_scale = model.meta.wcsinfo.cdelt1 * 3600.

    # Calculate the pixel corresponding to the RA, Dec value of the star
    world2det = model.meta.wcs.get_transform('world', 'detector')
    star_x, star_y = world2det(STAR_RA, STAR_DEC)

    # Check to see if the star is actually there
    sub_fwhm = get_fwhm(model.meta.instrument.filter)
    plot_file = '{}_source_map.png'.format(os.path.basename(filename))
    plot_file = os.path.join(out_dir, plot_file)
    found_sources = find_sources(model.data, threshold=50, fwhm=sub_fwhm, plot_name=plot_file)

    # Calculate the distance from each source to the target
    dx = found_sources['xcentroid'] - star_x
    dy = found_sources['ycentroid'] - star_y
    delta = np.sqrt(dx**2 + dy**2)

    # Add distances to the table and save
    found_sources['delta_from_target'] = delta
    basename = os.path.basename(filename)
    table_out = os.path.join(out_dir, '{}_sources.txt'.format(basename))
    ascii.write(found_sources, table_out, overwrite=True)

    # Closest source - assume this is our target star
    min_delta = np.nanmin(delta)
    print('Target expeted location: ({0:1.2f}, {1:1.2f})'.format(star_x, star_y))
    print(('{0}:\nDistance between calculated target location and measured target '
           'location: {1:1.3f} pixels'.format(basename, min_delta)))
    full_err = np.sqrt(ra_err**2 + dec_err**2)
    print(('Uncertainty in the source location from 2MASS: RA: {}", Dec: {}", '
           'Total: {}"'.format(ra_err, dec_err, full_err)))
    print(('                                             = RA: {0:1.2f} pix, Dec: {1:1.2f} '
           'pix, Total: {2:1.2f} pix\n\n'.format(ra_err / pix_scale, dec_err / pix_scale, full_err / pix_scale)))


def check_pointing_using_2mass_catalog(filename, out_dir='./'):
    """Check that stars from an external 2MASS catalog are present at the expected
    locations within filename.

    Parameters
    ----------
    filename : str
        Name of fits file containing the observation to check. This file must
        contain a valid GWCS (i.e. it must have gone through the assign_wcs
        pipeline step.)

    out_dir : str
        Name of directory into which source catalogs are written
    """
    model = ImageModel(filename)
    pix_scale = model.meta.wcsinfo.cdelt1 * 3600.

    # Read in catalog from 2MASS. We'll be using the positions in this
    # catalog as "truth"
    catalog_2mass = os.path.join(os.path.dirname(__file__), 'car19_2MASS_source_catalog.txt')
    cat_2mass = ascii.read(catalog_2mass)
    ra_unc = cat_2mass['err_maj'] * u.arcsec
    dec_unc = cat_2mass['err_min'] * u.arcsec
    total_unc = np.sqrt(ra_unc**2 + dec_unc**2)

    # Calculate the pixel corresponding to the RA, Dec value of the stars
    # in the 2MASS catalog
    world2det = model.meta.wcs.get_transform('world', 'detector')
    star_x, star_y = world2det(cat_2mass['ra'], cat_2mass['dec'])
    cat_2mass['x'] = star_x
    cat_2mass['y'] = star_y

    # Remove stars that are not on the detector
    ydim, xdim = model.data.shape
    good = ((star_x > 0) & (star_x < xdim) & (star_y > 0) & (star_y < ydim))
    print('{} 2MASS sources should be present on the detector.'.format(np.sum(good)))
    cat_2mass = cat_2mass[good]

    # Find sources in the data
    sub_fwhm = get_fwhm(model.meta.instrument.filter)
    plot_file = '{}_source_map.png'.format(os.path.basename(filename))
    plot_file = os.path.join(out_dir, plot_file)
    sources = find_sources(model.data, threshold=50, fwhm=sub_fwhm, plot_name=plot_file)

    # Calculate RA, Dec of the sources in the new catalog
    det2world = model.meta.wcs.get_transform('detector', 'world')
    found_ra, found_dec = det2world(sources['xcentroid'].data, sources['ycentroid'].data)
    sources['calc_RA'] = found_ra
    sources['calc_Dec'] = found_dec

    # Match catalogs
    found_sources = SkyCoord(ra=found_ra * u.degree, dec=found_dec * u.degree)
    from_2mass = SkyCoord(ra=cat_2mass['ra'] * u.degree, dec=cat_2mass['dec'] * u.degree)
    idx, d2d, d3d = from_2mass.match_to_catalog_3d(found_sources)

    # Filter out bad matches Remember, if you have sources
    # in the data other than 2MASS stars, good matches will not be
    # found for those stars and that will contaminate the median distance.
    #max_sep = 1.0 * u.arcsec
    #sep_constraint = d2d < max_sep
    #c_matches = found_sources[sep_constraint]
    #catalog_matches = from_2mass[idx[sep_constraint]]

    # Print table of sources
    d2d_arcsec = d2d.to(u.arcsec).value
    cat_2mass['image_x'] = sources[idx]['xcentroid']
    cat_2mass['image_y'] = sources[idx]['ycentroid']
    cat_2mass['image_ra'] = sources[idx]['calc_RA']
    cat_2mass['image_dec'] = sources[idx]['calc_Dec']
    cat_2mass['delta_sky'] = d2d_arcsec
    cat_2mass['delta_pix'] = d2d_arcsec / pix_scale

    for colname in ['x', 'y', 'image_x', 'image_y', 'delta_pix']:
        cat_2mass[colname].info.format = '7.3f'

    print(cat_2mass['x', 'y', 'image_x', 'image_y', 'delta_pix'])

    # Get info on median offsets
    med_dist = np.nanmedian(d2d_arcsec)
    dev_dist = np.nanstd(d2d_arcsec)
    print(("Median distance between sources in 2MASS catalog and those found in the "
           "data: {0:1.2f} arcsec = {1:1.2f} pixels".format(med_dist, med_dist / pix_scale)))

    mean_unc = np.mean(total_unc)
    mean_unc_pix = np.mean(total_unc.value) / pix_scale
    print(("Mean uncertainty in the source locations within the 2MASS catalog: {0:1.2f} = {1:1.2f} "
           "pixels\n\n".format(mean_unc, mean_unc_pix)))


def check_targ_ra_dec(hdu, expected_ra, expected_dec):
    """Confirm that the values of the TARG_RA and TARG_DEC are as expected
    The RA, Dec value at the reference location should match the target
    RA, Dec in the APT file.  This should only be true for the LW channel
    data though. It should not be true for each of the 4 SW detectors.

    Parameters
    ----------
    hdu : astropy.io.fits.header
        Header to be checked. Should be extension 0 in JWST files

    expected_ra : float
        Expected target RA in decimal degrees

    expected_dec : float
        Expected target Dec in decimal degrees
    """
    channel = hdu['CHANNEL']

    # This function won't work for SW data
    if channel.upper() == 'SW':
        print('check_targ_ra_dec not intended for SW data.')
        return 0

    assert np.isclose(float(hdu['TARG_RA']), expected_ra, rtol=0., atol=2.8e-7)
    assert np.isclose(float(hdu['TARG_DEC']), expected_dec, rtol=0, atol=2.8e-7)

    print("Pointing appears to be correct.")
