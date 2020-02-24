#! /usr/bin/env python

"""This module can be used to confirm that source countrates are the same in full
frame and subarray apertures

Do this in two places:

1) on rate.fits images
2) on the final combined mosaics from the level 3 pipeline

pipeline level3 output source catalog columns:
xcentroid ycentroid sky_centroid.ra sky_centroid.dec source_sum source_sum_er abmag abmag_err


Want to compare:

full frame image from one detector
subarray images from same detector

comparing level 3 outputs where data from multiple detectors have been combined complicates the issue

Probably best to simply work with the rate files and use groups of files that are all from the same pointing?

Could keep the level3 source catalog checks as a secondary check...


level 2a check:

input: rate.fits images
main function: compare_rate_images
find sources in subarrays
find sources in full frames
translate (x,y) from subarrays to full frame coords
compare rates

level3 check:

input: level3 catalog files
main function: compare_level3_catalogs
match catalogs
compare rates
"""
from astropy import units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from jwst import datamodels
import numpy as np
import os
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus

from nircam_calib.commissioning.utils.photometry import find_sources, get_fwhm


def compare_level3_catalogs(subarray_cat_file, fullframe_cat_file, output_name):
    """MAIN FUNCTION"""
    # Read in catalogs
    cat_sub = ascii.read(subarray_cat_file)
    cat_full = ascii.read(fullframe_cat_file)

    # Create output catalog and lists for results
    output_cat = Table()
    out_ra = []
    out_dec = []
    out_ff_flux = []
    out_sub_flux = []
    out_flux_err = []
    out_ff_mag = []
    out_sub_mag = []
    out_mag_err = []

    # Match catalogs
    c_sub = SkyCoord(ra=cat_sub['sky_centroid.ra'] * u.degree, dec=cat_sub['sky_centroid.dec'] * u.degree)
    c_full = SkyCoord(ra=cat_full['sky_centroid.ra'] * u.degree, dec=cat_full['sky_centroid.dec'] * u.degree)
    idx, d2d, d3d = c_sub.match_to_catalog_sky(c_full)

    # Remove bad matches
    max_sep = 0.06 * u.arcsec  # Match 1 LW pix or 2 SW pix
    sep_constraint = d2d < max_sep
    sub_catalog_matches = cat_sub[sep_constraint]
    ff_catalog_matches = cat_full[idx[sep_constraint]]
    print('Found {} matching sources.'.format(len(ff_catalog_matches)))

    delta_flux = ff_catalog_matches['source_sum'] - sub_catalog_matches['source_sum']
    delta_mag = ff_catalog_matches['abmag'] - sub_catalog_matches['abmag']
    flux_err = np.sqrt(ff_catalog_matches['source_sum_er']**2 + sub_catalog_matches['source_sum_er']**2)

    output_cat['fullframe_RA'] = ff_catalog_matches['sky_centroid.ra']
    output_cat['fullframe_Dec'] = ff_catalog_matches['sky_centroid.dec']
    output_cat['fullframe_flux'] = ff_catalog_matches['source_sum']
    output_cat['fullframe_fluxerr'] = ff_catalog_matches['source_sum_er']
    output_cat['fullframe_mag'] = ff_catalog_matches['abmag']
    output_cat['fullframe_magerr'] = ff_catalog_matches['abmag_err']
    output_cat['subarray_RA'] = sub_catalog_matches['sky_centroid.ra']
    output_cat['subarray_Dec'] = sub_catalog_matches['sky_centroid.dec']
    output_cat['subarray_flux'] = sub_catalog_matches['source_sum']
    output_cat['subarray_fluxerr'] = sub_catalog_matches['source_sum_er']
    output_cat['subarray_mag'] = sub_catalog_matches['abmag']
    output_cat['subarray_magerr'] = sub_catalog_matches['abmag_err']
    output_cat['delta_flux'] = delta_flux
    output_cat['delta_mag'] = delta_mag
    output_cat['fluxerr_full_minus_sub'] = flux_err
    ascii.write(output_cat, output_name)

    print("Median and stdev percentage flux difference between full frame and subarray: ")
    perc = delta_flux / ff_catalog_matches['source_sum'] * 100.
    print(np.median(perc), np.std(perc))


def compare_rate_images(subarray_file, fullframe_file):
    """MAIN FUNCTiON FOR COMPARING SOURCE RATES IN RATE FILES
    """
    # Read in data
    subarray = datamodels.open(subarray_file)
    fullframe = datamodels.open(fullframe_file)

    # Locate sources
    sub_fwhm = get_fwhm(subarray.meta.instrument.filter)
    sub_sources = find_sources(subarray.data, threshold=500, fwhm=sub_fwhm, plot_name='{}_subarray_source_map_datamodels.png'.format(os.path.basename(subarray_file)))
    ff_fwhm = get_fwhm(fullframe.meta.instrument.filter)
    ff_sources = find_sources(fullframe.data, threshold=500, fwhm=ff_fwhm, plot_name='{}_fullframe_map_datamodels.png'.format(os.path.basename(fullframe_file)))

    if sub_sources is None:
        print("No subarray sources to compare.")
        return 0

    if ff_sources is None:
        print("No full frame sources to compare")
        return 0

    # Put subarray sources in full frame detector coordinates
    # 1. Transform to RA, Dec
    # 2. Use WCS of full frame file to put coordinates in full frame coords

    # Define coord transform functions
    sub_xy_to_radec = subarray.meta.wcs.get_transform('detector', 'world')
    ff_radec_to_xy = fullframe.meta.wcs.get_transform('world', 'detector')
    ff_xy_to_radec = fullframe.meta.wcs.get_transform('detector', 'world')

    # Transform subarray x,y -> RA, Dec -> full frame x, y
    ra, dec = sub_xy_to_radec(sub_sources['xcentroid'], sub_sources['ycentroid'])
    ffx, ffy = ff_radec_to_xy(ra, dec)

    # Add RA, Dec and fullframe equivalent x,y to the subarray source catalog
    sub_sources['RA'] = ra
    sub_sources['Dec'] = dec
    sub_sources['fullframe_x'] = ffx
    sub_sources['fullframe_y'] = ffy

    # Find RA, Dec of sources in full frame catalog
    ffra, ffdec = ff_xy_to_radec(ff_sources['xcentroid'], ff_sources['ycentroid'])
    ff_sources['RA'] = ffra
    ff_sources['Dec'] = ffdec

    # Match catalogs
    ff_cat = SkyCoord(ra=ffra * u.degree, dec=ffdec * u.degree)
    sub_cat = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    idx, d2d, d3d = sub_cat.match_to_catalog_3d(ff_cat)

    # Remove bad matches
    pix_scale = subarray.meta.wcsinfo.cdelt1 * 3600.
    max_sep = 1.0 * pix_scale * u.arcsec  # Matches must be within a pixel
    sep_constraint = d2d < max_sep
    sub_catalog_matches = sub_sources[sep_constraint]
    ff_catalog_matches = ff_sources[idx[sep_constraint]]
    num_matched = len(ff_catalog_matches)
    print('Found {} matching sources.'.format(num_matched))
    if num_matched == 0:
        print("Quitting.")
        return 0

    # Now perform aperture photometry on the sources in the subarray
    # and full frame data. Keep the aperture small since the difference
    # in exposure time and SNR will be large
    sub_pos = [(m['xcentroid'], m['ycentroid']) for m in sub_catalog_matches]
    ff_pos = [(m['xcentroid'], m['ycentroid']) for m in ff_catalog_matches]

    sub_aperture = CircularAperture(sub_pos, r=3.)
    ff_aperture = CircularAperture(ff_pos, r=3.)
    sub_annulus = CircularAnnulus(sub_pos, r_in=10, r_out=15)
    full_annulus = CircularAnnulus(ff_pos, r_in=10, r_out=15)

    # Photometry
    sub_phot_table = aperture_photometry(subarray.data, sub_aperture)
    ff_phot_table = aperture_photometry(fullframe.data, ff_aperture)
    sub_annulus_masks = sub_annulus.to_mask(method='center')
    full_annulus_masks = full_annulus.to_mask(method='center')

    sub_bkg_median = median_background(sub_annulus_masks, subarray.data)
    sub_phot_table['annulus_median'] = sub_bkg_median
    sub_phot_table['aper_bkg'] = sub_bkg_median * sub_aperture.area
    sub_phot_table['aper_sum_bkgsub'] = sub_phot_table['aperture_sum'] - sub_phot_table['aper_bkg']

    full_bkg_median = median_background(full_annulus_masks, fullframe.data)
    ff_phot_table['annulus_median'] = full_bkg_median
    ff_phot_table['aper_bkg'] = full_bkg_median * ff_aperture.area
    ff_phot_table['aper_sum_bkgsub'] = ff_phot_table['aperture_sum'] - ff_phot_table['aper_bkg']

    # Compare photometry results
    delta_phot = ff_phot_table['aper_sum_bkgsub'].data - sub_phot_table['aper_sum_bkgsub'].data
    delta_phot_perc = delta_phot / ff_phot_table['aper_sum_bkgsub'].data * 100.
    sub_phot_table['delta_from_fullframe'] = delta_phot
    sub_phot_table['delta_from_fullframe_percent'] = delta_phot_perc

    print("Percent diff't from subarray: ")
    print(delta_phot_perc)
    print("Difference: full frame - subarray:")
    print(delta_phot)
    print('Median difference: {} MJy/sr'.format(np.median(delta_phot)))
    print('Median percentage difference: {}%'.format(np.median(delta_phot_perc)))


def median_background(annulus_masks, data):
    """Calculate the sigma clipped background value for a list of annulus
    apertures

    Parameters
    ----------
    annulus_masks : photutils.CircularAnnulus.mask
        Masks of the locations of the annuli within which to calculate background

    data : numpy.ndarray
        2D image within which to calculate backgrounds

    Returns
    -------
    bkg_median : numpy.ndarray
        1D array of sigma clipped median background values
    """
    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
    return np.array(bkg_median)



