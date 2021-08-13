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


def compare_level3_catalogs(cat_file_1, cat_file_2, output_name, out_dir='./', max_match_thresh=0.06,
                            ee_perc=70, limiting_vegamag=20.):
    """MAIN FUNCTION
    Compare source fluxes in level 3 catalogs

    Parameters
    ----------
    cat_file_1 : str
        Catalog name

    cat_file_2 : str
        Catalog name

    output_name : str
        Name of ascii file for outputs

    out_dir : str
        Directory to write the output file into

    max_match_thresh : float
        Maximum distance in arcsec that will constitute a match when matching sources between the two catalogs

    ee_perc : int
        Encircled energy percentage to use for the comparison. The default ee percentages used
        by the pipeline are 30%, 50%, and 70%

    limiting_vegamag : float
        Maximum vegamag (dimmest source) to consider when comparing catalogs
    """
    print('Comparing: {} and {}'.format(os.path.basename(cat_file_1), os.path.basename(cat_file_2)))
    # Read in catalogs
    cat1 = ascii.read(cat_file_1)
    cat2 = ascii.read(cat_file_2)

    # Remove sources that are dimmer than the magnitude threshold
    good1 = cat1['aper{}_vegamag'.format(ee_perc)] < limiting_vegamag
    cat1 = cat1[good1]
    print('Removing sources that are too dim. Keeping {} sources in {}'.format(len(cat1), cat_file_1))

    good2 = cat2['aper{}_vegamag'.format(ee_perc)] < limiting_vegamag
    cat2 = cat2[good2]
    print('Removing sources that are too dim. Keeping {} sources in {}'.format(len(cat2), cat_file_2))

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
    #c_sub = SkyCoord(ra=cat1['sky_centroid.ra'] * u.degree, dec=cat1['sky_centroid.dec'] * u.degree)
    #c_full = SkyCoord(ra=cat2['sky_centroid.ra'] * u.degree, dec=cat2['sky_centroid.dec'] * u.degree)
    idx, d2d, d3d = cat1['sky_centroid'].match_to_catalog_sky(cat2['sky_centroid'])

    # Remove bad matches
    max_sep = max_match_thresh * u.arcsec  # Match 1 LW pix or 2 SW pix
    sep_constraint = d2d < max_sep
    catalog_1_matches = cat1[sep_constraint]
    catalog_2_matches = cat2[idx[sep_constraint]]
    print('Found {} matching sources.'.format(len(catalog_2_matches)))

    aperstr = 'aper{}'.format(ee_perc)
    delta_flux = catalog_2_matches['{}_flux'.format(aperstr)] - catalog_1_matches['{}_flux'.format(aperstr)]
    delta_flux_perc = delta_flux / catalog_2_matches['{}_flux'.format(aperstr)] * 100.
    delta_mag = catalog_2_matches['{}_abmag'.format(aperstr)] - catalog_1_matches['{}_abmag'.format(aperstr)]
    flux_err = np.sqrt(catalog_2_matches['{}_flux_err'.format(aperstr)]**2 + catalog_2_matches['{}_flux_err'.format(aperstr)]**2)

    output_cat['cat1_pos'] = catalog_1_matches['sky_centroid']
    output_cat['cat1_flux'] = catalog_1_matches['{}_flux'.format(aperstr)]
    output_cat['cat1_fluxerr'] = catalog_1_matches['{}_flux_err'.format(aperstr)]
    output_cat['cat1_mag'] = catalog_1_matches['{}_abmag'.format(aperstr)]
    output_cat['cat1_magerr'] = catalog_1_matches['{}_abmag_err'.format(aperstr)]
    output_cat['cat2_pos'] = catalog_2_matches['sky_centroid']
    output_cat['cat2_flux'] = catalog_2_matches['{}_flux'.format(aperstr)]
    output_cat['cat2_fluxerr'] = catalog_2_matches['{}_flux_err'.format(aperstr)]
    output_cat['cat2_mag'] = catalog_2_matches['{}_abmag'.format(aperstr)]
    output_cat['cat2_magerr'] = catalog_2_matches['{}_abmag_err'.format(aperstr)]
    output_cat['delta_flux'] = delta_flux
    output_cat['delta_flux_perc'] = delta_flux_perc
    output_cat['delta_mag'] = delta_mag
    output_cat['fluxerr_full_minus_sub'] = flux_err
    output_cat.meta['comments'] = ['cat1 = {}'.format(cat_file_1), 'cat2 = {}'.format(cat_file_2)]
    ascii.write(output_cat, os.path.join(out_dir, output_name), overwrite=True)

    #print(output_cat['cat1_pos', 'cat1_flux', 'cat2_flux', 'delta_flux', 'delta_flux_perc'])

    #print("Median and stdev percentage flux difference between the matched sources in the two catalogs: ")
    perc = delta_flux / catalog_2_matches['{}_flux'.format(aperstr)] * 100.
    print("Median percentage flux difference between matched sources: {:.2f}%".format(np.median(perc)))
    print('Stdev of percentage flux differences between matched sources: {:.2f}%'.format(np.std(perc)))
    print('See {} for the list of matched sources\n\n'.format(os.path.join(out_dir, output_name)))
    return output_cat


def compare_rate_images(subarray_file, fullframe_file, out_dir='./'):
    """MAIN FUNCTiON FOR COMPARING SOURCE RATES IN RATE FILES

    Parameters
    ----------
    subarray_file : str
        Fits file containing the subarray data to be compared

    fullframe_file : sr
        Fits file containing the full frame data to be compared

    out_dir : str
        Output directory into which products are saved

    Returns
    -------
    sub_sources : astropy.table.Table
        Table of source positions, equivalent full frame positions, and photometry
        results from the data
    """
    # Read in data
    subarray = datamodels.open(subarray_file)
    fullframe = datamodels.open(fullframe_file)

    # Locate sources
    sub_fwhm = get_fwhm(subarray.meta.instrument.filter)
    sub_source_image = os.path.join(out_dir, '{}_subarray_source_map_countrate_compare.png'.format(os.path.basename(subarray_file)))
    sub_sources = find_sources(subarray.data, threshold=200, fwhm=sub_fwhm, show_sources=True, save_sources=True, plot_name=sub_source_image)
    ff_fwhm = get_fwhm(fullframe.meta.instrument.filter)
    full_source_image = os.path.join(out_dir, '{}_fullframe_map_datamodels.png'.format(os.path.basename(fullframe_file)))
    ff_sources = find_sources(fullframe.data, threshold=200, fwhm=ff_fwhm, show_sources=True, save_sources=True, plot_name=full_source_image)

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
        print("No matching sources found. Quitting.\n\n\n")
        return 0, 0

    # What if we just do photometry on the sources in the subarray and their
    # calculated positions in the full frame?
    sub_pos = []
    ff_pos = []
    good_indexes = []
    non_zero_sub_dq = []
    non_zero_full_dq = []
    for i, line in enumerate(sub_sources):
        if ((line['fullframe_x'] > 5) & (line['fullframe_x'] < 2039) & \
            (line['fullframe_y'] > 5) & (line['fullframe_y'] < 2039)):
            sub_pos.append((line['xcentroid'], line['ycentroid']))
            ff_pos.append((line['fullframe_x'], line['fullframe_y']))
            good_indexes.append(i)

            # Make a note if there are any non-zero DQ flags within the apertures
            # During development, found some cases with NO_LIN_CORR and NO_FLAT_FIELD
            # that were screwing up photometry
            yc_sub = np.int(np.round(line['ycentroid']))
            xc_sub = np.int(np.round(line['xcentroid']))
            sub_dq = subarray.dq[yc_sub-3:yc_sub+4, xc_sub-3:xc_sub+4]
            if np.sum(sub_dq) > 0:
                non_zero_sub_dq.append(True)
            else:
                non_zero_sub_dq.append(False)

            yc = np.int(np.round(line['fullframe_x']))
            xc = np.int(np.round(line['fullframe_y']))
            sub_dq = fullframe.dq[yc-1:yc+2, xc-1:xc+2]
            if np.sum(sub_dq) > 0:
                non_zero_full_dq.append(True)
            else:
                non_zero_full_dq.append(False)


    print('Performing photometry on a total of {} sources.'.format(len(sub_pos)))


    # Now perform aperture photometry on the sources in the subarray
    # and full frame data. Keep the aperture small since the difference
    # in exposure time and SNR will be large
    #sub_pos = [(m['xcentroid'], m['ycentroid']) for m in sub_catalog_matches]
    #ff_pos = [(m['xcentroid'], m['ycentroid']) for m in ff_catalog_matches]


    sub_aperture = CircularAperture(sub_pos, r=5.)
    ff_aperture = CircularAperture(ff_pos, r=5.)
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

    # Keep track of whether there are bad pixels in the apertures
    sub_dq = np.zeros(len(sub_sources), dtype=bool)
    sub_dq[good_indexes] = non_zero_sub_dq
    sub_sources['sub_dq'] = sub_dq
    full_dq = np.zeros(len(sub_sources), dtype=bool)
    full_dq[good_indexes] = non_zero_full_dq
    sub_sources['full_dq'] = full_dq

    # Add photometry to the table
    sub_phot_data = np.zeros(len(sub_sources))
    sub_phot_data[good_indexes] = sub_phot_table['aper_sum_bkgsub'].data
    ff_phot_data = np.zeros(len(sub_sources))
    ff_phot_data[good_indexes] = ff_phot_table['aper_sum_bkgsub'].data
    delta_phot_col = np.zeros(len(sub_sources))
    delta_phot_col[good_indexes] = delta_phot
    delta_phot_perc_col = np.zeros(len(sub_sources))
    delta_phot_perc_col[good_indexes] = delta_phot_perc

    sub_sources['sub_phot'] = sub_phot_data
    sub_sources['ff_phot'] = ff_phot_data
    sub_sources['d_phot'] = delta_phot_col
    sub_sources['d_phot_p'] = delta_phot_perc_col

    sub_sources['xcentroid'].info.format = '7.3f'
    sub_sources['ycentroid'].info.format = '7.3f'
    sub_sources['fullframe_x'].info.format = '7.3f'
    sub_sources['fullframe_y'].info.format = '7.3f'
    sub_sources['sub_phot'].info.format = '7.3f'
    sub_sources['ff_phot'].info.format = '7.3f'
    sub_sources['d_phot'].info.format = '7.3f'
    sub_sources['d_phot_p'].info.format = '7.3f'
    print(sub_sources['xcentroid', 'ycentroid', 'fullframe_x', 'fullframe_y', 'sub_phot', 'ff_phot', 'd_phot_p', 'sub_dq', 'full_dq'])
    print('')

    # Save the complete table
    sub_base = os.path.basename(subarray_file).replace('.fits', '')
    full_base = os.path.basename(fullframe_file).replace('.fits', '')
    table_name = os.path.join(out_dir, 'photometry_comparison_{}_{}.txt'.format(sub_base, full_base))
    ascii.write(sub_sources, table_name, overwrite=True)
    print('Photometry results saved to: {}'.format(table_name))

    # Try filtering out sources where there is a pixel flagged in the full
    # frame or subarray DQ mask at the source location
    clean = []
    for row in sub_sources:
        clean.append(row['sub_dq'] == False and row['full_dq'] == False)
    if np.sum(clean) > 0:
        clean_table = sub_sources[clean]
        print('Excluding sources with a pixel flagged in the DQ array:')
        med_clean_diff = np.around(np.median(clean_table['d_phot_p']), 1)
        print('Median photometry difference between subarray and full frame sources is: {}%'.format(med_clean_diff))
    else:
        clean_table = None
        print('No sources without a flagged pixel in the DQ arrays.')
    return sub_sources, clean_table


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


def ratio_images(file1, file2, out_dir='./'):
    """Ratio the two input images in order to compare signal rates.
    The input images should have the same pointing

    Parameters
    ----------
    file1 : str
        Fits file containing the first image

    file2 : str
        Fits file containing the second image

    out_dir : str
        Directory into which the ratio image will be saved
        as a fits file

    Returns
    -------
    ratio : numpy.ndarray
        2D ratio image

    med : float
        Median value of the ratio image
    """
    # Get the data
    fobj1 = fits.open(file1)
    fobj2 = fits.open(file2)

    # Crop arrays to match
    cropped_fobj1, cropped_fobj2 = crop_array(fobj1, fobj2)

    # Ratio image
    ratio = cropped_fobj1[1].data / cropped_fobj2[1].data
    med = np.median(ratio)

    # Save the cropped arrays and ratio image
    h0 = fits.PrimaryHDU(cropped_fobj1)
    h1 = fits.ImageHDU(cropped_fobj2)
    h2 = fits.ImageHDU(ratio)
    h0.header = fobj1[0].header
    h1.header = fobj2[0].header
    hdulist = fits.HDUList([h0, h1, h2])

    f1_base = os.path.basename(file1).replace('.fits', '')
    f2_base = os.path.basename(file2).replace('.fits', '')
    outfile = os.path.join(out_dir, 'ratio_image_{}_{}.fits'.format(f1_base, f2_base))
    hdulist.writeto(outfile, overwrite=True)
    print('Ratio image saved to: {}'.format(outfile))

    return ratio, med


def crop_array(hdu1, hdu2):
    """Crop the data from two hdu lists to matching sizes
    and locations on the detector

    Parameters
    ----------
    hdu1 : astropy.fits.HDUList

    hdu2 : astropy.fits.HDUList

    Returns
    -------
    hdu1
    hdu2
    """
    data1 = hdu1[1].data
    header1_0 = hdu1[0].header

    data2 = hdu2[1].data
    header2_0 = hdu2[0].header

    # Check the location of the images with respect to the
    # full detector. Remember that coordinates in the header
    # are indexed to 1, so adjust for python
    xstart1 = header1_0['SUBSTRT1'] - 1
    xsize1 = header1_0['SUBSIZE1']
    xend1 = xstart1 + xsize1
    ystart1 = header1_0['SUBSTRT2'] - 1
    ysize1 = header1_0['SUBSIZE2']
    yend1 = ystart1 + ysize1

    xstart2 = header2_0['SUBSTRT1'] - 1
    xsize2 = header2_0['SUBSIZE1']
    xend2 = xstart2 + xsize2
    ystart2 = header2_0['SUBSTRT2'] - 1
    ysize2 = header2_0['SUBSIZE2']
    yend2 = ystart2 + ysize2

    if xstart1 <= xstart2:
        dx = xstart2 - xstart1
    else:
        dx = xstart1 - xstart2

    if ystart1 <= ystart2:
        dy = ystart2 - ystart1
    else:
        dy = ystart1 - ystart2

    # Only need to compare ysizes to find which subarray is larger
    # and should be cropped. That's because all subarrays are square
    # with the exception of the substripe subarray for TSO data. That
    # subarray doesn't overlap with any of the other apertures besides
    # the full frame. Checking the ysize can distinguish all cases then.
    if ysize1 < ysize2:
        data2 = data2[dy: dy+ysize1, dx: dx+xsize1]
        hdu2[1].data = data2
    elif ysize2 < ysize1:
        data1 = data1[dy: dy+ysize2, dx: dx+xsize2]
        hdu1[1].data = data1

    return hdu1, hdu2
