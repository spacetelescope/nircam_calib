#! /usr/bin/env python

"""This module can be used to confirm that source countrates are the same in full
frame and subarray apertures

pipeline level3 output source catalog columns:
xcentroid ycentroid sky_centroid.ra sky_centroid.dec source_sum source_sum_er abmag abmag_err


Want to compare:

full frame image from one detector
subarray images from same detector

comparing level 3 outputs where data from multiple detectors have been combined complicates the issue

Probably best to simply work with the rate files and use groups of files that are all from the same pointing?

Could keep the level3 source catalog checks as a secondary check...
"""

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u


def matched_catalog_indexes(catalog1, catalog2):
    """Read in and compare two source catalogs produced by the level3 pipeline

    Parameters
    ----------
    catalog1 : astropy.table.Table
        Source catalog

    catalog2 : astropy.table.Table
        Source catalog

    Returns
    -------
    idx : numpy.ndarray
        Array of index values corresponding to the sources in catalog2 that match
        the sources in catalog1
    """
    # Match catalogs using astropy
    c1 = SkyCoord(ra=catalog1['sky_centroid.ra']*u.degree, dec=catalog1['sky_centroid.dec']*u.degree)
    c2 = SkyCoord(ra=catalog2['sky_centroid.ra']*u.degree, dec=catalog2['sky_centroid.dec']*u.degree)
    idx, d2d, d3d = c1.match_to_catalog_sky(c2)

    # c2[idx] matches up with sources in c1

    return idx

def compare_sources(subarray_cat_file, fullframe_cat_file, output_name):
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

    ff_indexes = matched_catalog_indexes(cat1, cat2)

    for index_sub, index_full in enumerate(ff_indexes):
        flux_sub = cat_sub[index_sub]['source_sum']
        fluxerr_sub = cat_sub[index_sub]['source_sum_er']
        mag_sub = cat_sub[index_sub]['abmag']
        magerr_sub = cat_sub[index_sub]['abmag_err']

        flux_full = cat_full[index_full]['source_sum']
        fluxerr_full = cat_full[index_full]['source_sum_er']
        mag_full = cat_full[index_full]['abmag']
        magerr_full = cat_full[index_full]['abmag_err']

        delta_flux = flux_sub - flux_full
        delta_mag = mag_sub - mag_full

        flux_err = np.sqrt(fluxerr_sub**2 + fluxerr_full**2)
        if delta_flux > flux_err:
            print("WARNING: for source {}, the measured fluxes are {} and {}, while the uncertainty is {}".format(index_sub, flux_sub, flux_full, flux_err))

        mag_err = np.sqrt(magerr_sub**2 + magerr_full**2)
        if delta_mag > mag_err:
            print("WARNING: for source {}, the measured abmags are {} and {}, while the uncertainty is {}".format(index_sub, mag_sub, mag_full, mag_err))

        out_ra.append(cat_sub[index_sub]['sky_centroid.ra'])
        out_dec.append(cat_sub[index_sub]['sky_centroid.dec'])
        out_ff_flux.append(flux_full)
        out_sub_flux.append(flux_sub)
        out_flux_err.append(flux_err)
        out_ff_mag.append(mag_full)
        out_sub_mag.append(mag_sub)
        out_mag_err.append(mag_err)

    output_cat['RA'] = out_ra
    output_cat['Dec'] = out_dec
    output_cat['Fullframe_flux'] = out_ff_flux
    output_cat['Subarray_flux'] = out_sub_flux
    output_cat['Flux_err'] = out_flux_err
    output_cat['Fullframe_ABMag'] = out_ff_mag
    output_cat['Subarray_ABMag'] = out_sub_mag
    output_cat['Mag_err'] = out_mag_err
    ascii.write(output_cat, output_name)












