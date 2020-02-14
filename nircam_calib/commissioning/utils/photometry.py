#! /usr/bin/env python

"""This module contains functions related to photometry that may be generally
useful to analysis of multiple CARs
"""
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
import numpy as np
from photutils import CircularAperture, DAOStarFinder


def find_sources(data, threshold=30, fwhm=3.0, show_sources=True, plot_name='sources.png'):
    """Use photutils' DAOFind to locate sources in the input image.

    Parameters
    ----------
    data : numpy.ndarray
        2D image

    threshold : float
        Number of sigma above the background used to determine the
        presence of a source

    show_sources : bool
        If True, create an image with overlain marks showing sources

    plot_name : str
        Name of file to save the image with marked sources

    Returns
    -------
    sources : astropy.table.Table
        Table of source positions
    """
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std)
    sources = daofind(data - median)
    print('{} sources found.'.format(len(sources)))

    if show_sources:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.)
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.savefig(plot_name)
        print('Plot saved to {}'.format(plot_name))
        plt.close()

    return sources


def fwhm(wavelength):
    """Calculate the approximate FWHM, in arcsec, for a given wavelength.
    This equation was calculated by hand using approximate values for the
    FWHM found on JDox:
    https://jwst-docs.stsci.edu/near-infrared-camera/nircam-predicted-performance/nircam-point-spread-functions#NIRCamPointSpreadFunctions-PSFFWHM

    Parameters
    ----------
    wavelength : float
        Wavelength, in microns

    Returns
    -------
    fwhm : float
        FWHM of the PSF, in units of arcsec
    """
    return 0.03 * wavelength + 0.004


def get_fwhm(filter_name):
    """Calculate the FWHM in units of pixels given a filter name

    Parameters
    ----------
    filter_name : str
        e.g. F200W

    Returns
    -------
    FWHM of the PSF, in units of pixels
    """
    filter_name = filter_name.upper()
    try:
        wave = float(filter_name[1:4]) / 100.
    except ValueError:
        raise ValueError("ERROR: Cannot get a wavelength from the filter name: {}".format(filter_name))

    fwhm_arcsec = fwhm(wave)
    if wave < 2.4:
        pix_scale = 0.031  # arcsec/pixel
    else:
        pix_scale = 0.062  # arcsec/pixel
    return fwhm_arcsec / pix_scale
