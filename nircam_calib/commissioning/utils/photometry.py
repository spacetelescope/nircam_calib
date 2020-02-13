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


def find_sources(data, threshold=30, show_sources=True, plot_name='sources.png'):
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
    daofind = DAOStarFinder(fwhm=3.0, threshold=threshold*std, peakmax=12000.)
    sources = daofind(data - median)

    if show_sources:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.)
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.savefig(plot_name)

    return sources
