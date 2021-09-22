#! /usr/bin/env python

"""This module contains functions related to photometry that may be generally
useful to analysis of multiple CARs
"""
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
import numpy as np
from photutils import CircularAnnulus, CircularAperture, DAOStarFinder, \
                      aperture_photometry, source_properties
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.utils import calc_total_error


def do_photometry(image, source_table, aperture_radius=3, subtract_background=False, annulus_radii=(6, 8)):
    """Perform aperture photometry on the input data using the source
    locations specified in the source_table

    Parameters
    ----------
    image : numpy.ndarray
        2D image

    source_table : astropy.table.Table
        Table of source locations (i.e. output from photutils's
        DAOStarFinder)

    subtract_background : bool
        If True, use photutils to estimate and subtract background

    Returns
    -------
    source_table : astropy.table.Table
        Updated source_table including photometry results
    """
    # Create list of tuples giving source locations
    positions = [(tab['xcentroid'], tab['ycentroid']) for tab in source_table]

    apertures = CircularAperture(positions, r=aperture_radius)

    # If background is to be subtracted, calculate the sigma-clipped median
    # value of the pixels in the annuli around the sources
    if subtract_background:
        annulus_apertures = CircularAnnulus(positions, r_in=annulus_radii[0], r_out=annulus_radii[1])
        annulus_masks = annulus_apertures.to_mask(method='center')

        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(image)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)

    # Do the photometry and add to the input table
    phot_table = aperture_photometry(image, apertures)
    source_table['aperture_sum'] = phot_table['aperture_sum']

    if subtract_background:
        # Add columns with that background subtraction info
        source_table['annulus_median'] = bkg_median
        source_table['aper_bkg'] = bkg_median * apertures.area
        source_table['aper_sum_bkgsub'] = source_table['aperture_sum'] - source_table['aper_bkg']
    return source_table


def find_sources(data, threshold=30, fwhm=3.0, show_sources=True, save_sources=False, plot_name='sources.png'):
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

    save_sources : bool
        If True, save the image of the source locations at ```plot_name```

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
    if sources is not None:
        print('{} sources found.'.format(len(sources)))

        if show_sources or save_sources:
            positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
            apertures = CircularAperture(positions, r=4.)
            norm = ImageNormalize(stretch=SqrtStretch())
            plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
            apertures.plot(color='blue', lw=1.5, alpha=0.5)
            if show_sources:
                plt.show()
            if save_sources:
                plt.savefig(plot_name)
                print('Plot saved to {}'.format(plot_name))
            plt.close()
    else:
        print('No sources found.')

    return sources


def find_sources_via_segmentation(data, sigma=3, fwhm=2.0, min_pix=5, make_plot=False):
    """
    """
    yd, xd = data.shape

    # Let's define the background using boxes of ~50x50 pixels
    nboxx = int(xd / 150)
    nboxy = int(yd / 150)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (nboxy, nboxx), filter_size=(3, 3),
                       bkg_estimator=bkg_estimator)

    data -= bkg.background  # subtract the background
    threshold = sigma * bkg.background_rms

    #threshold = detect_threshold(data, nsigma=sigma)

    gaussian_sigma = fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(gaussian_sigma, x_size=3, y_size=3)
    kernel.normalize(mode='integral')
    segm = detect_sources(data, threshold, npixels=min_pix, filter_kernel=kernel)
    segm_deblend = deblend_sources(data, segm, npixels=min_pix,
                                   filter_kernel=kernel, nlevels=32,
                                   contrast=0.001)

    if make_plot:
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12.5))
        ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        cmap = segm.make_cmap(seed=123)
        ax2.imshow(segm, origin='lower', cmap=cmap, interpolation='nearest')
        ax2.set_title('Segmentation Image')
        ax3.imshow(segm_deblend, origin='lower', cmap=cmap, interpolation='nearest')
        ax3.set_title('Deblended Segmentation Image')
        plt.show()
        #plt.save('testing.jpg')

    return data, segm_deblend, bkg


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


def photometry_via_segmentation_image(data, seg_image, background_image, effective_gain, wcs=None):
    """
    """
    error = calc_total_error(data, background_image.background_rms, effective_gain)
    cat = SourceCatalog(data, seg_image, background=background_image.background, error=error, apermask_method='mask', wcs=wcs)

    columns = ['area', 'background_mean', 'sky_centroid', 'cxx', 'cxy', 'cyy', 'ellipticity', 'fwhm',
               'local_background', 'xcentroid', 'ycentroid', 'kron_flux', 'kron_fluxerr', 'segment_flux', 'segment_fluxerr']

    tbl = cat.to_table(columns=columns)
    for colname in columns:
        if colname != 'sky_centroid':
            tbl[colname].info.format = '.4f'
    #tbl['xcentroid'].info.format = '.4f'  # optional format
    #tbl['ycentroid'].info.format = '.4f'
    #tbl['kron_flux'].info.format = '.4f'
    print(tbl)


    # Make background image




    #props = source_properties(data, seg_image, background=background_image.background)
    #tbl = properties_table(props)
    #print(tbl)

    #columns = ['id', 'background_at_centroid', 'background_mean',
    #           'background_sum']
    #tbl = properties_table(props, columns=columns)
    #print(tbl)
