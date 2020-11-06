#! /usr/bin/env python

"""
Search for latent sources in NIRCam data

Using _cal.fits files:
1. Locate and do aperture photometry on sources in exposure 1
2. Locate and do aperture photometry on sources in exposure 2
3. In exposure 2, do aperture photometry on locations from exposure 1, excluding locations where there is also a source in exposure 2.
4. Create plot of "empty" location photometry vs time since dither?
5. Repeat for subsequent exposures
"""
from astropy.io import ascii, fits
import copy
import matplotlib.pyplot as plt
import numpy as np
import os

from nircam_calib.commissioning.utils.fileio import integration_times
from nircam_calib.commissioning.utils.photometry import do_photometry, find_sources, fwhm, get_fwhm


def check(file_list, rate_sources_file=None):
    """file_list is a list of filenames. - rateints files. ALSO ASSUME THE RATE FILES ARE PRESENT
    Get the source locations from the first file in the list. Then do photometry on those
    source locations in all integrations of the first file, as well as all integrations
    in subsequent files, in order to try to create a plot of latent signal versus time.
    Note that in the first file the sources will be present for all integrations, so
    latent signal will show up only as "extra" counts on top of the expected counts.
    """
    # Put files in chronological order
    #filenames = put_in_order(file_list)

    # Read in rate version of first file. Find sources. Photometry.
    ratefile = file_list[0].replace('rateints.fits', 'rate.fits')
    with fits.open(ratefile) as hdulist:
        rate_data = hdulist[1].data
        filtername = hdulist[0].header['FILTER']
        rate_file_ints = hdulist[0].header['NINTS']

    # Find sources
    rate_fwhm = get_fwhm(filtername)
    rate_sources = find_sources(rate_data, threshold=50, fwhm=rate_fwhm, plot_name='{}_srcmap.png'.format(os.path.basename(ratefile)))
    rate_sources['xcentroid'].info.format = '%.5g'  # for consistent table output
    rate_sources['ycentroid'].info.format = '%.5g'  # for consistent table output

    print('Sources from rate file: ')
    print(rate_sources)

    # Average photometry of sources across the exposure
    rate_phot = do_photometry(rate_data, rate_sources, aperture_radius=3, subtract_background=True, annulus_radii=(25, 30))

    # Save the source catalog and photometry from the rate file
    if rate_sources_file is None:
        rate_sources_file = 'source_catalog_{}'.format(ratefile)
    print('Saving rate file source catalog to: {}'.format(rate_sources_file))
    ascii.write(rate_phot, rate_sources_file, overwrite=True)

    # Now work on the (user-provided) rateints files. Perform photometry
    # separately on each integration in the locations found above.
    # Remember that in the first file, the sources will be present in
    # all integrations
    all_times = []
    for index, filename in enumerate(file_list):

        print('Working on file: {}'.format(filename))

        # Read in file
        with fits.open(filename) as hdulist:
            data = hdulist[1].data
            filtername = hdulist[0].header['FILTER']

            # Get timing information
            print('FIX INT_TIMES bug in Mirage and then change over to using that extension')
            #int_starts, int_mids, int_ends = integration_times(hdulist)
            exstart = hdulist[0].header['EXPSTART']
            exend = hdulist[0].header['EXPEND']
            inttime = (exend - exstart) / data.shape[0]
            int_starts = exstart + inttime * np.arange(data.shape[0])
            int_mids = int_starts + inttime / 2.
            int_ends = int_starts + inttime

        datadims = data.shape
        if len(datadims) != 3:
            raise ValueError("Expecting a 3D data array in {}. Should be a rateints file.".format(filename))
        nint, ny, nx = datadims

        all_times = np.concatenate([all_times, int_starts])

        # For each file, do photometry on all the locations of the sources
        # in the subsequent exposures
        for integ in range(nint):
            phot = do_photometry(data[integ, :, :], rate_sources, aperture_radius=3, subtract_background=True, annulus_radii=(25, 30))

            if 'aper_sum_bkgsub' in phot.colnames:
                results = phot['aper_sum_bkgsub'].data
            else:
                results = phot['aperture_sum'].data
            if index == 0 and integ == 0:
                phot_results = copy.deepcopy(results)
            else:
                phot_results = np.vstack([phot_results, results])

    # Index times to the start time of the first file and put into seconds
    all_times = (all_times - all_times[0]) * 24. *3600.

    # Here we plot photometry of uncontaminated sources vs time?
    # phot_results[:, 0] is the photometry for source 0 at all times
    # If results are background subtracted, then we should look for the
    # time at which the signal is consistent with zero.
    ntimes, nsources = phot_results.shape
    for num in range(nsources):
        #for num in range(1):
        f, a = plt.subplots()
        a.scatter(all_times, phot_results[:, num], color='red')
        a.set_title('Source #{}: (x, y) = ({}, {})'.format(num, rate_sources['xcentroid'][num], rate_sources['ycentroid'][num]))
        a.set_xlabel('Time (sec)')
        a.set_ylabel('Source Signal (ADU/sec)')
        outname = 'signal_v_time_source_{}.png'.format(num)
        plt.tight_layout()
        f.savefig(outname)
        #plt.show()

