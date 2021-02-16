#! /usr/bin/env python

"""This module can be used to create a dark current refrence file that contains values
of zero for all nominal pixels. Any pixels with an elevated dark current signal will
have non-zero values.

For each hot/warm pixel, find a mean slope and then extrapolate to the RAPID frame times.
For this to be correct we'll also have to unlinearize the signal, since dark correction
is done before linearity correction. The other option would be to simply calculate averages
frame-by-frame and use those. That method may give more noisy results.
"""
from astropy.io import fits
from astropy.stats import sigma_clip
from copy import deepcopy
from glob import glob
from jwst.datamodels import DarkModel, dqflags
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
import numpy as np
import os
from scipy.stats import sigmaclip



def create_dark(input_dark_files, input_dark_slope_files=None, sigma_threshold_for_hot=5, total_groups=187, output_file=None,
                output_slope_file_dir='./', pedigree='GROUND', author='Hilbert', descrip=None, use_after=None, history=None):
    """Given a list of individual dark current ramps, create a dark current reference
    file with zeros for nominal pixels and non-zero values for hot/warm pixels.
    Assume that the individual input darks have been processed through the
    calibration pipeline to the appropriate level such that the data can be
    used without any further calibration (i.e. the dark pipeline)

    Parameters
    ----------
    input_dark_files : list
        List of calibrated dark exposure files

    input_dark_slope_files : list
        List of slope files corresponding to the files
        in ``input_dark_files``. If None, then the ramp_fit
        step of the pipeline is called on the the ``input_dark_files``
        and the resulting slope images are used.

    sigma_threshold_for_hot : int
        Number of sigma above the mean to use for the threshold that
        defines warm/hot pixels. Pixels with signal below this level
        are considered nominal and will have signal values of 0 used
        in the dark. Pixels above this threshold will have dark signals
        based on their mean dark rate.

    total_groups : int
        Number of groups in the resulting dark reference file. Since
        the signal values in the reference file will be created by
        extrapolating the slope out to the exposure time of each read,
        this can be any number. Default is 187, as that is enough
        reads to cover a 2000 second integration (assuming a frame time
        of 10.73677 seconds). This is double the recommended maximum
        exposure time, so it should be sufficient? Or is that too large,
        and therefore wasting time through I/O of data that will never
        be used?

    output_file : str
        Name of file to save new reference file data into. If None,
        a name will be constructed based on instrument/detector/etc

    output_slope_file_dir : str
        Name of directory where slope files will be saved if they are
        to be created.

    Returns
    -------

    """
    # Create slope images if needed
    if input_dark_slope_files is None:
        for filename in input_dark_files:
            run_pipeline(filename, outdir=output_slope_file_dir)
        input_dark_files = sorted(glob(os.path.join(output_slope_file_dir, '*linearity.fits')))
        input_dark_slope_files = sorted(glob(os.path.join(output_slope_file_dir, '*rate.fits')))
        print('Slope files created and saved.')
    else:
        print('Assuming that inputs are linearized ramps that can be analyzed directly, as well as slope images.')
        print('No pipeline steps will be run.')

    print('Ramp files and corresponding slope files are: ')
    for ramp_file, slope_file in zip(input_dark_files, input_dark_slope_files):
        print(ramp_file, slope_file)

    # Collect needed metadata from one of the input dark ramp files
    metadata = collect_metadata(input_dark_files[0])
    metadata['PEDIGREE'] = pedigree
    metadata['AUTHOR'] = author

    if descrip is None:
        raise ValueError('descrip cannot be None.')
    metadata['DESCRIP'] = descrip

    if use_after is None:
        raise ValueError('use_after cannot be None')
    metadata['USE_AFTER'] = use_after

    if history is None:
        raise ValueError('history cannot be None')
    metadata['HISTORY'] = history

    # Read in slope images and create a mean slope image as well as an uncertainty map
    mean_slope_image, unc_slope_image = create_mean_slope(input_dark_slope_files)

    # Read in the linearity files and create a mean zeroframe. We'll need this as
    # our offset value for hot pixels where the dark will be non-zero
    mean_zero_frame, zero_frame_unc = create_mean_zero_frame(input_dark_files)

    # Save the mean slope image
    mean_slope_filename = 'mean_slope_for_{}_{}_{}.fits'.format(metadata['INSTRUME'],
                                                                metadata['DETECTOR'],
                                                                metadata['SUBARRAY'])

    # Define the threshold for what constitues a hot/warm pixel
    warm_threshold = define_warm_threshold(mean_slope_image, threshold=sigma_threshold_for_hot)

    # Identify nominal pixels so we can zero them out. Also identify
    # warm/hot pixels, so that we can put those flags in the DQ array
    nominal = mean_slope_image < warm_threshold
    hot = mean_slope_image >= warm_threshold

    # Zero out the slopes of the nominal pixels
    orig_mean_slope_image = deepcopy(mean_slope_image)
    mean_slope_image[nominal] = 0.

    # Zero out the zero frame of the nominal pixels
    mean_zero_frame[nominal] = 0.

    # Save the original mean slope, as well as the updated mean
    # slope image with zeroed out pixels
    output_dir = os.path.dirname(output_file)
    mean_slope_filename = os.path.join(output_dir, mean_slope_filename)
    save_image([orig_mean_slope_image, mean_slope_image], mean_slope_filename)

    # Get a list of the exposure times associated with all groups
    tframe = 10.73677  # seconds
    frame_times = np.arange(total_groups) * tframe

    # Create empty RAPID ramp to hold the signals
    ramp = np.zeros((total_groups, 2048, 2048))
    unc_ramp = np.zeros((total_groups, 2048, 2048))
    for group, time in zip(range(total_groups), frame_times):
        ramp[group, :, :] = mean_slope_image * time + mean_zero_frame
        unc_ramp[group, :, :] = np.sqrt((unc_slope_image * time) + zero_frame_unc**2)

    # Create a DQ map where we note all the warm/hot pixels
    # Or should we leave the DQ map to be all zeros, assuming that
    # hot pixels are in the bad pixel mask reference file?
    # (Why in the world did we decide to do that in the working group??)
    ngrp, ny, nx = ramp.shape
    dq = np.zeros((ny, nx)).astype(np.int8)
    dq[hot] = dqflags.pixel['WARM']

    # Save everything in a reference file
    if output_file is None:
        output_file = '{}_{}_{}_dark_ref_file.fits'.format(metadata['INSTRUME'],
                                                           metadata['DETECTOR'],
                                                           metadata['SUBARRAY'])
    save_reffile(ramp, unc_ramp, dq, metadata, input_dark_files, output_file)


def collect_metadata(filename):
    """Read in the header of a fits file and collect metadata
    needed for the final reference file

    Parameters
    ----------
    filename : str
        Name of file to get metadata from

    Returns
    -------
    info : dict
        Dictionary of metadata
    """
    # List of header keywords to collect
    keywords = ['INSTRUME', 'DETECTOR', 'SUBARRAY', 'GROUPGAP', 'SUBSTRT1',
                'SUBSTRT2', 'SUBSIZE1', 'SUBSIZE2', 'FASTAXIS', 'SLOWAXIS']
    header = fits.getheader(filename)

    info = {}
    for key in keywords:
        info[key] = header[key]

    info['P_READPA'] = 'ANY'
    return info


def create_mean_slope(file_list):
    """Read in individual slope files and calculate a mean slope image

    Parameters
    ----------
    file_list : list
        List of slope image files

    Returns
    -------
    mean_image : numpy.ndarray
        2D mean image
    """
    img = fits.getdata(file_list[0])
    mean_stack = np.zeros((len(file_list), img.shape[0], img.shape[1]))
    mean_stack[0, :, :] = img

    for i, filename in enumerate(file_list[1:]):
        img = fits.getdata(filename)
        mean_stack[i+1, :, :] = img

    mean_image = np.median(mean_stack, axis=0)
    unc_mean_image = np.std(mean_stack, axis=0)

    # Set all reference pixels to zero
    mean_image[0:4, :] = 0.
    mean_image[2044:, :] = 0.
    mean_image[:, 0:4] = 0.
    mean_image[:, 2044:] = 0.
    unc_mean_image[0:4, :] = 0.
    unc_mean_image[2044:, :] = 0.
    unc_mean_image[:, 0:4] = 0.
    unc_mean_image[:, 2044:] = 0.

    return mean_image, unc_mean_image


def create_mean_zero_frame(file_list):
    """Create a mean zeroth frame image from a list of ramps.

    Parameters
    ----------
    file_list : list
        List of fits files

    Returns
    -------
    mean_zero : numpy.ndarray
        2D mean zeroth frame
    """
    zeros = np.zeros((len(file_list), 2048, 2048))
    for i, filename in enumerate(file_list):
        data = fits.getdata(filename)
        zeros[i, :, :] = data[0, 0, :, :]

    clipped = sigma_clip(zeros, sigma=3., cenfunc='median', axis=0, masked=False)
    mean_zero = np.nanmedian(clipped, axis=0)
    zero_unc = np.nanstd(clipped, axis=0)
    return mean_zero, zero_unc


def define_warm_threshold(image, threshold=5):
    """Pixels with signal rates above the threshold calculated here will
    be flagged as hot

    Parameters
    ----------
    image : numpy.ndarray
        2d image

    threshold : int
        Number of sigma above the mean to use for the threshold that
        defines warm/hot pixels. Pixels with signal below this level
        are considered nominal and will have signal values of 0 used
        in the dark. Pixels above this threshold will have dark signals
        based on their mean dark rate.

    Returns
    -------
    threshold : float
        Threshold value
    """
    clipped, lo, high = sigmaclip(image, low=3., high=3.)
    mean_val = np.mean(clipped)
    dev = np.std(clipped)
    threshold = mean_val + threshold * dev
    return threshold


def run_pipeline(filename, outdir='./'):
    """Run calwebb_detector1 on the input file

    Parameters
    ----------
    filename : str
        Name of input raw file

    Returns
    -------
    slope_filename : str
        Name of fits file with slope image
    """
    lin_file = filename.replace('_uncal', '_linearity')
    rate_file = filename.replace('_uncal', '_rate')

    if not os.path.isfile(lin_file) or not os.path.isfile(rate_file):
        cal = Detector1Pipeline()
        cal.dark_current.skip = True
        cal.persistence.skip = True

        cal.linearity.save_results = True

        cal.jump.rejection_threshold = 9
        cal.jump.maximum_cores = 'quarter'

        cal.save_results = True
        cal.output_dir = outdir
        cal.run(filename)


def save_image(images, filename):
    """Save the input image to a fits file

    Parameters
    ----------
    image : list
        List of numpy ndarrays

    filename : str
        Name of file to save ``image`` to
    """
    h0 = fits.PrimaryHDU(images[0])
    h1 = fits.ImageHDU(images[1])
    hlist = fits.HDUList([h0, h1])
    hlist.writeto(filename, overwrite=True)


def save_reffile(data, uncertainty, dq, metadata, file_list, outfile):
    """Save the given data in the official reference file format

    Parameters
    ----------
    data : numpy.ndarray
        3D array containing the mean dark current signals

    uncertainty : numpy.ndarray
        3D array containing associated uncertainties

    dq : numpy.ndarray
        2D array containing associated data quality flags

    metadata : dict
        Dictionary of metadata needed for the reference file
        header. Instrument, detector, subarray, etc as well
        as pedigree, author, descrip, and history

    file_list : list
        List of files used to construct ``data``

    outfile : str
        Name of file to save new reference file into
    """
    model = DarkModel()

    # Insert data
    model.data = data
    model.err = uncertainty
    model.dq = dq

    # Metadata copied from input files
    model.meta.reftype = 'DARK'
    model.meta.instrument.name = metadata['INSTRUME'].upper()
    model.meta.instrument.detector = metadata['DETECTOR'].upper()
    model.meta.exposure.type = 'NRC_DARK'
    model.meta.exposure.readpatt = 'RAPID'
    model.meta.exposure.nframes = 1
    model.meta.exposure.ngroups = data.shape[0]
    model.meta.exposure.groupgap = metadata['GROUPGAP']
    model.meta.subarray.name = metadata['SUBARRAY']
    model.meta.subarray.xstart = metadata['SUBSTRT1']
    model.meta.subarray.ystart = metadata['SUBSTRT2']
    model.meta.subarray.xsize = metadata['SUBSIZE1']
    model.meta.subarray.ysize = metadata['SUBSIZE2']
    model.meta.subarray.fastaxis = metadata['FASTAXIS']
    model.meta.subarray.slowaxis = metadata['SLOWAXIS']

    # Instrument-specific metadata
    model.meta.exposure.preadpatt = '{}|'.format(metadata['P_READPA'])

    # User-provided metadata
    model.meta.pedigree = metadata['PEDIGREE']
    model.meta.descrip = metadata['DESCRIP']
    model.meta.author = metadata['AUTHOR']
    model.meta.useafter = metadata['USE_AFTER']

    # Add HISTORY entries
    if metadata['HISTORY'] is not None:
        model.history.append(metadata['HISTORY'])
    model.history.append(("NIRCAM {} {} dark current file, created from data in files: "
                          .format(metadata['DETECTOR'].upper(), metadata['SUBARRAY'])))
    for filename in file_list:
        model.history.append(os.path.basename(filename))

    model.save(outfile)
    print('Dark current reference file saved to {}'.format(outfile))
