#! /usr/bin/env python

"""This module contains code to generate a skyflat from uncalibrated
NIRCam imaging data. It was developed to support CAP-10 analysis using 
data from PID1063, which will be taken during NIRCam commissioning.

Author
------
    - Ben Sunnquist, 2021
    - Brian Brooks, 2021

Use
---
    This module can be used from the command line as such:

        python make_skyflat.py

    To run on only certain files:

        python make_skyflat.py --detector nrcalong --filter f356w --pupil clear

    Or alternatively, from within python:

        make_skyflat.main(detector='nrcalong', fltr='f356w', pupil='clear')
"""

from glob import glob
import os

import argparse
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clip
from astropy.table import Table
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3 import Asn_Lv3Image
from jwst.datamodels import FlatModel, ImageModel
from jwst.outlier_detection.outlier_detection import gwcs_blot
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
from jwst.pipeline.calwebb_image2 import Image2Pipeline
from jwst.pipeline.calwebb_image3 import Image3Pipeline
import numpy as np
from photutils import detect_sources, detect_threshold
from stdatamodels import util

# -----------------------------------------------------------------------------

def blot_segmap(filename, files):
    """
    Blots the segmap values from the drizzled image onto the individual
    files than went into the drizzle, and write these out as new segmaps.

    Parameters
    ----------
    filename : str
        The drizzled file.

    files : list
        The files that went into creating the drizzled file.

    Outputs
    -------
    {file}_seg.fits
        A corresponding segmap for each file in files.
    """

    # Make an image model using the drizzled image WCS and the corresponding
    # segmentation maps data.
    drizzle_model = ImageModel(filename)
    drizzle_segmap = fits.getdata(filename.replace('.fits', '_seg.fits'))
    drizzle_segmap[drizzle_segmap != 0] = 1
    drizzle_model.data = drizzle_segmap

    # Blot the segmap data from the drizzled image onto each individual file,
    # and write out the corresponding segmap for each.
    for f in files:
        outfile = f.replace('.fits', '_seg.fits')
        model = ImageModel(f)
        blotted_data = gwcs_blot(drizzle_model, model, interp='nearest')
        fits.writeto(outfile, blotted_data, overwrite=True)

# -----------------------------------------------------------------------------

def create_file_table(files):
    """
    Creates a table of useful input file info that can be used for future
    selection steps.

    Parameters
    ----------
    files : list
        The list of input files.

    Returns
    -------
    t : astropy.table.table.Table
        The table containing the file info.
    """

    t = Table()
    t['filename'] = files
    t['detector'] = [fits.getheader(f)['DETECTOR'].lower() for f in files]
    t['filter'] = [fits.getheader(f)['FILTER'].lower() for f in files]
    t['pupil'] = [fits.getheader(f)['PUPIL'].lower() for f in files]
    t['target'] = [fits.getheader(f)['TARGNAME'].lower() for f in files]

    return t

# -----------------------------------------------------------------------------

def make_segmap(filename):
    """
    Generates a source segmentation map of the input file.
    
    Parameters
    ----------
    filename : str
        The file to make segmentation map for.
    
    Outputs
    -------
    {filename}_seg.fits
        The segmentation map.
    """

    # See if segmap already exists
    outfile = filename.replace('.fits', '_seg.fits')

    # Get the data
    data = fits.getdata(filename, 'SCI')
   
    # Detector sources; Make segmap
    threshold = detect_threshold(data, 3.0)
    sigma = 3.0 * gaussian_fwhm_to_sigma    # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold, npixels=3, filter_kernel=kernel)
    fits.writeto(outfile, segm.data, overwrite=True)

    return outfile

# -----------------------------------------------------------------------------

def make_skyflat(files, outfile=''):
    """
    Creates a skyflat from the stack of input images. For best results, each
    file should have a corresponding segmap in the same directory.

    The outputted skyflat is saved in CRDS format and can be used as a
    reference file in the flat_field step of the image2 JWST pipeline.

    Parameters
    ----------
    files : list
        A list of files to use to create the skyflat.

    outfile : str
        The filename of the outputted skyflat. Defaults to 
        {detector}_{filter}_{pupil}_skyflat.fits.

    Returns
    -------
    skyflat : numpy.ndarray
        The 2D skyflat image.

    skyflat_error : numpy.ndarray
        The 2D skyflat error image.

    skyflat_dq : numpy.ndarray
        The 2D skyflat data quality image.

    Outputs
    -------
    {detector}_{filter}_{pupil}_skyflat.fits
        The CRDS-formatted skyflat reference file.
    """

    # Get basic header info to use for the reference file creation
    h = fits.open(files[0])
    instrument = h[0].header['INSTRUME']
    detector = h[0].header['DETECTOR']
    fltr = h[0].header['FILTER']
    pupil = h[0].header['PUPIL']
    fastaxis = h[0].header['FASTAXIS']
    slowaxis = h[0].header['SLOWAXIS']
    substrt1 = h[0].header['SUBSTRT1']
    substrt2 = h[0].header['SUBSTRT2']
    h.close()
    if outfile == '':
        outfile = '{}_{}_{}_skyflat.fits'.format(detector.lower(), 
                                                 fltr.lower(), pupil.lower())

    # Create a stack with the source-masked image data for each file
    stack = np.zeros((len(files), 2048, 2048))
    for i,f in enumerate(files):
        data = fits.getdata(f, 'SCI')
        try:
            segmap = fits.getdata(f.replace('.fits', '_seg.fits'))
        except FileNotFoundError:
            print('WARNING: No segmap found for {}'.format(f))
            segmap = np.zeros(data.shape)
        data[segmap != 0] = np.nan
        stack[i] = data
    
    # Make the skyflat, along with its error and dq arrays
    skyflat = np.nanmedian(stack, axis=0)
    skyflat = skyflat / np.nanmedian(skyflat)
    skyflat_error = np.nanstd(stack, axis=0)
    skyflat_dq = np.zeros(skyflat.shape)

    # Fill missing pixel values
    n_empty = len(skyflat[~np.isfinite(skyflat)])
    print('{} empty flat pixels - setting these to 1.'.format(n_empty))
    skyflat_error[~np.isfinite(skyflat)] = 1
    skyflat_dq[~np.isfinite(skyflat)] = 5  # DO_NOT_USE + NO_FLAT_FIELD
    skyflat[~np.isfinite(skyflat)] = 1

    # Save the skyflat data in a CRDS-formmated reference file
    save_skyflat(skyflat, skyflat_error, skyflat_dq, instrument=instrument, 
                 detector=detector, fltr=fltr, pupil=pupil, fastaxis=fastaxis, 
                 slowaxis=slowaxis, substrt1=substrt1, substrt2=substrt2, 
                 filenames=files, outfile=outfile)
    print('Skyflat complete: {}'.format(outfile))

    return skyflat, skyflat_error, skyflat_dq

# -----------------------------------------------------------------------------

def run_pipeline(files, flat_field=True):
    """
    Runs all stages of the JWST Pipeline on uncalibrated imaging data.

    Parameters
    ----------
    files : list
        The files to run the pipeline on.

    flat_field : bool
        Option to run the flat field step in image2.

    Returns
    -------
    outname : str
        The name of the final i2d drizzled image.
    """

    # Create a name for image3 association and drizzled files
    detector = fits.getheader(files[0])['DETECTOR'].lower()
    fltr = fits.getheader(files[0])['FILTER'].lower()
    pupil = fits.getheader(files[0])['PUPIL'].lower()
    target = fits.getheader(files[0])['TARGNAME'].lower()
    name = '{}_{}_{}_{}'.format(detector, fltr, pupil, target)

    # Run detector1
    for f in files:
        m = Detector1Pipeline()
        #m.jump.override_readnoise = 'jwst_nircam_readnoise_0024_psub.fits'
        m.refpix.odd_even_rows = False  # only for MIRI
        m.ipc.skip = True
        m.persistence.skip = True
        m.save_results = True
        m.output_dir = os.getcwd()
        m.run(f)

    # Run image2
    files = [f.replace('_uncal.fits', '_rate.fits') for f in files]
    for f in files:
        m = Image2Pipeline()
        if not flat_field:
            m.flat_field.skip = True
        m.save_results = True
        m.output_dir = os.getcwd()
        m.run(f)

    # Run image3
    files = [f.replace('_rate.fits', '_cal.fits') for f in files]
    asn = asn_from_list(files, rule=Asn_Lv3Image, product_name=name)
    asn_file = '{}.json'.format(name)
    with open(asn_file, 'w') as f:
            f.write(asn.dump()[1])
    m = Image3Pipeline()
    m.save_results = True
    m.output_dir = os.getcwd()
    m.run(asn_file)

    return '{}_i2d.fits'.format(name)

# -----------------------------------------------------------------------------

def save_skyflat(skyflat, skyflat_error, skyflat_dq, instrument='', detector='', 
                 fltr='', pupil='', subarray='GENERIC', author='STScI', 
                 description='Pixel flat calibration file', pedigree='GROUND', 
                 useafter='2021-01-01T00:00:00', history='', fastaxis=-1, 
                 slowaxis=2, substrt1=1, substrt2=1, filenames=[], 
                 outfile='skyflat_reffile.fits'):
    """
    Saves skyflat data in a CRDS-formatted file that can be used
    as a reference file in the flat_field step of the image2 JWST pipeline.

    Parameters
    ----------
    skyflat : numpy.ndarray
        The 2D skyflat image.

    skyflat_error : numpy.ndarray
        The 2D skyflat error image.

    skyflat_dq : numpy.ndarray
        The 2D skyflat data quality image.

    instrument : str
        CRDS-required instrument for which to use this reference file for.

    detector : str
        CRDS-required detector for which to use this reference file for.

    fltr : str
        CRDS-required filter for which to use this reference file for.

    pupil : str
        CRDS-required pupil for which to use this reference file for.

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    author : str
        CRDS-required name of the reference file author, to be placed in the
        referece file header.

    description : str
        CRDS-required description of the reference file, to be placed in the
        reference file header.

    pedigree : str
        CRDS-required pedigree of the data used to create the reference file.

    useafter : str
        CRDS-required date of earliest data with which this referece file
        should be used. (e.g. '2019-04-01T00:00:00').

    history : str
        CRDS-required history section to place in the reference file header.

    fastaxis : int
        CRDS-required fastaxis of the reference file.

    slowaxis : int
        CRDS-required slowaxis of the reference file.

    substrt1 : int
        CRDS-required starting pixel in axis 1 direction.

    substrt2 : int
        CRDS-required starting pixel in axis 2 direction.

    filenames : list
        A list of files that went into the skyflat creation.

    outfile : str
        The filename of the outputted skyflat reference file. Defaults to 
        {detector}_{filter}_{pupil}_skyflat.fits.

    Outputs
    -------
    {outfile}.fits
        The CRDS-formatted skyflat reference file.
    """

    m = FlatModel()

    # Populate the data
    m.data = skyflat
    m.err = skyflat_error
    m.dq = skyflat_dq
    m.dq_def = [(0, 0, 'GOOD', ''),
                (0, 1, 'DO_NOT_USE', ''),
                (1, 2, 'UNRELIABLE_FLAT', ''),
                (2, 4, 'NO_FLAT_FIELD', '')]

    # Add CRDS-required keywords
    m.meta.instrument.name = instrument
    m.meta.instrument.detector = detector
    m.meta.instrument.filter = fltr
    m.meta.instrument.pupil = pupil
    m.meta.subarray.name = subarray
    m.meta.author = author
    m.meta.description = description
    m.meta.pedigree = pedigree
    m.meta.useafter = useafter
    m.meta.subarray.fastaxis = fastaxis
    m.meta.subarray.slowaxis = slowaxis
    m.meta.reftype = 'FLAT'
    yd, xd = skyflat.shape
    m.meta.subarray.xstart = substrt1
    m.meta.subarray.xsize = xd
    m.meta.subarray.ystart = substrt2
    m.meta.subarray.ysize = yd

    # Add the list of input files used to create the skyflat reference file
    m.history.append('DATA USED:')
    for f in filenames:
        f = os.path.basename(f)
        totlen = len(f)
        div = np.arange(0, totlen, 60)
        for val in div:
            if totlen > (val+60):
                m.history.append(util.create_history_entry(f[val:val+60]))
            else:
                m.history.append(util.create_history_entry(f[val:]))

    # Add any more history to the header
    if history != '':
        m.history.append(util.create_history_entry(history))

    # Save the flat reference file
    m.save(outfile, overwrite=True)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def main(detector='', fltr='', pupil='', flat_field=True):
    """
    The main processing function. See docstring for more details.

    Parameters
    ----------
    detector : str
        The detector of interest (nrcalong, nrcb1, etc).

    fltr : str
        The filter of interest (f356w, f070w, etc).

    pupil : str
        The pupil of interest (clear, maskbar, maskrnd, etc).

    flat_field : bool
        Option to run the flat field step in image2.
    """

    # Create a table of input files and select out relevent entries
    files = sorted(glob('*_uncal.fits'))
    t = create_file_table(files)
    if detector != '':
        t = t[t['detector']==detector.lower()]
    if detector != '':
        t = t[t['filter']==fltr.lower()]
    if detector != '':
        t = t[t['pupil']==pupil.lower()]
    print('Selected file table:')
    print(t)

    # Run the pipeline on the files belonging to each target, and create 
    # segmentation maps for all input files.
    for targ in np.unique(t['target']):
        files_to_drizzle = list(t['filename'][t['target']==targ])

        # Run all stages of the JWST pipeline
        print('Drizzling the following {} files: {}'.format(targ, files_to_drizzle))
        drizzle_filename = run_pipeline(files_to_drizzle, flat_field=flat_field)

        # Make a segmentation map of the drizzled image
        print('Creating segmap for {}...'.format(drizzle_filename))
        drizzle_segmap = make_segmap(drizzle_filename)

        # Blot the drizzled segmap back onto the individual cal files
        print('Blotting drizzled segmap back to the individual files...')
        cal_files = [f.replace('_uncal.fits', '_cal.fits') for f in files_to_drizzle]
        blot_segmap(drizzle_filename, cal_files)

    # Make the skyflat, using all available files regardless of target,
    # and save it as a CRDS-formatted reference file.
    print('Creating the skyflat...')
    files = [f.replace('_uncal.fits', '_cal.fits') for f in t['filename']]
    skyflat, skyflat_error, skyflat_dq = make_skyflat(files)
    print('make_skyflat.py complete')

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def parse_args():
    """
    Parses command line arguments.
    
    Returns
    -------
    args : object
        Contains the input arguments.
    """

    # Make the help strings
    detector_help = 'The detector of interest (nrcalong, nrcb1, etc).'
    filter_help = 'The filter of interest (f356w, f070w, etc).'
    pupil_help = 'The pupil of interest (clear, maskbar, maskrnd, etc).'
    no_flat_field_help = ('Option to omit the flat field step in image2 during'
                          'pipeline processing.')


    # Add the potential arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--detector', dest='detector', action='store', type=str, 
                        required=False, help=detector_help)
    parser.add_argument('--filter', dest='fltr', action='store', type=str, 
                        required=False, help=filter_help)
    parser.add_argument('--pupil', dest='pupil', action='store', type=str, 
                        required=False, help=pupil_help)
    parser.add_argument('--no_flat_field', dest='flat_field', action='store_false',
                        required=False, help=no_flat_field_help)
    
    # Set defaults
    parser.set_defaults(detector='', fltr='', pupil='', flat_field=True)

    # Get the arguments
    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the command line arguments
    args = parse_args()

    # Run the code
    main(detector=args.detector, fltr=args.fltr, pupil=args.pupil, 
         flat_field=args.flat_field)
