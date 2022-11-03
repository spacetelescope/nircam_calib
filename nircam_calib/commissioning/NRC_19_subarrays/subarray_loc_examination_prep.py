#! /usr/bin/env python

"""Given a two input files, make a quick check that the SUBSTRT and SUBSIZE keyword
values are correct, then crop the larger array down to contain only the pixels within the
smaller array. Scale the values in the larger array for the shorter integration time of
the shorter array, and save the result. The user can then blink these images as a
way of checking if the subarray is situated in the proper location.

This visual examination is the only way I can think of to confirm that the proper
pixels are being used, if you don't want to assume that the SUBSTRT and SUBSIZE
are correct. i.e. if SUBSTRT values say that the subarray starts at x,y = (1, 1)
on the detector, but the pixels actually read out started at (2, 2).

"""
from astropy.io import fits
from jwst import datamodels
import numpy as np
import os
import pysiaf


def check_location(fullframe_file, subarray_file):
    """Main function
    """
    # Open files with datamodels
    obj1 = open_file(fullframe_file)
    obj2 = open_file(subarray_file)

    # Check basic metadata
    metadata_check(obj1, obj2)

    # Get the location of the subarray on the detector from the header
    #ob1_loc_from_header = aperture_location(obj1)
    ob2_loc_from_header = aperture_location(obj2)

    # Crop the dq arrays in order to contain only overlapping pixels
    cropped_ff = match_arrays(obj1.data, obj2.data, ob2_loc_from_header)

    # Save the final group of each ramp in a fits file, so they can be blinked
    h0 = fits.PrimaryHDU(obj2.data[0, -1, :, :])
    h1 = fits.ImageHDU(cropped_ff[0, -1, :, :])
    hdulist = fits.HDUList([h0, h1])
    ff_base = os.path.basename(fullframe_file).replace('.fits', '')
    sub_base = os.path.basename(subarray_file).replace('.fits', '')
    outname = 'location_comparison_full_frame_vs_{}_{}_{}.fits'.format(obj2.meta.subarray.name.upper(), sub_base, ff_base)
    hdulist.writeto(outname, overwrite=True)
    print('Subarray group and cropped full frame saved to {} for visual inspection.'.format(outname))


def aperture_location(mod):
    """Get the location of the aperture in full frame coords
    """
    xstart = mod.meta.subarray.xstart - 1
    ystart = mod.meta.subarray.ystart - 1
    xend = xstart + mod.meta.subarray.xsize
    yend = ystart + mod.meta.subarray.ysize
    return xstart, ystart, xend, yend


def get_siaf_aperture_name(mod):
    """Translate the SUBARRAY name from the header to the matching
    SIAF aperture name
    """
    sub_name = mod.meta.subarray.name.upper()
    detector = mod.meta.instrument.detector.upper()
    filtername = mod.meta.instrument.filter.upper()
    if 'GRISM256' not in sub_name:
        if sub_name == 'SUB32TATS':
            sub_name = 'TAPSIMG32'
        elif sub_name == 'SUB32TATSGRISM':
            sub_name = 'TAGRISMTS_SCI_F322W2'
        siaf_aperture = '{}_{}'.format(detector, sub_name)
    else:
        if detector in ['NRCALONG', 'NRCA5']:
            siaf_aperture = '{}_GRISM256_{}'.format(detector, filtername)
        else:
            siaf_aperture = '{}_GRISMTS256'.format(detector)

    if 'LONG' in siaf_aperture:
        siaf_aperture = siaf_aperture.replace('LONG', '5')

    return siaf_aperture


def match_arrays(arr1, arr2, loc2):
    xstart, ystart, xend, yend = loc2
    arr1 = arr1[:, :, ystart:yend, xstart:xend]
    return arr1


def metadata_check(mod1, mod2):
    """Be sure that the detector used is the same
    """
    if mod1.meta.instrument.detector != mod2.meta.instrument.detector:
        raise ValueError("Files are from different detectors! Cannot compare.")

    if mod1.data.shape[2] != 2048:
        raise ValueError("First input file does not appear to be a full frame observation.")

    # Check the expected corner values of the subarray from SIAF
    siaf_aperture = get_siaf_aperture_name(mod2)
    x_siaf, y_siaf = sci_subarray_corners('nircam', siaf_aperture)

    if x_siaf[0]+1 != mod2.meta.subarray.xstart:
        print(('\n\nWARNING: SUBSTRT1 header keyword value ({}) does not match the corner index from SIAF ({})')
               .format(mod2.meta.subarray.xstart, x_siaf[0]+1))
    if y_siaf[0]+1 != mod2.meta.subarray.ystart:
        print(('\n\nWARNING: SUBSTRT2 header keyword value ({}) does not match the corner index from SIAF ({})')
               .format(mod2.meta.subarray.ystart, y_siaf[0]+1))


def open_file(filename):
    """Open with datamodels and check that we have ramp data
    """
    model = datamodels.open(filename)
    if not (isinstance(model, datamodels.RampModel) or isinstance(model, datamodels.Level1bModel)):
        raise ValueError("{} read in as {}. Expecting RampModel. Unable to continue.".format(filename, type(model)))

    return model


def sci_subarray_corners(instrument, aperture_name, siaf=None, verbose=False):
    """Return the two opposing aperture corners in the SIAF Science frame of the full-frame SCA.

    Parameters
    ----------
    instrument : str
        JWST instrument name with correct capitalization
    aperture_name : str
        SIAF aperture name
    siaf : pysiaf.Siaf
        SIAF instance for a single instrument
    verbose : bool
        Verbose output on/off
    Returns
    -------
    x_sci, y_sci : tuple of numpy arrays
        Subarray corner coordinates
    """
    # get SIAF
    if siaf is None:
        siaf = siaf = pysiaf.Siaf(instrument)

    # get master aperture names
    siaf_detector_layout = pysiaf.iando.read.read_siaf_detector_layout()
    master_aperture_names = siaf_detector_layout['AperName'].data

    # read pysiaf aperture definition file containing DetRef and SciRef values
    siaf_aperture_definitions = pysiaf.iando.read.read_siaf_aperture_definitions(instrument)

    # aperture object
    aperture = siaf[aperture_name]

    # aperture corners in SIAF detector coordinates
    x_det, y_det = aperture.corners('det', rederive=True)

    # determine parent aperture, i.e. the underlying full frame SCA aperture
    index = siaf_aperture_definitions['AperName'].tolist().index(aperture_name)
    aperture._parent_apertures = siaf_aperture_definitions['parent_apertures'][index]

    # If multiuple apertures are listed as parents keep only the first
    if ';' in aperture._parent_apertures:
        print('Multiple parent apertures: {}'.format(aperture._parent_apertures))
        aperture._parent_apertures = aperture._parent_apertures.split(';')[0]

    if aperture_name in master_aperture_names:
        # if master aperture, use it directly to transform to science frame
        x_sci, y_sci = aperture.det_to_sci(x_det, y_det)
    elif aperture._parent_apertures is not None:
        # use parent aperture for transformation
        if verbose:
            logger.info('Using parent {} for {}'.format(aperture._parent_apertures, aperture_name))
        x_sci, y_sci = siaf[aperture._parent_apertures].det_to_sci(x_det, y_det)
        aperture = siaf[aperture._parent_apertures]

    if instrument.lower() == 'nircam':
        if aperture.DetSciParity == 1:
            corner_index = np.array([1, 3])
        elif aperture.DetSciParity == -1:
            # NIRCam will always fall in here, except in case of non-dms orientation
            corner_index = np.array([0, 2])
        x_corner = x_sci[corner_index]
        y_corner = y_sci[corner_index]
    elif instrument.lower() == 'niriss':
        x_corner_index = np.array([0, 2])
        y_corner_index = np.array([0, 2])
        if aperture_name == 'NIS_CEN_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([3, 1])
        x_corner = x_sci[x_corner_index]
        y_corner = y_sci[y_corner_index]
        if aperture_name in ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256']:
            x_corner = [1, 2048]
            if aperture_name == 'NIS_SUBSTRIP96':
                y_corner = [1803, 1898]
            if aperture_name == 'NIS_SUBSTRIP256':
                y_corner = [1793, 2048]
    elif instrument.lower() == 'fgs':
        x_corner_index = np.array([0, 2])
        y_corner_index = np.array([0, 2])
        if aperture_name == 'FGS1_FULL_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([3, 1])
        if aperture_name == 'FGS2_FULL_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([1, 3])
        x_corner = x_sci[x_corner_index]
        y_corner = y_sci[y_corner_index]
    else:
        raise NotImplementedError(("Instrument {} not supported for SIAF subarray corners"
                                   .format(instrument)))

    # account for python conventions (e.g. 0-based indexing)
    # we also want integer values as these will be indexes
    x_corner = np.array([np.ceil(x_corner[0]) - 1, np.floor(x_corner[1]) - 1])
    y_corner = np.array([np.ceil(y_corner[0]) - 1, np.floor(y_corner[1]) - 1])
    return x_corner.astype(np.int), y_corner.astype(np.int)

