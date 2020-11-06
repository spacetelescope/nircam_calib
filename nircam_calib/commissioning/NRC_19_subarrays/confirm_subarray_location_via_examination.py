#! /usr/bin/env python

"""Confirm that the subarray aperture is in the correct location on the detector
using bias structure/bad pixels. This module contains some convenience functions
for comparing a subarray exposure to the expected corresponding pixels from
a full frame exposure. In the end, manual inspection of these two arrays
will be used to confirm subarray position
"""
import os

from astropy.io import fits


def compare(full_frame_file=None, subarray_file=None, output_dir='./'):
    """Main function for comparison

    Parameters
    ----------
    full_frame_file : str
        Name of fits file with full frame observation

    subarray_file : str
        Name of fits file with subarray observation

    output_dir : str
        Directory into which output files are saved

    Returns
    -------
    subarray : numpy.ndarray
        2D subarray image

    cropped : numpy.ndarray
        2D subarray cropped from full frame image

    subarray_loc : dict
        Information on the subarray location within the detector
    """
    # Get subarray location info from the subarray file
    subarray_loc, subarray = get_subarray_info(subarray_file)

    # Extract a subarray from the full frame file that
    # corresponds to the subarray location from the subarray
    # file
    cropped = extract_subarray(full_frame_file, subarray_loc)

    # Save the cropped image and subarray image together in 2
    # extensions of a fits file, to allow for later comparison
    comparison_file = save_comparison_file(full_frame_file, subarray_file, cropped, subarray, subarray_loc, out_dir=output_dir)

    # Return the subarray image and cropped image so they can
    # be displayed if desired (e.g. in a notbeook)
    return subarray, cropped, subarray_loc, comparison_file


def extract_subarray(filename, location):
    """Extract a subarray from a full frame image based on
    location data

    Parameters
    ----------
    filename : str
        Name of fits file with full frame image

    location : dict
        Subarray location information

    Returns
    -------
    cropped_array : numpy.ndarray
        2D array extracted from the full frame image
        (from the 0th group if the full frame exposure is
        3D or 4D)
    """
    # Open the full frame file
    with fits.open(filename) as fobj:
        data = fobj[1].data
        header0 = fobj[0].header

    # Quick consistency checks
    if header0['DETECTOR'] != location['DETECTOR']:
        raise ValueError('Subarray and full frame files have different detectors! Cannot continue.')

    # Coordinates in the headers are indexed to 1. Modify for python's
    # expectation that values are indexed to 0.
    xstart = location['SUBSTRT1'] - 1
    xend = xstart + location['SUBSIZE1']
    ystart = location['SUBSTRT2'] - 1
    yend = ystart + location['SUBSIZE2']

    # Extract the given subarray from the full frame data
    ndim = len(data.shape)
    if ndim == 4:
        cropped_array = data[0, 0, ystart:yend, xstart:xend]
    elif ndim == 3:
        cropped_array = data[0, ystart:yend, xstart:xend]
    elif ndim == 2:
        cropped_array = data[ystart:yend, xstart:xend]

    return cropped_array


def get_subarray_info(filename):
    """Retrieve information on the subarray location within the full
    frame from the header of the given file. Return this information
    in a dictionary. Also return the subarray image itself

    Parameters
    ----------
    filename : str
        Name of fits file containing the subarray exposure

    Returns
    -------
    location : dict
        Dictionary of subarray information

    image : numpy.ndarray
        2D array containing the subarray image. Will be the 0th
        group if the input file contains a 3D or 4D array
    """
    with fits.open(filename) as fobj:
        data = fobj[1].data
        header0 = fobj[0].header

    location = {}
    location['DETECTOR'] = header0['DETECTOR']
    location['SUBARRAY'] = header0['SUBARRAY']
    location['SUBSTRT1'] = header0['SUBSTRT1']
    location['SUBSTRT2'] = header0['SUBSTRT2']
    location['SUBSIZE1'] = header0['SUBSIZE1']
    location['SUBSIZE2'] = header0['SUBSIZE2']

    ndim = len(data.shape)
    if ndim == 4:
        image = data[0, 0, : ,:]
    elif ndim == 3:
        image = data[0, :, :]
    elif ndim == 2:
        image = data

    return location, image


def save_comparison_file(full_file, sub_file, cropped_from_full, sub_img, loc_data, out_dir='./'):
    """Save the subarray image and the extracted subarray from the full frame image
    in a fits file

    Parameters
    ----------
    full_file : str
        Filename of fits file with full frame data

    sub_file : str
        Filename of fits file with subarray data

    cropped_from_full : numpy.ndarray
        2D subarray image that was cropped from the full frame image

    sub_img : numpy.ndarray
        2D subarray image

    loc_data : dict
        Dictionary of subarray location information

    out_dir : str
        Directory in which the file is saved

    Returns
    -------
    out_filename : str
        Name of fits file containing cropped_from_full and sub_img
    """
    # Generate the filename of the output file
    sub_base = os.path.basename(sub_file).replace('.fits', '')
    full_dir, full_base = os.path.split(full_file)
    full_base = full_base.replace('.fits', '')
    comp_base = 'sub_loc_check_{}_{}.fits'.format(sub_base, full_base)
    out_filename = os.path.join(out_dir, comp_base)

    # Create and populate the fits object
    h0 = fits.PrimaryHDU(cropped_from_full)
    h1 = fits.ImageHDU(sub_img)
    for key in loc_data:
        h0.header[key] = loc_data[key]
    hlist = fits.HDUList([h0, h1])

    # Save the file
    hlist.writeto(out_filename, overwrite=True)
    return out_filename
