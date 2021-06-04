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

def check_location(fullframe_file, subaray_file):
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
    h1 = fits.ImageHDU(cropped_ff)
    hdulist = fits.HDUList([h0, h1])
    ff_base = fullframe_file.replace('.fits', '')
    sub_base = subarray_file.replace('.fits', '')
    outname = 'location_comparison_{}_{}.fits'.format(sub_base, ff_base)
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

def match_arrays(arr1, arr2, loc2):
    xstart, ystart, xend, yend = loc2
    arr1 = arr1[ystart:yend, xstart:xend]
    return arr1

def metadata_check(mod1, mod2):
    """Be sure that the detector used is the same
    """
    if mod1.meta.instrument.detector != mod2.meta.instrument.detector:
        raise ValueError("Files are from different detectors! Cannot compare.")

    if mod1.sci.shape[2] != 2048:
        raise ValueError("First input file does not appear to be a full frame observation.")


def open_file(filename):
    """Open with datamodels and check that we have ramp data
    """
    model = datamodels.open(filename)
    if not isinstance(model, datamodels.RampModel):
        raise ValueError("{} read in as {}. Expecting RampModel. Unable to continue.".format(filename, type(model)))

    return model




