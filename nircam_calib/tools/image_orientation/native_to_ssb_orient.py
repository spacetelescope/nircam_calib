#! /usr/bin/env python

"""
Re-orient a NIRCam SSB-format data file from the native orientation to the
orientation that comes through the DMS pipeline. This script is intended for 
data and reference files that have already been delivered or converted from
RAW format, but where we forgot to correctly orient the data and add
the FASTAXIS and SLOWAXIS keywords.
"""

from astropy.io import fits
import argparse
import sys

def native_to_science_image_flip(file):
    #flip the data to match the orientation expected by SSB
    
    #open file
    hdulist = fits.open(file)
    head0 = hdulist[0].header
    detector = head0['DETECTOR']

    if detector == 'NRCA5':
        detector = 'NRCALONG'
        head0['DETECTOR'] = detector
    if detector == 'NRCB5':
        detector = 'NRCBLONG'
        head0['DETECTOR'] = detector

    #first, add the correct values of fastaxis and slowaxis to the 
    #main header
    if detector in ['NRCA2','NRCA4','NRCB1','NRCB3','NRCBLONG']:
        head0['FASTAXIS'] = 1
        head0['SLOWAXIS'] = -2

    elif detector in ['NRCA1','NRCA3','NRCB2','NRCB4','NRCALONG']:
        head0['FASTAXIS'] = -1
        head0['SLOWAXIS'] = 2
    else:
        print("WARNING! I don't recognize {} as a valid detector!".format(detector))
        sys.exit(0)


    #We'll need to rotate all extensions except for
    #DQ_DEF and ASDF (for now).
    for i in xrange(1,len(hdulist)):
        hh=hdulist[i].header
        try:
            typ = hh['EXTNAME']
            if typ in ['SCI','ERR','PIXELDQ','GROUPDQ','DQ','COEFFS']:
                data = hdulist[i].data
                dim = len(data.shape)
                if detector in ['NRCA2','NRCA4','NRCB1','NRCB3','NRCBLONG']:
                    #vertical flip
                    #print("Vertical flip to get into DMS format")
                    if dim == 4:
                        flip = data[:,:,::-1]
                    if dim == 3:
                        flip = data[:,::-1]
                    if dim == 2:
                        flip = data[::-1]
                elif detector in ['NRCA1','NRCA3','NRCB2','NRCB4','NRCALONG']:
                    #horizontal flip
                    #print("Horizontal flip to get into DMS format")
                    if dim == 4:
                        flip = data[:,:,:,::-1]
                    if dim == 3:
                        flip = data[:,:,::-1]
                    if dim == 2:
                        flip = data[:,::-1]
                hdulist[i].data = flip
        except KeyError:
            print("Uh-oh. {} and {} are giving a problem.".format(detector,typ))
            sys.exit(0)

    #save the rotated and keyword-updated file with DMSorient appended to the input filename
    outname = file[0:-5] + '_DMSorient.fits'
    hdulist.writeto(outname,clobber=True)


def add_options(parser=None,usage=None):
    if parser==None:
        parser = argparse.ArgumentParser(usage=usage,description='For a file that is otherwise in SSB-format, rotate native (fitswriter) orientation data to the DMS orientation expected by SSB. Add appropriate fastaxis keyword.')
    parser.add_argument("file",help="Name of file to rotate")
    return parser



if __name__ == '__main__':
    usagestring = 'python native_to_ssb_orient.py file.fits'

    parser = add_options(usage=usagestring)
    args = parser.parse_args()
    
    native_to_science_image_flip(args.file)
