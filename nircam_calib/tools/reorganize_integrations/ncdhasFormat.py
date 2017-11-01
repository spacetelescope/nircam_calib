#!/usr/bin/env python



import argparse
from astropy.io import fits
import numpy as np
import sys



class ncdhasFormat:
    """Class to reformat files to match NCDHAS format. This script:

            1.) Changes orientation from SSB to FITSwriter format.
            2.) Changes data type to 'uint16' for NCDHAS pipeline.
            3.) Slices data to requested number of groups.
            4.) Pushes science data and headers into primary extension.

       Reformatted files can be calibrated with NCDHAS pipeline for
       comparison with SSB calibration pipeline. Reformatted file is saved with
       a different name. Changing data orientation is optional -- to just
       reformat headers and data extensions, leave default arguments.

       NOTE 1: To change the orientation, part of the 'native_to_ssb_orient.py'
       script was copied and added to this script (around line 78).

       NOTE 2: NCDHAS reference files seem to be in ADU and e-/sec depending on
       the file, need to be careful about units before comparing?


       2017-01-18 A. Canipe

    """



    def __init__(self):
        self.verbose = False



    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument("infile",help="File with format to change.")
        parser.add_argument("--max_groups",help="Number of groups to keep.",default=None,type=int)
        parser.add_argument("--outfile",help="Name of altered file.",default=None,type=str)
        parser.add_argument("--fits_orient",help="Option to change orientation to FITSwriter format.",default=False,type=bool)
        return parser



    # main function
    def run(self):


        # load the data that needs to be reformatted
        print('\nLoading data...\n')
        hdulist = fits.open(str(self.infile))
        indata = hdulist[1].data
        hdr = hdulist[0].header
        detector = hdr['DETECTOR']

        if detector == 'NRCA5':
            detector = 'NRCALONG'
            hdr['DETECTOR'] = detector
        if detector == 'NRCB5':
            detector = 'NRCBLONG'
            hdr['DETECTOR'] = detector


        # change data orientation by rotating all extensions except for DQ_DEF
        # and ASDF to get into FITSwriter format.
        if self.fits_orient == True:


            # find pixel to check for proper flipping at the end
            numypix = indata.shape[-2]
            numxpix = indata.shape[-1]
            y = np.int(numypix/3)
            x = np.int(numxpix/3)
            inpixel = indata[0,:,y,x]


            # start rotation for all extensions
            print('\nChanging orientation.')
            for i in range(1,len(hdulist)):
                print('Extension is:',i)
                hh=hdulist[i].header
                try:
                    typ = hh['EXTNAME']
                    if typ in ['SCI','ERR','PIXELDQ','GROUPDQ','DQ','COEFFS']:
                        data = hdulist[i].data
                        dim = len(data.shape)
                        if detector in ['NRCA2','NRCA4','NRCB1','NRCB3','NRCBLONG']:
                            #vertical flip
                            print("Vertical flip to get into FITSwriter format")
                            flipcoords = [2047-y,x]
                            if dim == 4:
                                flip = data[:,:,::-1]
                            if dim == 3:
                                flip = data[:,::-1]
                            if dim == 2:
                                flip = data[::-1]
                        elif detector in ['NRCA1','NRCA3','NRCB2','NRCB4','NRCALONG']:
                            #horizontal flip
                            print("Horizontal flip to get into FITSwriter format")
                            flipcoords = [y,2047-x]
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


        # check input data type
        data = hdulist[1].data
        print('\nOriginal data type:',data.dtype)
        print('Original data shape:',np.shape(data))


        # if option to change orientation is selected, check to make sure
        # flip worked properly using pixel from above by comparing coordinates
        # and pixel values before and after the flip
        if self.fits_orient == True:
            outpixel = data[0,:,flipcoords[0],flipcoords[1]]
            if np.array_equal(inpixel,outpixel) == True:
                print('\nFlip successful.')
                print('In coordinates [y,x]:',[y,x],'with values:\n',inpixel)
                print('Out coordinates [y,x]:',flipcoords,'with values:\n',outpixel)
            else:
                print('\nUh-oh. Problem with the flip. Quitting!')
                sys.exit(0)


        # make sure new data type is unsigned integer
        data = data.astype('uint16')
        print('\nNew data type:',data.dtype)


        # move SCI data into primary extension
        # slice to desired number of groups if requested
        if self.max_groups is not None:
            hdulist[0].data = data[0,:self.max_groups,:,:]
            print('New data shape:',np.shape(hdulist[0].data))
        else:
            hdulist[0].data = data[0,:,:,:]


        # move SCI extension headers into primary extension
        # (NCDHAS requires NAXIS keywords to be in primary extension)
        hdrSCI = hdulist[1].header
        for k,v in zip(hdrSCI.keys(), hdrSCI.values()):
            #print(k,v)
            hdulist[0].header[k]=v
        print('\nAll SCI data moved to primary extension.')


        # # delete data from all extensions to cut down on file size?
        # # (is this okay for rough comparison?)
        # for i in range(1,len(hdulist)):
        #     hdulist[i].data = None
        # print('Other extensions empty.')


        # write altered hdulist to new file with 'standard' name format
        if self.outfile == None:
            name, ext = self.infile.split('.')
            if self.max_groups is not None:
                outnamebase = name+'_ncdhasReformat_ngroups'+str(self.max_groups)
            elif self.max_groups is None:
                outnamebase = name+'_ncdhasReformat'

            if self.fits_orient == True:
                outnamebase = outnamebase+'_FITSorient'
            outname = outnamebase+'.'+ext


        # or save it with specified output filename
        else:
            outname = str(self.outfile)
        hdulist[0].writeto(outname,clobber=True)
        # hdulist.writeto(outname,clobber=True)
        print('\nFinished! Output written to:',outname,'\n')



if __name__ == '__main__':
    usagestring = 'USAGE: ncdhasFormat.py infile.fits --max_groups=10 --outfile=outfile.fits --fits_orient=True'

    reformat = ncdhasFormat()
    parser = reformat.add_options(usage=usagestring)
    args = parser.parse_args(namespace=reformat)

    reformat.run()
