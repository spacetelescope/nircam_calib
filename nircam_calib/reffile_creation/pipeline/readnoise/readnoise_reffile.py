#! /usr/bin/env python

from astropy.io import fits
import numpy as np
import argparse
from readnoise import Readnoise
from jwst.datamodels import ReadnoiseModel

class ReadnoiseReffile(Readnoise):
    def __init__(self):
        Readnoise.__init__(self)

        
    def create(self):
        '''run the readnoise calculator'''

        if self.outfile is None:
            self.final_outfile = self.outfile[0:-5] + '_final_Readnoise_reffile.fits'
        
        #get detector name and determine if there is a low-epoxy region
        detector = fits.getval(self.infile,'DETECTOR')

        if detector not in ['NRCA1','NRCA3','NRCA4','NRCALONG','NRCB1','NRCB4']:
            outside_file, outside_ron  = self.make_rdnoise()

            #read in via ReadnoiseModel
            ron_map = ReadnoiseModel(outside_file)

        else: #we need 2 runs of readnoise.py if there is an epoxy void
            final_outfile = self.outfile
            self.outfile = 'tmp_outside_epoxy_void_ronmap.fits'
            outside_file, outside_ron  = self.make_rdnoise()

            #read in via ReadnoiseModel
            ron_map = ReadnoiseModel(outside_file)
            
            self.invert_epoxy = True
            self.outfile = 'tmp_inside_epoxy_void_ronmap.fits'
            inside_file, inside_ron = self.make_rdnoise()

            #sum the two maps
            ron_map.data += inside_ron

            #save the total file
            self.outfile = final_outfile
            ron_map.save(self.outfile)
            
        return self.outfile,ron_map


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='calculate readnoise')
        parser.add_argument("infile",help="Name of file to process")
        parser.add_argument('-v','--verbose', default=False, action="store_true",help='Verbose output')
        parser.add_argument('-o','--outfile', default=None,help='file name of output file')
        parser.add_argument('--bpm'  , default=None ,help='file name of bad pixel mask')
        parser.add_argument('--bpmval'  , default=0xffff ,type=int,help='pixels with bpm&bpmval>0 are masked.')
        parser.add_argument('--xmin'  , default=None , type=int,help='xmin for stats.')
        parser.add_argument('--xmax'  , default=None , type=int,help='xmax for stats.')
        parser.add_argument('--ymin'  , default=None , type=int,help='ymin for stats.')
        parser.add_argument('--ymax'  , default=None , type=int,help='ymax for stats.')
        parser.add_argument('--forcexylimits', default=False, action="store_true",help='Do not use any pixels outside the xylimits, even if possible')
        parser.add_argument('--boxsizex'  , default=128 , type=int,help='boxsize for stats in x direction.')
        parser.add_argument('--boxsizey'  , default=128 , type=int,help='boxsize for stats in y direction.')
        parser.add_argument('--stepsizex'  , default=None , type=int,help='stepsize between positions in x direction.')
        parser.add_argument('--stepsizey'  , default=None , type=int,help='stepsize between positions in y direction.')
        parser.add_argument('--gmin'  , default=None , type=int,help='minimum group used. If None, then group 0.')
        parser.add_argument('--gmax'  , default=None , type=int,help='maximum group used. If None, then last group.')
        parser.add_argument('-f','--fill', default=False, action="store_true",help='if stepsize>1: fill out the output matrix')
        parser.add_argument('--Pclipmax'  , default=10.0 , type=float,help='maximum % (0-100) of pixels clipped for statistic.')
        parser.add_argument('--Npixmin'  , default=32 , type=int,help='minimum # of valid pixels in box required.')
        parser.add_argument('--invert_epoxy', default = False, action='store_true',help="Invert the epoxy void mask, in order to calculate readnoise values inside the void region.")
        parser.add_argument('--fill_iterations', default = 4, type=int, help="Number of times to run the filling function that calculates readnoise values for boxes with too few pixels by calculating the mean of surrounding good boxes.")
        parser.add_argument('--useafter',help='Use after date to put into readnoise reference file',default='')
        parser.add_argument('--description',help='String to place into the DESCRIP header keyword of output file',default='NIRCam Readnoise file')
        parser.add_argument('--author',help='Name of person creating the file',default='NIRCam Instrument Team')
        parser.add_argument('--pedigree',help='Source of the data used to make reffile',default='GROUND')
        return parser


if __name__ == '__main__':

    usagestring = 'readnoise_reffile.py file.fits'

    rnfile =  ReadnoiseReffile()
    parser = rnfile.add_options(usage=usagestring)
    args = parser.parse_args(namespace=rnfile)

    rnfile.create()
