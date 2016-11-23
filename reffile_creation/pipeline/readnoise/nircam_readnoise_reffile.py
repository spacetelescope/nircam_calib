#! /usr/bin/env python

import sys, os,re
from astropy.io import fits
import numpy as np
import argparse,subprocess
import subprocess
import datetime

# put the tools directory into the path
# add pythonmodules to PATH
if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
else:
    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
    sys.exit(0)

from multiprocessing import Pool
from jwst_lib.models import ReadnoiseModel
import mkimrdnoise_ssboutput_epoxy_faster as ron_calc


#DO NOT MAKE THIS A CLASS IF YOU WANT TO RUN IT IN PARALLEL!


def read_listfile(file):
    '''read in a list of files'''
    files = []
    with open(file) as f:
        for line in f:
            if len(line) > 2:
                files.append(line.strip())
    return files

def setup(args):
    '''create the input tuples that will be fed into the mutliprocesing calls'''

    #read in the list of files to run
    files = read_listfile(args.listfile)

    #loop over filenames and create the list of tuples
    inputs = []
    for file in files:
        inputs.append((file,args))

    return inputs


def run(in_tuple):
    '''run the readnoise calculator'''
    #in_tuple is (filenames,args)

    infile = in_tuple[0]
    args = in_tuple[1]

    ron = ron_calc.mkimrdnoiseclass()
    ron.verbose = args.verbose
    ron.debug = args.debug
    ron.outfile = args.outfile
    ron.bpm = args.bpm
    ron.bpmval = args.bpmval
    ron.xmin = args.xmin
    ron.xmax = args.xmax
    ron.ymin = args.ymin
    ron.ymax = args.ymax
    ron.forcexylimits = args.forcexylimits
    ron.boxsizex = args.boxsizex
    ron.boxsizey = args.boxsizey
    ron.stepsizex = args.stepsizex
    ron.stepsizey = args.stepsizey
    ron.gmin = args.gmin
    ron.gmax = args.gmax
    ron.fill = args.fill
    ron.Pclipmax = args.Pclipmax
    ron.Npixmin = args.Npixmin
    ron.invert_epoxy = False
    ron.fill_iterations = args.fill_iterations

    outfilebasename = args.outfilebasename + infile[0:-5] + '_outside_void'
    ron.mkimrdnoise(infile,outfilebasename)
    outsidevoid_filename = outfilebasename + '_ssbreadnoise.fits'

    outfilebasename = args.outfilebasename + infile[0:-5] + '_inside_void'
    ron.invert_epoxy = True
    ron.mkimrdnoise(infile,outfilebasename)
    insidevoid_filename = outfilebasename + '_ssbreadnoise.fits'

    #read in the results and combine
    outvoid = ReadnoiseModel(outsidevoid_filename)
    invoid = ReadnoiseModel(insidevoid_filename)

    outdata = outvoid.data
    indata = invoid.data

    #combine the two maps
    ronframe = outdata + indata
    
    #save to a new file using one of the inputs as a base
    outvoid.data = ronframe
    finalfile = args.outfilebasename + infile[0:-5] + '_final_ssbreadnoise.fits'
    outvoid.save(finalfile)

    #fits format checks
    check_ssb = fits.open(finalfile)
    print(check_ssb.info())
    print(check_ssb['SCI'].data[500,500])
    print(ronframe[500,500])

    #redcat team checks
    #subprocess.call(['fitsverify',finalfile])



def add_options(parser=None,usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage,description='calculate readnoise')
    parser.add_argument("listfile",help="Name of file containing the list of files to calculate readnoise for.")
    parser.add_argument("outfilebasename",help="Name of output file", default = 'NRC_readnoise.fits')
    parser.add_argument('-v','--verbose', default=False, action="store_true",help='Verbose output')
    parser.add_argument('-d','--debug', default=False, action="store_true",help='Debugging output')
    parser.add_argument('-o','--outfile', default='rdnoise.fits',help='file name of output file')
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
    parser.add_argument('--Npixmin'  , default=32 , type=int,help='minimum # of valid pixels in box required. If None, then .')
    parser.add_argument('--invert_epoxy', default = False, action='store_true',help="Invert the epoxy void mask, in order to calculate readnoise values inside the void region.")
    parser.add_argument('--fill_iterations', default = 4, type=int, help="Number of times to run the filling function that calculates readnoise values for boxes with too few pixels by calculating the mean of surrounding good boxes.")
    return parser



if __name__ == '__main__':

    starttime = datetime.datetime.now()
    print("Starting: {}".format(starttime))
    usagestring = 'nircam_readnoise_reffile.py files.list'

    parser = add_options(usage=usagestring)
    args = parser.parse_args()

    input_info = setup(args)
    print(input_info)
    n_cores = 8
    pool = Pool(n_cores)
    pool.map(run,input_info)
    endtime = datetime.datetime.now()
    print("Finished: {}".format(endtime))
    print("Duration: {}".format(endtime-starttime))
