#! /usr/bin/env python

'''
Run the dq_init and bias_drift (without sidepix) 
pipeline steps on the darks that will go into
creating the superbias
'''

from multiprocessing import Pool
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
import glob,argparse,os,sys
from astropy.io import fits
from .refpixcorr_g0 import refpix_g0


#bpmfiles = {'A1':'','A2':'','A3':'','A4':'','ALONG':'','B1':'','B2':'','B3':'','B4':'','BLONG':''}


def get_filenames(file):
    files = []
    with open(file) as f:
        for line in f:
            if len(line) > 2:
                files.append(line.strip())
    return files


def run_cal(file):
    #get list of bad pixel masks that can be used
    #bpmlist = glob.glob('/grp/jwst/wit/nircam/cv3_reffile_conversion/bpm/*ssbspmask.fits')
    bpmlist = glob.glob('/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/bpm/*DMSorient.fits')
    satlist = glob.glob('/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/welldepth/*ADU*DMSorient.fits')

    #run the DQ init and bias_drift steps on the file
    detstr = fits.getval(file,'DETECTOR')
    if 'LONG' in detstr:
        detstr = detstr[0:4]+'5'

    #find the appropriate bad pixel mask and saturation file for the detector
    bpm = [s for s in bpmlist if detstr in s]
    if len(bpm) > 1:
        print("More than one bad pixel mask found. Need a better list of possibilities.")
        stophere

    welldepth = [s for s in satlist if detstr in s]
    if len(welldepth) > 1:
        print("More than one saturation map found. Need a better list of possibilities.")
        stophere

    #run the dq_init step
    dqfile = file[0:-5] + '_dq_init.fits'
    dqstep = DQInitStep.call(file,config_file='dq_init.cfg',override_mask=bpm[0],output_file=dqfile)

    #run saturation flagging
    satfile = dqfile[0:-5] + '_saturation.fits'
    satstep = SaturationStep.call(dqstep,config_file='saturation.cfg',override_saturation=welldepth[0],output_file=satfile)

    #run the reference pixel subtraction step using the group 0 subtraction and re-addition
    #which is what the Build 4 pipeline used to use.
    refcor = refpix_g0()
    refcor.infile = satfile
    refcor.outfile = None
    refcor.run()


def add_options(parser=None,usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage,description='calibrate inputs')
    parser.add_argument("infile",help="File containing list of files to calibrate.")
    return parser


if __name__ == '__main__':
    usagestring = 'python calibrate_superbias_inputs.py listfile.list'

    parser = add_options(usage=usagestring)
    args = parser.parse_args()

    allfiles = get_filenames(args.infile)

    n_cores = 4
    pool = Pool(n_cores)
    pool.map(run_cal,allfiles)

