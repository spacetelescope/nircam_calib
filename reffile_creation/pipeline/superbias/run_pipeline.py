#! /usr/bin/env python

#generate a list of data to use, and run the appropriate pipeline steps
#on the files.

import glob
import numpy as np
#import sci2ssb
from jwst_pipeline.pipeline import SloperPipeline
from multiprocessing import Pool
from astropy.table import Table
import argparse
from itertools import izip
from datetime import datetime
from astropy.io import fits
from jwst_pipeline.flatfield import FlatFieldStep
from jwst_pipeline.dq_init import DQInitStep
from jwst_pipeline.bias_drift import BiasDriftStep
from jwst_pipeline.ipc import IPCStep

datadir = 'ifs/jwst/wit/nircam/isim_cv3_files_for_calibrations/'
config_dict = {'A1':'','A2':'','A3':'','A4':'','ALONG':'','B1':'','B2':'','B3':'','B4':'','BLONG':''}

def read_listfile(file):
    #read in and return the contents of a single column list file
    col = []
    with open(file) as f:
        for line in f:
            if len(line) > 2:
                col.append(line.strip())
    return col
               
         
def setup(listfile,do_ipc,ipc_only):
    ipcdir = '/grp/jwst/wit/nircam/reference_files/SSB/CV2/delivery_Dec_2015/IPC/'
    ipcfiles = glob.glob(ipcdir+'*DMSorient.fits')

    #maskdir = '/grp/jwst/wit/nircam/cv3_reffile_conversion/bpm/'
    maskdir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/cv3_reffile_conversion/bpm/'
    maskfiles = glob.glob(maskdir+'*DMSorient.fits')

    satfiles = glob.glob('/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/cv3_reffile_conversion/welldepth/*ADU*DMSorient.fits')

    files = read_listfile(listfile)
    inputs = []
    for file in files:

        det = fits.getval(file,'DETECTOR')
        if det == 'NRCBLONG':
            det = 'NRCB5'
        if det == 'NRCALONG':
            det = 'NRCA5'
    
        maskfile = [s for s in maskfiles if det in s][0]
        if len(maskfile) > 1:
            print("More than one bad pixel mask found. Need a better list of possibilities.")
            stophere

        welldepth = [s for s in satfiles if detstr in s]
        if len(welldepth) > 1:
            print("More than one saturation map found. Need a better list of possibilities.")
            stophere

        ipckernel = None
        if do_ipc == True:
            cfile = 'proc_superbias_inputs.cfg'
            outfile = file[0:-5] + '_dq_init_bias_drift_nosidepix_ipc.fits'
            if ipc_only:
                outfile = file[0:-5] + '_ipc.fits'
            ipckernel = [s for s in ipcfiles if det in s][0]
        else:
            cfile = 'proc_superbias_inputs_withIPC.cfg'
            outfile = file[0:-5] + '_dq_init_bias_drift_nosidepix.fits'

        inputs.append((file,cfile,outfile,maskfile,ipckernel,ipc_only))

    return inputs
        
        
def run(tup):
        
    #pipeline steps
    #dq_init - yes
    #saturation - no
    #ipc - yes/no
    #superbias - no
    #bias_drift - yes
    #linearity - no
    #dark - no
    #jump - no
    #ramp fit - no
            
    file = tup[0]
    cfile = tup[1]
    outfile = tup[2] 
    maskfile = tup[3]
    ipckernel = tup[4]
    ipconly = tup[5]
    
    if ipconly == False:
        input = DQInitStep.call(file,override_mask=maskfile)
        input = BiasDriftStep.call(input,config_file='bias_drift.cfg')
    
        if ipckernel != None:
            input = IPCStep.call(input,override_ipc=ipckernel)

    else:
        input = IPCStep.call(file,override_ipc=ipckernel)

    
    input.save(outfile)

    #cal = SloperPipeline.call(file,config_file=cfile)
    #cal.save(outfile)
            
                    

def add_options(parser=None,usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage,description='calibrate cv3 data')
    parser.add_argument('infile',help='list file that contains a list of files to process.')
    parser.add_argument('--do_ipc',help='include IPC correction step',action='store_true',default=False)
    parser.add_argument('--ipc_only',help="run only the IPC step. designed for files that dqinit and biasdrift are already done.",action='store_true',default=False)
    return parser

if __name__ == '__main__':
    usagestring = 'USAGE: run_pipeline.py'
    
    parser = add_options(usage=usagestring)
    args = parser.parse_args()

    input_info = setup(args.infile,args.do_ipc,args.ipc_only)
    n_cores = 8
    pool = Pool(n_cores)
    pool.map(run,input_info)
