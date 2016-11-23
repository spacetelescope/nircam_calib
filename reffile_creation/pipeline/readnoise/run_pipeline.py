#! /usr/bin/env python

#generate a list of data to use, and run the appropriate pipeline steps
#on the files.

import glob
import numpy as np
#import sci2ssb
from jwst_pipeline.pipeline import SloperPipeline
from astropy.table import Table
import argparse
from itertools import izip
from datetime import datetime
from astropy.io import fits
from jwst_pipeline.flatfield import FlatFieldStep

datadir = 'ifs/jwst/wit/nircam/isim_cv3_files_for_calibrations/'
config_dict = {'A1':'','A2':'','A3':'','A4':'','ALONG':'','B1':'','B2':'','B3':'','B4':'','BLONG':''}

class Caldata():
    def __init__(self):
        self.verbose = False

    
    def get_fname_list(self,file):
        #generate a list of fits files given a file that follows the format
        #of Karl's list of good calibraiton files.
        allf = []

        #open the listfile
        tab = Table.read(file,header_start=0,data_start=1,format='ascii',delimiter=' ')

        #keep only unique entries
        allobs = tab['Obsid'].data
        alljob = tab['Jobid'].data.astype(int)

        combo = []
        for obs,job in izip(allobs,alljob):
            combo.append(str(job)+obs)

        combo = np.array(combo)
        unique,unique_indices,unique_counts = np.unique(combo,return_index=True,return_counts=True)
        
        uniqueobs = allobs[unique_indices]
        uniquejob = alljob[unique_indices]
            
        #search for files listed in the table
        #for obs,job in izip(tab['Obsid'],tab['Jobid']):
        i=0
        for obs,job in izip(uniqueobs,uniquejob):
            jobstr = str(job)
            obsstr = obs[0:-2]
            searchstr = maindir+'*/*_'+jobstr+'_*/*'+obsstr+'*fits'
            if job < 30665:
                searchstr = maindir+'026498_cpt/*/*_'+jobstr+'_*/*'+obsstr+'*fits'
            f = glob.glob(searchstr)
            if len(f) > 0:
                allf = allf + f
                expect = len(f) - unique_counts[i]
                print("Searching: "+searchstr+"  Found {} files. Expecting {} files. Difference {}.".format(len(f),unique_counts[i],expect))
            else:
                print("Searching: "+searchstr+"Unable to find any files for ObsID {} and JobID {}".format(obs,job))

            i = i + 1

        allf = list(set(allf))
        return allf
            


    def get_fname_from_jobid(self,jobids):
        #generate a list of fits files given a list of job IDs.

        allf = []
        for job in jobids:
            f = glob.glob(maindir+'*/*_'+str(job)+'_*/*.fits')
            if len(f) > 0:
                allf.append(f)
        allf = np.ravel(allf)

        return allf


    def read_listfile(self,file):
        #read in and return the contents of a single column list file
        col = []
        with open(file) as f:
            for line in f:
                if len(line) > 2:
                    col.append(line.strip())

        return col

    def sort_by_detector(self,files):
        #sort a list of files by detector and return
        #dict = {}
        dets = np.empty(len(files), dtype=str)
        #dets = np.zeros(len(files))

        for i,file in enumerate(files):
            h = fits.getheader(file)
            dets[i] = h['DETECTOR']
            #dict[file] = det

        indices = np.argsort(dets)
        files = np.array(files)
        files = files[indices]
        return files

    def check_time(self,files):
        #remove files taken too early
        trimmed = []
        for file in files:
            h = fits.open(file)
            det = h[0].header['DETECTOR'][3:]
            timeobs = h[0].header['TIME-OBS'][0:5]
            dateobs = h[0].header['DATE-OBS']
            exptime = datetime.strptime(dateobs + ' ' + timeobs, '%Y-%m-%d %H:%M')
            print('{}, {}, Karltime: {}.'.format(det,exptime,mindates_obj[det]))
            if exptime >= mindates_obj[det]:
                print('Keeping file: ',file)
                trimmed.append(file)
            else:
                print('Removing file: ',file)
        return trimmed
            
            
    def run(self):
        #read in list of job IDs
        #jobs = self.read_listfile(listfile)

        #get corresponding file names
        #files = self.get_fname_list(listfile)

        #print('all files:',len(files))

        #go through files and get execution times.
        #comapre to Karl's earliest exposure times for each detector.
        #throw out files that were taken too early.
        #files = self.check_time(files)
        #print('After checking time: ',len(files))

        #sort files by detector, so that we can move on to 
        #the next steps for one detector while the others are still 
        #processing
        #files = self.sort_by_detector(files)


        #files = ['NRCN815A-LIN-5365090406_7_483_SE_2015-12-31T09h50m07_uncal.fits']

        files = self.read_listfile(self.infile)


        #loop through files. Convert to ssb format, and then run the pipeline
        for file in files:
            #using detector, determine which calibration files to use. Our good cal files are not
            #in CRDS yet, for some reason.
            
            #then call each pipeline step yourself, so that you can provide the correct cal file, and also
            #save the outputs that you want.

            det = fits.getval(file,'DETECTOR')
            det = det[3:]
            if 'DARK' in file:
                cfile = 'badpix_pipeline_dark_'+det+'.cfg'
            else:
                cfile = 'badpix_pipeline_flat_'+det+'.cfg'

            #pipeline steps
            #dq_init - yes
            #saturation - yes
            #ipc - yes
            #superbias - yes
            #bias_drift - yes
            #linearity - yes
            #dark - skip for darks (probably can also skip for flats)
            #jump - yes, save output ramp
            #ramp fit - yes
            

            #flt = FlatFieldStep.call(cal_model,override_flat='')
            
            #slash = file.rfind('/')
            outfile = file[0:-5] + '_slopeint.fits'

            #run pipeline steps
            cal = SloperPipeline.call(file,config_file=cfile)
            cal.save(outfile)
            
            #extract the jump-derived cosmic-ray masks and convert to a form the
            #bad pixel generator recognizes?
            #or re-write the bad pixel generator to recognize this format. 
            #The latter makes more sense long-term
            

    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='calibrate cv3 data')
        parser.add_argument('infile',help='list file that contains a list of files to process.')
        parser.add_argument('--cfg_file',help='configuration file to use for pipeline.',default=None)
        return parser

if __name__ == '__main__':
    usagestring = 'USAGE: run_pipeline.py'
    
    cal = Caldata()
    parser = cal.add_options(usage=usagestring)
    args = parser.parse_args()
    cal.infile = args.infile
    cal.cfg_file = args.cfg_file

    cal.run()
