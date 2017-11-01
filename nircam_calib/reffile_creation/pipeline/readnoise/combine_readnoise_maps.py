#! /usr/bin/env python

'''Combine a series of individual readnoise maps into a mean map
Assume the individual maps came our of Armin's mkimrdnoise.py script,
in which case the individual maps are for readnoise in a single frame,
in ADU. For Build 5, the pipeline assumes the readnoise reference file 
provides a map of CDS readnoise, in electrons. In this script we make
the corrections for those effects.
'''


from jwst_lib.models import ReadnoiseModel,GainModel
from astropy.io import fits
import numpy as np
import argparse,os,sys
import glob
import copy
from sigmacut import calcaverageclass

class mean_readnoise():
    def __init__(self):
        self.verbose = False


    def run(self):
        #read in listfile
        files = []
        with open(self.listfile) as f:
            for line in f:
                if len(line) > 2:
                    files.append(line.strip())
        
        #read in all of the individual readnoise maps
        numfiles = len(files)
        for i,file in enumerate(files):
            if self.verbose:
                print("Reading in "+file)
            map = ReadnoiseModel(file)
            mapdata = map.data
            
            if i == 0:
                all_maps = np.zeros((len(files),mapdata.shape[0],mapdata.shape[1]))

            all_maps[i,:,:] = mapdata

        print("Individual readnoise maps read in")
        print("Beginning averaging")

        #now calculate the sigma-clipped mean readnoise for each pixel
        #could make this smarter by having it look for boxes....might be
        #tough in the case of epoxy voids that are irregularly shaped
        readnoise_map = np.zeros(mapdata.shape)
        readnoise_err = np.zeros(mapdata.shape)
        for y in xrange(mapdata.shape[0]):
            for x in xrange(mapdata.shape[1]):
                pixdata = all_maps[:,y,x]
                sigmacut = calcaverageclass()
                sigmacut.calcaverage_sigmacutloop(pixdata)
                readnoise_map[y,x] = sigmacut.mean
                readnoise_err[y,x] = sigmacut.stdev


        #convert to readnoise for a CDS pair
        readnoise_map = readnoise_map * np.sqrt(2)

        #find the gain reference file if needed
        if self.gainfile == None:
            self.gainfile = self.find_gainfile(files[0])
            print("Using gain file: "+self.gainfile)

        #save the mean readnoise file in ADU, so we don't need to 
        #repeat the steps above when the gain reference file changes
        mapadu = map
        mapadu.data = readnoise_map

        mapadu = self.header_prep(mapadu,files)

        #save the readnoise reference file
        if self.outfile == None:
            self.outfile = self.listfile + '_MEAN_READNOISE_MAP.fits'
        dot = self.outfile.rfind('.') 
        adufile = self.outfile[0:dot] + '_ADU_FOR_CDS.fits'

        #if the specified output file exists, delete it.
        if os.path.isfile(adufile):
            os.remove(adufile)

        print("Saving mean CDS readnoise map in units of ADU as "+adufile)
        mapadu.save(adufile)



        print("Converting to electrons")

        #read in gain file, using SSB's GainModel
        gainmodel = GainModel(self.gainfile)
        gain = gainmodel.data

        #Now convert from ADU to electrons using the gain, and from single-frame to CDS
        #by multiplying by sqrt(2)
        readnoise_map = readnoise_map * gain 

        #put the readnoise_map back into the previously-opened ReadnoiseModel instance
        map.data = readnoise_map
        #take care of all the details for the reference file header keywords
        #Since these files have already gone through mkimreadnoise_ssboutput.py, all that should
        #need to be updated is the HISTORY section, with the names of all the input files.
        print("Preparing header")
        map = self.header_prep(map,files)

        #save the readnoise reference file
        if self.outfile == None:
            self.outfile = self.listfile + '_MEAN_READNOISE_MAP_electrons_FOR_CDS.fits'

        #if the specified output file exists, delete it.
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

        print("Saving mean CDS readnoise image in electrons as "+self.outfile)
        map.save(self.outfile)
    

    def get_datafiles(self,files):
        '''get the names of the data files from a list of individual readnoise files'''
        datafiles = []
        for file in files:
            m = ReadnoiseModel(file)

            data_used = -999
            diff = -999
            for i,lined in enumerate(m.history):
                line = lined['description']
                if 'DATA USED:' in line:
                    data_used = i
                if 'DIFFERENCES:' in line:
                    diff = i

            datafiles = datafiles + m.history[data_used+1:diff]
        return datafiles

    def header_prep(self,model,files):
        '''fine-tune the header of the mean readnoise reference file'''
        #start over with the HISTORY keyword?? how can I avoid duplicating the DATA USED/etc lists?
        new_history = copy.deepcopy(model.history)

        #collect the names of the input files used from the headers of the individual 
        #readnoise files
        datafiles = self.get_datafiles(files)

        #insert the name of the gain file used to convert to electrons
        new_history.append('GAIN FILE USED:')
        gainlen = len(self.gainfile)
        nlines = gainlen / 60
        for i in xrange(nlines+1):
            strt = 60*i
            fin = 60*(i+1)
            if fin > gainlen:
                fin = gainlen
            new_history.append({unicode('description'):self.gainfile[strt:fin]})

        #insert the names of the data files used to make the mean readnoise image
        software = -999
        du = -999
        for i,lined in enumerate(model.history):
            line = lined['description']
            if ('SOFTWARE:' in line):
                software = i
            if ('DATA USED:' in line):
                du = i

        #insert filenames first, since they are lower down in the HISTORY list
        #and won't affect the index number of the software
        for idx,file in enumerate(datafiles):
            #new_history.insert(du+1+idx,{unicode('description'):file})
            new_history.insert(du+1+idx,file)

        #insert software name last since it is towards the top
        #remove the single filename currently in the data used section, since it will 
        #be redundant.
        new_history.insert(software+1,{unicode('description'):'/grp/jwst/wit/nircam/nircam-tools/pythonmodules/'})
        new_history.insert(software+2,{unicode('description'):'combine_readnoise_maps.py'})

        model.history = new_history
        return model

    def find_gainfile(self,file):
        '''locate the appropriate gain file for the given readnoise maps'''
        #gaindir = '/grp/jwst/wit/nircam/cv3_reffile_conversion/gain/'
        gaindir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/cv3_reffile_conversion/gain/'
        gainlist = glob.glob(gaindir+'*DMSorient.fits')
        det = fits.getval(file,'DETECTOR')
        
        if 'LONG' in det:
            det = det[0:4]+'5'

        gainfile = [s for s in gainlist if det in s][0]
        return gainfile

    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='Create mean readnoise map.')
        parser.add_argument('listfile',help="Name of file containing a list of readnoise maps to be averaged.")
        parser.add_argument('--outfile',help="Name of file to output readnoise map to.",default=None)
        parser.add_argument('--gainfile',help="Name of the gain file to use to covert readnoise values to electrons",default=None)
        parser.add_argument('-v','--verbose',help='Verbose output',action='store_true',default=False)
        return parser

if __name__ == '__main__':
    usagestring = 'combine_readnoise_maps.py listfile.list'
    
    meanron = mean_readnoise()
    parser = meanron.add_options(usage=usagestring)
    args = parser.parse_args()
    meanron.listfile = args.listfile
    meanron.gainfile = args.gainfile
    meanron.outfile = args.outfile
    meanron.run()
    
    
