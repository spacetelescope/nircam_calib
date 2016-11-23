#! /usr/bin/env python

# Purpose:  This Python script re-groups a NIRCam exposure to simulate 
#the desired readout pattern.

#DOES NOT YET SUPPORT SUBARRAYS

import sys, os,re, copy
import numpy as np
import argparse
from itertools import izip
from jwst_lib.models import RampModel

# put the tools directory into the path
# add pythonmodules to PATH
if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
else:
    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
    sys.exit(0)


#dictionary of NIRCam readout patterns
deep8 = {}
deep8['tgroup'] = 212.
deep8['ngroup'] = 20.
deep8['nframe'] = 8
deep8['nskip'] = 12

deep2 = {}
deep2['tgroup'] = 212.
deep2['ngroup'] = 20.
deep2['nframe'] = 2
deep2['nskip'] = 18

medium8 = {}
medium8['tgroup'] = 106.
medium8['ngroup'] = 10.
medium8['nframe'] = 8
medium8['nskip'] = 2

medium2 = {}
medium2['tgroup'] = 106.
medium2['ngroup'] = 10.
medium2['nframe'] = 2
medium2['nskip'] = 8

shallow4 = {}
shallow4['tgroup'] = 53.
shallow4['ngroup'] = 10.
shallow4['nframe'] = 4
shallow4['nskip'] = 1

shallow2 = {}
shallow2['tgroup'] = 53.
shallow2['ngroup'] = 10.
shallow2['nframe'] = 2
shallow2['nskip'] = 3

bright2 = {}
bright2['tgroup'] = 21.2
bright2['ngroup'] = 10.
bright2['nframe'] = 2
bright2['nskip'] = 0

bright1 = {}
bright1['tgroup'] = 21.2
bright1['ngroup'] = 10.
bright1['nframe'] = 1
bright1['nskip'] = 1

rapid = {}
rapid['tgroup'] = 10.73676
rapid['ngroup'] = 10.
rapid['nframe'] = 1
rapid['nskip'] = 0

readpatts = {}
readpatts['deep8'] = deep8
readpatts['deep2'] = deep2
readpatts['medium8'] = medium8
readpatts['medium2'] = medium2
readpatts['shallow4'] = shallow4
readpatts['shallow2'] = shallow2
readpatts['bright2'] = bright2
readpatts['bright1'] = bright1
readpatts['rapid'] = rapid







class regroup:
    def __init__(self):
        self.verbose = False

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument("infile",help="Name of the fits file to re-group.")
        parser.add_argument("readpatt",help="Name of readpattern to re-group the data into.")
        parser.add_argument("--ngroup",help="Number of groups to have in the final re-grouped exposure.",default=None,type=int)
        parser.add_argument("--outfile",help="Name of file to output results into.",default=None)
        #parser.add_argument("--clobber",help="Overwrite an existing output file.",action="store_true")
        return parser

    def avg_frame(self,data,sigclip=False):
        '''Average the input frames into a single frame. On-board averaging will be simple averaging.
        Include an option to use sigma-clipped averaging here though'''
        shapes = data.shape
        if len(shapes) != 3:
            print("WARNING! MORE THAN ONE INTEGRATION PASSED TO AVG_FRAME")
            sys.exit(0)

        if sigclip == False:
            newdata = np.mean(data,axis=0)
        else:
            print("FIX ME!")
        return newdata


    def run(self):
        '''main function'''

        #check for the existance of the output file
        if self.outfile == None:
            self.outfile = self.infile[0:-5] + '_REGROUP_'+self.readpatt+'_ngroup'+str(self.ngroup)+'.fits'
            
        if (os.path.isfile(self.outfile)): # & self.clobber == False):
            print("WARNING: Proposed output file {} already exists. Removing.".format(self.outfile))
            os.remove(self.outfile)

        #read in the exposure to use. Read in with RampModel
        exposure = RampModel(self.infile)

        #assume that the readpattern of the input file is 'RAPID'. If not, throw an error.
        rp = exposure.meta.exposure.readpatt
        if rp != 'RAPID':
            print('WARNING! INPUT DATA WERE NOT COLLECTED USING THE RAPID READPATTERN. QUITTING.')
            sys.exit(0)

        #extract data
        data = exposure.data
        err = exposure.err
        groupdq = exposure.groupdq

        #sizes
        integrations = data.shape[0]
        ingroups = data.shape[1]
        ydim = data.shape[2]
        xdim = data.shape[3]

        #if the number of groups was not requested, use the maximum for the given readpattern
        if self.ngroup == None:
            self.ngroup = readpatts[self.readpatt.lower()]['ngroup']

        #group the input groups into collections of frames which will be averaged into the output groups
        #Only group as many input groups as you need to make the requested number of output groups
        frames_per_group = readpatts[self.readpatt.lower()]['nframe']
        frames_to_skip = readpatts[self.readpatt.lower()]['nskip']
        total_frames = (frames_per_group * self.ngroup) + (frames_to_skip * (self.ngroup-1))
        total_exposure_time = total_frames * readpatts['rapid']['tgroup']

        #if the total number of frames needed to make the requested integration don't exist
        #throw an error
        if total_frames > ingroups:
            print("WARNING: Requested regrouping requires more groups than are contained in the input file {}. Quitting.".format(self.infile))
            sys.exit(0)

        #starting and ending indexes of the input groups to be averaged to create the new groups
        groupstart_index = np.arange(0,total_frames,frames_per_group+frames_to_skip)
        groupend_index = groupstart_index + frames_per_group

        #prepare for averaging
        newdata = np.zeros((integrations,self.ngroup,ydim,xdim))
        newerrs = np.zeros((integrations,self.ngroup,ydim,xdim))
        newgroupdq = np.zeros((integrations,self.ngroup,ydim,xdim))

        #average the input data to create the output data
        for integration in xrange(integrations):
            newgp = 0
            for gs,ge in izip(groupstart_index,groupend_index):

                #average the data frames
                print("Averaging groups {} to {}.".format(gs,ge-1))
                newframe = self.avg_frame(data[integration,gs:ge,:,:])

                newdata[integration,newgp,:,:] = newframe

                #reduce the error in the new frames by sqrt(number of frames) for now
                newerrs[integration,newgp,:,:] = err[integration,gs+frames_per_group/2,:,:] / np.sqrt(frames_per_group)

                #just keep the DQ array from the final frame of the group
                newgroupdq[integration,newgp,:,:] = groupdq[integration,ge,:,:]

                #increment the counter for the new group number
                newgp += 1

        #place the updated data back into the model instance
        exposure.data = newdata
        exposure.err = newerrs
        exposure.groupdq = newgroupdq

        #update header 
        exposure.meta.exposure.ngroups = self.ngroup
        exposure.meta.exposure.nframes = frames_per_group
        exposure.meta.exposure.groupgap = frames_to_skip
        exposure.meta.exposure.group_time = readpatts[self.readpatt.lower()]['tgroup']
        exposure.meta.exposure.exptime = total_exposure_time
        exposure.meta.exposure.readpatt = self.readpatt.upper()
        
        #write the regrouped file out to a new file
        exposure.save(self.outfile)
        



if __name__ == '__main__':
    usagestring = 'USAGE: regroup_exposure.py filename.fits deep8'

    rearrange = regroup()
    parser = rearrange.add_options(usage=usagestring)
    args = parser.parse_args(namespace=rearrange)

    rearrange.run()
