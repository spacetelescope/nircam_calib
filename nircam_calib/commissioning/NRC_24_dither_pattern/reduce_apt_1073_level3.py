#!/usr/bin/env python
import sys
import os
import glob
import re
import photutils
import shutil
import json
import asdf
#
from astropy.io import fits
#
# The entire calwebb_image3 pipeline
import jwst
from jwst.pipeline import calwebb_image3

# Individual steps that make up calwebb_image3
from jwst.tweakreg import TweakRegStep
from jwst.skymatch import SkyMatchStep
from jwst.outlier_detection import OutlierDetectionStep
from jwst.resample import ResampleStep
from jwst.source_catalog import SourceCatalogStep
from jwst import datamodels
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
#
#-------------------------------------------------------------------------------
# This can be expanded to include other parameters
#
def read_fits_image(file, verbose):
    hdulist = fits.open(file)
    if(verbose >0 ) :
        hdulist.info()
    header = hdulist[0].header
    filter = header['FILTER']
    observation =  header['OBSERVTN']
    visit       =  header['VISIT']
    visit_group =  header['VISITGRP']
    seq_id      =  header['SEQ_ID']
    act_id      =  header['ACT_ID']
    expnumber   =  header['EXPOSURE']
    detector    =  header['DETECTOR']
    module      =  header['MODULE']
    
#    naxis1 = header['NAXIS1']
#    naxis2 = header['NAXIS2']
#    image = hdulist[0].data
#    return naxis1, naxis2, image
    hdulist.close()
    return filter, observation, visit, visit_group, seq_id, act_id, expnumber, detector, module
#
#-------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------
#
def get_file_list(obs_number, filter, detectors, channel_module, file_list):

    swa = ['NRCA1','NRCA2','NRCA3','NRCA4']
    swb = ['NRCB1','NRCB2','NRCB3','NRCB4']
    lwa = ['NRCALONG']
    lwb = ['NRCBLONG']
    sw  = swa + swb
    lw  = lwa + lwb
    
    detector_list = detectors
    
    list = []

    if(channel_module == 'swa'):
        detector_list = swa
    if(channel_module == 'swb'):
        detector_list = swb
    if(channel_module == 'sw'):
        detector_list = sw
    if(channel_module == 'lwa'):
        detector_list = lwa
    if(channel_module == 'lwb'):
        detector_list = lwb
    if(channel_module == 'lw'):
        detector_list = lw
    
    if(channel_module == 'swa' or \
       channel_module == 'swb' or \
       channel_module == 'sw' or \
       channel_module == 'lwa' or \
       channel_module == 'lwb' or \
       channel_module == 'lw'):

        for index in range(0, len(detector_list)):
            sca = detector_list[index]
            for file in  sorted(file_list):
                if(file_filter[file] == filter and \
                   file_detector[file] == sca and \
                   file_observation[file] == obs_number):
                    list.append(file)
        return list
    else:
# for single detectors    
        for index in range(0, len(detector_list)):
            sca = detector_list[index]
            for file in  sorted(file_list):
                if(file_filter[file] == filter and \
                   file_detector[file] == sca and \
                   file_observation[file] == obs_number):
                    list.append(file)
                    print("filter: sca ", filter, sca, len(list))
        return list
#
#-------------------------------------------------------------------------------
#
def run_image3(file_list, prefix, obs_number, suffix, filter, output_dir, verbose):
               
    mosaic = prefix+obs_number+'_'+suffix+'_'+filter

    print("run_image3 : mosaic ", mosaic)
    
    asn = asn_from_list(file_list,rule=DMS_Level3_Base, product_name=mosaic)
    
    json_file = output_dir+prefix+obs_number+'_'+suffix+'_'+filter+'_level-3_association.json'
    if(verbose > 0) :
        print ("def run_image3: json_file is ", json_file)

    with open(json_file, 'w') as fh:
        fh.write(asn.dump()[1])

    asn_file = os.path.join(output_dir, json_file)
    with open(asn_file) as f_obj:
        asn_data = json.load(f_obj)

    if(verbose > 0) :
        print("def run_image3: asn_data", asn_data)
    asn_data
#    return
# Create an instance of the pipeline class
    image3 = calwebb_image3.Image3Pipeline()

# Set some parameters that pertain to the
# entire pipeline
    image3.output_dir = output_dir
    image3.save_results = True

# Set some parameters that pertain to some of
# the individual steps
#    image3.tweakreg.snr_threshold = 10.0  # 5.0 is the default
#    image3.tweakreg.kernel_fwhm = 2.302  # 2.5 is the default
#    image3.tweakreg.brightest = 20  # 100 is the default
    image3.tweakreg.skip = True
    image3.source_catalog.kernel_fwhm = 2.302  # pixels
    image3.source_catalog.snr_threshold = 10.
# Call the run() method
    image3.run(asn_file)
    return
    #
#-------------------------------------------------------------------------------
#

verbose    = 0 
car_number ='24'
apt_number = '1073'
nvisits    = 7
#
# Path on orange for mirage simulations
#
version = 'mirage'
version = 'guitarra'

channel_module =['lwa', 'lwb', 'swa','swb']

sw_filters = ('F070W', 'F115W', 'F150W')
lw_filters = ('F277W')

filters    = ('F070W', 'F115W', 'F150W', 'F277W')
detectors  = ('NRCA1','NRCA2','NRCA3','NRCA4','NRCB1','NRCB2','NRCB3','NRCB4','NRCALONG','NRCBLONG')

if(version == 'mirage'):
    path    = '/data1/car_'+car_number+'_'+'apt_0'+apt_number+'/mirage/'
    reduced = path+'reduced/'
    sca = 'nrca5'
    prefix  = 'jw0'+apt_number
    swa = ('nrca1','nrca2','nrca3','nrca4')
    swb = ('nrcb1','nrcb2','nrcb3','nrcb4')
    lwa = ('nrca5')
    lwb = ('nrcb5')
    sw  = swa + swb
    lw  = lwa + lwb
    detector_names= ('nrca1','nrca2','nrca3','nrca4','nrcb1','nrcb2','nrcb3','nrcb4','nrca5','nrcb5')

if(version == 'guitarra'):
    path    = '/data1/car_'+car_number+'_'+'apt_0'+apt_number+'/guitarra/'
    reduced = path+'st_reduced/'
    prefix  = 'jw0'+apt_number+'_'
#    swa = ('481','482','483','484')
#    swb = ('486','487','488','489')
#    lwa = ('485')
#    lwb = ('490')
#    detector_names = ('481','482','483','484','486','487','488','489','485','490')
    swa = ('nrca1','nrca2','nrca3','nrca4')
    swb = ('nrcb1','nrcb2','nrcb3','nrcb4')
    lwa = ('nrca5')
    lwb = ('nrcb5')
    sw  = swa + swb
    lw  = lwa + lwb
    detector_names= ('nrca1','nrca2','nrca3','nrca4','nrcb1','nrcb2','nrcb3','nrcb4','nrca5','nrcb5')

#detectors = detector_names
output_dir = reduced
#
# Create lists where images are sorted according to filter
# channel (SW/LW), possibly module (A/B) and SCA name
#
# DMS name contains
# observation+visit+*+sca_id.fits
# it does not store the filter in file name, thus one has
# to read the file header to recover the filter name
# Guitarra does contain the filter name
#
file_list        = sorted(glob.glob(reduced+'*cal.fits'))
file_filter      = {}
file_observation = {}
file_detector    = {}


for file in sorted(file_list):
    (filter, observation, visit, visit_group, seq_id, act_id, expnumber, detector, module)\
        = read_fits_image(file, verbose)
    file_filter[file]      = filter
    file_observation[file] = observation
    file_detector[file]    = detector
#    print (file, observation, filter, detector)


for observation in range(1, nvisits+1):
    obs_number = "%03d" % (observation)

    for filter_index in range(0, len(filters)):
        filter = filters[filter_index]

        for field  in  range(0, len(channel_module)):
            list = get_file_list(obs_number, filter, detectors, channel_module[field], file_list)
            suffix     = channel_module[field]
            if(len(list) == 0):
                if(verbose > 0) :
                    print("no files in list: ", list, ' for ', obs_number, filter, channel_module[field])
                continue
            else:
                print ("\nobservation", version, prefix, obs_number, suffix, filter, ' nfiles: ',len(list))
                if(verbose > 0) :
                    print ("observation", obs_number, list)
                run_image3(list, prefix, obs_number, suffix, filter, output_dir, verbose)
            
