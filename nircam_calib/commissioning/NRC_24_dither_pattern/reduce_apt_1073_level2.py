#!/usr/bin/env python
import sys
import os
import glob
import re
#import photutils
import shutil
import json
import asdf
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase

# The entire calwebb_image2 pipeline
from jwst.pipeline import calwebb_image2

# Individual steps that make up calwebb_image2
from jwst.background import BackgroundStep
from jwst.assign_wcs import AssignWcsStep
from jwst.flatfield import FlatFieldStep
from jwst.photom import PhotomStep
from jwst.resample import ResampleStep
from jwst import datamodels

# This is how it is done via the command line:
# asn_from_list -o mosaic_asn.json --product-name  junk ./reduced/jw01073007*nrcb5*cal.fits
#
# python version
#
car_number ='24'
apt_number = '1073'
nvisits    = 7
#
# Path on orange
#
path    = '/data1/car_'+car_number+'_'+'apt_0'+apt_number+'/mirage/'
reduced = path+'reduced/'

version = 'mirage'
version = 'guitarra'

if(version == 'mirage'):
    path    = '/data1/car_'+car_number+'_'+'apt_0'+apt_number+'/mirage/'
    reduced = path+'reduced/'
#    sca = 'nrca5'
    sca = ''
    prefix  = 'jw0'+apt_number

if(version == 'guitarra'):
    path    = '/data1/car_'+car_number+'_'+'apt_0'+apt_number+'/guitarra/'
    reduced = path+'st_reduced/'
#    sca = '*_485_*'
    sca = ''
    prefix  = 'jw0'+apt_number

output_dir = reduced
print("output_dir ", output_dir)
# exit(0)
for visit in range(1, nvisits+1):
    visit_number = "%03d" % (visit)
    visit_prefix = reduced+prefix+visit_number+sca+'*rate.fits'
    print ("visit_prefix is", visit_prefix)
    file_list = sorted(glob.glob(visit_prefix))
    if(len(file_list) == 0) :
        print('no files for ',visit_prefix)
        continue
    
    print(file_list)

    asn = asn_from_list(file_list,rule=DMSLevel2bBase, products='calibrated_image')
    json_file = reduced+'apt_0'+apt_number+visit_number+'_level-2b_association.json'

    with open(json_file, 'w') as fh:
        fh.write(asn.dump()[1])

    asn_file = os.path.join(output_dir, json_file)
    with open(asn_file) as f_obj:
        asn_data = json.load(f_obj)

    print("asn_data")
    asn_data

# Create an instance of the pipeline class
    image2 = calwebb_image2.Image2Pipeline()

# Set some parameters that pertain to the
# entire pipeline
    image2.output_dir = output_dir
    image2.save_results = True

# Set some parameters that pertain to some of
# the individual steps
    image2.resample.pixfrac = 1.0    # this is the default. Set here as an example

# Call the run() method
    image2.run(asn_file)
