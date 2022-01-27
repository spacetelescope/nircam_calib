#!/usr/bin/env python
import sys
import os
import glob
import re
from pathlib import Path

print("********************************")
print(" ")
print("use    conda  activate jwst_env ")
print(" ")
print("********************************")

car_number ='24'
apt_number = '1073'
nvisits    = 7
version = 'mirage'
version = 'guitarra'
overwrite = False
#overwrite = True

if(version == 'mirage'):
    path    = '/data1/car_'+car_number+'_'+'apt_0'+apt_number+'/mirage/'
    raw     = path+'raw/'
    reduced = path+'reduced/'
    sca = 'nrca5'
    sca = ''
    prefix  = 'jw0'+apt_number
    suffix = '*_uncal.fits'

if(version == 'guitarra'):
    path    = '/data1/car_'+car_number+'_'+'apt_0'+apt_number+'/guitarra/'
    raw     = path+'st_raw/'
    reduced = path+'st_reduced/'
    sca = '_nrca5_'
    sca = ''
    prefix  = 'jw0'+apt_number
    suffix = '*.fits'

output_dir = reduced
print("output_dir ", output_dir)
os.chdir(reduced)
nreduced = 0
for visit in range(1, nvisits + 1):
    visit_number = "%03d" % (visit)
    visit_prefix = raw+prefix+visit_number+sca+suffix
    print ("visit_prefix is", visit_prefix)
    file_list = sorted(glob.glob(visit_prefix))
    if(len(file_list) == 0) :
        print("no files satisfying  ",visit_prefix)
        continue
    
#    print(file_list)
#    exit(0)
    for file in sorted(file_list):
        if(version == 'mirage'):
            if(re.search('_uncal.fits', file)):
                rate_file =  re.sub('_uncal.fits','_rate.fits',file)
#                rate_file =  re.sub(path,'',rate_file)
                rate_file = re.sub('raw','reduced',rate_file)
                file_exists = Path(rate_file)
                print ("rate file is ", rate_file)
                if(file_exists.is_file() and overwrite == False) :
                    nreduced = nreduced +1
                    print("rate file exists: ",rate_file, " ", nreduced)
                else:
                    print(file)
                    command = 'strun jwst.pipeline.Detector1Pipeline '+ file+' --output_dir '+reduced
                    print(command)
                    os.system(command)

        else:
            print("file is ",file," version is ", version)
            slope_file  = file
            slope_file = re.sub('_uncal.fits','_rate.fits', file)
            slope_file = re.sub('st_raw','st_reduced',slope_file)
            file_exists = Path(slope_file)
            if(file_exists.is_file() and overwrite == False) :
                nreduced = nreduced +1
                print("slope file exists: ",slope_file, " nreduced is ", nreduced)
#                exit(0)
            else:
                print("uncal is ", file)
                command = 'strun jwst.pipeline.Detector1Pipeline '+ file+' --output_dir '+reduced
                print(command)
                os.system(command)
#                exit(0)
