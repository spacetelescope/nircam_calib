#! /usr/bin/env python

"""Run the pipeline on the CAR-19 data
"""
from glob import glob
import os

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
from jwst.pipeline.calwebb_image2 import Image2Pipeline
from jwst.pipeline.calwebb_image3 import Image3Pipeline
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.asn_from_list import asn_from_list


def run_calwebb_detector1(filename):
    m = Detector1Pipeline(config_file='pipeline_config_files/calwebb_detector1.cfg')

    # make changes to the parameters/reference files used
    m.refpix.odd_even_rows = False

    # jump step is way too sensitive
    m.jump.rejection_threshold = 91

    # skip steps you don't want to run
    m.group_scale.skip = True
    m.ipc.skip = True
    m.rscd.skip = True
    m.firstframe.skip = True
    m.lastframe.skip = True

    # name your output file
    m.save_results = True
    m.output_dir = '/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/'
    m.output_file = os.path.basename(filename.replace('_uncal', '_rate'))

    # run the pipeline with these paramters
    m.run(filename)
    print('')
    print("Done running CALDETECTOR1 on {}".format(filename))
    print("Output saved to {}".format(os.path.join(m.output_dir, m.output_file)))
    print('')


def run_calwebb_image2(filename, tso=False):
    if not tso:
        result2 = Image2Pipeline(config_file='pipeline_config_files/calwebb_image2.cfg')
    else:
        result2 = Image2Pipeline(config_file='pipeline_config_files/calwebb_tso-image2.cfg')
    result2.save_results = True
    result2.output_dir = '/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/'
    result2.run(filename)


def run_calwebb_image3(filename, tso=False):
    if not tso:
        result = Image3Pipeline(config_file='pipeline_config_files/calwebb_image3.cfg')
    else:
        result = Image3Pipeline(config_file='pipeline_config_files/calwebb_tso3.cfg')
    result.save_results = True
    result.source_catalog.save_results = True
    result.source_catalog.output_dir = '/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level3/'
    result.output_dir = '/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level3/'
    result.run(filename)


def make_level2_association(file_list, asn_filename):
    idx = file_list[0].find('nrca1')
    prod_name = file_list[0][0: idx+5]
    asn = asn_from_list(file_list, rule=DMSLevel2bBase, product_name=prod_name)
    outfile = os.path.join('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/', asn_filename)
    with open(outfile, 'w') as fh:
        fh.write(asn.dump()[1])


def make_level3_association(file_list, asn_filename):
    asn = asn_from_list(file_list, product_name='lw_imaging')
    outfile = os.path.join('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level3/', asn_filename)
    with open(outfile, 'w') as fh:
        fh.write(asn.dump()[1])


# Run Detector1
uncal_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Mirage_Output/j*uncal.fits')
#for filename in uncal_files:
#    run_calwebb_detector1(filename)

# Create association files. For level2 it's not as important, but let's make one association
# file for each subarray size in the extended subarray data
sub160_asn = 'Pipeline_Level2/level2_sub160_files_asn.json'
sub160_rate_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/jw01068001001*rate.fits')
#make_level2_association(sub160_rate_files, sub160_asn)

sub320_asn = 'Pipeline_Level2/level2_sub320_files_asn.json'
sub320_rate_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/jw01068002001*rate.fits')
#make_level2_association(sub320_rate_files, sub320_asn)

sub640_asn = 'Pipeline_Level2/level2_sub640_files_asn.json'
sub640_rate_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/jw01068003001*rate.fits')
#make_level2_association(sub640_rate_files, sub640_asn)

full_asn = 'Pipeline_Level2/level2_full_files_asn.json'
full_rate_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/jw01068004001*rate.fits')
#make_level2_association(full_rate_files, full_asn)

sub400p_asn = 'Pipeline_Level2/level2_sub400p_files_asn.json'
sub400p_rate_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/jw01068005001*rate.fits')
#make_level2_association(sub400p_rate_files, sub400p_asn)

sub64p_asn = 'Pipeline_Level2/level2_sub64p_files_asn.json'
sub64p_rate_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/jw01068006001*rate.fits')
#make_level2_association(sub64p_rate_files, sub64p_asn)

substripe256_asn = 'Pipeline_Level2/level2_substrip256_files_asn.json'
substripe256_rate_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level1/jw01068007001*rate.fits')
#make_level2_association(substrip256_rate_files, substrip256_asn)

#print('Manually add asn_pool to the level2 association files')

#for s160file in full_rate_files:
#    run_calwebb_image2(s160file)
#stop

# Run Image2
#association_files = [sub160_asn, sub320_asn, sub640_asn, full_asn, sub400p_asn, sub64p_asn]
association_files = [sub640_asn, full_asn, sub400p_asn, sub64p_asn]
#for asn in association_files:
#    if asn != substripe256_asn:
#        run_calwebb_image2(asn)
#    else:
#        run_calwebb_image2(asn, tso=True)



# Create Image3 association files
sub160_sw_asn_3 = 'level3_sub160_sw_files_asn.json'
sub160_sw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068001001*nrcb[1234]_cal.fits'))
make_level3_association(sub160_sw_cal_files, sub160_sw_asn_3)

sub160_lw_asn_3 = 'level3_sub160_lw_files_asn.json'
sub160_lw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068001001*nrcb5_cal.fits'))
make_level3_association(sub160_lw_cal_files, sub160_lw_asn_3)

sub320_sw_asn_3 = 'level3_sub320_sw_files_asn.json'
sub320_sw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068002001*nrcb[1234]_cal.fits'))
make_level3_association(sub320_sw_cal_files, sub320_sw_asn_3)

sub320_lw_asn_3 = 'level3_sub320_lw_files_asn.json'
sub320_lw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068002001*nrcb5_cal.fits'))
make_level3_association(sub320_lw_cal_files, sub320_lw_asn_3)

sub640_sw_asn_3 = 'level3_sub640_sw_files_asn.json'
sub640_sw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068003001*nrcb[1234]_cal.fits'))
make_level3_association(sub640_sw_cal_files, sub640_sw_asn_3)

sub640_lw_asn_3 = 'level3_sub640_lw_files_asn.json'
sub640_lw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068003001*nrcb5_cal.fits'))
make_level3_association(sub640_lw_cal_files, sub640_lw_asn_3)

full_sw_asn_3 = 'level3_full_sw_files_asn.json'
full_sw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068004001*nrcb[1234]_cal.fits'))
make_level3_association(full_sw_cal_files, full_sw_asn_3)

full_lw_asn_3 = 'level3_full_lw_files_asn.json'
full_lw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068004001*nrcb5_cal.fits'))
make_level3_association(full_lw_cal_files, full_lw_asn_3)

sub400p_sw_asn_3 = 'level3_sub400p_sw_files_asn.json'
sub400p_sw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068005001*nrcb[1234]_cal.fits'))
make_level3_association(sub400p_sw_cal_files, sub400p_sw_asn_3)

sub400p_lw_asn_3 = 'level3_sub400p_lw_files_asn.json'
sub400p_lw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068005001*nrcb5_cal.fits'))
make_level3_association(sub400p_lw_cal_files, sub400p_lw_asn_3)

sub64p_sw_asn_3 = 'level3_sub64p_sw_files_asn.json'
sub64p_sw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068006001*nrcb[1234]_cal.fits'))
make_level3_association(sub64p_sw_cal_files, sub64p_sw_asn_3)

sub64p_lw_asn_3 = 'level3_sub64p_lw_files_asn.json'
sub64p_lw_cal_files = sorted(glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068006001*nrcb5_cal.fits'))
make_level3_association(sub64p_lw_cal_files, sub64p_lw_asn_3)

substripe256_sw_asn_3 = 'level3_substripe256_sw_files_asn.json'
#substripe256_sw_cal_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068007001*nrc?[1234]_cal.fits')
#make_level3_association(substripe256_sw_cal_files, substripe256_sw_asn_3)

substripe256_lw_asn_3 = 'level3_substripe256_lw_files_asn.json'
#substripe256_lw_cal_files = glob('/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Pipeline_Level2/jw01068007001*nrc?5_cal.fits')
#make_level3_association(substripe256_lw_cal_files, substripe256_lw_asn_3)

# Run Image3
association_files_3 = [sub160_sw_asn_3, sub160_lw_asn_3, sub320_sw_asn_3, sub320_lw_asn_3, sub640_sw_asn_3, sub640_lw_asn_3,
                       full_sw_asn_3, full_lw_asn_3, sub400p_sw_asn_3, sub400p_lw_asn_3, sub64p_sw_asn_3, sub64p_lw_asn_3]
for asn in association_files_3:
    if asn not in [substripe256_sw_asn_3, substripe256_lw_asn_3]:
        run_calwebb_image3(os.path.join('Pipeline_Level3/', asn))
    else:
        run_calwebb_image3(os.path.join('Pipeline_Level3/', asn), tso=True)

