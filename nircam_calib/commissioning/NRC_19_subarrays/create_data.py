#! /usr/bin/env python

"""Re-create the NIRCam 001 simulated dataset, which is one of the first datasets
given to DMS. We need to recreate because JWST pipeline formats have changed since
the original data were made. This will put the data into a format that is acceptable
for the pipeline at the moment.
"""

import os
# For examining outputs
from glob import glob
from scipy.stats import sigmaclip
import numpy as np
from astropy.io import fits
from mirage.apt import apt_inputs
from mirage.yaml import yaml_generator
from mirage import imaging_simulator
from jwst.pipeline import calwebb_detector1


def create_data():
    """MAIN FUNCTION"""
    # String the three steps together
    #print('\n\nCREATING YAML FILES\n\n')
    #yam_object = create_yaml_files()



    #Partial run. extended source subarrays are done. Just do point source subarrays and time series
    #yam_object.yaml_files = [e for e in yam_object.yaml_files if 'jw01068005001' in e or 'jw01068006001' in e or 'jw01068007001' in e]


    extsrc_files = sorted(glob('yaml_files/jw0106800[1234]*yaml'))
    tso_files = sorted(glob('yaml_files/jw0106800[567]*yaml'))

    # Remove the TA files from the file list
    ta_files = ['yaml_files/jw01068005001_01101_00001_nrcb5.yaml',
                'yaml_files/jw01068006001_01101_00001_nrcb5.yaml',
                'yaml_files/jw01068007001_01101_00001_nrca5.yaml']
    for ta_file in ta_files:
        tso_files.remove(ta_file)
    #    yam_object.yaml_files.remove(ta_file)


    print('\n\n')
    #for yaml_file in yam_object.yaml_files:
    #for yaml_file in extsrc_files:
    for yaml_file in tso_files:
        print('\nRUNNING MIRAGE: {}\n'.format(yaml_file))
        run_simulator(yaml_file)

    # Find uncal files here
    #uncals = sorted(glob('Mirage_Output/*uncal.fits'))
    #uncals = sorted(glob(os.path.join(yam_object.simdata_output_dir, 'jw01068*_uncal.fits')))

    #print('\n\n')
    # Run CALDETECTOR1 on the resulting files
    #for uncal_file in uncals:
    #    print('\nRUNNING PIPELINE: {}\n'.format(uncal_file))
    #    run_calwebb_detector1(uncal_file)


def create_yaml_files():
    """Create a series of data simulator input yaml files
    from APT files
    """
    xml = 'APT_files/NCam019.xml'
    pointing_file = 'APT_files/NCam019.pointing'
    output_dir = 'yaml_files/'
    simdata_output_dir = 'Mirage_Output'
    catalogs = {'LMC-ASTROMETRIC-FIELD': {'point_source': 'Input_Catalogs/car19_ptsrc.cat'},
                'LMC-ASTROM-CENTER-JMAG14-STAR': {'point_source': 'Input_Catalogs/car19_ptsrc.cat'}
                }
    background = 'low'
    yam = yaml_generator.SimInput(xml, pointing_file, output_dir=output_dir, simdata_output_dir=simdata_output_dir, catalogs=catalogs,
                                  datatype='raw', dateobs_for_background=False, reffile_defaults='crds', background=background)
    yam.use_linearized_darks = True
    yam.create_inputs()
    return yam


def run_calwebb_detector1(filename):
    m = calwebb_detector1.Detector1Pipeline(config_file='pipeline_config_files/calwebb_detector1.cfg')

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
    m.output_dir, m.output_file = os.path.split(filename.replace('_uncal', '_rate'))

    # run the pipeline with these paramters
    m.run(filename)
    print('')
    print("Done running CALDETECTOR1 on {}".format(filename))
    print("Output saved to {}".format(m.output_file))
    print('')


def run_simulator(filename):
    img_sim_custom_yamls = imaging_simulator.ImgSim()
    img_sim_custom_yamls.paramfile = filename
    img_sim_custom_yamls.create()
    print('')
    print("Done creating file from {}".format(filename))
    print('')


if __name__ == '__main__':
    create_data()
