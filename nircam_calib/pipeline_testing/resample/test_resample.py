import pytest
import numpy as np
from glob import glob
import json
from astropy.io import fits, ascii
from astropy.coordinates import Angle
from astropy.table import Table, vstack, unique
from jwst.pipeline import calwebb_image3
from jwst.datamodels import RampModel
from jwst.resample import ResampleStep
from jwst.assign_wcs import AssignWcsStep

def get_input_table(num_source_files, sourcelists):
    '''Function to read in and access the simulator source input files.'''

    all_source_table = Table()

    # point source and galaxy source tables have different headers
    # change column headers to match for filtering later
    for i in np.arange(0, num_source_files):
        if "point" in sourcelists[i]:
            col_names = ["RA", "Dec", "RA_degrees", "Dec_degrees",
                         "PixelX", "PixelY", "Magnitude",
                         "counts_sec", "counts_frame"]
        elif "galaxy" in sourcelists[i]:
            col_names = ["PixelX", "PixelY", "RA", "Dec",
                         "RA_degrees", "Dec_degrees", "V2", "V3", "radius",
                         "ellipticity", "pos_angle", "sersic_index",
                         "Magnitude", "countrate_e_s", "counts_per_frame_e"]
        else:
            print('Error! Source list column names need to be defined.')
            sys.exit(0)

        # read in the tables
        input_source_table = ascii.read(sourcelists[i])
        orig_colnames = input_source_table.colnames

        for col, n in zip(orig_colnames, np.arange(0,len(orig_colnames))):
            input_source_table[col].name = col_names[n]

        # only grab values for source catalog analysis
        short_source_table = Table({'In_RA': input_source_table['RA_degrees'],
                                  'In_Dec': input_source_table['Dec_degrees']},
                                  names=['In_RA', 'In_Dec'])

        # combine source lists into one master list
        all_source_table = vstack([all_source_table, short_source_table])

    # set up columns to track which sources were detected by Photutils
    all_source_table['Out_RA'] = np.nan
    all_source_table['Out_Dec'] = np.nan
    all_source_table['Detected'] = 'N'
    all_source_table['RA_Diff'] = np.nan
    all_source_table['Dec_Diff'] = np.nan

    # filter by RA, Dec (for now)
    no_duplicates = unique(all_source_table,keys=['In_RA','In_Dec'])

    return no_duplicates


def  get_output_table(input_name_base):
    '''Function to read in and access the output source catalog.'''

    # use association file to read in output catalog from Image3 pipeline
    print('\nWorking on output table: '+input_name_base+'_cat.ecsv')
    cat_file = glob(input_name_base+'*cat.ecsv')[0]
    output_source_table = ascii.read(cat_file)
    try:
        short_output_table = Table({'RA': output_source_table['sky_centroid'].ra,
                                  'Dec': output_source_table['sky_centroid'].dec,
                                  'PixelX': output_source_table['xcentroid'],
                                  'PixelY': output_source_table['ycentroid']},
                                  names=['RA', 'Dec', 'PixelX', 'PixelY'])
    except:
        short_output_table = Table({'RA': output_source_table['sky_centroid.ra'],
                                  'Dec': output_source_table['sky_centroid.dec']},
                                  names=['RA', 'Dec'])

    return short_output_table


def test_resample(asnfile, RAtol, Dectol):
    '''Test to check the defaults.'''

    # run the pipeline
    img3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image3.cfg')
    img3.tweakreg.skip = True
    img3.output_file = asnfile[:5]+'_default_resample.fits'
    img3.run(asnfile)

    # get association file and source file names
    # create complete list of input source files
    sourcelists = []
    with open(asnfile) as json_data:
         d = json.load(json_data)
    members = d['products'][0]['members']
    for row in np.arange(0,len(members)):
        idx = members[row]['expname'].find("W_")+1
        input_name_base = members[row]['expname'][:idx]
        s1 = glob(input_name_base+'*sources.list')
        s2 = glob(input_name_base+'*Sources.list')
        sourcelists = sourcelists + s1 + s2
    num_source_files = len(sourcelists)

    # get master table of input sources
    all_source_table = get_input_table(num_source_files,sourcelists)

    # read in output source catalog
    output_name_base = d['products'][0]['name']
    if '_cal' in output_name_base:
        cal_idx = output_name_base.find('_cal')
        output_name_base = output_name_base[:cal_idx]
    short_output_table = get_output_table(output_name_base)

    # flag input sources with RA,Dec close to output source RA, Decs
    difference = []
    for i in np.arange(0,len(all_source_table['In_RA'])):
        for j in np.arange(0,len(short_output_table['RA'])):
            RA_diff = np.absolute(all_source_table['In_RA'][i]-short_output_table['RA'][j])
            Dec_diff = np.absolute(all_source_table['In_Dec'][i]-short_output_table['Dec'][j])
            if RA_diff < RAtol and Dec_diff < Dectol:
                all_source_table['Detected'][i] = 'Y'
                all_source_table['Out_RA'][i] = short_output_table['RA'][j]
                all_source_table['Out_Dec'][i] = short_output_table['Dec'][j]
                all_source_table['RA_Diff'][i] = RA_diff
                all_source_table['Dec_Diff'][i] = Dec_diff

    # count number of sources found by Photutils
    total_found = (all_source_table['Detected']=='Y').sum()
    total_missed = (all_source_table['Detected']=='N').sum()
    total_percent_found = (total_found/len(all_source_table))*100
    mask = (all_source_table['Detected']=='Y')

    print('\n')
    print('Total found:',total_found)
    print('Total number of sources:',len(all_source_table))
    print('Total percent found:',total_percent_found)
    print('Number of false detections: ',len(short_output_table)-len(all_source_table))
    print('\n')

    # generate tables with coordinates for comparison
    ascii.write(all_source_table, output_name_base+'_resampled_coords_comparison.dat',format='fixed_width_two_line',overwrite=True)
    ascii.write(all_source_table, output_name_base+'_resampled_coords_comparison.csv',format='csv',overwrite=True)


def test_good_bits(asnfile,good_bits):
    '''Test to check good_bits values.'''

    img3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image2.cfg')
    img3.tweakreg.skip = True
    img3.resample.good_bits = good_bits
    img3.output_file = asnfile[:5]+'_good_bits'+str(good_bits)+'_resample.fits'
    img3.run(asnfile)


def test_kernel(asnfile,kernel):
    '''Test to check kernel.'''

    calib = ResampleStep.call(asnfile,kernel=kernel,output_file=asnfile[:5]+'_kernel'+str(kernel)+'_resample.fits')
    img3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image2.cfg')
    img3.tweakreg.skip = True
    img3.resample.kernel = kernel
    img3.output_file = asnfile[:5]+'_kernel'+str(kernel)+'_resample.fits'
    img3.run(asnfile)


def test_image(asnfile):
    '''Test to single FITS input.'''

    with open(asnfile) as json_data:
         d = json.load(json_data)
    members = d['products'][0]['members'][0]['expname']

    wcs = AssignWcsStep.call(members)
    resample = ResampleStep.call(wcs)
    resample.save(asnfile[:5]+'_singleimage_resample.fits')
