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


percent_found_tolerance = 80.0
RAtol = 0.002
Dectol = 0.002


def check_catalog(input_sourcelists,out_file):

    # get new source catalog generated from Image3 pipeline
    cat_file_list = glob('*cat.ecsv')
    newest_catalog = max(cat_file_list, key=os.path.getctime)
    print('\nGetting output table: '+newest_catalog)

    output_source_table = ascii.read(newest_catalog)

    try:
        short_output_table = Table({'RA': output_source_table['sky_centroid'].ra,
                                  'Dec': output_source_table['sky_centroid'].dec},
                                  names=['RA', 'Dec'])
    except:
        short_output_table = Table({'RA': output_source_table['sky_centroid.ra'],
                                  'Dec': output_source_table['sky_centroid.dec']},
                                  names=['RA', 'Dec'])


    # flag input sources with RA,Dec close to output source RA, Decs
    difference = []
    for i in np.arange(0,len(input_sourcelists['In_RA'])):
        for j in np.arange(0,len(short_output_table['RA'])):
            RA_diff = np.absolute(input_sourcelists['In_RA'][i]-short_output_table['RA'][j])
            Dec_diff = np.absolute(input_sourcelists['In_Dec'][i]-short_output_table['Dec'][j])
            if RA_diff < RAtol and Dec_diff < Dectol:
                all_source_table['Detected'][i] = 'Y'
                all_source_table['Out_RA'][i] = short_output_table['RA'][j]
                all_source_table['Out_Dec'][i] = short_output_table['Dec'][j]
                all_source_table['RA_Diff'][i] = RA_diff
                all_source_table['Dec_Diff'][i] = Dec_diff

    # count number of sources found by Photutils
    total_found = (input_sourcelists['Detected']=='Y').sum()
    total_missed = (input_sourcelists['Detected']=='N').sum()
    total_percent_found = (total_found/len(input_sourcelists))*100
    mask = (input_sourcelists['Detected']=='Y')

    print('\n')
    print('Total found:',total_found)
    print('Total number of sources:',len(input_sourcelists))
    print('Total percent found:',total_percent_found)
    print('Number of false detections: ',len(short_output_table)-len(input_sourcelists))
    print('\n')

    # generate tables with coordinates for comparison
    ascii.write(input_sourcelists, out_file[:-5]+'_resampled_coords_comparison.dat',format='fixed_width_two_line',overwrite=True)
    ascii.write(input_sourcelists, out_file[:-5]+'_resampled_coords_comparison.csv',format='csv',overwrite=True)

    return total_percent_found


# @pytest.fixture(scope='module')
def fits_input(input_file):

    # open the input_file defined above once for each module
    if 'fits' in input_file:
        return fits.open(input_file)
    elif 'json' in input_file:
        return input_file


# @pytest.fixture(scope='module')
def in_datamodel(fits_input):
    '''open input file as a datamodel or leave as ASN table.'''

    if 'json' in fits_input:

        with open(fits_input) as json_data:
             d = json.load(json_data)
        return d
    else:
        return fits_input[0].header['FILENAME']


# @pytest.fixture(scope='module')
def input_sourcelists(in_datamodel):
    '''Function to read in and access the simulator source input files.'''

    sourcelists = []

    if type(in_datamodel) is dict:
        members = in_datamodel['products'][0]['members']

        for row in np.arange(0,len(members)):
            idx = members[row]['expname'].find('_')
            input_name_base = members[row]['expname'][:idx]
            s1 = glob(input_name_base+'*sources.list')
            s2 = glob(input_name_base+'*Sources.list')
            sourcelists = sourcelists + s1 + s2

    else:
        # information to find simulator source lists
        under = in_datamodel.find('_')
        input_name_base = in_datamodel[:under]

        # get a list of source files to generate master input source list
        s1 = glob(input_name_base+'*sources.list')
        s2 = glob(input_name_base+'*Sources.list')
        sourcelists = sourcelists + s1 + s2

    all_source_table = Table()

    for i in np.arange(0, len(sourcelists)):

        # point source and galaxy source tables have different headers
        # change column headers so they match
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

        # read in the source list tables
        input_source_table = ascii.read(sourcelists[i])
        orig_colnames = input_source_table.colnames

        for col, n in zip(orig_colnames, np.arange(0,len(orig_colnames))):
            input_source_table[col].name = col_names[n]

        # only grab columns for source catalog analysis
        short_source_table = Table({'In_RA': input_source_table['RA_degrees'],
                                  'In_Dec': input_source_table['Dec_degrees']},
                                  names=['In_RA', 'In_Dec'])

        # combine source lists into final master list
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



def test_resample(in_datamodel,input_file,input_sourcelists):
    '''Test to check the defaults.'''

    # run image3 pipeline to generate source catalog (1 step after resample)
    img3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image3.cfg')
    img3.output_file = input_file[:-5]+'_default_resample.fits'
    img3.run(input_file)

    results = check_catalog(input_file,img3.output_file,RAtol,Dectol)
    print('Percent of input sources found with defaults: ',results)

    assert (results > 80.0) == True



def test_gaussian_kernel(in_datamodel,input_file,input_sourcelists):
    '''Test to check gaussian kernel.'''

    # run image3 pipeline to generate source catalog (1 step after resample)
    img3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image3.cfg')
    img3.resample.kernel = 'gaussian'
    img3.output_file = input_file[:-5]+'_kernel'+str(kernel)+'_resample.fits')
    img3.run(input_file)

    results = check_catalog(input_file,img3.output_file,RAtol,Dectol)
    print('Percent of input sources found with defaults: ',results)

    assert (results > 80.0) == True
