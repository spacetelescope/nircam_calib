"""
This package can be used for basic validation of GSEG data.
Beginning with the APT file, this software will construct a
table of expected files and their properties (subarray, number of groups, etc)

It will then check the actual data to see if it matches the expected
properties.

Checks to perform:
1. file exists?
2. data matches the expected shape
3. confirm header keyword values:
    a. subarray
    b. detector
    c. starting coordinates on detector
    d. num ints, num groups, ny, nx
    e. exposure time



Summary
-------

Read in APT file using Mirage. Create table of exposure parameters
Loop through files
  perform checks above
Any checks on pixel values? These would need to be based on the
signal value pattern that John figured out (see his email)

"""

import os
from glob import glob
from astropy.io import fits
from mirage.yaml import yaml_generator
from mirage.apt import apt_inputs
from mirage.utils.siaf_interface import sci_subarray_corners
from mirage.utils.utils import calc_frame_time
from mirage.yaml.generate_observationlist import get_observation_dict
import numpy as np
import pysiaf

UNCAL_KEYWORDS = ['SUBARRAY', 'DETECTOR', 'NINTS', 'NGROUPS', 'NAXIS', 'EFFEXPTM',
                  'LONGFILTER', 'LONGPUPIL', 'SHORTFILTER', 'SHORTPUPIL', 'READPATT',
                  'OBSLABEL', 'EXP_TYPE', 'TITLE', 'OBSERVTN', 'TEMPLATE',
                  'EXPRIPAR', 'SUBSTRT1', 'SUBSTRT2', 'SUBSIZE1', 'SUBSIZE2',
                  'FASTAXIS', 'SLOWAXIS', 'PATTTYPE']

# These must correspond one-to-one with UNCAL_KEYWORDS ABOVE
UNCAL_TABLE_KEYWORDS = ['Subarray', None, 'Integrations', 'Groups', None, None,
                        'LongFilter', 'LongPupil', 'ShortFilter', 'ShortPupil', 'ReadoutPattern',
                        'ObservationName', 'Mode', 'Title', 'ObservationID', 'APTTemplate',
                        'ParallelInstrument', None, None, None, None,
                        None, None, 'PrimaryDitherType']

INTEGER_KEYWORDS = ['Integrations', 'Groups']
FLOAT_KEYWORDS = ['EFFEXPTM']
FILTER_KEYWORDS = ['LONGFILTER', 'LONGPUPIL', 'SHORTFILTER', 'SHORTPUPIL']

HORIZONTAL_FLIP = ['NRCA1', 'NRCA3', 'NRCALONG', 'NRCB2', 'NRCB4']
VERTICAL_FLIP = ['NRCA2', 'NRCA4', 'NRCB1', 'NRCB3', 'NRCBLONG']


def adjust_exptype(value):
    """Modify the exposure type as listed in the exposure table
    to match one of the strings as used in the fits files.
    e.g. 'imaging' becomes 'NRC_IMAGE'
    Remember that currently, Mirage only knows imaging and wfss
    """
    if value == 'imaging':
        return 'NRC_IMAGE'
    elif value == 'wfss':
        return 'NRC_GRISM'


def calculate_total_files(exp_dict, index):
    """Calculate the total number of files expected for an
    observation based on the number of dithers and the module.
    ASSUME that all detectors in a given module are used.
    """
    module = exp_dict['Module'][index]
    number_of_dithers = exp_dict['number_of_dithers'][index]
    if module in ['A', 'B']:
        dets = 5
    else:
        dets = 10
    total = number_of_dithers * dets
    return total


def equalize_file_lists(uncal, rate):
    """Given lists of uncal and rate files corresponding to a single
    observation, adjust the lists to be the same length, adding in
    None for any files that are missing in a given list
    """
    udict = {}
    rdict = {}
    expanded_rate = []
    expanded_uncal = []

    # Loop through uncal files and look for matching rate files
    for ufile in uncal:
        dirname, filename = os.path.split(ufile)
        base = filename.strip('_uncal.fits')
        fullbase = os.path.join(dirname, base)
        found = False
        for rfile in rate:
            if fullbase in rfile:
                found = True
                break
        udict[base] = found

    # Loop through rate files and look for matching uncal files
    for rfile in rate:
        dirname, filename = os.path.split(rfile)
        base = filename.strip('_rate.fits')
        fullbase = os.path.join(dirname, base)
        found = False
        for ufile in uncal:
            if fullbase in ufile:
                found = True
                break
        rdict[base] = found

    # Fill in missing files, in either uncal or rate lists,
    # with None
    for ukey in udict:
        expanded_uncal.append(ukey + '_uncal.fits')
        if udict[key]:
            expanded_rate.append(ukey + '_rate.fits')
        else:
            expanded_rate.append(None)
    for rkey in rdict:
        if not rdict[key]:
            expanded_rate.append(rkey + '_rate.fits')
            expanded_uncal.append(None)
    return expanded_uncal, expanded_rate


def find_fastaxis(detector):
    """Identify the values of FASTAXIS and SLOWAXIS based on the detector
    name
    """
    if detector in HORIZONTAL_FLIP:
        fast = -1
        slow = 2
    elif detector in VERTICAL_FLIP:
        fast = 1
        slow = -2
    return fast, slow


def get_data(filename):
    """Read in the given fits file and return the data and header
    """
    with fits.open(filename) as h:
        signals = h['SCI'].data
        header0 = h[0].header
        header1 = h[1].header
    return signals, header0, header1


def uncal_header_keywords(head):
    """Extract values for the desired keywords from the given header
    """
    file_info = {}
    for keyword in UNCAL_KEYWORDS:
        try:
            info = header[keyword]
        except KeyError:
            if 'FILTER' in keyword:
                info = header['FILTER']
            elif 'PUPIL' in keyword:
                info = header['PUPIL']
            else:
                info = None

        file_info[keyword] = info
    return file_info


def uncal_table_info(values, index):
    """Extract information from the exposure table that matches the
    header keyword values in uncal_header_keyword
    """
    values_dict = {}
    for table_keyword, file_keyword in zip(UNCAL_TABLE_KEYWORDS, UNCAL_KEYWORDS):
        if table_keyword is not None:
            if table_keyword in INTEGER_KEYWORDS:
                value = int(values[table_keyword][index])
            else:
                value = values[table_keyword][index]
            values_dict[file_keyword] = value
        else:
            values_dict[file_keyword] = None
    return values_dict


#xml_file = '../dry_run_may_2_2019/apt_file/00617.xml'
#pointing_file = xml_file.replace('.xml', '.pointing')
#output_dir = '../dry_run_may_2_2019/expected_parameters/'
#gseg_uncal_files = glob('../dry_run_may_2_2019/MAST_2019-05-02T2112/JWST/jw*/*uncal.fits')
#gseg_rate_files = [f.replace('uncal', 'rate') for f in gseg_uncal_files]


def validate(xml_file, output_dir, gseg_uncal_files):
    """MAIN FUNCTION"""
    pointing_file = xml_file.replace('.xml', '.pointing')
    gseg_rate_files = [f.replace('uncal', 'rate') for f in gseg_uncal_files]

    catalogs = {'nircam': {'sw': 'nothing.cat', 'lw': 'nothing.cat'}}

    observation_list_file = os.path.join(output_dir, 'observation_list.yaml')
    apt_xml_dict = get_observation_dict(xml_file, observation_list_file, catalogs,
                                    verbose=True)

    observation_list = set(apt_xml_dict['ObservationID'])
    int_obs = sorted([int(o) for o in observation_list])
    str_obs_list = [str(o).zfill(3) for o in int_obs]

    for observation_to_check in str_obs_list:
        print('')
        print('')
        print('OBSERVATION: {}'.format(observation_to_check))
        print('')

        good = np.where(np.array(apt_xml_dict['ObservationID']) == observation_to_check)

        try:
            total_expected_files = calculate_total_files(apt_xml_dict, good[0][0])
            print('Total number of expected files: {}'.format(total_expected_files))
        except IndexError:
            print("No files found.")
            continue

        # The complication here is that the table created by Mirage does not have a filename
        # attached to each entry. So we need a way to connect an actual filename
        # to each entry
        subdir_start = 'jw' + apt_xml_dict['ProposalID'][good[0][0]] + observation_to_check.zfill(3)
        matching_uncal_files = sorted([filename for filename in gseg_uncal_files if subdir_start in filename])
        matching_rate_files = sorted([filename for filename in gseg_rate_files if subdir_start in filename])
        print('Found uncal files:')
        for i in range(len(matching_uncal_files)):
            print(matching_uncal_files[i])
        print('')
        print('Found rate files:')
        for i in range(len(matching_rate_files)):
            print(matching_rate_files[i])
        print('')

        # Check to see if any files are missing
        if len(matching_uncal_files) != total_expected_files:
            print("WARNING: Missing uncal files for observation {}. Expected {} files, found {}.".format(observation_to_check, total_expected_files, len(matching_uncal_files)))
        if len(matching_rate_files) != total_expected_files:
            print("WARNING: Missing rate files for observation {}. Expected {} files, found {}.".format(observation_to_check, total_expected_files, len(matching_rate_files)))

        # Deal with the case of matching_uncal_files and matching_rate_files having
        # different lengths here. In order to loop over them they must have the same length
        if len(matching_uncal_files) != len(matching_rate_files):
            (matching_uncal_files, matching_rate_files) = equalize_file_lists(matching_uncal_files, matching_rate_files)
            print('Equalized file lists (should have a 1:1 correspondence):')
            for idx in range(len(matching_uncal_files)):
                print(matching_uncal_files[idx], matching_rate_files[idx])

        # Create siaf instance for later calculations
        siaf = pysiaf.Siaf('NIRCam')

        for uncal, rate in zip(matching_uncal_files, matching_rate_files):
            good_uncal = uncal != None
            good_rate = rate != None

            if good_uncal:
                print("Checking {}".format(os.path.split(uncal)[1]))
                print('-----------------------------------------------')
            elif good_rate:
                print("Checking {}".format(os.path.split(rate)[1]))
                print('-----------------------------------------------')


            if good_uncal:
                data, header, sci_header = get_data(uncal)
                detector_from_filename = uncal.split('_')[-2].upper()
                header_detector = header['DETECTOR']
                if 'LONG' in header_detector:
                    header_detector = header_detector.replace('LONG', '5')
                if header_detector not in header['APERNAME']:
                    print(("WARNING: Detector name and aperture name in file header appear to be incompatible: {}, {}"
                          .format(header['DETECTOR'], header['APERNAME'])))
                    print("Detector listed in filename: {}".format(detector_from_filename))
                    print('If the aperture is incorrect then the calculated subarray location from pysiaf will also be incorrect.')
                data_shape = data.shape

                # Get info from header to be compared
                header_vals = uncal_header_keywords(header)

                # Get matching data from the exposure table
                table_vals = uncal_table_info(apt_xml_dict, good[0][0])

                # Make some adjustments to the exposure table info

                # Calucate the exposure time
                aperture = header['APERNAME']  # could also try APERNAME, PPS_APER

                print('Aperture listed in header is: {}'.format(aperture))

                num_amps = 1
                frametime = calc_frame_time('NIRCam', aperture, data_shape[-1], data_shape[-2], num_amps)
                table_vals['EFFEXPTM'] = frametime * int(table_vals['NGROUPS'])

                # NAXIS
                table_vals['NAXIS'] = len(data.shape)
                header_vals['NAXIS'] = sci_header['NAXIS']

                # Use pysiaf to calculate subarray locations
                try:
                    xc, yc = sci_subarray_corners('NIRCam', aperture, siaf=siaf)
                    table_vals['SUBSTRT1'] = xc[0] + 1
                    table_vals['SUBSTRT2'] = yc[0] + 1
                    table_vals['SUBSIZE1'] = siaf[aperture].XSciSize
                    table_vals['SUBSIZE2'] = siaf[aperture].YSciSize
                except KeyError:
                    print("ERROR: Aperture {} is not a valid aperture in pysiaf".format(aperture))
                    xc = [-2, -2]
                    yc = [-2, -2]
                    table_vals['SUBSTRT1'] = xc[0] + 1
                    table_vals['SUBSTRT2'] = yc[0] + 1
                    table_vals['SUBSIZE1'] = 9999
                    table_vals['SUBSIZE2'] = 9999

                # Create FASTAXIS and SLOWAXIS values based on the detector name
                fast, slow = find_fastaxis(header_vals['DETECTOR'])
                table_vals['FASTAXIS'] = fast
                table_vals['SLOWAXIS'] = slow

                # Remove whitespace from observing template in file
                header_vals['TEMPLATE'] = header_vals['TEMPLATE'].replace(' ', '').lower()
                table_vals['TEMPLATE'] = table_vals['TEMPLATE'].lower()

                # Adjust prime/parallel boolean from table to be a string
                if not table_vals['EXPRIPAR']:
                    table_vals['EXPRIPAR'] = 'PRIME'
                else:
                    table_vals['EXPRIPAR'] = 'PARALLEL'

                # Change exposure type from table to match up with
                # types of strings in the file
                table_vals['EXP_TYPE'] = adjust_exptype(table_vals['EXP_TYPE'])

                # Set the DETECTOR field to be identical. This info is not in the
                # exposure table, so we can't actually check it
                table_vals['DETECTOR'] = header_vals['DETECTOR']

                # Compare the actual data shape to the shape given in the header
                header_shape = (header_vals['NINTS'], header_vals['NGROUPS'], header_vals['SUBSIZE2'], header_vals['SUBSIZE1'])
                if header_shape != data_shape:
                    print("WARNING: Shape of data in the file does not match that specified in the header.")
                    print('Data shape: {}'.format(data_shape))
                    print('Header shape: {}'.format(header_shape))

                # Now compare the data in the dictionary from the file versus that
                # from the exposure table created from the APT file
                err = False
                for key in header_vals:
                    if header_vals[key] != table_vals[key]:
                        if key not in FLOAT_KEYWORDS and key not in FILTER_KEYWORDS:
                            err = True
                            print('MISMATCH: {}, in exp table: {}, in file: {}'.format(key, table_vals[key], header_vals[key]))
                        elif key in FLOAT_KEYWORDS:
                            if not np.isclose(header_vals[key], table_vals[key], rtol=0.01, atol=0.):
                                err = True
                                print('MISMATCH: {}, in exp table: {}, in file: {}'.format(key, table_vals[key], header_vals[key]))

                        if key in ['LONGFILTER', 'LONGPUPIL'] and 'LONG' in header_vals['DETECTOR']:
                            err = True
                            print('MISMATCH: {}, in exp table: {}, in file: {}'.format(key, table_vals[key], header_vals[key]))
                        if key in ['SHORTFILTER', 'SHORTPUPIL'] and 'LONG' not in header_vals['DETECTOR']:
                            err = True
                            print('MISMATCH: {}, in exp table: {}, in file: {}'.format(key, table_vals[key], header_vals[key]))

                if not err:
                    print('No inconsistencies. File header info correct.')

            print('')
            print('')

