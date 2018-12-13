
#! /usr/bin/env python
"""

Convert DMS (flight-like) to FITSWriter (ground-like) format.
This code will perform the transformations from flight-like to
ground-like format, which includes taking the data out of the
SCIENCE extension and putting it in the PRIMARY extension.

Authors:
--------
    - Matt Hill
    - Alicia Canipe (updates for NIRCam)

Dependencies:
-------------

    - yaml

Use:
----
    python flight2ground.py config.yml input_file.fits --output_name_base

"""

# Required packages
import argparse
import requests
import numpy as np
from astropy.time import Time
from astropy.io import fits
import yaml


def dms_to_detector(data, detector):
    """Transformations for different instruments from create_data code here:
       https://github.com/spacetelescope/jwst/tree/master/jwst/fits_generator

    Parameters
    ----------
    data : numpy.ndarray
        Data that needs to be flipped

    detector : string
        Detector name

    Returns
    -------
    data : numpy.ndarray
        Flipped exposure in FITSwriter orientation

    """

    if detector == 'NRS1':
        # NRS1 is flipped over the line X=Y
        data = np.swapaxes(data, 2, 3)

    if detector == 'NRS2':
        # NRS2 is flipped over the line Y=X, then rotated 180 degrees
        data = np.swapaxes(data, 2, 3)[:, :, ::-1, ::-1]

    if detector in ['NRCA1', 'NRCA3', 'NRCALONG', 'NRCB2', 'NRCB4']:
        # NRCA1, NRCA3, NRCALONG, NRCB2, NRCB4 are just flipped in X
        data = data[:, :, :, ::-1]

    if detector in ['NRCA2', 'NRCA4', 'NRCB1', 'NRCB3', 'NRCBLONG']:
        # NRCA2, NRCA4, NRCB1, NRCB3, NRCBLONG are just flipped in Y
        data = data[:, :, ::-1]

    if detector == 'NIS':
        # NIRISS has a 180 degree rotation followed by a flip across the line
        # X=Y
        data = np.swapaxes(data[:, :, ::-1, ::-1], 2, 3)

    if detector == 'GUIDER1':
        # GUIDER1 is flipped in X and Y
        data = data[:, :, ::-1, ::-1]

    if detector == 'GUIDER2':
        # GUIDER2 is just flipped in X
        data = data[:, :, :, ::-1]

    # MIRI data doesn't need transforming

    return data


def extract_from_engdb(old_hdulist, new_hdulist, url_base, config):
    """Provide a link to the EngDB and a config file containing the desired
       mnemonics to extract values for headers to be added to the output file.

       e.g.
       http://iwjwdmsbemweb.stsci.edu/JWDMSEngSpAccB71/TlmMnemonicDataSrv.svc/MetaData/TlmMnemonics/


    Parameters
    ----------
    old_hdulist : FITS HDUList
        Old exposure headers

    new_hdulist : FITS HDUList
        New exposure headers

    url_base : string
        Engineering database URL

    config : string
        YAML format configuration file that includes "mnemonics" section

    Returns
    -------
    new_hdulist : FITS HDUList
        Updated FITS HDUList with mnemonic headers added

    """

    # Get the exposure start and end times
    start_time = Time(old_hdulist[0].header['EXPSTART'], format='mjd').isot
    end_time = Time(old_hdulist[0].header['EXPEND'], format='mjd').isot
    params = {'sTime' : start_time, 'eTime' : end_time}

    # Start HTTP request session
    s = requests.Session()

    for keyword, mnemonic in config['mnemonics'].items():

        # Get request to server.
        url = url_base + mnemonic
        r = s.get(url, params=params, verify=False)

        # Parse json
        parsed_json = r.json()

        # json ObsTime has format '/Date(1358619814230+0000)/' which is 1358619814.230 in UNIX time
        # isotime = Time(float(parsed_json['Data'][0]['ObsTime'][6:-7])/1000., format='unix').isot

        # Take the first value of the series (there are no values right now,
        # so use EUType instead of EUValue just to test functionality)
        # new_hdulist[0].header[keyword] = (parsed_json['TlmMnemonics'][0]['EUValue'], mnemonic.upper())
        new_hdulist[0].header[keyword] = (parsed_json['TlmMnemonics'][0]['EUType'], mnemonic.upper())

    # Add the Engineering Mnemonics section heading
    *_, last = old_hdulist[0].header
    new_hdulist[0].header.set('', '', before=list(config['mnemonics'].keys())[0])
    new_hdulist[0].header.set('', 'Engineering Mnemonics', before=list(config['mnemonics'].keys())[0])
    new_hdulist[0].header.set('', '', before=list(config['mnemonics'].keys())[0])

    return new_hdulist


def add_missing_headers(old_hdulist, new_hdulist, config):
    """Add missing keywords for FITSwriter data (e.g., COLCORNR or COLSTART)

    Parameters
    ----------
    old_hdulist : FITS HDUList
        Old exposure headers

    new_hdulist : FITS HDUList
        New exposure headers

    config : string
        YAML format configuration file that includes "headers" section

    Returns
    -------
    new_hdulist : FITS HDUList
        Updated FITS HDUList with new headers added

    """

    # Make sure all strings are capitalized
    for keyword, header in config['headers'].items():
        if isinstance(header, str):
            new_hdulist[0].header[keyword] = header.upper()
        else:
            new_hdulist[0].header[keyword] = header

    # Add the Added Headers section heading
    *_, last = old_hdulist[0].header
    new_hdulist[0].header.set('', '', before=list(config['headers'].keys())[0])
    new_hdulist[0].header.set('', 'Added Headers', before=list(config['headers'].keys())[0])
    new_hdulist[0].header.set('', '', before=list(config['headers'].keys())[0])

    return new_hdulist


def main(args):
    """Main function.

    """

    # Load the configuration file
    with open(args.config_file, 'r') as cfgfile:
        config = yaml.load(cfgfile)

    # Load the HDUList for the input file
    old_hdulist = fits.open(args.input_file)

    # JWST ENGINEERING DATABASE
    # If there are mnemonics headers to add, then add them
    if 'mnemonics' in config:

        # Create a new HDUList for the output file and add the old headers
        new_hdulist = fits.HDUList()
        new_hdulist.append(fits.PrimaryHDU())
        new_hdulist[0].header = old_hdulist[0].header

        # Get the exposure start and end times
        start_time = Time(old_hdulist[0].header['EXPSTART'], format='mjd').isot
        end_time = Time(old_hdulist[0].header['EXPEND'], format='mjd').isot
        params = {'sTime' : start_time, 'eTime' : end_time}

        # Link to the JWST Engineering Database
        eng_db_url = 'http://iwjwdmsbemweb.stsci.edu/JWDMSEngSpAccB71/TlmMnemonicDataSrv.svc/MetaData/TlmMnemonics/'

        # Extract mnemonics from the database and add them to the HDUList
        new_hdulist = extract_from_engdb(old_hdulist, new_hdulist, eng_db_url, config)

    # NEW HEADERS
    # If there are new headers to add, then add them
    if 'headers' in config:
        new_hdulist = add_missing_headers(old_hdulist, new_hdulist, config)

    # Transform data from DMS to detector orientation
    pixel_data = dms_to_detector(old_hdulist['SCI'].data, old_hdulist['PRIMARY'].header['DETECTOR'])

    # Collapse data shape from 4D to 3D
    nints, ngroups, nx, ny = pixel_data.shape

    # Add reference output for MIRI if necessary
    if old_hdulist['PRIMARY'].header['INSTRUME'] == 'MIRI':
        new_hdulist[0].data = np.append(old_hdulist['SCI'].data.reshape((nints*ngroups, nx, ny)),
            old_hdulist['REFOUT'].data.reshape((nints*ngroups, 256, 1032)), axis=1)
    else:
        new_hdulist[0].data = pixel_data.reshape((nints*ngroups, nx, ny))

    # Remove the NEXTEND keyword since there is only one extension now
    try:
        new_hdulist[0].header.remove('NEXTEND')
    except:
        print('No NEXTEND header, skipping.')

    # Write out the reformatted file
    if args.output_name is not 'None':
        print('Saving file: ', args.output_name+".fits")
        new_hdulist.writeto(args.output_name+".fits", overwrite=True)
    else:
        print('Saving file: ', args.input_file[:-5]+"_fitswriter.fits")
        new_hdulist.writeto(args.input_file[:-5]+"_fitswriter.fits", overwrite=True)


if __name__ == '__main__':

    # Command line argument handler.
    parser = argparse.ArgumentParser(
        description='Convert JWST data from DMS format to FITSWriter format',
        epilog='example: flight2ground.py config.yml input_file.fits --output_name_base')

    parser.add_argument('config_file',
                        help='config file with Telemetry FITS keyword/mnemonic pairs and headers to add.')
    parser.add_argument('input_file',
                        help='level 1B data file to reformat')
    parser.add_argument('--output_name',
                        help='level 1B data file to reformat', default='None')

    args = parser.parse_args()
    main(args)
