#!/usr/bin/env python

"""Create a CRDS-formatted saturation reference files for persistence step
from a list of exposures.

This script is constructed with the help of B. Hilbert's
make_NIRCAM_dark_reffile.py script to make dark reference files. Some functions
to read in the data and save output as a fits file were taken directly from the
script:

    read_listfile.py
    input_consistency_check.py
    save_reffile.py
    make_proper_header.py

The steps in this script are:

    1.) read in a list of uncalibrated exposures to use
    2.) use signal difference between consecutive reads to find saturation
    3.) take average of saturated reads as hard saturation
    4.) when final saturation levels are found, average together values
        for each input exposure (masking out bad pixels, weird pixels, etc.)
        and find standard deviation
    5.) save saturation values to a FITS file using the Saturation model
    6.) save errors and other values to FITS files for checking later

The output files are:

    1.) SSB formatted FITS file with averaged saturation values in the SCI
        extension, along with corresponding DQ and DQ_DEF extensions
    2.) FITS file with errors on saturation values (PersistenceSatModel doesn't have
        ERR extension yet?).
    3.) FITS file with saturation values for each input exposure in the listfile.
    4.) FITS file containing first saturated group values for each exposure.

Version 1.0 (Created 2017) - A. Canipe

"""

import argparse
import sys
import itertools
import datetime
import numpy as np
from scipy import stats
from astropy.io import fits
from astropy.stats import sigma_clip
from jwst.dq_init import DQInitStep
from jwst.superbias import SuperBiasStep
from jwst.refpix import RefPixStep
from jwst.datamodels import PersistenceSatModel

DET_NUM = {'NRCA1': 17004,
           'NRCA2': 17006,
           'NRCA3': 17012,
           'NRCA4': 17048,
           'NRCALONG': 17158,
           'NRCB1': 16991,
           'NRCB2': 17005,
           'NRCB3': 17011,
           'NRCB4': 17047,
           'NRCBLONG': 17161}

class MakeSatRef:
    '''Class to create the reference files.'''


    def __init__(self):
        self.verbose = False


    def add_options(self, parser=None, usage=None):
        '''Adds in command line arguments.'''

        if parser is None:
            parser = argparse.ArgumentParser(usage=usage)
            parser.add_argument("listfile",
                                help="List of uncalibrated exposures to use.")
            parser.add_argument("--intermediates",
                                default=False,
                                action="store_true",
                                help="Save output for pipeline steps.")

        return parser


    def read_listfile(self, myfile):
        '''Read in a listfile and return list of files.'''

        files = []
        with open(myfile) as fil:
            for line in fil:
                if len(line) > 2:
                    files.append(line.strip())

        return files


    def input_consistency_check(self, files):
        '''Check basic input file attributes.'''

        readpatt = set()
        ngroups = set()
        n_x = set()
        n_y = set()
        pupil = set()
        detector = set()
        instrument = set()

        for file in files:
            with fits.open(file) as fil:
                head0 = fil[0].header

            pupil.add(head0['PUPIL'].strip())
            instrument.add(head0['INSTRUME'].strip())
            detector.add(head0['DETECTOR'].strip())
            readpatt.add(head0['READPATT'].strip())
            ngroups.add(head0['NGROUPS'])
            n_x.add(head0['SUBSIZE1'])
            n_y.add(head0['SUBSIZE2'])

        if len(pupil) > 1:
            print("WARNING! MORE THAN ONE PUPIL SETTING USED FOR INPUT DATA!!")
            sys.exit(0)

        if len(instrument) > 1:
            print("WARNING! MORE THAN ONE INSTRUMENT USED FOR INPUT DATA!!")
            sys.exit(0)

        if len(detector) > 1:
            print("WARNING! MORE THAN ONE DETECTOR USED FOR INPUT DATA!!")
            sys.exit(0)

        if len(readpatt) > 1:
            print("WARNING! MORE THAN ONE READ PATTERN USED FOR INPUT DATA!!")
            sys.exit(0)

        if list(readpatt)[0] != 'RAPID':
            print("WARNING! Input data don't use RAPID read pattern! Quitting")
            sys.exit(0)

        if len(ngroups) > 1:
            print("WARNING! MORE THAN ONE NUMBER OF GROUPS IN INPUT DATA!!")
            print("Okay because saturation step discards saturated frames.")
            # sys.exit(0)

        if len(n_x) > 1:
            print("WARNING! MORE THAN ONE ARRAY X-SIZE IN INPUT DATA!!")
            sys.exit(0)

        if len(n_y) > 1:
            print("WARNING! MORE THAN ONE ARRAY Y-SIZE IN INPUT DATA!!")
            sys.exit(0)

        return head0['DETECTOR'], head0['SUBSIZE1'], head0['SUBSIZE2']


    def init_arrays(self, files, xshape, yshape):
        '''Initialize arrays to collect output from steps.'''

        sat_vals = np.zeros((len(files), yshape, xshape), dtype='float32')
        grp_arr = np.zeros((len(files), yshape, xshape), dtype='float32')
        new_dq = np.zeros((yshape, xshape), dtype='int32')
        big_dq = np.zeros((len(files), yshape, xshape), dtype='int32')
        tmp_mask = np.zeros((len(files), yshape, xshape), dtype='bool')

        return sat_vals, grp_arr, new_dq, big_dq, tmp_mask


    def get_lin_regime(self, sig, method):
        '''Find approximate linear regime.'''

        max_signal = np.nanmax(sig)

        # Method 1: fit a line to the first part of the ramp to determine
        # where the ramp is roughly linear (slower?)
        if method == "method1":
            x = np.arange(0, len(sig))
            early = np.argmin(sig < np.nanmax(sig)*(3**-1))
            if early == 0:
                early = int(20)
            slope, intercept, r_value, p_value, std_err = \
                                        stats.linregress(x[:early], sig[:early])
            true = intercept + slope * x
            ratio = true * sig**-1
            regime = np.argwhere((ratio > 0.995) & (ratio < 1.005))

        # Method 2: grab a chunk of the ramp before significant nonlinearity
        elif method == "method2":
            early = float(max_signal/5)
            mid = float(early*2)
            regime = np.argwhere((sig > early) & (sig < mid))

        return regime


    def get_saturation(self, signal, signal_range):
        '''Find hard saturation using signal differences.'''

        # Determine saturation by making array of signal differences between
        # consecutive frames. Signal differences after saturation will be
        # close to zero.
        sig_diff = np.diff(signal)

        # Choose 3-sigma region around mean of zero to find saturation limit.
        three_sigma = 3 * np.std(sig_diff[int(signal_range[0]):int(signal_range[-1])])
        upperbound = three_sigma
        lowerbound = -1 * three_sigma
        first_sat_grp = np.argmax((lowerbound < sig_diff) &
                                  (sig_diff < upperbound))

        # Hard saturation limit = average signal level of saturated frames.
        hard_sat = np.average(signal[first_sat_grp:])

        return hard_sat, first_sat_grp


    def calc_stats(self, sat_arr, mask):
        '''Calculate statistics for values.'''

        # use astropy's sigma-clipping function to mask bad pixels
        clipped = sigma_clip(sat_arr, axis=0)
        clipped.data[clipped.mask] = np.nan
        clipped.data[(mask == 1)] = np.nan

        # pixels have nans where saturation couldn't be determined
        averages = np.nanmean(clipped.data, axis=0)
        errs = np.nanstd(clipped.data, axis=0)

        return averages, errs


    def make_proper_header(self, model, files, output):
        '''Make sure the reference file has the correct headers.'''

        # Read in headers from one of the input files.
        with fits.open(files[0]) as fil:
            header0 = fil[0].header

        # Change the header values.
        model.meta.author = 'Canipe'
        meta.date = str(datetime.datetime.now())
        model.meta.description = 'Saturation reference file from CV3 data'
        model.meta.filename = str(output)
        model.meta.filetype = 'REFERENCE'
        model.meta.reftype = 'PERSAT'
        model.meta.model_type = 'PersistenceSatModel'
        model.meta.origin = 'STSCI'
        model.meta.pedigree = 'DUMMY'
        model.meta.useafter = '2017-11-29'
        model.meta.telescope = 'JWST'
        model.meta.time_sys = 'UTC'
        model.meta.subarray.xsize = header0['SUBSIZE1']
        model.meta.subarray.ysize = header0['SUBSIZE2']
        model.meta.subarray.xstart = header0['SUBSTRT1']
        model.meta.subarray.ystart = header0['SUBSTRT2']
        model.meta.subarray.name = header0['SUBARRAY']
        model.meta.instrument.name = 'NIRCAM'
        detector = header0['DETECTOR']
        if detector == 'NRCA5':
            detector = 'NRCALONG'
        if detector == 'NRCB5':
            detector = 'NRCBLONG'
        model.meta.instrument.detector = detector

        try:
            model.meta.subarray.fastaxis = header0['FASTAXIS']
            model.meta.subarray.slowaxis = header0['SLOWAXIS']
        except KeyError:
            print('===============================================')
            print("FASTAXIS and SLOWAXIS header keys not found in input data.")
            print("Assuming they are in native (fitswriter) orientation, and")
            print("adding native orientation values for those keywords.")
            print('===============================================')
            model.meta.subarray.fastaxis = 1
            model.meta.subarray.slowaxis = 2


        # Add HISTORY (this needs to be edited still).
        model.history.append('Description of Reference File Creation')
        model.history.append('Hard saturation was constructed from long ramps')
        model.history.append('that reach hard saturation.')
        model.history.append('Steps:')
        model.history.append('1. For each file containing long data ramps:')
        model.history.append('     DQ Initialization')
        model.history.append('     Superbias subtraction')
        model.history.append('     Reference pixel correction')
        model.history.append('2. Use delta(signal) to find hard saturation,')
        model.history.append('   since delta(signal) ~ 0 when saturated.')
        model.history.append('3. Average samples at hard saturation.')
        model.history.append('4. Calculate sigma-clipped mean saturation level')
        model.history.append('   of input ramps.')
        model.history.append('SOFTWARE:')
        model.history.append('https://github.com/spacetelescope/nircam_calib')

        # Put the list of input files into the HISTORY keyword.
        model.history.append('DATA USED:')
        for file in files:
            totlen = len(file)
            div = np.arange(0, totlen, 60)
            for val in div:
                if totlen > (val+60):
                    model.history.append(file[val:val+60])
                else:
                    model.history.append(file[val:])

        model.history.append('DIFFERENCES:')
        model.history.append('N/A. No previous version.')

        return model


    def save_reffile(self, sat, err, passedmap, files, output):
        '''Save the reference file using JWST data models.'''

        # Use jwst.datamodels PersistenceSatModel to put data in correct FITS format.
        finalsaturation = PersistenceSatModel()
        finalsaturation.data = sat

        # Need to add errors but PersistenceSatModel doesn't have ERR extension yet.
        finalsaturation.err = err

        # Create DQ flag definition table.
        dqdef = []
        dqssb_lin = {'DO_NOT_USE':np.uint8(1), 'NO_SAT_CHECK':np.uint8(2),
                     'WEIRD_PIXEL':np.uint8(3)}
        dqssb_desc = {'Bad pixel, do not use.', 'Unable to determine saturation.',
                      'Values may not be trustable.'}
        for bitname, bitdescription in zip(dqssb_lin, dqssb_desc):
            bitvalue = dqssb_lin[bitname]
            bitnumber = int(np.log(bitvalue)/np.log(2))
            newrow = (bitnumber, bitvalue, bitname, bitdescription)
            dqdef.append(newrow)
            dqmap = passedmap

        # Insert dq array into saturation reffile.
        finalsaturation.dq = dqmap
        finalsaturation.dq_def = dqdef

        # Create the proper headers.
        finalsaturation = self.make_proper_header(finalsaturation, files, output)

        # Save the reference file.
        outfile = output
        finalsaturation.save(outfile)

        return outfile


    def run(self):
        '''Main function.'''

        # Read in list of files to use to calculate saturation.
        files = self.read_listfile(self.listfile)

        # Check that input files have the same readpattern, array size, etc.
        detector, xshape, yshape = self.input_consistency_check(files)
        det_short = str(detector[3:])
        if 'long' in det_short:
            det_short = det_short[0] + '5'

        # Set up arrays to hold output data.
        sat_arr, grp_arr, new_dq, big_dq, tmp_mask = \
                                        self.init_arrays(files, xshape, yshape)

        # Create mask to isolate reference pixels
        tmp_mask[:, 4:-4, 4:-4] = True

        # Set reference pixel values so they don't get used in calculations.
        sat_arr[~tmp_mask] = 1e6
        grp_arr[~tmp_mask] = np.nan
        big_dq[~tmp_mask] = np.uint8(0)

        # Loop over files.
        for file, n in zip(files, np.arange(0, len(files))):

            # Run dq, superbias, refpix steps and save outputs (optional)
            bpm = DQInitStep.call(file)
            sup = SuperBiasStep.call(bpm)
            ref = RefPixStep.call(sup, odd_even_rows=False)

            if self.intermediates:
                sup.save(file[:-5]+'_dq_superbias.fits')
                ref.save(file[:-5]+'_dq_superbias_refpix.fits')

            # Grab the name of the mask file used from the headers
            bpmcalfile = ref.meta.ref_file.mask.name
            if 'crds' in bpmcalfile:
                jwst = bpmcalfile.find('jwst')
                bpmfile = '/grp/crds/cache/references/jwst/'+bpmcalfile[jwst:]
            else:
                bpmfile = bpmcalfile

            # Get data values
            mask = fits.getdata(bpmfile, 1)
            data = ref.data
            xstart = ref.meta.subarray.xstart
            ystart = ref.meta.subarray.ystart
            xend = ref.meta.subarray.xsize
            yend = ref.meta.subarray.ysize

            # Loop over pixel combinations for given array (no ref pixels).
            for i, j in itertools.product(np.arange(xstart+3, xend-4),
                                          np.arange(ystart+3, yend-4)):

                # Set values for bad pixels so they don't get used in calculations.
                if mask[j, i] == np.uint8(1):

                    sat_arr[n, j, i] = np.nan
                    grp_arr[n, j, i] = np.nan
                    big_dq[n, j, i] = np.uint8(1)

                else:

                    # Get signal values for each pixel.
                    signal = data[0, :, j, i].astype('float32')

                    # Get linear region early in the ramp with method 1
                    signal_range = self.get_lin_regime(signal, "method1")

                    # If signal_range can't be determined, must be weird ramp
                    if np.shape(signal_range)[0] == 0:

                        # Try again to get linear region with different method
                        signal_range = self.get_lin_regime(signal, "method2")

                        # If that still doesn't work, quit.
                        if np.shape(signal_range)[0] == 0:

                            sat_arr[n, j, i] = np.nan
                            grp_arr[n, j, i] = np.nan
                            big_dq[n, j, i] = np.uint8(2)

                        else:
                            hard_sat, first_sat_grp = \
                                self.get_saturation(signal, signal_range)

                            # Save all values.
                            sat_arr[n, j, i] = hard_sat.astype('float32')
                            grp_arr[n, j, i] = first_sat_grp.astype('float32')
                            big_dq[n, j, i] = np.uint8(3)

                    # Otherwise, must be good ramp?
                    elif np.shape(signal_range)[0] > 0:

                        # Get hard saturation.
                        hard_sat, first_sat_grp = \
                            self.get_saturation(signal, signal_range)

                        # Save all saturation values.
                        sat_arr[n, j, i] = hard_sat.astype('float32')
                        grp_arr[n, j, i] = first_sat_grp.astype('float32')

                    # Catch errors
                    else:
                        print('ERROR for pixel ', i, j)
                        sys.exit(0)

        # If each file gave same pixel DQs, make sure output DQ matches
        locs = np.all(big_dq == big_dq[0, :], axis=0)
        new_dq[locs] = big_dq[0][locs]

        # Get statistics for saturation values, averaging over exposures.
        avg, err = self.calc_stats(sat_arr, big_dq)

        # Save saturation values for each exposure to a FITS file
        newhdu = fits.PrimaryHDU(sat_arr)
        newhdulist = fits.HDUList([newhdu])
        grpname = detector + '_'  \
                      + str(DET_NUM[detector])  \
                      + '_WellDepthADU_'  \
                      + str(datetime.date.today())  \
                      + '_beforeAverage.fits'
        newhdulist.writeto(grpname, overwrite=True)

        # Save first saturated groups array to a FITS file
        newhdu = fits.PrimaryHDU(grp_arr)
        newhdulist = fits.HDUList([newhdu])
        grpname = detector + '_'  \
                      + str(DET_NUM[detector])  \
                      + '_WellDepthADU_'  \
                      + str(datetime.date.today())  \
                      + '_firstSatGroup.fits'
        newhdulist.writeto(grpname, overwrite=True)

        # Save averaged saturation values to a FITS file.
        outfilename = detector + '_'  \
                          + str(DET_NUM[detector])  \
                          + '_WellDepthADU_'  \
                          + str(datetime.date.today())  \
                          + '_ssbsaturation_DMSorient.fits'
        outfile = self.save_reffile(avg, err, new_dq, files, outfilename)

        # Save saturation errors, since the pipeline doesn't currently use them
        errhdu = fits.PrimaryHDU(err)
        errhdulist = fits.HDUList([errhdu])
        errname = detector + '_'  \
                          + str(DET_NUM[detector])  \
                          + '_WellDepthADU_'  \
                          + str(datetime.date.today())  \
                          + '_saturationErrors.fits'
        errhdulist.writeto(errname, overwrite=True)

        z = fits.open(outfile)
        z0 = z[0]
        z1 = z[1]
        z2 = z[2]
        z3 = z[3]

        # Add other things that the pipeline doesn't use, but are helpful.
        z0.header['S_DQINIT'] = ('COMPLETE', 'Data Quality Initialization')
        z0.header['S_SUPERB'] = ('COMPLETE', 'Superbias Subtraction')
        z0.header['S_REFPIX'] = ('COMPLETE', 'Reference Pixel Correction')
        newhdu = fits.HDUList([z0, z1, z2, z3])
        newhdu.writeto(outfile, overwrite=True)



if __name__ == '__main__':
    usagestring = 'USAGE: make_NIRCam_saturation_reffile.py listfile.list --intermediates'

    starting = datetime.datetime.now()

    satcalc = MakeSatRef()
    parse = satcalc.add_options(usage=usagestring)
    args = parse.parse_args(namespace=satcalc)

    satcalc.run()

    ending = datetime.datetime.now()
    difftime = ending - starting
    print("DONE. Elapsed time: {}".format(difftime))
