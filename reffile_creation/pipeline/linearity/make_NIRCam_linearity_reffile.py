#!/usr/bin/env python

"""Create a CRDS-formatted linearity correction reference file
from a list of exposures.

This script is a work-in-progress. It is constructed with the help of B.
Hilbert's make_NIRCAM_dark_reffile.py script to make dark reference files.
Some functions to read in the data and save output as a fits file were
taken directly from the script:

    read_listfile.py
    input_consistency_check.py
    save_reffile.py
    make_proper_header.py
    multi_step_sigma_clipped_mean.py

NOTE: This script takes a while to run because it performs several
processes. Improvements to speed are being made.

The steps are:

    1.) read in a list of exposures to use and a parameter file
    2.) determine hard saturation levels to figure out how much of the
        ramps to use
    3.) fit the ramps to find ideal slope and intercept (bias) for
        correction curves used to find coefficients
    4.) fit correction curves up to saturation to find coefficients -- if
        residuals between corrected signals and ideal linear signals are
        large, subtract a frame and fit again until soft saturation limit
        is reached
    5.) when final coefficients are found, average together coefficients
        for each input exposure (masking out bad pixels, weird pixels, fit
        error pixels, etc) and find standard deviation
    6.) save coefficients to a FITS file using the Linearity model and save
        saturation values to separate FITS file for analysis


Version 1.0 (Created 2016) - A. Canipe

"""

import argparse
import warnings
import sys
import yaml
import itertools
import datetime
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
from jwst.datamodels import LinearityModel


class MakeNonlinRef:
    '''Class to create the reference files.'''


    def __init__(self):
        self.verbose = False


    def add_options(self, parser=None, usage=None):
        '''Adds in command line arguments.'''

        if parser is None:
            parser = argparse.ArgumentParser(usage=usage)
            parser.add_argument("paramfile",
                                help="Parameter file for function steps.")
        return parser


    def load_params(self, paramfile):
        '''Read in yaml file containing parameters.'''

        try:
            with open(paramfile, 'r') as infile:
                params = yaml.load(infile)
        except:
            print("WARNING: unable to open {}".format(paramfile))
            sys.exit()

        return params


    def read_listfile(self, myfile):
        '''Read in a listfile and return list of files.'''

        files = []
        with open(myfile) as fil:
            for line in fil:
                if len(line) > 2:
                    files.append(line.strip())

        return files


    def input_consistency_check(self, files):
        '''Check basic input file attributes to make sure everything is
        consistent.'''

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


    def init_arrays(self, files, num_coefs, coords):
        '''Initialize arrays to collect output from steps.'''

        coeffs_arr = np.zeros((len(files), num_coefs, coords[3], coords[2]))
        dq_arr = np.zeros((coords[3], coords[2]))
        sat_arr = np.zeros((len(files), 4, coords[3], coords[2]))
        mask_arr = np.zeros(np.shape(coeffs_arr), dtype='bool')

        return coeffs_arr, dq_arr, sat_arr, mask_arr


    def load_data(self, data):
        '''Load the data and define variables for other functions to use.'''

        with fits.open(data) as fil:
            hdr = fil[0].header
            data = fil[1].data
            tgroup = hdr['TGROUP']
            ngroups = hdr['NGROUP']
            skip = hdr['DRPFRMS3']
            nframe = hdr['NFRAME']

        return data, tgroup, ngroups, skip, nframe


    def get_saturation(self, signal, sigma_range, limit):
        '''Find hard saturation for each ramp.'''

        # Determine saturation by making array of signal differences between
        # consecutive frames. Signal differences after saturation will be
        # close to zero.
        sig_diff = np.diff(signal)

        # Choose 3-sigma region around mean of zero to find saturation limit.
        three_sigma = limit * np.std(sig_diff[sigma_range[0]:sigma_range[1]])
        upperbound = three_sigma
        lowerbound = -1 * three_sigma
        first_sat_grp = np.argmax((lowerbound < sig_diff) &
                                  (sig_diff < upperbound))

        # Hard saturation limit = average signal level of saturated frames.
        hard_sat = np.average(signal[first_sat_grp:])

        return hard_sat, first_sat_grp


    def ramp_fit(self, pixel, zeroframe, tgroup, ramp, poly_order, unc_ron):
        '''Function to fit the ramp and determine ideal linear signal.'''

        # Define minimum frame to include in ramp fitting (either 0 or 1, 0
        # gives a larger bias value). Make sure arrays are floats.
        frame_cut = ramp[0][zeroframe:]
        new_ramp = [np.asfarray(frame_cut*tgroup),
                    np.asfarray(ramp[1][zeroframe:])]

        # Define polynomial function of the true signal rate for different
        # polynomial orders. Here, 'b' will be constant signal rate for
        # perfect detector and 'c', 'd', 'e', etc. describe magnitudes of
        # deviations from ideal signal.
        if poly_order == 5:
            def func(t, a, b, c, d, e, f):
                bt = b*t
                # Using Horner's method to optimize speed, so this polynomial:
                # a+ bt+ c*(bt)**2+ d*(bt)**3+ e*(bt)**4+ f*(bt)**5
                # becomes:
                return a + bt*(1 + bt*(c + bt*(d + bt*(e + bt*f))))

        elif poly_order == 4:
            def func(t, a, b, c, d, e):
                bt = b*t
                return a + bt*(1 + bt*(c + bt*(d + bt*e)))

        else:
            print('Ramp fit polynomial order not available. Quitting.')
            sys.exit(0)

        # Initial values for curve-fitting.
        initvals = np.array([100, 4.e+01, 7.e-03, -6.e-06, 3.e-08, -8.e-11])

        # Provide uncertainties in the signal, either an estimate or an
        # array of signal uncertainties.
        if isinstance(unc_ron, int) is True or isinstance(unc_ron, float) is True:
            unc = np.sqrt(unc_ron + new_ramp[1])
        elif isinstance(unc_ron, np.ndarray) is True:
            unc = np.asfarray(unc_ron[zeroframe:len(frame_cut), pixel[1],
                              pixel[0]])

        # Fit data with scipy's built in curve_fit function, a curve-fitting
        # wrapper around scipy's optimize.leastsq function that minimizes the
        # sum of squares of a set of equations (optimization with
        # Levenberg-Marquardt algorithm).
        # inputs:  ramp-fitting function defined above
        #          frame and signal values
        #          sigma is array of signal uncertainties
        #          absolute_sigma = False, sigma denotes relative weights of
        #                           data, covariance matrix is based on
        #                           estimated errors in the data and is not
        #                           affected by the overall magnitude of sigma.
        #                           Only relative magnitudes of sigma matter.
        #                         = True, sigma describes 1 standard deviation
        #                           errors of input data Estimated covariance
        #                           is based on these values.
        #          p0 are initial values
        #          maxfev is max number of iterations to find convergence
        # outputs: beta will be coefficients for hybrid polynomial, so
        #               beta[1] is our "ideal ramp rate" and b[0] is bias
        #          pcov is covariance matrix
        ramp_coefs, ramp_cov = curve_fit(func, new_ramp[0], new_ramp[1],
                                         sigma=unc, absolute_sigma=False,
                                         p0=initvals[:poly_order+1],
                                         maxfev=20000)

        # Compute 1 standard deviation errors on parameters as
        # perr = np.sqrt(np.diag(ramp_cov)).
        fit_unc = [np.sqrt(ramp_cov[j,j]) for j in
                   range(ramp_coefs.size)]
        fit_cov = ramp_cov[1,0]

        # Get ideal ramp slope and intercept.
        new_ideal_rate = ramp_coefs[1]
        intercept = ramp_coefs[0]

        # Ideal signal is y = b + a * t.
        ideal = intercept + new_ideal_rate * new_ramp[0]
        ideal_uncsq = fit_unc[0]**2 + (fit_unc[1]*new_ramp[0])**2+ 2.*fit_cov**2
        ideal_unc = np.sqrt(ideal_uncsq)

        # Get fractional nonlinearity to produce the correction curve.
        nonlin = (ideal * new_ramp[1]**(-1))
        non_uncsq = nonlin**2*(unc**2*new_ramp[1]**(-2)+ideal_unc**2*ideal**(-2))
        nonlin_unc = np.sqrt(non_uncsq)

        return new_ramp, new_ideal_rate, intercept, nonlin, nonlin_unc


    def correct_signal(self, new_ramp, nonlin, nonlin_unc, intercept,
                      new_ideal_rate, first_sat_grp, nonlin_accuracy, signal,
                      tgroup, corr_order, ramp_percent):
        '''Function to correct signal using nonlinearity curve from
            ramp_fit().'''

        # Define new function to fit to nonlinearity correction curve.
        # (curve = nonlinearity vs. signal)
        # aa, bb, cc, dd, ee, ff are correction coefficients.
        if corr_order == 5:

            def poly_fit_func(t, bb, cc, dd, ee, ff):
                # with Horner's method again,
                # 1 + bb*t + cc*t**2 + dd*t**3 + ee*t**4 + ff*t**5 is:
                return 1 + t*(bb + t*(cc + t*(dd + t*(ee + t*ff))))

        elif corr_order == 4:
            def poly_fit_func(t, bb, cc, dd, ee):
                return 1 + t*(bb + t*(cc + t*(dd + t*ee)))

        else:
            print('Correction polynomial order not available. Quitting.')
            sys.exit(0)

        # Initial values for curve-fitting.
        nonlin_init = np.array([3e-6, -1e-10, 2e-15, 1e-18, 1e-21])

        # Fit the correction curve with curve_fit as before.
        # beta contains correction coefficients
        # pcov contains uncertainties in coefficients
        beta, pcov = curve_fit(poly_fit_func, new_ramp[1], nonlin,
                               sigma=nonlin_unc, p0=nonlin_init[:corr_order],
                               absolute_sigma=False, maxfev=20000)
        # poly_fit_unc = [np.sqrt(pcov[j,j]) for j in range(beta.size)]
        beta = np.asfarray(beta, dtype='float')
        nonlin_accuracy = np.float(nonlin_accuracy)

        # # Get curve fit residuals.
        # curve_fit_resid = (nonlin-poly_fit_func(new_ramp[1],*beta))*nonlin**(-1)

        # Correct signal using coefficients.
        if corr_order == 5:
            # c = y * (1 + b1*y + b2*y**2 + b3*y**3 + b4*y**5 + b5*y**6)
            corrected_signal = new_ramp[1] * (1 + new_ramp[1] * (beta[0]
                                              + new_ramp[1] * (beta[1]
                                              + new_ramp[1] * (beta[2]
                                              + new_ramp[1] * (beta[3]
                                              + new_ramp[1] * beta[4])))))

        elif corr_order == 4:
            corrected_signal = new_ramp[1] * (1 + new_ramp[1] * (beta[0]
                                              + new_ramp[1] * (beta[1]
                                              + new_ramp[1] * (beta[2]
                                              + new_ramp[1] * beta[3]))))

        else:
            print('Correction polynomial order not available. Quitting.')
            sys.exit(0)


        # Find correction accuracy, stdev(corrected - ideal / measured).
        corr_ideal = intercept + new_ideal_rate * new_ramp[0]
        nonnorm_resid = np.subtract(corrected_signal,corr_ideal)
        correction_resid =  nonnorm_resid*new_ramp[1]**(-1)
        correction_resid_std = np.std(correction_resid)

        # Place limit on how much of ramp is cut off to find soft saturation
        # (e.g. not more than 20%)
        cut_limit = np.int(ramp_percent*first_sat_grp)

        # If accuracy is worse than desired soft saturation limit (accuracy),
        # this is the loop to subtract last frame and try again.
        for i in np.arange(0, 100):

            if correction_resid_std > nonlin_accuracy and i < cut_limit:

                # Subtract a group from saturation until desired accuracy
                # level is reached. Cut data from previous steps for new limit.
                test_lim = first_sat_grp-(i+1)
                test_frame = new_ramp[0][:test_lim]
                test_sig = new_ramp[1][:test_lim]
                test_nonlin = nonlin[:test_lim]
                test_nonlin_unc = nonlin_unc[:test_lim]

                # Make sure data is in correct format for curve_fit function.
                test_sig = np.asfarray(test_sig, dtype='float')
                test_nonlin = np.asfarray(test_nonlin, dtype='float')
                test_nonlin_unc = np.asfarray(test_nonlin_unc, dtype='float')

                # Fit the nonlinearity correction curve again.
                beta, pcov = curve_fit(poly_fit_func, test_sig, test_nonlin,
                                       sigma=test_nonlin_unc, p0=nonlin_init,
                                       absolute_sigma=False, maxfev=10000)
                # poly_fit_unc = [np.sqrt(pcov[j,j]) for j in range(beta.size)]
                beta = np.asfarray(beta, dtype='float')

                # # Get new curve fit residuals.
                # curve_fit_resid = (test_nonlin-poly_fit_func(test_sig,*beta)
                #                     )/test_nonlin

                # Crrect signal using new coefficients.
                if corr_order == 5:
                    corrected_signal = test_sig * (1 + test_sig * (beta[0]
                                                      + test_sig * (beta[1]
                                                      + test_sig * (beta[2]
                                                      + test_sig * (beta[3]
                                                      + test_sig * beta[4])))))
                elif corr_order == 4:
                    corrected_signal = test_sig * (1 + test_sig * (beta[0]
                                                      + test_sig * (beta[1]
                                                      + test_sig * (beta[2]
                                                      + test_sig * beta[3]))))

                else:
                    print('Correction polynomial order not valid. Quitting.')
                    sys.exit(0)

                # Find new correction accuracy.
                corr_ideal = intercept + new_ideal_rate * test_frame
                nonnorm_resid = np.subtract(corrected_signal,corr_ideal)
                correction_resid = nonnorm_resid*test_sig**(-1)
                correction_resid_std = np.std(correction_resid)

            # Prevent loop from cutting off too much of the ramp to find soft
            # saturation using cut_limit.
            # If no soft saturation can be found, flag pixel with NaN value.
            elif correction_resid_std > nonlin_accuracy and i >= cut_limit:
                upperlimit = first_sat_grp
                soft = np.nan
                break

            # Stop loop if soft saturation is reached.
            elif correction_resid_std < nonlin_accuracy and i <= cut_limit:

                # If upperlimit was changed to improve accuracy, get new limit.
                if 'test_lim' in locals():
                    upperlimit = test_lim
                    soft = signal[upperlimit]

                # If no new limit needs to be found, just get the max frame #.
                else:
                    upperlimit = np.int(np.max(new_ramp[0]*tgroup**(-1)))
                    soft = signal[upperlimit]
                break

            else:
                upperlimit = np.nan
                soft = np.nan
                break

        return beta, correction_resid_std, soft, upperlimit


    def calc_stats(self, coeffs_arr):
        '''Calculate statistics for correction coefficients.'''

        # coeffs_arr probably has nan values where coefficients couldn't be
        # determined, so use "np.nanmean" to take average while ignoring nans.
        averages = np.nanmean(coeffs_arr, axis=0)
        errs = np.nanstd(coeffs_arr, axis=0)

        return averages, errs


    def make_proper_header(model, files):
        '''Make sure the reference file has the correct headers.'''

        # Read in headers from one of the input files.
        with fits.open(files[0]) as fil:
            header0 = fil[0].header
            header1 = fil[1].header

        # Change the header values.
        model.meta.reffile.type = 'LINEARITY'
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

        model.meta.reffile.author = 'Canipe'
        model.meta.reffile.description = 'Linearity reffile from CV3 data'
        model.meta.reffile.pedigree = 'DUMMY'
        model.meta.reffile.useafter = '2016-1-15'

        # Add HISTORY (this needs to be edited still).
        model.history.append('Description of Reference File Creation')
        model.history.append('DOCUMENT:')
        model.history.append('JWST-STScI-TR-XXXX')
        model.history.append('SOFTWARE:')
        model.history.append('/ifs/jwst/wit/witserv/data7/nrc/')
        model.history.append('/linearity/make_NIRCAM_linearity_reffile.py')

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


    def save_reffile(self, lin, err, passedmap, files, output):
        '''Save the reference file using JWST data models.'''

        # Use jwst.datamodels LinearityModel to put data in correct FITS format.
        finallinearity = LinearityModel()
        finallinearity.coeffs = lin

        # Need to add errors but LinearityModel doesn't have ERR extension yet.
        finallinearity.err = err

        # Create DQ flag definition table.
        dqdef = []
        dqssb_lin={'DO_NOT_USE':np.uint8(1),'NONLINEAR':np.uint8(2),
                  'NO_LIN_CORR':np.uint8(4)}
        dqssb_desc={'Bad pixel, do not use.','Highly nonlinear pixel.',
                  'No linearity correction available.'}
        for bitname, bitdescription in zip(dqssb_lin, dqssb_desc):
            bitvalue = dqssb_lin[bitname]
            bitnumber = int(np.log(bitvalue)/np.log(2))
            newrow = (bitnumber, bitvalue, bitname, bitdescription)
            dqdef.append(newrow)
            dqmap = passedmap

        # Insert dq array into linearity reffile.
        finallinearity.dq = dqmap
        finallinearity.dq_def = dqdef

        # Create the proper headers.
        finallinearity = make_proper_header(finallinearity, files)

        # Save the reference file.
        outfile = output
        finallinearity.save(outfile)

        return outfile


    def run(self):
        '''Main function.'''

        # Read in parameter file and get parameters for steps.
        params = self.load_params(self.paramfile)
        coords = params['Inputs']['array_bounds']
        zeroframe = params['Inputs']['first_frame']
        sigma_range = params['Saturation']['sigma_range']
        limit = params['Saturation']['limit']
        poly_order = params['Rampfit']['poly_order']
        corr_order = params['CorrectSignal']['poly_order']
        accuracy = params['Softsat']['accuracy']
        ramp_percent = params['Softsat']['ramp_percent']

        if params['Uncertainty']['infile'] == 'None':
            unc = params['Uncertainty']['estimate']
        elif params['Uncertainty']['infile'] is not 'None':
            unc = fits.getdata(params['Uncertainty']['infile'], 2)
            unc = unc[0,:,:,:]

        # Read in list of files.
        files = self.read_listfile(params['Inputs']['infile'])

        # Check that input files have the same readpattern, array size, etc.
        if params['Inputs']['input_check'] is True:
            self.input_consistency_check(files)

        # Set up arrays to hold output data.
        num_coefs = params['CorrectSignal']['pipeline_coefs'] +  \
                    params['CorrectSignal']['poly_order']
        coeffs_arr, dq_arr, sat_arr, mask_arr = \
            self.init_arrays(files, num_coefs, coords)

        # Loop over files.
        for file, n in zip(files, np.arange(0, len(files))):

            # Load the data from each file.
            data, tgroup, ngroups, nskip, nframe = self.load_data(file)

            # Get frame values.
            frame = np.arange(nframe, ngroups, np.int(nframe+nskip))

            # Loop over pixel combinations for given array.
            for j,i in itertools.product(np.arange(coords[0], coords[2]),
                                         np.arange(coords[1], coords[3])):

                # Get signal values for each pixel.
                pixel = [j,i]
                signal = data[0,:,pixel[1],pixel[0]]

                # Catch warnings to flag them in the dq_array.
                with warnings.catch_warnings():
                    warnings.simplefilter("error", OptimizeWarning)
                    warnings.simplefilter("error", RuntimeWarning)
                    try:

                        # Get hard saturation and ramp with no saturated frames.
                        hard_sat, first_sat_grp = \
                            self.get_saturation(signal, sigma_range, limit)
                        nosat_ramp = [frame[:first_sat_grp],
                                      signal[:first_sat_grp]]

                        # Fit the ramp below hard saturation.
                        new_ramp, new_ideal_rate, intercept, nonlin, \
                                nonlin_unc = \
                            self.ramp_fit(pixel, zeroframe, tgroup, nosat_ramp,
                                         poly_order, unc)

                        # Do the non-linearity correction and get coefficients.
                        beta, correction_resid_std, soft, upperlimit = \
                            self.correct_signal(new_ramp, nonlin,
                                               nonlin_unc, intercept,
                                               new_ideal_rate, first_sat_grp,
                                               accuracy, signal, tgroup,
                                               corr_order, ramp_percent)

                        # Save all correction coefficients.
                        coeffs_arr[n,0,pixel[1],pixel[0]] = \
                                        params['CorrectSignal']['first_coeff']
                        coeffs_arr[n,1,pixel[1],pixel[0]] = \
                                        params['CorrectSignal']['second_coeff']
                        coeffs_arr[n,2:,pixel[1],pixel[0]] = beta[:]

                        # Find highly non-linear ramps, flag them in DQ array.
                        if correction_resid_std > \
                                        params['Softsat']['corr_limit']:
                            dq_arr[pixel[1],pixel[0]] = np.uint8(2)
                            mask_arr[n,:,pixel[1],pixel[0]] = True

                        # Load results for saturation into separate array.
                        sat_arr[n,0,pixel[1],pixel[0]] = hard_sat
                        sat_arr[n,1,pixel[1],pixel[0]] = soft
                        sat_arr[n,2,pixel[1],pixel[0]] = upperlimit
                        sat_arr[n,3,pixel[1],pixel[0]] = correction_resid_std

                        # Create a mask to use for averaging later that blocks
                        # out pixels for which soft saturation could not be
                        # determined (weird ramps, possibly). Update DQ array.
                        if np.isnan(soft) == True:
                            mask_arr[n,:,pixel[1],pixel[0]] = True
                            dq_arr[pixel[1],pixel[0]] = np.uint8(4)

                    # Flag pixels fitting or correction problems.
                    except (RuntimeError, ValueError, OptimizeWarning,
                            TypeError, RuntimeWarning):
                        dq_arr[pixel[1],pixel[0]] = np.uint8(4)
                        mask_arr[n,:,pixel[1],pixel[0]] = True
                        coeffs_arr[n,:,pixel[1],pixel[0]] = np.nan

        # Get statistics for correction coefficients, averaging over exposures.
        # Average is taken for masked array, with bad pixels masked out.
        # Trust numpy to do the averaging if bad pixels are filtered?
        coeffs_arr[mask_arr] = np.nan
        averages, errs = self.calc_stats(coeffs_arr)

        # Save saturation results and mask array to a fits file?
        if params['Output']['save_intermediates'] == True:
            newhdu = fits.PrimaryHDU(sat_arr)
            newhdulist = fits.HDUList([newhdu])
            newhdulist.writeto(params['Output']['saturation'],clobber=True)

            newhdu2 = fits.PrimaryHDU(mask_arr.astype('int'))
            newhdulist2 = fits.HDUList([newhdu2])
            newhdulist2.writeto('mask_for_filtering.fits',clobber=True)


        # Save averaged coefficients (and errors) to a FITS file.
        if params['Output']['coeffs'] is not 'None':
            outfilename = params['Output']['coeffs']
        else:
            outfilename = hdr['DETECTOR'] + '_'  \
                          + str(params['Output'][hdr['DETECTOR']])  \
                          + '_LinearityCoeff_ADU0_'  \
                          + str(datetime.date.today())  \
                          + '_ssblinearity_DMSorient.fits'
        outfile = self.save_reffile(averages,errs,dq_arr,files,outfilename)

        z = fits.open(outfile)
        z0=z[0]
        z1=z[1]
        z2=z[2]
        z3=z[3]

        #Add other things that the pipeline doesn't use, but are helpful.
        z0.header['S_DQINIT']= ('COMPLETE','Data Quality Initialization')
        #z0.header['S_SATURA']= ('COMPLETE','Saturation Checking')
        #z0.header['S_IPC']   = ('COMPLETE','Interpixel Capacitance Correction')
        z0.header['S_SUPERB']= ('COMPLETE','Superbias Subtraction')
        z0.header['S_REFPIX']= ('COMPLETE','Reference Pixel Correction')
        #z0.header['S_LINEAR']= ('COMPLETE','Linearity Correction')
        newhdu = fits.HDUList([z0,z1,z2,z3])
        newhdu.writeto(outfile,clobber=True)



if __name__ == '__main__':
    usagestring = 'USAGE: make_NIRCam_linearity_reffile.py paramfile.yaml'

    starting = datetime.datetime.now()

    nonlinearity = MakeNonlinRef()
    parser = nonlinearity.add_options(usage=usagestring)
    args = parser.parse_args(namespace=nonlinearity)

    nonlinearity.run()

    ending = datetime.datetime.now()
    difftime = ending - starting
    print("DONE. Elapsed time: {}".format(difftime))
