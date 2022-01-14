#!/usr/bin/env python
#
# Calculate average PSFs from images obtained as part of CAR 24
# (dither pattern confirmation)
# multi-processing code by Jarron Leisenring (UofAz)
#
# cnaw 2022-01-13
#
# For passing argument to python program such as:
#   $ ./psf_template_mp.py --ncpu 8
# Try:
#   $ ./psf_template_mp.py --help
import argparse

import os, re, sys
from glob import glob
import multiprocessing as mp
import numpy as np
from os.path import exists
from pathlib import Path
import photutils
import scipy

# Progress bar
from tqdm.auto import tqdm
import traceback
#
#-----------------------------------------------------------------------
#
def rand_val():

    # Must call random seed to ensure unique values during multiprocessing
    np.random.seed()
    num = np.random.random()
    print(num)

#
#-----------------------------------------------------------------------
#
def reduce(args):
    """This might work in place of `reduce`, but I don't know how `stellar_photometry` is coded"""

    import stellar_photometry

    file_to_read, target_dir, plot_results, overwrite, test = args
    try:
        if test:
            # Pause for 1 second
            import time
            time.sleep(1)
        else:
            stellar_photometry.main(file_to_read, target_dir, plot_results=plot_results, overwrite=overwrite)

        return True
    except Exception as e:
        print('Caught exception in worker thread ({}):'.format(file_to_read))
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        print('')
        return False
#
#-----------------------------------------------------------------------
#
def run_multiprocess(func, worker_args, nsplit):
    """Run multiple instances of some function"""

    # Select the lesser of number of files that exist or requested CPUs
    nfiles = len(worker_args)
    nsplit = np.min([nsplit, nfiles])

    results = []
    if nsplit>1:
        try:
            with mp.Pool(nsplit) as pool:
                for res in tqdm(pool.imap_unordered(func, worker_args), total=nfiles):
                    results.append(res)
                pool.close()
            if results[0] is None:
                raise RuntimeError('Returned None values. Issue with multiprocess??')
        except Exception as e:
            print('Caught an exception during multiprocess.')
            print('Closing multiprocess pool with error.')
            pool.terminate()
            pool.close()
            raise e
        else:
            print('Closing multiprocess pool without error.')
    else:
        results = [func(args) for args in worker_args]

    return np.asarray(results)
#
#-----------------------------------------------------------------------
#
def main(ncpu, target_dir, plot_results=False, overwrite=False, test=False):
    """Main function call"""

    ncpu_total = mp.cpu_count()
    print('Total CPUs: ', ncpu_total)
    print('Requested CPUs: ', ncpu)

    # Ensure we aren't requesting more CPUs than actually exist!
    if ncpu > ncpu_total:
        ncpu = ncpu_total

    # Get files but first, what computer are we on ?
    host = os.environ.get('HOST')
    root_dir = './data/'
    if(host == 'ema.as.arizona.edu'):
        root_dir   = '/home/cnaw/commissioning/car_24_apt_01073/'
        target_dir = '/home/cnaw/commissioning/car_24_apt_01073/analysis/'
    if(host == 'orange.as.arizona.edu'):
        root_dir   = '/data1/car_24_apt_01073/mirage/reduced/'
        target_dir = '/data1/car_24_apt_01073/mirage/analysis/'
    status = exists(target_dir)
    
# create analysis directory if it does not exist
    if(status == False) :
        os.makedirs(target_dir)
        
    print("root_dir is   ", root_dir)
    print("target_dir is ", target_dir)
    if test:
        files_to_read = ['example_file_{}'.format(i) for i in range(150)]
    else:
        files_to_read = sorted(glob(root_dir+'*i2d.fits'))
        
    # files_to_read = ['one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight']
    nfiles = len(files_to_read)
    print("there are ",nfiles," files to process")
    # Set all the work arguments to pass to reduce()
    worker_args = [(files_to_read[ii], target_dir, plot_results, overwrite, test) for ii in range(nfiles)]

    
    # How many simultaneous instances to we want to split over?
    # Select the lesser of number of files that exist or requested CPUs
    nsplit = np.min([nfiles, ncpu])

    results = run_multiprocess(reduce, worker_args, nsplit)
    if np.alltrue(results):
        print('All files completed without issue.')
    else:
        print('The following files failed:')
        ind_bad = np.where(results==False)[0]
        for ii in ind_bad:
            print('  {}'.format(files_to_read[ii]))
#
#=======================================================================
#
if __name__ == "__main__":

    # Use argument parser to pass keywords to program
    parser = argparse.ArgumentParser(description="Run multiple instances of stellar_photometry.py")
    parser.add_argument('--ncpu', help='Number of requested instances. Default = 4.', 
                        type=int, default=4)
    parser.add_argument("--target_dir", help="directory where results will be stored", type=str, default="./analysis/")

    parser.add_argument("--plot", help="Plot results.", action="store_true")
    parser.add_argument("--overwrite", help="Set to overwrite results.", action="store_true")
    parser.add_argument("--test", help="Testing. Skips run os.system() command. Fake file names.", action="store_true")
    args = parser.parse_args()

    main(args.ncpu, args.target_dir, plot_results=args.plot, overwrite=args.overwrite, test=args.test)
