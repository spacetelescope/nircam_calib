#! /usr/bin/env python

"""Compute basic statistics on input slope image

Detector-average and median slope in each integration for science pixels
Average slope in each integration for reference pixels
Create CDS images for all groups
Calculate median value for each CDS image.
Plot CDS medians vs time and FPA temp
Count rate histogram for the detector slope image
"""

from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from jwst import datamodels




def detector_avg_and_median_per_integration(filename):
    """Calculate the average and median values for each integration and
    for the science pixels and reference pixels

    Parameters
    ----------
    data : numpy.ndarray
        3D array of data (e.g. from rateints.fits file)

    Returns
    -------
    lists
        Lists of science pixel means, medians, and sigma-clipped means and medians, as
        well as reference pixel means, medians, and sigma-clipped means and medians for
        each integration
    """
    model = datamodels.open(filename)

    nints = model.data.shape[0]
    sci_medians = []
    ref_medians = []
    sci_means = []
    ref_means = []
    sci_clip_mean = []
    sci_clip_median = []
    sci_clip_dev = []
    ref_clip_mean = []
    ref_clip_median = []
    ref_clip_dev = []

    #refpix = (model.dq & datamodels.dqflags.pixel['REFERENCE_PIXEL'] > 0)
    #scipix = (model.dq & datamodels.dqflags.pixel['REFERENCE_PIXEL'] == 0)

    for integ in range(nints):
        dq = model.dq#[integ, :, :]
        data = model.data#[integ, :, :]
        refpix = (dq[integ,: ,:] & datamodels.dqflags.pixel['REFERENCE_PIXEL'] > 0)
        scipix = (dq[integ,: ,:] & datamodels.dqflags.pixel['REFERENCE_PIXEL'] == 0)

        frame = data[integ, :, :]

        sci_medians.append(np.nanmedian(frame[scipix]))
        ref_medians.append(np.nanmedian(frame[refpix]))
        sci_means.append(np.nanmean(frame[scipix]))
        ref_means.append(np.nanmean(frame[refpix]))

        mn, med, stdev = sigma_clipped_stats(frame[scipix], sigma=3)
        sci_clip_mean.append(mn)
        sci_clip_median.append(med)
        sci_clip_dev.append(stdev)

        mn, med, stdev = sigma_clipped_stats(frame[refpix], sigma=3)
        ref_clip_mean.append(mn)
        ref_clip_median.append(med)
        ref_clip_dev.append(stdev)

    return sci_medians, sci_means, sci_clip_mean, sci_clip_median, sci_clip_dev, ref_medians, ref_means, ref_clip_mean, ref_clip_median, ref_clip_dev



def cds_all_groups(filename, save=True, outfile=None):
    """Create CDS images for all groups and all integrations in the input file.
    Save to a fits file.
    """
    # If no output filename is given, save the results to the same directory
    # as the input data.
    if outfile is None:
        in_dir, in_base = os.path.split(filename)
        outfile = os.path.join(in_dir, f'CDS_images_for_{in_base}')

    with fits.open(filename) as fileobj:
        data = fileobj[1].data
        int_times = fileobj[4].data
        grptime = h[0].header['TGROUP']

    if len(data.shape) != 4:
        raise ValueError('Expecting 4D input file.')

    nints, ngroups, ny, nx = data.shape
    cds = np.zeros((nints, ngroups-1, ny, nx))
    for integ in range(data.shape[0]):
        # Convert to floats before subtracting
        cds[integ, :, :, :] = data[integ, 1:, :, :]*1. - data[integ, 0:-1, :, :]
        starttime = int_times['int_start_MJD_UTC'][integ]
        #cds_times = [starttime + grptime/3600./24.*(i+1) for i in range(ngroups-1)]  absolute time
        cds_times = [grptime * (i+1) for i in range(ngroups-1)]  # relative time

    if save:
        h0 = fits.PrimaryHDU(cds)
        hlist = fits.HDUList([h0])
        hlist.writeto(outfile, overwrite=True)

    return cds, cds_times


def cds_medians(data, dq):
    """Calculate the median value of the science pixels in a series of CDS images

    Paramters
    ---------
    data : numpy.ndarray
        4D array of data - output from cds_all_groups

    dq : numpy.ndarray
        4D array of DQ values
    """
    if len(data.shape) != 4:
        raise ValueError("Expecting 4D input array.")

    scipix = (dq & datamodels.dqflags.pixel['REFERENCE_PIXEL'] == 0)

    nints, ngrps, ny, nx = data.shape
    meds = []
    for integ in range(nints):
        for grp in range(ngrps):
            frame = data[integ, grp, :, :]
            meds.append(np.median(frame[scipix[integ, grp, :, :]]))
    return meds

def plot_cds_means_vs_time_temp(cds_means, times, temperatures, outfile):
    f, a = plt.subplots(ncols=2)
    a[0].plot(times, cds_means)
    #a.plot(times-times[0], cds_means)
    a[1].plot(temperatures, cds_means)
    plt.show()
    f.savefig(outfile)


def create_histogram(filename):
    """Create a histogram of the data values in each group/integration of the input data

    Parameters
    ----------
    filename : str
        Name of fits file with data to be examined. For NRC-03, rate.fits file? rateints.fits file?

    """
    data = fits.open(filename)
    header = data[0].header
    filtername = header['FILTER']
    detector = header['DETECTOR']
    subarr = header['SUBARRAY']

    in_dir, in_base = os.path.split(filename)
    outfile = os.path.join(in_dir, f'histogram_{in_base}.png')

    obs = int(in_base[7:10])

    if len(data[1].data.shape) == 2:
        # rate file
        hist, bin_edges = np.histogram(data)

        # convert bin edges to bin centers
        bins = bin_edges[0:-1] + (bin_edges[1:] - bin_edges[0:-1]) / 2.

        # Plot
        f, a = plt.subplots()
        a.bar(bin_edges[:-1], hist, width=1)

        f.savefig(outfile)
        #plt.savefit(outfile)

    elif len(data[1].data.shape) == 3:
        # rateints file
        for integ in range(data[1].shape[0]):
            integration = data[1].data[integ, :, :]
            finite = np.isfinite(integration)
            hist, bin_edges = np.histogram(integration[finite], bins=np.linspace(-10,30,400))


            #print(bin_edges)
            #print(hist)
            # convert bin edges to bin centers
            #bins = bin_edges[0:-1] + (bin_edges[1:] - bin_edges[0:-1]) / 2.

            # Plot
            """
            f, a = plt.subplots(figsize=(8, 6))
            a.bar(bin_edges[:-1], hist, color='blue', width=0.1)
            a.set_xlim(-0.1, 4)
            a.set_ylabel('Num. Pixels')
            a.set_xlabel('Signal Rate (DN/sec)')
            a.set_title(f'Obs 4 Signal Rates, {detector}, {filtername}, Integration {integ+1}')

            outfile = os.path.join(in_dir, f'histogram_obs4_{detector}_{filtername}_integration{integ+1}.jpg')
            #plt.show()
            f.savefig(outfile)
            plt.close()
            """

            n, bins, patches = plt.hist(integration[finite], bins=np.linspace(-0.1,4,41), alpha=0.25)
            plt.ylabel('Num. Pixels')
            plt.xlabel('Signal Rate (DN/sec)')
            plt.title(f'Obs {obs} Signal Rates, {detector}, {subarr}, {filtername} per Integration')
            plt.subplots_adjust(left=0.15)
        outfile = os.path.join(in_dir, f'histogram_obs{obs}_{subarr}_{detector}_{filtername}.png')
        plt.savefig(outfile)
        plt.close()
        #save_multi_image(outfile)

def save_multi_image(filename):
    """Save multiple matplotlib plots to a single PDF
    """
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()

