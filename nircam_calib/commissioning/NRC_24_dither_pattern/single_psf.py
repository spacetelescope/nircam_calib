#!/usr/bin/env python
#
import argparse
from astropy.io import fits
from astropy.table import Table,Column, MaskedColumn
from astropy.stats import sigma_clipped_stats
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.nddata import NDData
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.visualization import simple_norm

from astropy.stats import sigma_clipped_stats
#
import inspect
#
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse
#
import numpy as np
#
import os
from os.path import exists
#
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
#
import photutils
from photutils.detection import DAOStarFinder
from photutils.detection import IRAFStarFinder
from photutils.detection import find_peaks
#
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.psf import BasicPSFPhotometry
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
# for regular expressions
import re
import sep
import sys
sys.path.append("/home/cnaw/python/commissioning/")
from sub_profile import *
#=======================================================================
#-----------------------------------------------------------------------
def single_psf(file, plot_results=False, use_sep=False, debug=False):
    norm_profile = True
        
    # Centroid search parameters
    print("\nsingle_psf.py: file is ",file)
    if(debug == 1 or debug == True):
        print("single_psf.py: plot_results is ", plot_results)
        print("single_psf.py: use SEP      is ", use_sep)

    fwhm = 1.3
    if(re.search('PSF_NIRCam',file) != None):
        nsigma = 200.0
    else:
        nsigma = 1.5
        
    # define output file names for plots
    output1 = re.sub('fits','_diff_profile.png',file)
    output2 = re.sub('fits','_int_profile.png',file)
    output1 = re.sub('/home/cnaw/guitarra/data/WebbPSF_NIRCam_PSFs/','',output1)
    output2 = re.sub('/home/cnaw/guitarra/data/WebbPSF_NIRCam_PSFs/','',output2)

    hdunum = 0
    image, header, oversampling_rate,filter, pixel_scale = read_image(file, debug)
    if('FILTER' in header):
        filter = header['FILTER']
        pixel_scale = scale_from_filter(filter)
        # Try to figure out the pixel scale from file name
    else:
        pixel_scale, filter = scale_from_filename(file)

    if('PIXELSCL' in header):
        pixel_scale = float(header['PIXELSCL'])

    fwhm        = fwhm * oversampling_rate

    print("single_psf.py: pixel_scale, filter, oversampling_rate : ", pixel_scale, filter, oversampling_rate)
    output = re.sub('.fits','_prof.fits',file)
        
    print("single_psf.py: output file is ", output)
    print("single_psf.py: png    file is ", output1)

    nan_mask = np.zeros(image.shape, dtype=bool)
    nan_mask = np.where(np.isnan(image),True,False)
    naxis1 = int(header['NAXIS1'])
    naxis2 = int(header['NAXIS2'])

    xc = naxis1/2.0 -0.5
    yc = naxis2/2.0 -0.5
    if(use_sep == 0 or use_sep == False):
        (xc, yc) = find_centroid(image, fwhm, nsigma, debug)
    else :
        ext = 0
        (xc, yc) = run_sep(file,ext, nsigma, debug,plot_results=plot_results)

    rmax = 80.0
    nsky = 10
    dr = 1.0

    xmax = min(xc, naxis1-xc)
    ymax = min(yc, naxis2-yc)
#    rmax = min(xmax, ymax)
    rmax = int(rmax)
    rsky  = [rmax-nsky, rmax]
    objmax = rsky[0] - 1
    #    if(debug == 1 or debug == True):

    nr = int((objmax/dr))-1
    fcenter = image[int(yc),int(xc)]
#    print("naxis1, naxis2 ",naxis1, naxis2)
#    print("at line ", lineno(), 'xc, yc, rmax, nr ',xc, yc, rmax, nr,fcenter)
    radii = []
    for index in range(0,nr):
        radii.append((index+1) * dr)

    nradii = len(radii)

    (apflux, area) = aper_phot(image, nan_mask, xc, yc, radii, rsky, debug)
    diff,encircled0 = differential(apflux, area, norm_profile, True)
    diff,encircled1 = differential(apflux, area, norm_profile, False)

    # print(lineno()," radii ", radii)

    rt        = []
    for index in range(0,len(radii)):
        rt.append(radii[index] * pixel_scale)
        #        if(debug == 1) :
#        print("line : %d %10.5f %10.5f %10.5f %11.5e %10.5f" % (lineno(), rt[index], area[index],diff[index], apflux[index], encircled[index]))

    if(debug == 1 or debug == True):
        print("single_psf.py: naxis1, naxis2, xc, yc ", naxis1, naxis2, xc, yc)
        print("single_psf.py: xmax, ymax, rmax ",xmax, ymax, rmax)
        print("single_psf.py: rsky ", rsky)
        print("single_psf.py: nr ", nr, len(rt))
        print("single_psf.py: len(rt), len(encircled0)",len(rt), len(encircled0))
        # print("single_psf.py: rat", len(rt),rt)

    a1 = np.array(rt)
    a2 = np.array(area)
    a3 = np.array(diff)
    a4 = np.array(apflux)
    a5 = np.array(encircled0)
    a6 = np.array(radii)
    a7 = np.array(encircled1)
    
    col1 = fits.Column(name='r',format='E', unit='arcsec',array=a1)
    col2 = fits.Column(name='area',format='E',unit='pixel**2', array=a2)
    col3 = fits.Column(name='flux_r',format='E',unit='counts/s', array=a3)
    col4 = fits.Column(name='flux_t',format='E',unit='counts/s', array=a4)
    col5 = fits.Column(name='encircled0',format='E', unit='flux/arcsec^2',array=a5)
    col6 = fits.Column(name='rpix',format='E', unit='pixel',array=a6)
    col7 = fits.Column(name='encircled1',format='E', unit='flux/arcsec^2',array=a7)
    cols = fits.ColDefs([col1, col2, col3, col4, col5,col6, col7])
    hdu  = fits.BinTableHDU.from_columns(cols)
    hdu.writeto(output,overwrite=True)

    
    png_plot = filter+'_diff_profile.png'
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(left=-0.01, right=1.0)
    ax.set_ylim(bottom=1.e-05, top=1.5)
    plt.yscale("log")
    plt.xlabel("radius (arc sec)")
    plt.ylabel("Normalized encircled flux (==dflux(i)/flux(i)")
    plt.title(filter+" "+file)
    plt.plot(rt, encircled0,'o', color='brown')
    plt.plot(rt, encircled1,'x', color='blue')
    png_plot = re.sub('.fits','.png',output1)
    plt.savefig(png_plot,bbox_inches='tight')

    if(plot_results == True):
        plt.show()
    return

#    png_plot = filter+'_int_profile.png'
#    fig = plt.figure(figsize=(8,6))
#    plt.xlabel("radius (pixels)")
#    plt.ylabel("totalintensity")
#    plt.title(filter)
#    plt.plot(rt, apflux,'+', color='brown')
#    png_plot = re.sub('.fits','.png',output2)
#    plt.savefig(png_plot,bbox_inches='tight')
#    plt.show()
#
#=======================================================================
#
if __name__ == "__main__":

    # Use argument parser to pass keywords to program
    parser = argparse.ArgumentParser(description="Calculate 1D PSF profile")
#    parser.add_argument('--file', help='file to process', type=str,default="/home/cnaw/commissioning/car_24_apt_01073/analysis/jw01073006001_01101_00006_nrcb3_cal_F070W_psf.fits")
    parser.add_argument('--file', help='file to process', type=str)
    
    parser.add_argument("--plot", help="show plot results",
                        action="store_true")
    parser.add_argument("--use_sep", help="use SEP for object detection",
                        action="store_true")
    
    parser.add_argument("--debug", help="debug mode", action="store_true")
    args = parser.parse_args()

    single_psf(args.file, plot_results=args.plot, use_sep=args.use_sep, debug=args.debug)
