#!/usr/bin/env python
# 
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import NDData
#
from astropy.stats import sigma_clipped_stats
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import biweight_location
from astropy.stats import mad_std
#
from astropy.table import Table,Column, MaskedColumn
from astropy.visualization import simple_norm
#
from glob import glob
#
import matplotlib.pyplot as plt
import numpy
import numpy as np
#
from os.path import exists
#
import pathlib
#
import photutils
import photutils.background

from photutils.centroids import centroid_com, centroid_quadratic
from photutils.centroids import centroid_1dg, centroid_2dg

from photutils.detection import DAOStarFinder
from photutils.detection import IRAFStarFinder
from photutils.detection import find_peaks
#
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
#
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.psf import BasicPSFPhotometry
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
#
import re
import scipy
import sys
sys.path.append("/home/cnaw/python/commissioning/")
from sub_profile import *

import argparse
#--------------------------------------------------------------------------
#
# 

def aper_phot_old(image, xc, yc, radii, plot_results,profile_png, mask=None,verbose=False):

    positions = [(xc, yc)]

    apertures = [CircularAperture(positions, r=r) for r in radii ]
    phot_table = aperture_photometry(image, apertures,mask=mask)
#
    junk = []
#    print(phot_table)
#    print(type(phot_table))
    for index  in phot_table.colnames:
        if('aperture_sum' in index):
            array = phot_table[index].data[0]
            temp  = array.tolist()
            junk.append(temp)
    values = np.array(junk)
    diff   = np.zeros(len(values))
    diff[0] = values[0]
    for ii in range(1,len(values)):
        diff[ii]= values[ii]-values[ii-1]
        if(verbose == True or verbose == 1):
            print(values[ii], diff[ii])

    if(verbose == True or verbose == 1):
        print("positions: " ,positions)
        print("Radii  :", radii)
        print("values :",values)

    if(plot_results == True) :
        fig = plt.figure(figsize=(10,8))
        plt.subplot(2, 1, 1)
        plt.plot(radii, values,'bo',markersize=4)
        plt.xlabel("radius (pixels)")
        plt.ylabel("totalintensity")

        plt.subplot(2, 1, 2)
        diff = diff/diff[0]
        plt.plot(radii, diff,'bo',markersize=4)
        plt.yscale("log")
        plt.xlabel("radius (pixels)")
        plt.ylabel("Normalised intensity")
        plt.savefig(profile_png,bbox_inches='tight')
        plt.show()
        return
#
#--------------------------------------------------------------------------
# This follows ipsis litteris (well, almost)
# https://photutils.readthedocs.io/en/stable/epsf.html?highlight=EPSFBuilder
#
def calculate_psf_template(image, xx, yy, size, oversampling_rate,
                           maxiters=11, mask=None, plot_candidates=False) :
    hsize = (size+2)/2
# remove stars close to the sample edges
    pos_mask  = ((xx > hsize) & (xx < (image.shape[1] - 1 - hsize)) &
                 (yy > hsize) & (yy < (image.shape[0] - 1 - hsize)))

# Table of good star positions
    stars_tbl      = Table()
    stars_tbl['x'] = xx[pos_mask]
    stars_tbl['y'] = yy[pos_mask]


# Subtract background
#    print("image.shape ", image.shape)
    mean_val, median_val, std_val = sigma_clipped_stats(image, mask=mask, sigma=2.)
    image -= median_val
    nddata = NDData(data=image)
    
    stars = extract_stars(nddata, stars_tbl, size=size)
    print("number of stars ", len(stars))
#
# Plot candidates (use with care!)
#
    if(plot_candidates == True) :
        nrows = 1+int(np.sqrt(float(len(stars))))
        ncols = 1+ int(float(len(stars))/float(nrows))
#        print("nrows, ncols, len(stars) ",nrows, ncols, len(stars))
#        nrows = 7
#        ncols = 8
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20),squeeze=True)
        ax = ax.ravel()
        for i in range(nrows*ncols):
            norm = simple_norm(stars[i], 'log', percent=99.)
            ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')
        plt.show()
#
# build average PSF
#
    epsf_builder = EPSFBuilder(oversampling=oversampling_rate, maxiters=3,progress_bar=False)
    epsf, fitted_stars = epsf_builder(stars)
    return epsf, fitted_stars

#--------------------------------------------------------------------------

def plot_psf_fit(image, residual_image):
    plt.subplot(1, 2, 1)
    norm = simple_norm(image, 'log', percent=99.)
    plt.imshow(image, norm=norm, origin='lower', cmap='viridis')
#    plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest',origin='lower')
    plt.title('Data')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    plt.subplot(1, 2, 2)
    norm = simple_norm(residual_image, 'log', percent=99.)
    plt.imshow(residual_image, norm=norm, origin='lower', cmap='viridis')
#    plt.imshow(residual_image, cmap='viridis', aspect=1,interpolation='nearest', origin='lower')
    plt.title('Residual')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    plt.show()
    return

#----------------------------------------------------------------------
#======================================================================


def main(file, target_dir, plot_results=False, overwrite=False):
    fwhm    = 1.70
    nsigma  = 1.0
    size  = 64
    verbose = False
    psf_resamp = 2
    debug = 0

    if(re.search('psf.fits',file) != None):
       print("**** skipping PSF file:", file)
       return
       
    print("file is ", file," target_dir is ", target_dir, " plot_results is ", plot_results," overwrite is ", overwrite)

    (image, header, samp_rate, filter, pixel_scale) = read_image(file, debug)
    print("at line ", lineno(), " Filter is ", filter," pixel_scale ", pixel_scale)
    if(filter != "null" and pixel_scale == 0.0):
        pixel_scale = scale_from_filter(filter)

    if(filter == "null" and pixel_scale == 0.0):
        pixel_scale, filter = scale_from_filename(file)
        print("at line ", lineno(), " Filter is ", filter," pixel_scale ", pixel_scale)

    if(filter == "null"):
        print("at line ", lineno(), " no filter keyword for file ", file)
        return


    junk = re.split('/', file)
    if(target_dir != None):
        new_file = target_dir + junk[len(junk)-1]
    else:
        new_file = file
    print("Filter is ", filter)
    print("oversampling rate  is ", samp_rate)
    print("pixel_scale is " , pixel_scale)
#
    psf_output = re.sub('.fits','_'+filter+'_psf.fits',new_file)
    output_exists = exists(psf_output)
    
    if(output_exists and overwrite == False) :
        print ("file has already been reduced ", file, psf_output)
        return
    else:
        print("reducing      ", file)
        print("psf_output is ", psf_output)

    #
    #-----------------------------------------------------------------------#
    # create a pixel  mask for NaNs and noughts
    mask_file = re.sub('.fits','_mask.fits',new_file)
    zero_mask = numpy.where(image == 0,0,1)
    nan_mask  = numpy.where(np.isnan(image),0,1)
    zero_mask = nan_mask * zero_mask
    fits.writeto(mask_file,zero_mask,header=header,overwrite=True)

    # photutils uses boolean masks
    nan_mask = numpy.where(zero_mask == 0,True,False)

    # detect sources using the photutils algorithms
    sources = find_objects(image, fwhm, nsigma, debug, mask=nan_mask)

    if(verbose == True) :
        for col in sources.colnames:
            sources[col].info.format = "%.8g"
        print(sources)

    #
    # Write regions file for detected stars
    #
    np.xc = sources['xcentroid']
    np.yc = sources['ycentroid']
    #
    output = open('test1.reg','w')
    for ii in range(0,len(np.xc)):
        line = 'point '+str(np.xc[ii]+1)+' '+str(np.yc[ii]+1) +' # point=diamond color=blue\n'
        output.write(line)
    #
    output.close()
    #
    # Fit average PSF, following:
    # https://photutils.readthedocs.io/en/stable/epsf.html
    #
    # Make cuttouts
    #
    xx    = sources['xcentroid']
    yy    = sources['ycentroid']

    epsf, fitted_stars = calculate_psf_template(image, xx, yy, size, psf_resamp,
                                                maxiters=3, mask=nan_mask,plot_candidates=False) 
    comment = 'Filter'
    header.set('FILTER',filter, comment)
    comment = 'oversampling rate == 1 unless WEBBPSF template'
    header.set('OVERSAMP',samp_rate, comment)
    comment = 'arcsec/pixel, corrected for DET_SAMP '
    header.set('PIXELSCL',pixel_scale/psf_resamp, comment)
    comment = 'oversampling rate used in EPSFBuilder'
    header.set('DET_SAMP',psf_resamp, comment)
    fits.writeto(psf_output,epsf.data,header=header,overwrite=True)
    #
    # PSF profile using the average PSF (may not be the best idea)
    #
    (xc, yc) =  centroid_com(epsf.data)
    # print ("centroid_com ", (xc, yc))
    print("ouput is ", psf_output)
    print("PSF resampling rate ", psf_resamp, "shape", epsf.data.shape, "centroid", xc, yc)
    radii = np.arange(1, 15, dtype=float)
    radii = radii * samp_rate
    print("array ", radii)

    radii = radii.tolist() 


    if(plot_results == True or plot_results == 1) :
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        png_plot = re.sub('.fits','.png',psf_output)
        plt.savefig(png_plot,bbox_inches='tight')
        plt.show()

        profile_png = re.sub('psf.png','profile.png',png_plot)
        aper_phot_old(epsf.data, xc, yc, radii, plot_results, profile_png)
        print("average PSF saved as ", output)
        print("PNG plot    saved as ", png_plot)
#
    print("file analysed ", file,"\n\n")
    return
#-------------------------------------------------------------------------------
#    # PSF photometry
#def psf_photometry(image, mask, sources, epsf):
#
#    # It is unclear how to add masks to PSF photometry;
#    print ("Stellar  photometry:")
#    sigma_psf = 2.0
#    daogroup  = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
#    fitter    = LevMarLSQFitter()
#    psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
#    mmm_bkg   = MMMBackground()
#    psf_model.x_0.fixed = True
#    psf_model.y_0.fixed = True
##
#    pos = Table(names=['x_0', 'y_0'],
#                data=[sources['xcentroid'],sources['ycentroid']])
#    print("entering BasicPSFPhotometry")
#    z_photometry = BasicPSFPhotometry(group_maker=daogroup,bkg_estimator=mmm_bkg,
#                                    psf_model=psf_model,fitter=LevMarLSQFitter(),
#                                    fitshape=(11, 11))
##
#    result_tab = z_photometry(image=image, init_guesses=pos)
#    print(result_tab)
##
#    xc   = result_tab['x_0']
#    yc   = result_tab['y_0']
#    flux = result_tab['flux_0']
#    fit  = numpy.where(np.isnan(flux),0,1)
##
#    output = open('test3.reg','w')
#    for ii in range(0,len(xc)):
#        if(fit[ii] == 0):
#            line = 'point '+str(xc[ii]+1)+' '+str(yc[ii]+1) +' # point=box color=red\n'
#        else:
#            line = 'point '+str(xc[ii]+1)+' '+str(yc[ii]+1) +' # point=circle color=green\n'
#        output.write(line)
##
#    output.close()

if __name__ == "__main__":

    # Use argument parser to pass keywords to program
    parser = argparse.ArgumentParser(description="Run multiple instances of average_psf.py")
    parser.add_argument('--file', help='Input file path.', type=str)
    parser.add_argument('--target_dir', help='Output directory', type=str,default="./")
    parser.add_argument("--plot_results", help="Plot results.", action="store_true")
    parser.add_argument("--overwrite", help="Set to overwrite results.", action="store_true")
    args = parser.parse_args()

#    root_dir = '/home/cnaw/commissioning/car_24_apt_01073/'
#    file    = 'jw01073001_lwb_F277W_i2d.fits'
#    file    = 'jw01073002_lwb_F277W_i2d.fits'
#    file    = 'jw01073003_lwb_F277W_i2d.fits'
#    file    = 'jw01073007_swa_F115W_i2d.fits'
#    file    = 'jw01073006001_01101_00001_nrcb3_cal.fits'

    file         = args.file
    target_dir   = args.target_dir
    plot_results = args.plot_results
    overwrite    = args.overwrite
    print("target_dir is ", target_dir)
    main(file, target_dir,plot_results=plot_results, overwrite=overwrite)
