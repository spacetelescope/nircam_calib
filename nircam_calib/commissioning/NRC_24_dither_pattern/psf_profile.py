#!/usr/bin/env python
#
import inspect
import re
#
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse
import numpy as np
#
from astropy.io import fits
from astropy.table import Table,Column, MaskedColumn
from astropy.stats import sigma_clipped_stats
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.nddata import NDData
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.visualization import simple_norm

from astropy.stats import sigma_clipped_stats
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
#
import sep
#

import sys
sys.path.append("/home/cnaw/python/commissioning/")
from sub_profile import *
#------------------------------------------------------------------------

print("\nPSF_PROFILE.PY\n")

debug = 0
plot_results = True
norm_profile = True
det_samp     = 1

fwhm = 1.3
nsigma = 200.0
filter = 'F115W'
pixel_scale  = 0.0311
dr = 1.0

radii = [1., 2., 3.,4., 5., 6.,7., 8., 9., 10., 11., 12., 13., 14., 15.]
rsky  = [20., 30.]

template_r = []
for ii in range(0, 40):
    template_r.append(float(ii)+1.0)
# print(lineno()," template_r ", template_r)
# Centroid search parameters

psf= dict([
    ('F070W', 'PSF_NIRCam_F150W_fov_671_os_2.fits'),
    ('F115W', 'PSF_NIRCam_F115W_fov_508_os_2.fits'),
    ('F150W', 'PSF_NIRCam_F150W_fov_671_os_2.fits'),
    ('F277W', 'PSF_NIRCam_F277W_fov_595_os_2.fits')
])

template_dir = './templates/'
tempfile = template_dir+psf[filter]
#
# 1. analyse webbpsf template with oversampling = 2
#
ext = 0
radiit, apmagt, areat, difft, encircledt, temp_samp = read_webbpsf_ext(tempfile, ext, template_r, norm_profile, debug)
a1 = np.array(radiit)* pixel_scale
a2 = np.array(areat)
a3 = np.array(difft)
a4 = np.array(apmagt)
a5 = np.array(encircledt)

output = './templates/test_webbpsf_f115w_os2_prof.fits'
col1 = fits.Column(name='r',format='E', unit='arcsec',array=a1)
col2 = fits.Column(name='area',format='E',unit='pixel**2', array=a2)
col3 = fits.Column(name='flux_r',format='E',unit='counts/s', array=a3)
col4 = fits.Column(name='flux_t',format='E',unit='counts/s', array=a4)
col5 = fits.Column(name='encircled',format='E', unit='counts/arcsec**2',array=a5)
cols = fits.ColDefs([col1, col2, col3, col4, col5])
hdu  = fits.BinTableHDU.from_columns(cols)
hdu.writeto(output,overwrite=True)

#
# 2. analyse webbpsf template with oversampling = 1
#
ext = 1
radiit1, apmagt1, areat1, difft1, encircledt1, temp_samp1 = read_webbpsf_ext(tempfile, ext, template_r, norm_profile, debug)

#
# 3. analyse ISIM_CV3 images with photutils
#
file= []
file.append("/home/cnaw/isim_cv3/34314/NRCN812A-F115-5361113733_1_483_SE_2015-12-27T11h50m32_I000.wfs.fits")

file.append("/home/cnaw/isim_cv3/34314/NRCN812A-F115-5361113733_1_483_SE_2015-12-27T11h50m32_I001.wfs.fits")

#print("read_image")
(image, header, osamp, filter_name, scale) = read_image(file[0], debug)
# create a mask identifying NaNs and infs
#print("make mask")
mask = np.zeros(image.shape, dtype=bool)
bad  = np.where(image != image+0.0)
mask[bad] = True
#print("entering find_centroid")
#
#xc0, yc0 = find_centroid(image, fwhm, nsigma, debug)
ext = 0
(xc0, yc0) = run_sep(file[0] ,ext, nsigma, debug,plot_results=False)
print("entering aper_phot")
(apmag0, area0) = aper_phot(image, mask, xc0, yc0, radii, rsky, debug)
diff0,encircled0 = differential(apmag0,area0, norm_profile)
radii0 =[]
for ii in range(0,len(radii)):
    radii0.append(radii[ii])
    
#if(debug == 1):
print("xc, yc, image ", xc0, yc0, file[0])

(image, header, osamp, filter_name, scale) = read_image(file[1], debug)
# create a mask identifying NaNs and infs
mask = np.zeros(image.shape, dtype=bool)
bad  = np.where(image != image+0.0)
mask[bad] = True

#print("entering find_centroid")
#xc1, yc1 = find_centroid(image, fwhm, nsigma, debug)
ext = 0
(xc1, yc1) = run_sep(file[1] ,ext, nsigma, debug,plot_results=False)
(apmag1, area1) = aper_phot(image, mask, xc1, yc1, radii, rsky, debug)
diff1,encircled1 = differential(apmag1, area1, norm_profile)
radii1 =[]
for ii in range(0,len(radii)):
    radii1.append(radii[ii])
if(debug == 1):
    print("xc, yc, image ", xc1[1], yc1[1], file[1])

if(plot_results == True):

    # 
    # 4. webbpsf profile calculated using ~/src/photometry/test_pixwt
    #
    file = '/home/cnaw/src/photometry/profile_f115w_oversampled_2.dat'
    (r_webbpsf, f_webbpsf, f_area) = read_ascii(file, debug)
    r_webbpsf = r_webbpsf/temp_samp
    webbpsf_diff, webbpsf_encircled = differential(f_webbpsf, f_area, norm_profile)    
    #
    # 5. ISIM_CV3 image profiles using ~/src/photometry/test_pixwt
    #
    file = '/home/cnaw/python/commissioning/isim_cv3/profile1.dat'
    (r_for1, f_for1, area_for1) = read_ascii(file, debug)
    fortran_diff1, fortran_encircled1 = differential(f_for1, area_for1, norm_profile)
##
    file = '/home/cnaw/python/commissioning/isim_cv3/profile2.dat'
    (r_for2, f_for2, area_for2) = read_ascii(file, debug)
    fortran_diff2, fortran_encircled2 = differential(f_for2, area_for2, norm_profile)
    #
    # 6. ISIM_CV3 profiles calculated using IRAF.PHOT
    #
    file = '/home/cnaw/isim_cv3/isim_cv3_f115_a.dat'
    (r_iraf1, f_iraf1, area_iraf1) = read_iraf_isim_cv3(file,debug)
    iraf_diff1, iraf_encircled1 = differential(f_iraf1, area_iraf1,norm_profile)
    #
    file = '/home/cnaw/isim_cv3/isim_cv3_f115_b.dat'
    (r_iraf2, f_iraf2, area_iraf2) = read_iraf_isim_cv3(file, debug)
    iraf_diff2, iraf_encircled2 = differential(f_iraf2, area_iraf2,norm_profile)
    #
    # convert radii from pixels in to arc seconds
    #
    radii0    = pixel_scale * np.array(radii0)
    radii1    = pixel_scale * np.array(radii1)
    r_iraf1   = pixel_scale * r_iraf1
    r_iraf2   = pixel_scale * r_iraf2
    r_for1    = pixel_scale * r_for1
    r_for2    = pixel_scale * r_for2
    r_webbpsf = pixel_scale * r_webbpsf
    radiit    = pixel_scale * np.array(radiit)
    radiit1   = pixel_scale * np.array(radiit1)

    png_plot = filter+'_comp_profile.png'
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1, 1, 1)
##    ax.set_xlim(left=-0.01, right=20.0)
    ax.set_xlim(left=-0.01, right=1.0)
    ax.set_ylim(bottom=1e-05, top=1.8)
    plt.yscale("log")
#    plt.xlabel("radius (pixels)")
    plt.xlabel("radius (arc sec)")
    plt.ylabel("Normalised surface brightess profile (flux/arcsec^2)")
    plt.title(filter)
    plt.plot(radii0, encircled0,'-',marker=(5,2),color='grey')
    plt.plot(radii1, encircled1,'-',marker=(5,2),color='blue')
    plt.plot(r_for1, fortran_encircled1,'o',color='red')
    plt.plot(r_for2, fortran_encircled2,'o',color='orange')
    plt.plot(r_iraf1, iraf_encircled1,'^',color='green')
    plt.plot(r_iraf2, iraf_encircled2,'^',color='cyan')
#    plt.plot(r_webbpsf, webbpsf_encircled,'d',color='lime')
    plt.scatter(r_webbpsf, webbpsf_encircled,s=50,facecolors='none',edgecolors='lime')
    plt.plot(radiit, encircledt,'+',color='brown')
    plt.plot(radiit1, encircledt1,'x',color='gold')
##    png_plot = re.sub('.fits','.png',output)
    plt.text(0.5,1,           'photutils ISIM_CV3 exposure 1',color='grey')
    plt.text(0.5,10.**(-0.20),'IRAF      ISIM_CV3 exposure 1',color='green')
    plt.text(0.5,10.**(-0.40),'Fortran    ISIM_CV3 exposure 1',color='red')
    plt.text(0.5,10.**(-0.70),'photutils ISIM_CV3 exposure 2',color='blue')
    plt.text(0.5,10.**(-0.90),'IRAF       ISIM_CV3 exposure 2',color='cyan')
    plt.text(0.5,10.**(-1.10),'Fortran   ISIM_CV3 exposure 2',color='orange')
    plt.text(0.5,10.**(-1.40),'Fortran webbpsf os2',color='lime')
    plt.text(0.5,10.**(-1.60),'photutils webbpsf os2',color='brown')
    plt.text(0.5,10.**(-1.80),'photutils webbpsf os1',color='gold')
#   
    plt.savefig(png_plot,bbox_inches='tight')
    plt.show()
    exit(0)
    png_plot = filter+'_int_profile.png'
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1, 1, 1)
#    ax.set_xlim(left=-0.01, right=20.0)
    ax.set_xlim(left=-0.01, right=0.65)
    ax.set_ylim(bottom=310, top=1050)
    plt.plot(radii0, apmag0,'o-',markersize=4,color='grey')
    plt.plot(radii, apmag1,'o-',markersize=4,color='blue')
#    plt.xlabel("radius (pixels)")
    plt.xlabel("radius (arc sec)")
    plt.ylabel("totalintensity")
    plt.title(filter)
    plt.plot(r_for1, f_for1,'o',color='red')
    plt.plot(r_for2, f_for2,'o',color='orange')
    plt.plot(r_iraf1, f_iraf1,'^',color='green')
    plt.plot(r_iraf2, f_iraf2,'^',color='cyan')
    plt.text(8,680,'Photutils  ISIM_CV3 exposure 1',color='grey')
    plt.text(8,650,'IRAF    ISIM_CV3 exposure 1',color='green')
    plt.text(8,620,'Fortran ISIM_CV3 exposure 1',color='red')

    plt.text(8,580,'Photutils  ISIM_CV3 exposure 2',color='blue')
    plt.text(8,550,'IRAF    ISIM_CV3 exposure 2',color='cyan')
    plt.text(8,520,'Fortran ISIM_CV3 exposure 2',color='orange')
    plt.savefig(png_plot,bbox_inches='tight')
    plt.show()
