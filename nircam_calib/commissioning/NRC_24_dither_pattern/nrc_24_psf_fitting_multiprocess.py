import sys
import numpy as np
import pysiaf
from datetime import datetime
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import imshow_norm, PercentileInterval, SqrtStretch
from photutils import BasicPSFPhotometry
from photutils import DBSCANGroup, MMMBackground
from photutils.detection import DAOStarFinder
import webbpsf
from webbpsf.gridded_library import display_psf_grid

import jwst
from jwst.datamodels import ImageModel

import multiprocessing


# This is the function that actually runs the PSF fitting on the data files from the same SCA
def psf_fit(data_array, data_file, psf_grid, fheader, imagemodel):

	# Let's run a quick fitting using DAOStarFinder
	mean, median, std = sigma_clipped_stats(data_array, sigma=3.0) 
	daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)  
	sources = daofind(data_array - median) 
	
	# And now, let's make some cuts to remove objects that are too faint or too bright
	flux_min = 5#5 #0
	flux_max = 50#50 #1000
	flux_range = np.where((sources['flux'] > flux_min) & (sources['flux'] < flux_max))[0]
	
	init_tbl = Table()
	init_tbl['x_0'] = sources['xcentroid'][flux_range]
	init_tbl['y_0'] =  sources['ycentroid'][flux_range]
	init_tbl['flux_0'] =  sources['flux'][flux_range]
	
	# And now, let's make a plot of the original image.  
	plt.figure(figsize=(9, 9))
	imshow_norm(data_array, interval=PercentileInterval(99.), stretch=SqrtStretch())
	plt.colorbar()
	plt.savefig(data_file+'.png', dpi=300)
	plt.clf()

	# And now, let's make a plot of the image showing the positions of the objects. 
	plt.figure(figsize=(9, 9))
	imshow_norm(data_array, interval=PercentileInterval(99.), stretch=SqrtStretch())
	plt.scatter(init_tbl['x_0'], init_tbl['y_0'], s=10, color = 'black')
	plt.colorbar()
	plt.savefig(data_file+'_with_daostarfinder_objects.png', dpi=300)
	plt.clf()

	eval_xshape = int(np.ceil(psf_grid.data.shape[2] / psf_grid.oversampling))
	eval_yshape = int(np.ceil(psf_grid.data.shape[1] / psf_grid.oversampling))

	# And now, let's run the PSF Photometry
	sigma_psf = 3.
	daogroup = DBSCANGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
	mmm_bkg = MMMBackground()
	fit_shape = (eval_yshape, eval_xshape)
	phot = BasicPSFPhotometry(daogroup, mmm_bkg, psf_grid, fit_shape, finder=None, aperture_radius=3.)
	
	# This is the part that takes the longest, so I print out the date/time before and after.
	now = datetime.now()
	dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	print("Starting the fit: date and time = ", dt_string)
	
	tbl = phot(data_array, init_guesses=init_tbl)
	
	now = datetime.now()
	dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	print("Ending the fit: date and time = ", dt_string)
	
	# Now I format the output
	tbl['x_fit'].format = '%.1f'
	tbl['y_fit'].format = '%.1f'
	tbl['flux_fit'].format = '%.4e'
	tbl['flux_unc'].format = '%.4e'
	tbl['x_0_unc'].format = '%.4e'
	tbl['y_0_unc'].format = '%.4e'
	
	diff = phot.get_residual_image()
	hdu_out = fits.PrimaryHDU(diff, header = fheader)
	hdul_out = fits.HDUList([hdu_out])
	hdul_out.writeto(data_file+'_residual.fits')
	
	# And create a residual image from the fit 
	plt.figure(figsize=(9, 9))
	imshow_norm(diff, interval=PercentileInterval(99.), stretch=SqrtStretch())
	#plt.scatter(tbl['x_fit'], tbl['y_fit'], s=80, facecolors='none', edgecolors='r')
	plt.colorbar()
	plt.savefig(data_file+'_residual.png', dpi=300)
	plt.clf()
	
	# Calculate the RA and DEC values from the x_fit and y_fit values.
	RA_fit = np.zeros(len(tbl['x_fit']))
	DEC_fit = np.zeros(len(tbl['x_fit']))
	RA_fit, DEC_fit = imagemodel.meta.wcs(tbl['x_fit'], tbl['y_fit'])

	tbl.add_column(DEC_fit, index=0, name='DEC_fit')
	tbl.add_column(RA_fit, index=0, name='RA_fit')
	
	# And write out the table to a file. 
	tbl.write(data_file+'_psf_fit_output.fits')

def main():
	
	# Here's the path to the data 
	data_path = '/data1/car_24_apt_01073/mirage/reduced/'

	nrc_sca = 'nrca1'

	data_files = ['jw01073001001_01101_00001_'+nrc_sca+'_cal',
		'jw01073001001_01101_00002_'+nrc_sca+'_cal',
		'jw01073001002_01101_00003_'+nrc_sca+'_cal',
		'jw01073001002_01101_00004_'+nrc_sca+'_cal',
		'jw01073001003_01101_00005_'+nrc_sca+'_cal',
		'jw01073001003_01101_00006_'+nrc_sca+'_cal',
		'jw01073001004_01101_00007_'+nrc_sca+'_cal',
		'jw01073001004_01101_00008_'+nrc_sca+'_cal']

	# Let's load in the first file 
	hdu = fits.open(data_path+data_files[0]+'.fits')

	# Let's start by loading up the psf grid, which will be 3x3 for now. 
	nrc = webbpsf.NIRCam()
	nrc.filter = hdu[0].header['FILTER']#"F150W"
	if (hdu[0].header['DETECTOR'] == 'NRCALONG'):
		nrc.detector = 'NRCA5'
	elif (hdu[0].header['DETECTOR'] == 'NRCBLONG'):
		nrc.detector = 'NRCB5'
	else:
		nrc.detector = hdu[0].header['DETECTOR']#'NRCA3'
	#nrc_grid = nrc.psf_grid(num_psfs=9, all_detectors=False)
	nrc_grid = nrc.psf_grid(num_psfs=25, all_detectors=False, fov_pixels=33, oversample=2)

	# And let's close the first data file. 
	hdu.close()
	
	# Now, let's spawn processes that will run on the various images.

	number_of_processes_at_the_same_time = 8

	process_list = []
	for i in range(number_of_processes_at_the_same_time):
		data_file_name = data_files[i]
		hdu = fits.open(data_path+data_files[i]+'.fits')
		image = ImageModel(data_path+data_files[i]+'.fits')
		data_array = hdu[1].data
		pheader = hdu[0].header
		fheader = hdu[1].header
		hdu.close()
		#print(data_file_name)
		p = multiprocessing.Process(target=psf_fit, args=(data_array, data_file_name, nrc_grid, fheader, image))
		p.start()
		process_list.append(p)
	
	for process in process_list:
		process.join()


if __name__ == '__main__':

    print('starting main')
    main()
    print('finishing main')
