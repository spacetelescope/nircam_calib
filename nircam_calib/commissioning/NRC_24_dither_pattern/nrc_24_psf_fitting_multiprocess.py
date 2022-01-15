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

import multiprocessing

def convert_sci_to_sky(xsci, ysci, pheader, fheader, siaf):
    """
    Convert science pixel coordinates to RA/Dec

    Uses header info from an input file to generate siaf aperture object,
    and uses the pointing information to then create attitude matrix,
    which combines with pysiaf distortion information to convert from 
    'sci' pixel coordinates to 'sky' coorinates.
    
    Parameters
    ==========
    xsci : float or ndarray
        pixel positions along x-axis ('sci' orientation of image)
    ysci : float or ndarray
        pixel positions along y-axis ('sci' orientation of image)
    file : string
        Location of DMS file containing necessary header info

    Returns
    =======
    RA and Dec arrays.
    """

    #hdul = fits.open(file)
    #pheader = hdul[0].header
    #fheader = hdul[1].header
    #hdul.close

    # Get SIAF info
    apername = pheader.get('APERNAME')
    apsiaf = siaf[apername]

    # V3 PA
    pa_v3 = fheader.get('V3I_YANG')
    # RA/Dec located at aperture reference position
    ra_ref = fheader.get('RA_REF')
    dec_ref = fheader.get('DEC_REF')
    # V2/V3 aperture reference
    v2_ref, v3_ref = apsiaf.reference_point('tel')

    # Create attitude matrix
    att = pysiaf.utils.rotations.attitude(v2_ref, v3_ref, ra_ref, dec_ref, pa_v3)
    apsiaf.set_attitude_matrix(att)

    # Conver to sky coordinates
    ra_deg, dec_deg = apsiaf.convert(xsci, ysci, 'sci', 'sky')

    return (ra_deg, dec_deg)

# This is the function that actually runs the PSF fitting on the data files from the same SCA
def psf_fit(data_array, pheader, fheader, data_file, psf_grid, siaf):

	# Let's run a quick fitting using DAOStarFinder
	mean, median, std = sigma_clipped_stats(data_array, sigma=3.0) 
	daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)  
	sources = daofind(data_array - median) 
	
	# And now, let's make some cuts to remove objects that are too faint or too bright
	flux_min = 5
	flux_max = 50
	flux_range = np.where((sources['flux'] > flux_min) & (sources['flux'] < flux_max))[0]
	
	init_tbl = Table()
	init_tbl['x_0'] = sources['xcentroid'][flux_range]
	init_tbl['y_0'] =  sources['ycentroid'][flux_range]
	init_tbl['flux_0'] =  sources['flux'][flux_range]
	
	# And now, let's make a plot of the image showing the positions of the objects. 
	plt.figure(figsize=(9, 9))
	imshow_norm(data_array, interval=PercentileInterval(99.), stretch=SqrtStretch())
	plt.scatter(init_tbl['x_0'], init_tbl['y_0'], s=10, color = 'black')
	plt.savefig(data_file+'_with_daostarfinder_objects.png', dpi=300)

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
	
	# And create a residual image from the fit 
	plt.figure(figsize=(9, 9))
	imshow_norm(diff, interval=PercentileInterval(99.), stretch=SqrtStretch())
	#plt.scatter(tbl['x_fit'], tbl['y_fit'], s=80, facecolors='none', edgecolors='r')
	plt.colorbar()
	plt.savefig(data_file+'_residual.png', dpi=300)
	
	RA_fit = np.zeros(len(tbl['x_fit']))
	DEC_fit = np.zeros(len(tbl['x_fit']))
	RA_fit, DEC_fit = convert_sci_to_sky(tbl['x_fit'], tbl['y_fit'], pheader, fheader, siaf)

	tbl.add_column(DEC_fit, index=0, name='DEC_fit')
	tbl.add_column(RA_fit, index=0, name='RA_fit')
	
	# And write out the table to a file. 
	tbl.write(data_file+'_psf_fit_output.fits')

def main():
	
	# Here's the path to the data 
	data_path = '/data1/car_24_apt_01073/mirage/reduced/'

	data_files = ['jw01073001001_01101_00001_nrca1_cal',
		'jw01073001001_01101_00002_nrca1_cal',
		'jw01073001002_01101_00003_nrca1_cal',
		'jw01073001002_01101_00004_nrca1_cal',
		'jw01073001003_01101_00005_nrca1_cal',
		'jw01073001003_01101_00006_nrca1_cal',
		'jw01073001004_01101_00007_nrca1_cal',
		'jw01073001004_01101_00008_nrca1_cal']

	# Let's load in the first file 
	#hdu = fits.open('jw01073001001_01101_00001_nrca1_cal.fits')
	hdu = fits.open(data_path+data_files[0]+'.fits')

	# Let's start by loading up the psf grid, which will be 3x3 for now. 
	nrc = webbpsf.NIRCam()
	nrc.filter = hdu[0].header['FILTER']#"F150W"
	nrc.detector = hdu[0].header['DETECTOR']#'NRCA3'
	#nrc_grid = nrc.psf_grid(num_psfs=9, all_detectors=False)
	nrc_grid = nrc.psf_grid(num_psfs=25, all_detectors=False, fov_pixels=33, oversample=2)

	# And let's close the first data file. 
	hdu.close()
	
	# Now, let's spawn processes that will run on the various images.

	number_of_processes_at_the_same_time = 8

	# Create an instance of siaf to pass through along to the psf_fit function. 
	siaf = pysiaf.Siaf('NIRCam')

	process_list = []
	for i in range(number_of_processes_at_the_same_time):
		data_file_name = data_files[i]
		hdu = fits.open(data_path+data_files[i]+'.fits')
		data_array = hdu[1].data
		pheader = hdu[0].header
		fheader = hdu[1].header
		hdu.close()
		#print(data_file_name)
		p = multiprocessing.Process(target=psf_fit, args=(data_array, pheader, fheader, data_file_name, nrc_grid, siaf))
		p.start()
		process_list.append(p)
	
	for process in process_list:
		process.join()


if __name__ == '__main__':

    print('starting main')
    main()
    print('finishing main')

