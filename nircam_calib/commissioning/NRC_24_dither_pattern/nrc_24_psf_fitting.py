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

data_path = '/data1/car_24_apt_01073/mirage/reduced/'
data_file = 'jw01073001001_01101_00001_nrca1_cal' # sys.argv[1]#

# Jarron's Convert_sci_to_sky function! 
def convert_sci_to_sky(xsci, ysci, hdul):
    """
    Convert science pixel coordinates to RA/Dec

    Uses header info from an input file to generate siaf aperture object,
    and uses the pointing information to then create attitude matrix,
    which combines with pysiaf distortion information to convert from 
    'sci' pixel coordinates to 'sky' coordinates.
    
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
    pheader = hdul[0].header
    fheader = hdul[1].header
    #hdul.close

    # Get SIAF info
    siaf = pysiaf.Siaf('NIRCam')
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

# Let's load in the file 
hdu = fits.open(data_path+data_file+'.fits')

#hdu.info()
#Filename: jw01073001001_01101_00001_nrca1_cal.fits
#No.    Name      Ver    Type      Cards   Dimensions   Format
#  0  PRIMARY       1 PrimaryHDU     258   ()
#  1  SCI           1 ImageHDU        85   (2048, 2048)   float32
#  2  ERR           1 ImageHDU        10   (2048, 2048)   float32
#  3  DQ            1 ImageHDU        11   (2048, 2048)   int32 (rescales to uint32)
#  4  AREA          1 ImageHDU         9   (2048, 2048)   float32
#  5  VAR_POISSON    1 ImageHDU         9   (2048, 2048)   float32
#  6  VAR_RNOISE    1 ImageHDU         9   (2048, 2048)   float32
#  7  VAR_FLAT      1 ImageHDU         9   (2048, 2048)   float32
#  8  ASDF          1 BinTableHDU     11   1R x 1C   [17197B]

# The data is the hdu[1]
data = hdu[1].data

# Let's run a quick fitting using DAOStarFinder
mean, median, std = sigma_clipped_stats(data, sigma=3.0) 
daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)  
sources = daofind(data - median) 

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
imshow_norm(data, interval=PercentileInterval(99.), stretch=SqrtStretch())
plt.scatter(init_tbl['x_0'], init_tbl['y_0'], s=10, color = 'black')
plt.savefig(data_file+'_with_daostarfinder_objects.png', dpi=300)

# Let's go and load up the psf grid, which will be 3x3 for now. 
nrc = webbpsf.NIRCam()
nrc.filter = hdu[0].header['FILTER']#"F150W"
nrc.detector = hdu[0].header['DETECTOR']#'NRCA3'
nrc_grid = nrc.psf_grid(num_psfs=9, all_detectors=False)

eval_xshape = int(np.ceil(nrc_grid.data.shape[2] / nrc_grid.oversampling))
eval_yshape = int(np.ceil(nrc_grid.data.shape[1] / nrc_grid.oversampling))

# And now, let's run the PSF Photometry
sigma_psf = 3.
daogroup = DBSCANGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fit_shape = (eval_yshape, eval_xshape)
phot = BasicPSFPhotometry(daogroup, mmm_bkg, nrc_grid, fit_shape, finder=None, aperture_radius=3.)

# This is the part that takes the longest, so I print out the date/time before and after.
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print("Starting the fit: date and time = ", dt_string)

tbl = phot(data, init_guesses=init_tbl)

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

# Let's convert the PSF photometry fit X and Y values to RA and DEC values
RA_fit = np.zeros(len(tbl['x_fit']))
DEC_fit = np.zeros(len(tbl['x_fit']))
for star in range(len(tbl['x_fit'])):
	RA_fit[star], DEC_fit[star] = convert_sci_to_sky(tbl['x_fit'][star], tbl['y_fit'][star], hdu)
	
tbl.add_column(DEC_fit, index=0, name='DEC_fit')
tbl.add_column(RA_fit, index=0, name='RA_fit')

# And write out the table to a file. 
tbl.write(data_file+'_psf_fit_output.fits')
