#! /usr/bin/env python

"""Use source locations to confirm that the correct pixels are being extracted from the full detector when
subarray observations are made.

Input:
subarray slope images
full frame slope images

Method:
Detect sources and create source catalog for each exposure
Get expected offset between full frame and subarray coordinates from pysiaf
Compare subarray source locations + coordinate offsets (+ any dithers) to source locations on full frame data

Note that the sub320, sub640, and full data all have the same set of 4 dithers, so we should be able to
compare exposure 1 of each to each other, and exposure 2, etc.
sub160 has more dithers, but the dither positions of the full frame observations are a subset of those
in the sub160 observations. So we should be able to pull our the appropriate sub160 exposures and compare
those directly as well. The observations with sub64p and sub400p have no dithers, and should be comparable
to the initial exposure in the full data.

For the GRISM256 data, we would need to back out the location of the source in the direct data in order to compare,
since the grism shifts the source in y by 1.449", and we'd need to calculate location in x using the dispersion solution

========================================
* sub160 RAPID (Obs 1)
** Visit 1:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -87.081 -491.730   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.160   +0.144  -87.241 -491.586   +0.160   +0.144 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1    3 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.080   +0.224  -87.161 -491.506   +0.080   +0.224 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1    4 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.240   +0.048  -87.321 -491.682   +0.240   +0.048 SUBDITHER  SCIENCE        0      0  0.238
  1    1   1    5 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -95.274 -483.543   +8.190   +8.190 DITHER     SCIENCE        0      0 11.380
  1    1   1    6 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.350   +8.334  -95.434 -483.399   +8.350   +8.334 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1    7 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.270   +8.414  -95.354 -483.319   +8.270   +8.414 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1    8 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.430   +8.238  -95.514 -483.495   +8.430   +8.238 SUBDITHER  SCIENCE        0      0  0.238
  1    1   1    9 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -95.269 -499.923   +8.190   -8.190 DITHER     SCIENCE        0      0 16.430
  1    1   1   10 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.350   -8.046  -95.429 -499.779   +8.350   -8.046 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1   11 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.270   -7.966  -95.349 -499.699   +8.270   -7.966 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1   12 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.430   -8.142  -95.509 -499.875   +8.430   -8.142 SUBDITHER  SCIENCE        0      0  0.238
  1    1   1   13 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -78.889 -499.918   -8.190   -8.190 DITHER     SCIENCE        0      0 16.620
  1    1   1   14 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.030   -8.046  -79.049 -499.774   -8.030   -8.046 SUBDITHER  SCIENCE        0      0  0.215
  1    1   1   15 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.110   -7.966  -78.969 -499.694   -8.110   -7.966 SUBDITHER  SCIENCE        0      0  0.113
  1    1   1   16 NRCB5_SUB160              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -7.950   -8.142  -79.129 -499.870   -7.950   -8.142 SUBDITHER  SCIENCE        0      0  0.238

========================================
* sub320 RAPID (Obs 2)
** Visit 2:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -87.081 -491.730   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -95.274 -483.543   +8.190   +8.190 DITHER     SCIENCE        0      0 11.582
  1    1   1    3 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -95.269 -499.923   +8.190   -8.190 DITHER     SCIENCE        0      0 16.380
  1    1   1    4 NRCB5_SUB320              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -78.889 -499.918   -8.190   -8.190 DITHER     SCIENCE        0      0 16.380

========================================
* sub640 RAPID (Obs 3)
** Visit 3:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -87.081 -491.730   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -95.274 -483.543   +8.190   +8.190 DITHER     SCIENCE        0      0 11.582
  1    1   1    3 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -95.269 -499.923   +8.190   -8.190 DITHER     SCIENCE        0      0 16.380
  1    1   1    4 NRCB5_SUB640              1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -78.889 -499.918   -8.190   -8.190 DITHER     SCIENCE        0      0 16.380

========================================
* full (Obs 4)
** Visit 4:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  1    1   1    1 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +0.000   +0.000  -82.287 -496.206   +0.000   +0.000 TARGET     SCIENCE        0      0  0.000 (base)
  1    1   1    2 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   +8.190  -90.486 -488.025   +8.190   +8.190 DITHER     SCIENCE        0      0 11.582
  1    1   1    3 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   +8.190   -8.190  -90.468 -504.405   +8.190   -8.190 DITHER     SCIENCE        0      0 16.380
  1    1   1    4 NRCBS_FULL                1 LMC-ASTR  +80.48750  -69.49750   +0.000   +0.000   -8.190   -8.190  -74.088 -504.387   -8.190   -8.190 DITHER     SCIENCE        0      0 16.380

========================================
* sub400p (Obs 5)
** Visit 5:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  0    0   0    0 NRCB5_TAPSIMG32           2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -152.003 -440.340   +0.000   +0.000 TARGET     T_ACQ          0      0  0.000 (base)
  1    1   1    1 NRCB1_SUB400P             2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -146.051 -432.206   +0.000   +0.000 TARGET     SCIENCE        0      0 10.079

========================================
* sub64p (Obs 6)
** Visit 6:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  0    0   0    0 NRCB5_TAPSIMG32           2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -152.003 -440.340   +0.000   +0.000 TARGET     T_ACQ          0      0  0.000 (base)
  1    1   1    1 NRCB1_SUB64P              2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000 -151.093 -427.403   +0.000   +0.000 TARGET     SCIENCE        0      0 12.969

========================================
* subgrism256 (Obs 7)
** Visit 7:1
Tar Tile Exp Dith Aperture Name             Target       RA         Dec        BaseX    BaseY    DithX    DithY    V2       V3       IdlX     IdlY   Level      Type      ExPar  DkPar   dDist
  0    0   0    0 NRCA5_TAGRISMTS32         2 LMC-ASTR  +80.48258  -69.49396   +0.000   +0.000   +0.000   +0.000  +73.262 -551.572   +0.000   +0.000 TARGET     T_ACQ          0      0  0.000 (base)
  1    1   1    1 NRCA5_GRISM256_F322W2     2 LMC-ASTR  +80.48258  -69.49396   +0.000   +1.449   +0.000   +0.000  +50.773 -555.028   +0.000   +1.449 TARGET     SCIENCE        0      0 22.752


"""
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
import pysiaf


def run():
    """MAIN FUNCTION
    """
    # Determine which full frame image to compare with which subarray images
    # Detect sources in full-frame image
    # Detect sources in subarray image
    # Get subarray location from pysiaf
    # Add pysiaf result to subarray catalog
    # Compare full frame catalog to shifted subarray catalog to see how well they match.
    # For subarrays not centered on the detector (sub64p/sub400p), this will fail because the reference locations are different!!!
        # In that case, you'll have to use knowledge of the refrence locations to predict source locations
    # Compare


def find_sources(data, threshold=30, show_sources=True, plot_name='sources.png'):
    """
    Parameters
    ----------
    data : numpy.ndarray
        2D image

    threshold : float
        Number of sigma above the background used to determine the
        presence of a source

    show_sources : bool
        If True, create an image with overlain marks showing sources

    plot_name : str
        Name of file to save the image with marked sources

    Returns
    -------
    sources : astropy.table.Table
        Table of source positions
    """
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=threshold*std)
    sources = daofind(data - median)

    if show_sources:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.)
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.savefig(plot_name)

    return sources

def subarray_corners(aperture_name):
    """Find the offset in x, y between a full frame and the given subarray
    This function has been stolen from Mirage

    Parameters
    ----------
    subarray_name : str
        Name of aperture

    Returns
    -------
    dx : int
        X-location on the full frame of the subarray lower left corner

    dy : int
        Y-location on the full frame of the subarray lower left corner
    """
    # Get SIAF instance
    siaf = pysiaf.Siaf('nircam')

    # get master aperture names
    siaf_detector_layout = pysiaf.iando.read.read_siaf_detector_layout()
    master_aperture_names = siaf_detector_layout['AperName'].data

    # read pysiaf aperture definition file containing DetRef and SciRef values
    siaf_aperture_definitions = pysiaf.iando.read.read_siaf_aperture_definitions('nircam')

    # aperture object
    aperture = siaf[aperture_name]

    # aperture corners in SIAF detector coordinates
    x_det, y_det = aperture.corners('det', rederive=True)

    # determine parent aperture, i.e. the underlying full frame SCA aperture
    index = siaf_aperture_definitions['AperName'].tolist().index(aperture_name)
    aperture._parent_apertures = siaf_aperture_definitions['parent_apertures'][index]

    # If multiuple apertures are listed as parents keep only the first
    if ';' in aperture._parent_apertures:
        print('Multiple parent apertures: {}'.format(aperture._parent_apertures))
        aperture._parent_apertures = aperture._parent_apertures.split(';')[0]

    if aperture_name in master_aperture_names:
        # if master aperture, use it directly to transform to science frame
        x_sci, y_sci = aperture.det_to_sci(x_det, y_det)
    elif aperture._parent_apertures is not None:
        # use parent aperture for transformation
        if verbose:
            print('Using parent {} for {}'.format(aperture._parent_apertures, aperture_name))
        x_sci, y_sci = siaf[aperture._parent_apertures].det_to_sci(x_det, y_det)
        aperture = siaf[aperture._parent_apertures]

    if instrument.lower() == 'nircam':
        if aperture.DetSciParity == 1:
            corner_index = np.array([1, 3])
        elif aperture.DetSciParity == -1:
            # NIRCam will always fall in here, except in case of non-dms orientation
            corner_index = np.array([0, 2])
        x_corner = x_sci[corner_index]
        y_corner = y_sci[corner_index]
    elif instrument.lower() == 'niriss':
        x_corner_index = np.array([0, 2])
        y_corner_index = np.array([0, 2])
        if aperture_name == 'NIS_CEN_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([3, 1])
        x_corner = x_sci[x_corner_index]
        y_corner = y_sci[y_corner_index]
        if aperture_name in ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256']:
            x_corner = [1, 2048]
            y_corner = [1, 2048]
    elif instrument.lower() == 'fgs':
        x_corner_index = np.array([0, 2])
        y_corner_index = np.array([0, 2])
        if aperture_name == 'FGS1_FULL_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([3, 1])
        if aperture_name == 'FGS2_FULL_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([1, 3])
        x_corner = x_sci[x_corner_index]
        y_corner = y_sci[y_corner_index]
    else:
        raise NotImplementedError(("Instrument {} not supported for SIAF subarray corners"
                                   .format(instrument)))

    # account for mirage conventions (e.g. 0-based indexing)
    # we also want integer values as these will be indexes
    x_corner = np.array([np.ceil(x_corner[0]) - 1, np.floor(x_corner[1]) - 1])
    y_corner = np.array([np.ceil(y_corner[0]) - 1, np.floor(y_corner[1]) - 1])
    return x_corner.astype(np.int), y_corner.astype(np.int)

