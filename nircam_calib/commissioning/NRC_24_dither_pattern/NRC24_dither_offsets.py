#!/usr/bin/env python
#
# NRC 24 Analysis Script
#
#	#################################
#
#	Author: Anton M. Koekemoer, STScI
#
#	#################################
#
#	v1.1	2022-01-31	Added input args and env variables, tidied some up.
#
# Usage:
#	If run without input args, as:
#
#		python NRC24_dither_offsets.py
#
#	then it will read all *_cal.fits files from path "pipeline_outputs_stage2"
#	and place all output files, plots etc in path "analysis_dir"




import os, glob, sys, argparse
from glob import glob
import shutil
import urllib


# Imports for pysiaf, astropy, matplotlib, math, etc.
import pysiaf
from astropy.io import ascii as asc
from astropy.io import fits
from matplotlib import cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator

import numpy as np

from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import sigma_clipped_stats
from astropy.table import Table

from photutils import BasicPSFPhotometry
from photutils import DBSCANGroup, MMMBackground
from photutils.detection import DAOStarFinder

from math import *
import numpy as np

# Mirage Imports (nircam_simulator)
from mirage import imaging_simulator
from mirage.catalogs import create_catalog
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator



# Imports for JWST pipeline
import jwst

from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import calwebb_image2
from jwst.pipeline import calwebb_image3
from jwst import assign_wcs


# List of possible data quality flags
from jwst.datamodels import dqflags

# The entire calwebb_detector1 pipeline
from jwst.pipeline import calwebb_detector1

# Individual steps that make up calwebb_detector1
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.ipc import IPCStep
from jwst.refpix import RefPixStep
from jwst.linearity import LinearityStep
from jwst.persistence import PersistenceStep
from jwst.dark_current import DarkCurrentStep
from jwst.jump import JumpStep
from jwst.ramp_fitting import RampFitStep
from jwst import datamodels

# Individual steps that make up calwebb_image2
from jwst.background import BackgroundStep
from jwst.assign_wcs import AssignWcsStep
from jwst.flatfield import FlatFieldStep
from jwst.photom import PhotomStep
from jwst.resample import ResampleStep
from jwst import datamodels

# Individual steps that make up calwebb_image3
from jwst.tweakreg import TweakRegStep
from jwst.skymatch import SkyMatchStep
from jwst.outlier_detection import OutlierDetectionStep
from jwst.resample import ResampleStep
from jwst.source_catalog import SourceCatalogStep
from jwst import datamodels
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

# Needed for multiprocessing
import multiprocessing










#========================++++++++++++
def setup_dither_pattern(prop_obsid):
#========================++++++++++++
  #
  # Set up the dither offset patterns, and x,y boundaries for the final plots
  #
  # This is highly specific to Program ID 1073.
  #
  #
  #		Primary			Secondary
  #
  #	Obs 1	FULLBOX 8NIRSPEC
  #
  #	Obs 2	FULLBOX 4
  #
  #	Obs 3	INTRAMODULEBOX 4
  #
  #	Obs 4	--			STANDARD 2
  #
  #	Obs 5	--			STANDARD 4
  #
  #	Obs 6	--			STANDARD 9
  #
  #	Obs 7	INTRAMODULEBOX 4	STANDARD 2
  #
  #
  # Can read these from the pointing file produced by APT, or from the header for each FITS file
  #	*cal.fits[0]
  #			XOFFSET
  #			YOFFSET
  #
  # Hardcode them here for now, just to be robust against possible changes in APT.
  #
  #
  if (prop_obsid == 'jw01073001'):		# FULLBOX 8NIRSPEC
    #
    xlim = [-100.,100.]
    ylim = [-100.,100.]
    xmajor = 25
    ymajor = 25
    xminor = 5
    yminor = 5
    #
    dither_offsets = [ [-24.6,  -64.1],
                       [-24.4,  -89.0],
                       [24.6,   -88.8],
                       [24.4,   -63.9],
                       [24.6,    64.1],
                       [24.4,    89.0],
                       [-24.6,   88.8],
                       [-24.4,   63.9] ]
  #
  if (prop_obsid == 'jw01073002'):		# FULLBOX 4
    #
    xlim = [-30.,30.]
    ylim = [-30.,30.]
    xmajor = 10
    ymajor = 10
    xminor = 1
    yminor = 1
    #
    dither_offsets = [ [-24.35, -12.58],
                       [-24.55,  12.38],
                       [24.35,   12.58],
                       [24.55,  -12.38] ]
  #
  if (prop_obsid == 'jw01073003'):		# INTRAMODULEBOX 4
    #
    xlim = [-4.,4.]
    ylim = [-4.,4.]
    xmajor = 1
    ymajor = 1
    xminor = 0.1
    yminor = 0.1
    #
    dither_offsets = [ [-2.9,  -3.1],
                       [-3.1,   2.9],
                       [2.9,    3.1],
                       [3.1,   -2.9] ]
  #
  if (prop_obsid == 'jw01073004'):		# STANDARD 2
    #
    xlim = [-0.1,0.4]
    ylim = [-0.1,0.4]
    xmajor = 0.1
    ymajor = 0.1
    xminor = 0.01
    yminor = 0.01
    #
    dither_offsets = [ [0.,      0.],
                       [0.2262,  0.1659] ]
  #
  if (prop_obsid == 'jw01073005'):		# STANDARD 4
    #
    xlim = [-0.1,0.4]
    ylim = [-0.1,0.4]
    xmajor = 0.1
    ymajor = 0.1
    xminor = 0.01
    yminor = 0.01
    #
    dither_offsets = [ [0.,      0.],
		       [0.1509,  0.1357],
		       [0.0754,  0.2112],
		       [0.2262,  0.0452] ]
  #
  if (prop_obsid == 'jw01073006'):		# STANDARD 9
    #
    xlim = [-0.1,0.4]
    ylim = [-0.1,0.4]
    xmajor = 0.1
    ymajor = 0.1
    xminor = 0.01
    yminor = 0.01
    #
    dither_offsets = [ [0.,     0.],
                       [0.0302, 0.1308],
                       [0.0603, 0.2916],
                       [0.1006, 0.0603],
                       [0.1609, 0.1609],
                       [0.1911, 0.3218],
                       [0.2618, 0.0302],
                       [0.2916, 0.1308],
                       [0.3219, 0.2916] ]
  #
  if (prop_obsid == 'jw01073007'):		# INTRAMODULEBOX 4 + STANDARD 2
    #
    xlim = [-4.,4.]
    ylim = [-4.,4.]
    xmajor = 1
    ymajor = 1
    xminor = 0.1
    yminor = 0.1
    #
    dither_offsets = [ [-2.9,           -3.1],
                       [-2.9 + 0.2262,  -3.1 + 0.1659],
                       [-3.1,            2.9],
                       [-3.1 + 0.2262,   2.9 + 0.1659],
                       [2.9,             3.1],
                       [2.9 + 0.2262,    3.1 + 0.1659],
                       [3.1,            -2.9],
                       [3.1 + 0.2262,   -2.9 + 0.1659] ]
  #
  #
  return dither_offsets, xlim, ylim, xmajor, ymajor, xminor, yminor










#=================================
def plot_title(prop_obsid_filter):
#=================================

  propid = prop_obsid_filter[3:7]
  obs = prop_obsid_filter[7:10]
  filtername = prop_obsid_filter[11:16]
  #
  if (obs == '001'): pattern = 'FULLBOX 8NIRSPEC'
  if (obs == '002'): pattern = 'FULLBOX 4'
  if (obs == '003'): pattern = 'INTRAMODULEBOX 4'
  if (obs == '004'): pattern = 'STANDARD 2'
  if (obs == '005'): pattern = 'STANDARD 4'
  if (obs == '006'): pattern = 'STANDARD 9'
  if (obs == '007'): pattern = 'INTRAMODULEBOX 4 + STANDARD 2'
  #
  title = 'NRC-24  '+propid+'-'+obs+'  '+filtername+'     '+pattern

  fig.suptitle(title, weight='bold')










#=========================
if __name__ == '__main__':
#=========================

  # There are two choices for the type of analysis to be done:
  #
  #	'absolute'	- aligns all exposures to external catalog, currently LMC
  #			- best for displaying the actual dither offsets (eg for patterns that have no (0,0) dither)
  #			- this is the option most redcently tested
  #
  #	'relative'	- aligns all exposures to the first one
  #			- best for cases where an external catalog might not be available (eg other areas of the sky)
  #			- limited in that all dithers are relative to first exposure, ie won't match dithers for patterns that have bo (0,0) dither
  #			- this was implemented earlier but hasn't been fully tested in a while.


  parser = argparse.ArgumentParser(description='Run NRC24_dither_offsets script.')
  parser.add_argument('-x',  '--xmlfile',        default='nrc24-1073_same_PA_141deg31.xml',      type=str, help='Input xml file from APT.')
  parser.add_argument('-p',  '--pointing',       default='nrc24-1073_same_PA_141deg31.pointing', type=str, help='Input pointing file from APT.') 
  parser.add_argument('-a',  '--analysis_type',  default='absolute', type=str, help='Type of analysis, either "absolute" or "relative.')
  parser.add_argument('-r',  '--refcat',         default='lmc_catalog_flag1.cat', type=str, help='Reference astrometric catalog.')

  options = parser.parse_args()

  xml_file       = options.xmlfile
  pointing_file  = options.pointing
  analysis_type  = options.analysis_type
  catfile_lmc    = options.refcat

  # Environment variables needed
  #
  pipeline_outputs_stage2 = os.getenv('pipeline_outputs_stage2')
  analysis_dir            = os.getenv('analysis_dir')
  #
  if (pipeline_outputs_stage2 == None):  pipeline_outputs_stage2 = './'
  if (analysis_dir == None):  analysis_dir = './'

  if (not os.path.exists(analysis_dir)):  os.mkdir(analysis_dir)




  # Read in the LMC catalog
  # -----------------------
  #
  if (analysis_type == 'absolute') and (not os.path.exists(catfile_lmc))):
    #
    print('For analysis_type="absolute", need to have a catalog! Did not find catalog file: ',catfile_lmc)
    sys.exit()
  #
  catalog_lmc = Table.read(catfile_lmc, format='ascii')
  #
  nsources_lmc = len(catalog_lmc)



 
  # Convert the LMC catalog to the Ideal frame, and also to V2,V3, for the nominal pointing with dither (0,0)
  # ---------------------------------------------------------------------------------------------------------
  #
  catfile_lmc_radec_xyidl_v2v3 = analysis_dir+'lmc_catalog_radec_xyidl_v2v3.csv'
  #
  if (not os.path.exists(catfile_lmc_radec_xyidl_v2v3)):
    #
    # Determine spacecraft attitude for "baseline" exposure, with dither (0,0)
    # -------------------------------------------------------------------------
    #
    #filename_dither_x0y0 = 'jw01073004001_01101_00001_nrca5_cal.fits'	# for LW
    filename_dither_x0y0 = 'jw01073004001_01101_00001_nrca1_cal.fits'	# for SW
    #
    hdr0 = fits.getheader(pipeline_outputs_stage2+filename_dither_x0y0, 0)
    hdr  = fits.getheader(pipeline_outputs_stage2+filename_dither_x0y0, 1)
    #
    instrument      = hdr0['INSTRUME']
    apername        = hdr0['APERNAME']
    #
    ra_ref_orig     = hdr['RA_REF']
    dec_ref_orig    = hdr['DEC_REF']
    roll_ref_orig   = hdr['ROLL_REF']
    #
    siaf = pysiaf.Siaf(instrument)
    #
    aper_orig = siaf[apername]
    #
    v2_ref_orig = aper_orig.V2Ref     # same as hdr['V2_REF']
    v3_ref_orig = aper_orig.V3Ref     # same as hdr['V2_REF']
    #
    attitude = pysiaf.utils.rotations.attitude(v2_ref_orig, v3_ref_orig, ra_ref_orig, dec_ref_orig, roll_ref_orig)
    #
    aper_orig.set_attitude_matrix(attitude)			# apply "attitude matrix" to this aperture
    #
    print(aper_orig)
    print(attitude)
    #
    aper_NRCALL_FULL = siaf['NRCALL_FULL']
    aper_NRCALL_FULL.set_attitude_matrix(attitude)	# apply "attitude matrix" to this aperture
    #
    catalog_lmc_radec_xyidl_v2v3 = Table(names=('id', 'ra',  'dec', 'x_idl', 'y_idl', 'v2',    'v3',    'mag_f150w'),
                                          dtype=(int,  float, float, float,   float,   float,   float,   float))
    for i in range(nsources_lmc):
      #
      idnumber, ra, dec, mag = catalog_lmc['index'][i], catalog_lmc['x_or_RA'][i], catalog_lmc['y_or_Dec'][i], catalog_lmc['nircam_f150w_magnitude'][i]
      #
      x_idl, y_idl      = aper_NRCALL_FULL.convert(ra, dec, 'sky', 'idl')		# identical to aper_NRCALL_FULL.sky_to_idl(ra, dec)
      v2, v3            = aper_NRCALL_FULL.convert(ra, dec, 'sky', 'tel')		# identical to aper_NRCALL_FULL.sky_to_tel(ra, dec)
      #
      catalog_lmc_radec_xyidl_v2v3.add_row([idnumber, ra, dec, x_idl, y_idl, v2, v3, mag])
    #
    print('')
    print('Writing out  ',catfile_lmc_radec_xyidl_v2v3)
    print('')
    a = os.system('/bin/rm -f '+catfile_lmc_radec_xyidl_v2v3)
    #
    catalog_lmc_radec_xyidl_v2v3.write(catfile_lmc_radec_xyidl_v2v3)
    #
  else:
    #
    print('')
    print('Reading in  ',catfile_lmc_radec_xyidl_v2v3)
    print('')
    catalog_lmc_radec_xyidl_v2v3 = Table.read(catfile_lmc_radec_xyidl_v2v3)
  #
  print('nsources_lmc      = ',nsources_lmc)
  print('')





  # Create a trimmed version of the LMC catalog:
  # --------------------------------------------
  #
  # Fairly quick and simple, just do an 8 arcmin radius around TARG coords.
  # Just hardocde these values for now.
  #
  ra_targ   = 80.4875
  dec_targ = -69.4975
  #
  catfile_lmc_radec_xyidl_v2v3_trim = analysis_dir+'lmc_catalog_radec_xyidl_v2v3_trim.csv'
  #
  if (not os.path.exists(catfile_lmc_radec_xyidl_v2v3_trim)):
    #
    catalog_lmc_radec_xyidl_v2v3_trim = Table(names=('id', 'ra',  'dec', 'x_idl', 'y_idl', 'v2',    'v3',    'mag_f150w'),
                                               dtype=(int,  float, float, float,   float,   float,   float,   float))
    #
    cosdec = cos(np.deg2rad(dec_targ))
    #
    for i in range(nsources_lmc):
      #
      ra, dec, mag = catalog_lmc_radec_xyidl_v2v3['ra'][i], catalog_lmc_radec_xyidl_v2v3['dec'][i], catalog_lmc_radec_xyidl_v2v3['mag_f150w'][i]
      #
      dra  = (ra - ra_targ) * cosdec
      ddec = dec - dec_targ
      #
      dr = sqrt(dra**2 + ddec**2) * 60.	# convert to arcmin
      #
      if ((dr < 8.) and (mag < 20.)):
        #
        catalog_lmc_radec_xyidl_v2v3_trim.add_row(catalog_lmc_radec_xyidl_v2v3[i])
    #
    print('Writing out  ',catfile_lmc_radec_xyidl_v2v3_trim)
    print('')
    a = os.system('/bin/rm -f '+catfile_lmc_radec_xyidl_v2v3_trim)
    #
    catalog_lmc_radec_xyidl_v2v3_trim.write(catfile_lmc_radec_xyidl_v2v3_trim)
    #
  else:
    #
    print('Reading in  ',catfile_lmc_radec_xyidl_v2v3_trim)
    print('')
    catalog_lmc_radec_xyidl_v2v3_trim = Table.read(catfile_lmc_radec_xyidl_v2v3_trim)
  #
  nsources_lmc_trim = len(catalog_lmc_radec_xyidl_v2v3_trim)
  #
  print('nsources_lmc_trim = ',nsources_lmc_trim)
  print('')









  # Determine the list of observations and filters:
  # ===============================================
  #
  # First set up the dictionary "proposal_obsid":
  # 	- each entry will be a list of "proposal_obsid_filter"
  #	- each "proposal_obsid_filter" will be a list of "proposal_obsid_filter_expnum"
  #
  filenames = sorted(glob(pipeline_outputs_stage2+'*_cal.fits'))
  #
  filename_dict = {}
  #
  for filename_full in filenames:
    #
    filename       = str.split(filename_full,'/')[-1]
    filename_split = str.split(filename,'_')
    prop_obsid     = filename_split[0][:10]		# eg jw01073001
    expnum         = filename_split[2]			# eg 00001
    hdr            = fits.getheader(filename_full)
    filtername     = hdr['FILTER']
    rootname       = filename[:-5]
    #
    filename_current_dict = {}
    filename_current_dict['prop_obsid'] = prop_obsid
    filename_current_dict['filename']   = filename
    filename_current_dict['rootname']   = rootname
    filename_current_dict['expnum']     = expnum
    filename_current_dict['filtername'] = filtername
    #
    filename_dict[filename_full] = filename_current_dict
  #
  #
  # Create the list of prop_obdsid's
  # --------------------------------
  #
  prop_obsid_list = []
  #
  for filename_full in filenames:
    #
    prop_obsid = filename_dict[filename_full]['prop_obsid']
    #
    if (not prop_obsid in prop_obsid_list):
      #
      prop_obsid_list.append(prop_obsid)






  # tolerance for matching detected sources to reference catalog
  #
  tol = 0.10 / 3600.	# Set maximum tolerance for matching detected sources to the existing astrometric catalog




  # MAIN LOOP OVER "prop_obsid"
  # ===========================
  #
  for prop_obsid in prop_obsid_list:
    #
    dither_offsets, xlim, ylim, xmajor, ymajor, xminor, yminor = setup_dither_pattern(prop_obsid)
    #
    print('')
    print('')
    print(prop_obsid)
    print('==========')
    #
    filenames = sorted(glob(pipeline_outputs_stage2+prop_obsid+'*_cal.fits'))
    #
    # For each prop_obsid, build the list of filters
    # ----------------------------------------------
    #
    prop_obsid_filter_list = []
    #
    for filename_full in filenames:
      #
      filtername = filename_dict[filename_full]['filtername']
      #
      prop_obsid_filter = '_'.join([prop_obsid,filtername])
      #
      if (not prop_obsid_filter in prop_obsid_filter_list):
        #
        prop_obsid_filter_list.append(prop_obsid_filter)
    #
    #
    # For each filter, build the list of exposures
    # --------------------------------------------
    #
    prop_obsid_filter_expnum_list = []
    #
    for prop_obsid_filter in prop_obsid_filter_list:
      #
      for filename_full in filenames:
        #
        filtername = filename_dict[filename_full]['filtername']
        #
        if (filtername == prop_obsid_filter[-5:]):
          #
          expnum   = filename_dict[filename_full]['expnum']
          #
          prop_obsid_filter_expnum = prop_obsid_filter + '_' + expnum
          #
          if (not prop_obsid_filter_expnum in prop_obsid_filter_expnum_list):
            #
            prop_obsid_filter_expnum_list.append(prop_obsid_filter_expnum)
    #
    #
    # For each exposure, build the list of fits files (one for each SCA)
    # ------------------------------------------------------------------
    #
    prop_obsid_filter_dict = {}
    #
    prop_obsid_filter_expnum_dict = {}
    #
    for prop_obsid_filter_expnum in prop_obsid_filter_expnum_list:
      #
      filtername = prop_obsid_filter_expnum[-11:-6]
      expnum     = prop_obsid_filter_expnum[-5:]
      #
      for filename_full in filenames:
        #
        filename_filter = filename_dict[filename_full]['filtername']
        #
        if (filename_filter == filtername):
          #
          filename        = filename_dict[filename_full]['filename']
          filename_expnum = filename_dict[filename_full]['expnum']
          #
          if (filename_expnum == expnum):
            #
            # Now, finally, fill up the dictonary of exposures
            #
            if (not prop_obsid_filter_expnum in prop_obsid_filter_expnum_dict.keys()):
              #
              prop_obsid_filter_expnum_dict[prop_obsid_filter_expnum] = [filename_full]
              #
            else:
              #
              prop_obsid_filter_expnum_dict[prop_obsid_filter_expnum].append(filename_full)




    # LOOP over each prop_obsid_filter
    # ================================
    #
    # For each one, loop over all the exposures for that filter and generate catalogs.
    #
    for prop_obsid_filter in prop_obsid_filter_list:
      #
      print('')
      print(prop_obsid_filter)
      print('----------------')
      #
      prop_obsid_filter_expnum_list = prop_obsid_filter_expnum_dict.keys()
      #
      first_exposure_found = False
      first_exposure_prop_obsid_filter_expnum = ''
      #
      # Initialize the "reference" catalog:
      #	- This will be the list of RA,Dec positions of all sources detected in the first exposure
      #
      catalog_reference = Table(names=('id', 'ra',  'dec', 'x_idl', 'y_idl', 'v2',  'v3'),
                                dtype=(str,   float, float, float,   float,   float, float))
      #
      catfile_all_matched = analysis_dir + prop_obsid_filter + '_all_matched_cat.csv'
      #
      if (os.path.exists(catfile_all_matched)):
        #
        catalog_all_matched = Table.read(catfile_all_matched)
        #
        dither_offsets_dict = {}
        #
        for prop_obsid_filter_expnum in prop_obsid_filter_expnum_list:
          #
          if (prop_obsid_filter in prop_obsid_filter_expnum):	# be sure to only select the current filter
            #
            #print(prop_obsid_filter_expnum)
            #
            current_prop_obsid_filter_expnum_list = prop_obsid_filter_expnum_dict[prop_obsid_filter_expnum]
            #
            filename_full = current_prop_obsid_filter_expnum_list[0]
            #
            hdr0 = fits.getheader(filename_full, 0)
            #
            commanded_xdither = hdr0['XOFFSET']
            commanded_ydither = hdr0['YOFFSET']
            #
            indices_catalog_all_matched_current = np.where(catalog_all_matched['prop_obsid_filter_expnum'] == prop_obsid_filter_expnum)[0].tolist()
            #
            xdither_array = catalog_all_matched['dx_idl'][indices_catalog_all_matched_current]
            ydither_array = catalog_all_matched['dy_idl'][indices_catalog_all_matched_current]
            #
            measured_xdither = np.mean(xdither_array)
            measured_ydither = np.mean(ydither_array)
            #
            difference_xdither_array = xdither_array - commanded_xdither
            difference_ydither_array = ydither_array - commanded_ydither
            #
            difference_xdither = np.mean(difference_xdither_array)
            difference_ydither = np.mean(difference_ydither_array)
            #
            difference_xdither_rms = np.std(difference_xdither_array)
            difference_ydither_rms = np.std(difference_ydither_array)
            #
            dither_offsets_dict[prop_obsid_filter_expnum] = {'commanded_xdither':      commanded_xdither,
                                                             'commanded_ydither':      commanded_ydither,
                                                             'measured_xdither':       measured_xdither,
                                                             'measured_ydither':       measured_ydither,
                                                             'difference_xdither':     difference_xdither,
                                                             'difference_ydither':     difference_ydither,
                                                             'difference_xdither_rms': difference_xdither_rms,
                                                             'difference_ydither_rms': difference_ydither_rms}
        #
      else:
        #
        catalog_all_matched = Table(names=('prop_obsid_filter_expnum', 'fitsfilename',  'x_pix',  'y_pix',  'ra',   'dec',  'ra_ref', 'dec_ref',  'x_idl', 'y_idl', 'v2',  'v3',  'dra',  'ddec',  'dx_idl', 'dy_idl', 'dv2', ' dv3',  'dither_dx', 'dither_dy'),
                                    dtype=(str,                         str,             float,    float,    float,  float,  float,    float,      float,   float,   float, float, float,  float,   float,    float,    float,  float,  float,       float))
        #
        #
        dither_offsets_dict = {}
        #
        for prop_obsid_filter_expnum in prop_obsid_filter_expnum_list:
          #
          if (prop_obsid_filter in prop_obsid_filter_expnum):	# be sure to only select the current filter
            #
            catfile_current_exposure = analysis_dir + prop_obsid_filter_expnum + '_allSCAs_cat.csv'
            #
            first_exposure = False
            #
            if (not first_exposure_found):
              #
              first_exposure       = True
              first_exposure_found = True
              first_exposure_prop_obsid_filter_expnum = prop_obsid_filter_expnum
            #
            print('')
            if ( (analysis_type == 'relative') and (first_exposure)):
              print(prop_obsid_filter_expnum, '  First exposure ***')
            else:
              print(prop_obsid_filter_expnum)
            #
            current_prop_obsid_filter_expnum_list = prop_obsid_filter_expnum_dict[prop_obsid_filter_expnum]
            #
            filename_full = current_prop_obsid_filter_expnum_list[0]
            #
            hdr0 = fits.getheader(filename_full, 0)
            #
            dither_xoffset = hdr0['XOFFSET']
            dither_yoffset = hdr0['YOFFSET']
            #
            dither_offsets_dict[prop_obsid_filter_expnum]['commanded_xdither'] = dither_xoffset
            dither_offsets_dict[prop_obsid_filter_expnum]['commanded_ydither'] = dither_yoffset
            #
            if (os.path.exists(catfile_current_exposure)):
              #
              catalog_current_exposure= Table.read(catfile_current_exposure)
              #
            else:
              #
              catalog_current_exposure = Table(names=('prop_obsid_filter_expnum', 'fitsfilename',  'x_pix', 'y_pix', 'ra',   'dec',  'x_idl', 'y_idl', 'v2',  'v3'),
                                               dtype=(str,                         str,             float,   float,   float,  float,  float,   float,   float, float))
              #
              # Loop through all the files for this "prop_obsid_filter_expnum":
              # ==============================================================
              #
              rootnames = []
              #
              current_prop_obsid_filter_expnum_list = prop_obsid_filter_expnum_dict[prop_obsid_filter_expnum]
              #
              for filename_full in current_prop_obsid_filter_expnum_list:
                #
                filename = filename_dict[filename_full]['filename']
                rootname = filename_dict[filename_full]['rootname']
                #
                rootnames.append(rootname)
                #
                catfile_daofind = rootname + '_daofind_cat.csv'
                #
                if (not os.path.exists(analysis_dir+catfile_daofind)):
                  #
                  print('Creating DAOfind catalog:  ',catfile_daofind)
                  #
                  img_data = fits.getdata(filename_full)
                  #
                  #	SNR	nsources
                  #	100	3000
                  #	200	1750
                  #	500	1000
                  #	1500	488
                  #	2000	273
                  #
                  mean, med, rms = sigma_clipped_stats(img_data, sigma=3.0) 
                  daofind = DAOStarFinder(fwhm=3.0, threshold=1500.*rms)  
                  tbl_sources = daofind.find_stars(img_data - med)
                  tbl_sources.write(analysis_dir+catfile_daofind, overwrite=True)
              #
              #
              # Loop through all the SCA files for this exposure and create the necessary catalogs
              # ----------------------------------------------------------------------------------
              #
              for rootname in rootnames:
                #
                filename_full     = pipeline_outputs_stage2 + rootname + '.fits'
                catfile_daofind   = analysis_dir            + rootname + '_daofind_cat.csv'
                #
                print('\nrootname = ',rootname)
                #
                # Read in the x, y pixel coords for this SCA
                # ------------------------------------------
                #
                if (os.path.exists(catfile_daofind)):
                  #
                  sources_tbl = Table.read(catfile_daofind)
                  x_pix_catfile = sources_tbl['xcentroid']  +  1.0
                  y_pix_catfile = sources_tbl['ycentroid']  +  1.0
                #
                nsources = len(x_pix_catfile)
                print('nsources detected = ',nsources)
                print('')


                # Set up SIAF config for this exposure
                # ------------------------------------
                #
                hdr0 = fits.getheader(filename_full, 0)
                hdr  = fits.getheader(filename_full, 1)
                #
                instrument     = hdr0['INSTRUME']
                apername       = hdr0['APERNAME']
                #
                ra_ref         = hdr['RA_REF']
                dec_ref        = hdr['DEC_REF']
                roll_ref       = hdr['ROLL_REF']
                #
                print(instrument)
                print(apername)
                #
                siaf = pysiaf.Siaf(instrument)
                #
                aper = siaf[apername]
                #
                v2_ref = aper.V2Ref     # same as hdr['V2_REF']
                v3_ref = aper.V3Ref     # same as hdr['V3_REF']
                #
                print('ra_ref, dec_ref, roll_ref, v2_ref, v3_ref:  ',ra_ref, dec_ref, roll_ref, v2_ref, v3_ref)
                print('')
                print(aper)
                #
                attitude = pysiaf.utils.rotations.attitude(v2_ref, v3_ref, ra_ref, dec_ref, roll_ref)
                aper.set_attitude_matrix(attitude)
                #
                aper_NRCALL_FULL   = siaf['NRCALL_FULL']
                aper_NRCALL_FULL.set_attitude_matrix(attitude)	# apply "attitude matrix" to this aperture
                #
                # Convert the x,y pixels from detector science frame ('sci') to ideal frame ('idl') and to RA,Dec ('sky')
                # -------------------------------------------------------------------------------------------------------
                #
                for i in range(nsources):
                  #
                  x_pix = x_pix_catfile[i]
                  y_pix = y_pix_catfile[i]
                  #
                  ra, dec       = aper.convert(x_pix, y_pix, 'sci', 'sky')
                  #
                  v2, v3        = aper.convert(x_pix, y_pix, 'sci', 'tel')
                  #
                  x_idl, y_idl  = aper_NRCALL_FULL.convert(v2, v3, 'tel', 'idl')
                  #
                  catalog_current_exposure.add_row([prop_obsid_filter_expnum, rootname+'.fits', x_pix, y_pix, ra, dec, x_idl, y_idl, v2, v3])
              #
              print('')
              print('Writing out:  ',catfile_current_exposure)
              print('')
              catalog_current_exposure.write(catfile_current_exposure, overwrite=True)
          


            # Also store this as the "reference catalog" if this is the firet exposure, and if we're doing "relative" analysis
            #
            if (analysis_type == 'relative'):
              #
              if (first_exposure):
                #
                print(prop_obsid_filter_expnum, '  First exposure  ***')
                print('Copying "catalog_current_exposure" to "catalog_reference"')
                #
                catalog_reference = catalog_current_exposure
            #
            #
            if (analysis_type == 'absolute'):
              #
              catalog_reference = catalog_lmc_radec_xyidl_v2v3_trim



            # Now match sources
            # -----------------
            #
            if ( ((analysis_type == 'relative') and (not first_exposure)) or
                  (analysis_type == 'absolute') ):
              #
              nsources_reference               = len(catalog_reference)
              print('nsources_reference        = ', nsources_reference)
              #
              nsources_current_exposure        = len(catalog_current_exposure)
              print('nsources_current_exposure = ',nsources_current_exposure)
              #
              n = 0
              #
              for i in range(nsources_current_exposure):
                #
                ra_current    = catalog_current_exposure['ra'][i]
                dec_current   = catalog_current_exposure['dec'][i]
                #
                dr_min = 1.e10
                #
                match_j = -1
                #
                for j in range(nsources_reference):
                  #
                  dec_ref = catalog_reference['dec'][j]
                  #
                  ddec = dec_current - dec_ref
                  #
                  if (abs(ddec) < tol):
                    #
                    cosdec = cos(np.deg2rad(dec_current))
                    #
                    ra_ref = catalog_reference['ra'][j]
                    #
                    dra = (ra_current - ra_ref) * cosdec
                    #
                    if (abs(dra) < tol):
                      #
                      dr = ddec**2 + dra**2
                      #
                      if (dr < dr_min):
                        #
                        dr_min   = dr
                        dra_min  = dra
                        ddec_min = ddec
                        match_j  = j
                #
                #
                # Save the best-matching source:
                #
                if (match_j >= 0):
                  #
                  prop_obsid_filter_expnum_current = catalog_current_exposure['prop_obsid_filter_expnum'][i]
                  fitsfilename_current             = catalog_current_exposure['fitsfilename'][i]
                  #
                  x_pix_current = catalog_current_exposure['x_pix'][i]
                  y_pix_current = catalog_current_exposure['y_pix'][i]
                  x_idl_current = catalog_current_exposure['x_idl'][i]
                  y_idl_current = catalog_current_exposure['y_idl'][i]
                  v2_current    = catalog_current_exposure['v2'][i]
                  v3_current    = catalog_current_exposure['v3'][i]
                  #
                  ra_ref    = catalog_reference['ra'][match_j]
                  dec_ref   = catalog_reference['dec'][match_j]
                  x_idl_ref = catalog_reference['x_idl'][match_j]
                  y_idl_ref = catalog_reference['y_idl'][match_j]
                  v2_ref    = catalog_reference['v2'][match_j]
                  v3_ref    = catalog_reference['v3'][match_j]
                  #
                  # Difference between RA,Dec of source in current exposure vs RA,Dec of matched source in reference catalog (LMC)
                  # This was already calculated and stored above, so just copy/ save it here.
                  #
                  dra       = dra_min
                  ddec      = ddec_min
                  #
                  # Measured dither offset of current exposure w.r.t reference exposure which is at dither (0,0)
                  #
                  dx_idl    = x_idl_current - x_idl_ref
                  dy_idl    = y_idl_current - y_idl_ref
                  #
                  # Difference between measured and commanded dither offsets.
                  #
                  dither_dx = dx_idl - dither_xoffset
                  dither_dy = dy_idl - dither_yoffset
                  #
                  dv2       = -(v2_current  - v2_ref)	# Recall that sign of v2 needs to be swapped w.r.t. sign of y_idl
                  dv3       =   v3_current  - v3_ref
                  #
                  catalog_all_matched.add_row([prop_obsid_filter_expnum_current, fitsfilename_current, x_pix_current, y_pix_current, ra_current, dec_current, ra_ref, dec_ref, x_idl_current, y_idl_current, v2_current, v3_current, dra, ddec, dx_idl, dy_idl, dv2, dv3, dither_dx, dither_dy])
              #
              nsources_matched_total = len(catalog_all_matched)
              print('nsources matched (total)  = ',nsources_matched_total)
              #
              hdr0 = fits.getheader(current_prop_obsid_filter_expnum_list[0], 0)
              #
              commanded_xdither = hdr0['XOFFSET']
              commanded_ydither = hdr0['YOFFSET']
              #
              indices_catalog_all_matched_current = np.where(catalog_all_matched['prop_obsid_filter_expnum'] == prop_obsid_filter_expnum)[0].tolist()
              #
              xdither_array = catalog_all_matched['dx_idl'][indices_catalog_all_matched_current]
              ydither_array = catalog_all_matched['dy_idl'][indices_catalog_all_matched_current]
              #
              measured_xdither = np.mean(xdither_array)
              measured_ydither = np.mean(ydither_array)
              #
              difference_xdither_array = xdither_array - commanded_xdither
              difference_ydither_array = ydither_array - commanded_ydither
              #
              difference_xdither = np.mean(difference_xdither_array)
              difference_ydither = np.mean(difference_ydither_array)
              #
              difference_xdither_rms = np.std(difference_xdither_array)
              difference_ydither_rms = np.std(difference_ydither_array)
              #
              dither_offsets_dict[prop_obsid_filter_expnum] = {'commanded_xdither':      commanded_xdither,
                                                               'commanded_ydither':      commanded_ydither,
                                                               'measured_xdither':       measured_xdither,
                                                               'measured_ydither':       measured_ydither,
                                                               'difference_xdither':     difference_xdither,
                                                               'difference_ydither':     difference_ydither,
                                                               'difference_xdither_rms': difference_xdither_rms,
                                                               'difference_ydither_rms': difference_ydither_rms}
        #
        print('')
        print('Writing out:  ',catfile_all_matched)
        print('')
        catalog_all_matched.write(catfile_all_matched, overwrite=True)






      # Make all the plots
      # ==================

      plottype = 'png'
      #plottype = 'pdf'


      # For the dithers - just need to first reshape from [x,y] to separate [x] and [y[
      #
      dither_x = []
      dither_y = []
      #
      for dither_xy in dither_offsets:
        #
        dither_x.append(dither_xy[0])
        dither_y.append(dither_xy[1])



      # Plot all panels on the same page.
      # ---------------------------------
      #
      plotfile = analysis_dir+'plot_'+prop_obsid_filter+'.'+plottype

      mpl.rcParams.update({
                'font.family': 'sans-serif',
                'font.weight': 'normal',
                'figure.titlesize': 15,
                'axes.titlesize': 12,
                'axes.labelsize': 10,
                'xtick.labelsize': 8,
                'ytick.labelsize': 8})


      fig = plt.figure(constrained_layout=True, figsize=(10, 7))
      subfigs = fig.subfigures(2, 1, wspace=0., height_ratios=[1.8, 1])
      fig.tight_layout(pad=1.0)

      plot_title(prop_obsid_filter)

      axsTop = subfigs[0].subplots(1, 2)


      ax = axsTop[0]
      #
      ax.plot(catalog_all_matched['dx_idl'], catalog_all_matched['dy_idl'], marker='o', fillstyle='full', markerfacecolor='blue', markeredgecolor='blue', markersize=0.5, markeredgewidth=0.1, linestyle="None")
      ax.plot(dither_x, dither_y, marker='+', fillstyle='full', markerfacecolor='red', markeredgecolor='red', markersize=10., markeredgewidth=1.0, linestyle='None')
      ax.set_xlim(xmin=xlim[0], xmax=xlim[1])
      ax.set_ylim(ymin=ylim[0], ymax=ylim[1])
      ax.xaxis.set_major_locator(MultipleLocator(xmajor))
      ax.xaxis.set_minor_locator(MultipleLocator(xminor))
      ax.yaxis.set_major_locator(MultipleLocator(ymajor))
      ax.yaxis.set_minor_locator(MultipleLocator(yminor))
      ax.set_xlabel('dx_IDL (")')
      ax.set_ylabel('dy_IDL (")')
      ax.set_title('Dither offsets (Ideal coord. system)')
      ax.set_aspect(1.)


      ax = axsTop[1]
      #
      x = catalog_all_matched['dra']*3600.
      y = catalog_all_matched['ddec']*3600.
      #
      ax.plot(x, y, marker='o', fillstyle='full', markerfacecolor='blue', markeredgecolor='blue', markersize=0.5, markeredgewidth=0.1, linestyle="None")
      ax.plot([0.], [-0.], marker='+', fillstyle='full', markerfacecolor='red', markeredgecolor='red', markersize=10., markeredgewidth=1.0, linestyle='None')
      ax.set_xlim(xmin=-0.02, xmax=0.02)
      ax.set_ylim(ymin=-0.02, ymax=0.02)
      ax.xaxis.set_major_locator(MultipleLocator(0.01))
      ax.xaxis.set_minor_locator(MultipleLocator(0.001))
      ax.yaxis.set_major_locator(MultipleLocator(0.01))
      ax.yaxis.set_minor_locator(MultipleLocator(0.001))
      ax.set_xlabel('dR.A. (")')
      ax.set_ylabel('dDec. (")')
      ax.set_aspect(1.)
      #
      divider = make_axes_locatable(ax)
      ax_histx = divider.append_axes("top",   size=0.8, pad=0.1, sharex=ax)	# size, pad are in inches by default
      ax_histy = divider.append_axes("right", size=0.8, pad=0.1, sharey=ax)	# size, pad are in inches by default
      #
      ax_histx.xaxis.set_tick_params(labelbottom=False)
      ax_histy.yaxis.set_tick_params(labelleft=False)
      #
      binwidth = 0.001
      xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
      lim = (int(xymax/binwidth) + 1)*binwidth
      #
      bins = np.arange(-lim, lim + binwidth, binwidth)
      ax_histx.hist(x, bins=bins, color='blue')
      ax_histy.hist(y, bins=bins, color='blue', orientation='horizontal')
      #
      ax_histx.set_yticks([])
      ax_histy.set_xticks([])
      ax_histx.set_title('dRA, dDec')
      #
      axsBottom = subfigs[1].subplots(1, 1)


      ax = axsBottom
      #
      ax.plot()
      ax.axis('off')
      ax.set_xlim(xmin=0,  xmax=10)
      ax.set_ylim(ymin=10, ymax=0)
      #
      yinc = 0.8
      #
      mpl.rcParams.update({'font.weight': 'bold'})
      y = 0.0;  mpl.pyplot.text(1, y, 'Commanded dithers');       mpl.pyplot.text(3, y, 'Measured dithers');        mpl.pyplot.text(5, y, 'Diff. (meas - comm.)');    mpl.pyplot.text(7, y, 'Measurement r.m.s.')
      y = yinc; mpl.pyplot.text(1, y, 'dx_IDL (")   dy_IDL (")'); mpl.pyplot.text(3, y, 'dx_IDL (")   dy_IDL (")'); mpl.pyplot.text(5, y, 'dx_IDL (")   dy_IDL (")'); mpl.pyplot.text(7, y, 'dx_IDL (")   dy_IDL (")')
      mpl.rcParams.update({'font.weight': 'normal'})
      #
      if (prop_obsid_filter == 'jw01073007_F277W'):
        explist = list(dither_offsets_dict.keys())[:8]	# just list the first 8 exposures, the next set are duplicates.
      else:
        explist = list(dither_offsets_dict.keys())
      #
      y += 0.2
      for prop_obsid_filter_expnum in explist:
        #
        if (prop_obsid_filter_expnum[0:10] in ['jw01073001', 'jw01073002']):
          formatstr = '%8.3f  %10.3f'
        else:
          formatstr = '%8.4f  %10.4f'
        y += yinc
        mpl.pyplot.text(1, y, formatstr % (dither_offsets_dict[prop_obsid_filter_expnum]['commanded_xdither'],      dither_offsets_dict[prop_obsid_filter_expnum]['commanded_ydither']))
        mpl.pyplot.text(3, y, formatstr % (dither_offsets_dict[prop_obsid_filter_expnum]['measured_xdither'],       dither_offsets_dict[prop_obsid_filter_expnum]['measured_ydither']))
        mpl.pyplot.text(5, y, formatstr % (dither_offsets_dict[prop_obsid_filter_expnum]['difference_xdither'],     dither_offsets_dict[prop_obsid_filter_expnum]['difference_ydither']))
        mpl.pyplot.text(7, y, formatstr % (dither_offsets_dict[prop_obsid_filter_expnum]['difference_xdither_rms'], dither_offsets_dict[prop_obsid_filter_expnum]['difference_ydither_rms']))


      plt.show()

      a = os.system('/bin/rm -f '+plotfile)
      if (plottype == 'pdf'):
        plt.savefig(plotfile)
      else:
        plt.savefig(plotfile, dpi=300)



