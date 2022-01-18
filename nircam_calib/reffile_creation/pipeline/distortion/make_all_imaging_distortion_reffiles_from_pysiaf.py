"""This script creates distortion reference files for NIRCam using the coefficients from
pysiaf.
"""
import nircam_distortion_reffiles_from_pysiaf as ref
import numpy as np
import os

# To use a version of the SIAF other than the default in pysiaf, enter the
# name of a SIAF xml file here.
#siaf_xml_file = None
siaf_xml_file = '/Users/hilbert/Documents/PRD/PRDOPSSOC-041/SIAFXML/SIAFXML/NIRCam_SIAF.xml'

nrc_a_apertures = ['NRCA{}_FULL'.format(i+1) for i in range(5)]
nrc_b_apertures = ['NRCB{}_FULL'.format(i+1) for i in range(5)]
nrc_apertures = nrc_a_apertures + nrc_b_apertures

hist = ("This reference file was created from the distortion coefficients contained in pysiaf "
        "(and therefore the PRD) as of 24 Oct 2019. This includes the fix for the round trip (V2,V3 -> x,y -> V2, V3) error "
        "present in previous versions of the coefficients. This update is described in JIRA issues: "
        "JWSTSIAF-161, JWSTSIAF-123. pysiaf version 0.6.1 was used to access the appropriate "
        "distortion coefficients. The translation model converts from units of pixels on the detector to "
        "V2,V3 in units of arcseconds, as well as the inverse.")

#basedir = 'reference_files/distortion/24Oct2019_round_trip_error_fixed'

hist = ("This reference file was created from the distortion coefficients contained in pysiaf "
        "version 0.13.0, which uses "
        "version 39 of the PRD. This version of the PRD makes use of the new <detector>_FULL_WEDGE_RND"
        "and <detector>_FULL_WEDGE_BAR parent apertures for coronagraphic observations. All "
        "distortion coefficients are calculated from ground based data.")


basedir = '/grp/jwst/wit/nircam/reference_files/distortion/coron_distortion_files_2021_Sept'

# IMAGING metadata-----------------------------------------------
sw_imaging_pupil = ['CLEAR', 'F162M', 'F164N', 'GDHS0', 'GDHS60', 'WLM8', 'WLP8', 'PINHOLES', 'MASKIPR', 'FLAT']
lw_imaging_pupil = ['CLEAR', 'F323N', 'F405N', 'F466N', 'F470N', 'PINHOLES', 'MASKIPR', 'GRISMR', 'GRISMC', 'FLAT']

# Put grism/wfss exp_types in SW list of exposure types because simultaneous SW data were
# tagged as e.g. NRC_GRISM when grism observations were made with the LW detector
# early on.
sw_exptype = ['NRC_IMAGE', 'NRC_TSIMAGE', 'NRC_FLAT', 'NRC_LED',
              'NRC_WFSC', 'NRC_TACQ', 'NRC_TACONFIRM', 'NRC_FOCUS',
              'NRC_DARK', 'NRC_WFSS', 'NRC_TSGRISM', 'NRC_GRISM']
lw_exptype = sw_exptype #+ ['NRC_WFSS', 'NRC_TSGRISM', 'NRC_GRISM']

for aperture in nrc_apertures:
    detector, apname = aperture.split('_')
    outname = os.path.join(basedir, '{}_distortion.asdf'.format(aperture))

    if int(detector[-1]) < 5:
        pupil = sw_imaging_pupil
        #subarr = sw_imaging_subarr
        # Leaving subarrays as any means that new subarrays can be added later and the file won't need to be updated
        subarr = ['GENERIC']
        exp_type = sw_exptype
    else:
        pupil = lw_imaging_pupil
        #subarr = lw_imaging_subarr
        subarr = ['GENERIC']
        exp_type = lw_exptype

    #ref.create_nircam_distortion(detector, apname, outname, pupil, subarr, exp_type, hist, siaf_xml_file=siaf_xml_file)

# CORONAGRAPHY metadata------------------------------------------
# For coronagraphy, we expect the bar and round wedges to have different astrometric solutions,
# in the future if not now, so treat these separately. Also, there is currently no aperture in SIAF for
# A1 and A3 that include the wedge, and so to make reference files for them, I would need to manually
# combine the NRCA1_FULL coefficients with the wedge offsets in the SIAF file. John is putting in a
# ticket to add these apertures to SIAF. It might be worth waiting for that to happen before creating
# these reference files?


# These apertures use a distortion model that is a shift + the imaging mode coefficients. A good approximation,
# but not exact. These are intended as a first attempt at getting this distortion model correct. We now have a
# real distortion model that includes the wedge for the A module, so the e.g. NRCA1_FULL_WEDGE_RND apertures
# should be used instead.
#a2_round_aperture = ['NRCA2_FULL_MASK210R']
#a4_bar_aperture = ['NRCA4_FULL_MASKSWB']
#a5_round_aperture = ['NRCA5_FULL_MASK335R']
#a5_bar_aperture = ['NRCA5_FULL_MASKLWB']

# These apertures use a distortion model that is a shift + the imaging mode coefficients. A good approximation,
# but not exact. These are intended as a first attempt at getting this distortion model correct. Later we will
# hopefully have e.g. NRCB1_FULL_WEDGE_RND apertures that contain the real distortion model, just as we currently
# have for the A module.
#b1_round_aperture = ['NRCB1_MASK210R']
#b3_bar_aperture = ['NRCB3_MASKSWB']
#b5_round_aperture = ['NRCB5_MASK335R']
#b5_bar_aperture = ['NRCB5_MASKLWB']
#
#nrc_coron_apertures = a2_round_aperture + a4_bar_aperture + a5_round_aperture + a5_bar_aperture + b1_round_aperture + b3_bar_aperture + b5_round_aperture + b5_bar_aperture


hist = ("This reference file was created using the distortion coefficients in PRD version PRDOPSSOC-041. This "
        "version of the PRD contains new FULL_WEDGE_RND and FULL_WEDGE_BAR apertures for all A module detectors. "
        "The distortion model for these new apertures is a complete model. Previously, the distortion model for "
        "these apertures, which was only available for the A2, A4, and A5 detectors, was an approximation created "
        "by shifting the imaging mode coefficients to account for the effects of the wedge. With this delivery, "
        "all A module detectors will contain the complete distortion model. See CRDS-506 for more information.")

basedir = 'reference_files/distortion/coron_distortion_files_amod_wedge_apertures_2022_Jan'

nrc_coron_apertures = ['NRCA1_FULL_WEDGE_RND', 'NRCA2_FULL_WEDGE_RND', 'NRCA3_FULL_WEDGE_RND',
                       'NRCA4_FULL_WEDGE_RND', 'NRCA5_FULL_WEDGE_RND', 'NRCA1_FULL_WEDGE_BAR',
                       'NRCA2_FULL_WEDGE_BAR', 'NRCA3_FULL_WEDGE_BAR', 'NRCA4_FULL_WEDGE_BAR',
                       'NRCA5_FULL_WEDGE_BAR']

coron_exptype = ['NRC_CORON', 'NRC_TACQ', 'NRC_TACONFIRM']

for aperture in nrc_coron_apertures:
    parts = aperture.split('_')
    detector = parts[0]
    if len(parts) == 2:
        apname = parts[1]
    elif len(parts) == 3:
        apname = '{}_{}'.format(parts[1], parts[2])
    elif len(parts) == 4:
        apname = '{}_{}_{}'.format(parts[1], parts[2], parts[3])

    if 'RND' in aperture:
        pupil = ['MASKRND']
    elif 'BAR' in aperture:
        pupil = ['MASKBAR']
    else:
        raise ValueError('Bad pupil value: {}'.format(aperture))

    outname = os.path.join(basedir, '{}_{}_distortion.asdf'.format(detector, pupil[0]))

    # Leaving subarrays as any means that new subarrays can be added later and the file won't need to be updated
    subarr = ['GENERIC']
    exp_type = coron_exptype

    ref.create_nircam_distortion(detector, apname, outname, pupil, subarr, exp_type, hist, siaf_xml_file=siaf_xml_file)
