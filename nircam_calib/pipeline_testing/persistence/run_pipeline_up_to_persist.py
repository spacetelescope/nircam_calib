#! /usr/bin/env python

'''
Run calwebb_detector1 up to the persist step,
to prepare a file for persistence testing
'''
from glob import glob
from astropy.io import fits
from jwst.pipeline import calwebb_detector1
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.refpix import RefPixStep
from jwst.linearity import LinearityStep


files = glob('V*uncal.fits')


sbdir = '/ifs/jwst/wit/witserv/data4/nrc/hilbert/superbias/cv3/'
lindir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/linearity/'
satdir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/welldepth/'
gaindir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/gain/'
rondir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/readnoise/'
wcsdir = '/ifs/jwst/wit/witserv/data4/nrc/hilbert/distortion_reference_file/jwreftools/nircam/'


reffiles = {}
reffiles['a1']={'superbias':'A1/A1_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCA1_17004_LinearityCoeff_ADU0_2016-05-14_ssblinearity_v2_DMSorient.fits','saturation':'NRCA1_17004_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCA1_17004_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits','readnoise':'NRCA1_17004_CDSNoise_ISIMCV3_ADU_2016-06-24_ssbreadnoise_DMSorient.fits','wcs':'NRCA1_FULL_distortion.asdf'}
reffiles['a2'] = {'superbias':'A2/A2_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCA2_17006_LinearityCoeff_ADU0_2016-05-14_ssblinearity_v2_DMSorient.fits','saturation':'NRCA2_17006_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCA2_17006_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits','readnoise':'NRCA2_17006_CDSNoise_ISIMCV3_ADU_2016-06-24_ssbreadnoise_DMSorient.fits','wcs':'NRCA2_FULL_distortion.asdf'}
reffiles['a3'] = {'superbias':'A3/A3_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCA3_17012_LinearityCoeff_ADU0_2016-05-14_ssblinearity_v2_DMSorient.fits','saturation':'NRCA3_17012_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCA3_17012_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits','readnoise':'NRCA3_17012_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCA3_FULL_distortion.asdf'}
reffiles['a4'] = {'superbias':'A4/A4_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCA4_17048_LinearityCoeff_ADU0_2016-05-15_ssblinearity_v2_DMSorient.fits','saturation':'NRCA4_17048_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCA4_17048_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits','readnoise':'NRCA4_17048_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCA4_FULL_distortion.asdf'}
reffiles['a5'] = {'superbias':'ALONG/A5_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCALONG_17158_LinearityCoeff_ADU0_2016-05-16_ssblinearity_v2_DMSorient.fits','saturation':'NRCA5_17158_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCA5_17158_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits','readnoise':'NRCA5_17158_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCA5_FULL_distortion.asdf'}
reffiles['b1'] = {'superbias':'B1/B1_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCB1_16991_LinearityCoeff_ADU0_2016-05-17_ssblinearity_v2_DMSorient.fits','saturation':'NRCB1_16991_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCB1_16991_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits','readnoise':'NRCB1_16991_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCB1_FULL_distortion.asdf'}
reffiles['b2'] = {'superbias':'B2/B2_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCB2_17005_LinearityCoeff_ADU0_2016-05-18_ssblinearity_v2_DMSorient.fits','saturation':'NRCB2_17005_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCB2_17005_Gain_ISIMCV3_2016-02-25_ssbgain_DMSorient.fits','readnoise':'NRCB2_17005_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCB2_FULL_distortion.asdf'}
reffiles['b3'] = {'superbias':'B3/B3_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCB3_17011_LinearityCoeff_ADU0_2016-05-20_ssblinearity_v2_DMSorient.fits','saturation':'NRCB3_17011_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCB3_17011_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits','readnoise':'NRCB3_17011_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCB3_FULL_distortion.asdf'}
reffiles['b4'] = {'superbias':'B4/B4_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCB4_17047_LinearityCoeff_ADU0_2016-05-20_ssblinearity_v2_DMSorient.fits','saturation':'NRCB4_17047_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCB4_17047_Gain_ISIMCV3_2016-02-25_ssbgain_DMSorient.fits','readnoise':'NRCB4_17047_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCB4_FULL_distortion.asdf'}
reffiles['b5'] = {'superbias':'BLONG/B5_superbias_from_list_of_biasfiles.list.fits','linearity':'NRCBLONG_17161_LinearityCoeff_ADU0_2016-05-22_ssblinearity_v2_DMSorient.fits','saturation':'NRCB5_17161_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits','gain':'NRCB5_17161_Gain_ISIMCV3_2016-02-25_ssbgain_DMSorient.fits','readnoise':'NRCB5_17161_CDSNoise_ISIMCV3_ADU_2016-06-25_ssbreadnoise_DMSorient.fits','wcs':'NRCB5_FULL_distortion.asdf'}

for file in files:
    loc = file.rfind('uncal')
    outfile = file.replace('uncal','lin')
    det = fits.getval(file,'DETECTOR')
    det = det.lower()[3:]
    if 'long' in det:
        det = det[0] + '5'
    refdict = reffiles[det]

    #m = calwebb_detector1.Detector1Pipeline(config_file='calwebb_detector1.cfg')
    #m.saturation.override_saturation = satdir+refdict['saturation']
    #m.superbias.override_superbias = sbdir+refdict['superbias']
    #m.refpix.odd_even_rows = False
    #m.group_scale.skip = True
    #m.ipc.skip = True
    #m.rscd.skip = True
    #m.lastframe.skip = True
    #m.dark_current.skip = True
    #m.persistence.skip = True
    #m.jump.skip = True
    #m.ramp_fit.skip = False #bug in pipeline means this must
    #be run. 
    #m.linearity.override_linearity = lindir+refdict['linearity']
    #m.output_file = outfile
    #m.run(file)

    m = DQInitStep.call(file,config_file = 'dq_init.cfg')
    m = SaturationStep.call(m,config_file = 'saturation.cfg')
    m = SuperBiasStep.call(m,config_file = 'superbias.cfg')
    m = RefPixStep.call(m,config_file='refpix.cfg')
    m = LinearityStep.call(m,config_file='linearity.cfg',output_file = outfile)
    


