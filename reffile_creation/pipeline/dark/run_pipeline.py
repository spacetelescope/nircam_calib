#! /usr/bin/env python

from jwst.datamodels import RampModel
from jwst.superbias import SuperBiasStep
from jwst.refpix import RefPixStep
from jwst.linearity import LinearityStep
from jwst.jump import JumpStep



listfile = 'alldarks_thru_saturation.list'

#read in list of dark files to process
files=[]
with open(listfile) as f:
    for line in f:
        if len(line) > 3:
            files.append(line.strip())


#assume they've been run thru bpm and saturation
#so now do superbias subtraction, refpix corr,
#linearity corr, jump step

#hardwire for A1 at the moment:
reffile_dir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/'
sbfile = '/ifs/jwst/wit/witserv/data4/nrc/hilbert/superbias/cv3/A1/A1_superbias_from_list_of_biasfiles.list.fits'
linfile = reffile_dir+'linearity/NRCA1_17004_LinearityCoeff_ADU0_2016-05-14_ssblinearity_DMSorient.fits'
gainfile = reffile_dir+'gain/NRCA1_17004_Gain_ISIMCV3_2016-01-23_ssbgain_DMSorient.fits'
ronfile = '/grp/jwst/wit/nircam/reference_files/SSB/CV2/delivery_Dec_2015/Read_Noise/NRCA1_16989_CDSNoise_2014-10-24_ssbreadnoise_DMSorient.fits'


for file in files:
    sbout = file[0:-5] + '_superbias.fits'
    data = SuperBiasStep.call(file,override_superbias=sbfile,output_file=sbout)
    refout = sbout[0:-5] + '_refpix.fits'
    data = RefPixStep.call(data,output_file=refout)
    linout = refout[0:-5] + '_linearity.fits'
    data = LinearityStep.call(data,override_linearity=linfile,output_file=linout)
    jumpout = linout[0:-5] + '_jump.fits'
    data = JumpStep.call(data,override_gain=gainfile,override_readnoise=ronfile,output_file=jumpout)

