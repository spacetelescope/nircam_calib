#! /usr/bin/env python

'''
test readnoise_reffile.py from within python. Command line testing was done 
using test_readnoise_reffile.bash.
'''

import readnoise_reffile as ron

rnfile = ron.ReadnoiseReffile()
rnfile.infile = 'NRCNRCA1-DARK-60011933591_1_481_SE_2016-01-01T20h04m44_uncal_dq_init_bias_drift_no1overfcorr_NN2corrected.fits'
rnfile.outfile = 'test_readnoise_reffile_withinpython.fits'
rnfile.boxsizex = 32
rnfile.boxsizey = 32
rnfile.fill = True
rnfile.bpm = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/bpm/NRCA1_17004_BPM_ISIMCV3_2016-01-21_ssbspmask_DMSorient.fits'
rnfile.Pclipmax = 50
rnfile.Npixmin = 300
rnfile.fill_iterations = 4
outfilename,ron_map_model = rnfile.create()
print(ron_map_model.data.shape)
