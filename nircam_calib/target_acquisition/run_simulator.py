#! /usr/bin/env python

'''
Run the simulator and create all of the TA image to use for 
this study
'''
import os
from glob import glob
from nircam_simulator.scripts import imaging_simulator

basestr = 'nrca5_TA_timeseries_NRCFLATA5GRTS'
endstr = '_uncal.fits'
endstr2 = '_linear.fits'

pfiles = glob('*/ta_sim_timeseries_A5_nrcflata5grts*yaml')

for paramfile in pfiles:
    p = imaging_simulator.ImgSim()
    p.paramfile = paramfile
    p.override_dark = 'nrca5_TA_timeseries_NRCFLATA5GRTS_mag04_uncal_linear_dark_prep_object.fits'
    p.create()
    print(paramfile)
    magloc = paramfile.find('_mag')
    noiseloc = paramfile.find('liz')
    midstr = paramfile[magloc:noiseloc+5]
    outfile = basestr + midstr + endstr
    linfile = basestr + midstr + endstr2
    yamin = os.path.split(paramfile)
    dirout = yamin[0]
    os.rename(outfile,os.path.join(dirout,outfile))
    allfiles = os.listdir('./')
    os.remove(linfile)
    for file in allfiles:
        if file.endswith('seed_image.fits'):
            os.remove(file)
        if file.endswith('cosmicrays.list'):
            os.remove(file)
        if file.endswith('pointsources.list'):
            os.remove(file)

