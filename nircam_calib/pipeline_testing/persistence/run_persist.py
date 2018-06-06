#! /usr/bin/env python

'''
Run the persistence step of the pipeline 

Current persistence saturation (persat), and
trap density (trapdensity) files in CRDS are
nonsense. Need more realistic files in order 
to test. Using the current dummy files, the 
proper outputs are constructed, but they are
all zeros.
'''

files = ['dither1_B5_F250M_linear.fits',
         'dither2_B5_F250M_linear.fits',
         'dither3_B5_F250M_linear.fits']


from jwst.persistence import PersistenceStep

m = PersistenceStep(files[0])
m.save_persistence = True
#m.override_persat = ''
#m.override_trapdensity = ''
m.run()

m2 = PersistenceStep(files[1])
m2.input_trapsfilled = 'dither1_B5_F250M_linear_trapsfilled.fits'
m2.save_persistence = True
m2.run()

m3 = PersistenceStep(files[2])
m3.input_trapsfilled = 'dither2_B5_F250M_linear_trapsfilled.fits'
m3.save_persistence = True
m3.run()

