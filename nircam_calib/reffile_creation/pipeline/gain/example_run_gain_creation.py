#! /usr/bin/env python

'''
Test the gain creation scripts from beginning to end
(from a pile of flats to a final CRDS reffile)
'''

from glob import glob
from nircam_calib.reffile_creation.pipeline.gain import gain
from nircam_calib.reffile_creation.pipeline.gain import final_gain_map

files = glob('NRCN815B-LIN*saturation.fits')
files.sort()
print("Creating gain reference file from: {}".format(files))

# Create individual gain files
for i in range(0,len(files),2):
    g = gain.Gainimclass()
    g.flatfile1 = files[i]
    g.flatfile2 = files[i+1]
    g.boxsize = 128
    g.gmax = 10
    g.gainim(files[i],files[i+1])

# Calculate the mean and make a reffile
indfiles = glob('*dsubij.fits')
gref = final_gain_map.Map()
gref.gainlist = indfiles
gref.author = 'Bryan Hilbert'
gref.descrip = 'This is a test reference file'
gref.pedigree = 'GROUND then until now'
gref.useafter = '2014-01-01T00:00:00.0'
gref.save_output = True
gref.outfile = 'NIRCam_TESTTEST_B3_CV3_128x128_gain.fits'
gref.history = ("This is a long history entry "
                "talking about how much blood, "
                "sweat, and tears went into making "
                "this reference file. But what used "
                "to be manual is now automated, so "
                "let's see how it comes out.")
gref.create_map()

