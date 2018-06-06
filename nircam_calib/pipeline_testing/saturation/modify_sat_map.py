#! /usr/bin/env python

'''
Modify an existing saturation map reference file in order
to more completely test the pipeline step

Flag a couple pixels with no_sat_check, one with the saturation
value set to 0 (as it usually is in the reference file), and
the other with some non-zero value, just to make sure that 
it is being ignored.
'''

from astropy.io import fits

file = '../reffiles/NRCA1_17004_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits'
outfile = 'nrca1_modified_saturation_reffile.fits' 

#pixel to set saturation value to 0 and flag with no_sat_check
x0 = 100
y0 = 100

#pixel to set non-zero sat value and flag with no_sat_check
x1 = 101
y1 = 101
sat1 = 100.


h = fits.open(file)

#set pixel to no_sat_check and
#set saturation value to zero
h[1].data[y0,x0] = 0.
h[2].data[y0,x0] = 2

#set pixel to no_sat_check and
#set saturation value to non-zero
h[1].data[y1,x1] = sat1
h[2].data[y1,x1] = 2

#save
h.writeto(outfile,overwrite=True)
