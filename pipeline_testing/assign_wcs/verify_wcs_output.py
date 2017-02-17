#! /usr/bin/env python

'''
Verify that the assign_wcs step correctly copied the appropriate
header keywords from the reference file into the data file.
'''

filenamefull = 'NRCNRCB1-DARK-60090405201_1_486_SE_2016-01-09T05h33m56_uncal_header_update_sliced_slp_assign_wcs.fits' #full frame
filenamesub = 'NRCV82600012001P0000000002102_1_486_SE_2016-01-18T04h25m10_uncal_sliced_slp_assign_wcs.fits'  #subarray
#rawsubfile = 'NRCV82600012001P0000000002102_1_486_SE_2016-01-18T04h25m10.fits'

from jwst import datamodels as models
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#read in full frame and subarray
imagefull = models.open(filenamefull)
imagesub = models.open(filenamesub)

yd,xd = imagesub.data.shape

#colcorner and rowcorner which indicate the location of the subarray on the detector
#are in the header, but not the model
#so we need to use astropy to pull them out.
with fits.open(filenamesub) as h:
    colcorner = h[0].header['COLCORNR'] - 1
    rowcorner = h[0].header['ROWCORNR'] - 1

#colcorner = imagesub.meta.subarray.colcorner - 1
#rowcorner = imagesub.meta.subarray.rowcorner - 1
print("subarray size and colcorner, rowcorner are {}x{} and {},{}".format(xd,yd,colcorner,rowcorner))

print('Full frame: CRPIX1,2, CRVAL1,2, CDELT1,2')
print(imagefull.meta.wcsinfo.crpix1,imagefull.meta.wcsinfo.crpix2,imagefull.meta.wcsinfo.crval2,imagefull.meta.wcsinfo.crval2,imagefull.meta.wcsinfo.cdelt1,imagefull.meta.wcsinfo.cdelt2)

print('Subarray: CRPIX1,2, CRVAL1,2, CDELT1,2')
print(imagesub.meta.wcsinfo.crpix1,imagesub.meta.wcsinfo.crpix2,imagesub.meta.wcsinfo.crval2,imagesub.meta.wcsinfo.crval2,imagesub.meta.wcsinfo.cdelt1,imagesub.meta.wcsinfo.cdelt2)

x1full = colcorner
y1full = rowcorner
x2full = x1full+40
y2full = y1full+40

x1sub = x1full - colcorner
y1sub = y1full - rowcorner
x2sub = x2full - colcorner
y2sub = y2full - rowcorner

#RA,Dec of reference pixel in the full frame image
x1fullrefx = imagefull.meta.wcsinfo.crpix1
y1fullrefy = imagefull.meta.wcsinfo.crpix2

fullra,fulldec = imagefull.meta.wcs(x1full,y1full)
subra,subdec = imagesub.meta.wcs(x1sub,y1sub)
fullrarefpix,fulldecrefpix = imagefull.meta.wcs(x1fullrefx,y1fullrefy)

print("RA,DEC in full frame pixel {},{} which is rowcorner,colcorner: {},{}".format(x1full,y1full,fullra,fulldec))
print("RA,DEC in same pixel {},{} in subarray: {},{}".format(x1sub,y1sub,subra,subdec))
print("RA,DEC difference between full frame and subarray: {},{}".format(fullra-subra,fulldec-subdec))
print("RA,DEC at reference pixel {},{} in full frame image {},{}".format(x1fullrefx,y1fullrefy,fullrarefpix,fulldecrefpix))

#plot the three pixels according to RA,Dec
f,a = plt.subplots()
a.plot([fullra],[fulldec],marker='8',color='red',label='FF rowc,colc')
a.plot([subra],[subdec],marker='8',color='blue',label='Sub 0,0')
a.plot([fullrarefpix],[fulldecrefpix],marker='8',color='green',label='FF refpix')
#a.legend(loc='best')
#f.savefig('radec_compare.pdf')
#plt.show()

print('RA and Dec of two pixels that are separated by 40 pixels:')
print(imagefull.meta.wcs(80,80))
print(imagefull.meta.wcs(80,120))
diffra = imagefull.meta.wcs(80,80)[0] - imagefull.meta.wcs(80,120)[0]
diffdec = imagefull.meta.wcs(80,80)[1] - imagefull.meta.wcs(80,120)[1]
print("difference: {},{}".format(diffra,diffdec))

#look at the mean distance (in arcsec) between pixels, and see if 
#it agrees with the values given in the CDELT keywords

#get the RA and Dec of two pixels that are 40 pixels apart
(subra1,subdec1) = imagesub.meta.wcs(80,80)
(subra2,subdec2) = imagesub.meta.wcs(80,120)
(fullra1,fulldec1) = imagefull.meta.wcs(80,1967)
(fullra2,fulldec2) = imagefull.meta.wcs(80,2007)
print(subra1-fullra1)
print(subdec1-fulldec1)
print(subra1-subra2)
print(subdec1-subdec2)

ydim,xdim = imagefull.data.shape
y,x = np.mgrid[:ydim,:xdim]
ra,dec = imagefull.meta.wcs(x,y)

print(ra[450:550:30])
print(dec[450:550:30])
