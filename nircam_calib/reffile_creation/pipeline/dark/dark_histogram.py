#! /usr/bin/env python

'''create a histogram of a dark current reference file'''

file = 'A1_CV3_dark_reffile_mylinefit.fits'


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

with fits.open(file) as h:
    data = h[1].data

ny,nx = data.shape

qstart = [4,512,1024,1536,2040]
colors = ['red','blue','black','plum']

f,a = plt.subplots()
for i in xrange(4):
    amp = data[4:2044,qstart[i]:qstart[i+1]]
    mn = np.mean(amp)
    dv = np.std(amp)
    hmin = mn - 7*dv
    hmax = mn + 7*dv
    nn,bb,pp = plt.hist(np.ravel(amp),bins=(np.arange(hmin,hmax,dv/29)),facecolor=colors[i],alpha=0.25,label='Sector '+str(i+1))

a.set_yscale('log')
a.set_xlim(-0.04,0.07)
a.set_ylabel('Number of Pixels')
a.set_xlabel('Dark Current Rate (DN/sec)')
a.legend(loc='upper right')
f.savefig('A1_CV3_dark_reffile_linefit_histogram_superzoom.pdf')
