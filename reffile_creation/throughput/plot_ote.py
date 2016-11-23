#! /usr/bin/env/ python

'''plot the OTE...for the TR'''

from astropy.table import Table
import matplotlib.pyplot as plt

otefile = '/itar/jwst/tel/share/Mirror_Reflectivity/jwst_telescope_ote_thruput.fits'
ote = Table.read(otefile)
ote_wave = ote['wavelength'].data
ote_thru = ote['throughput'].data

f,a = plt.subplots()
a.plot(ote_wave,ote_thru,color='blue')
a.set_xlim(0.5,5.5)
a.set_xlabel('Wavelength (microns)')
a.set_ylabel('OTE Throughput')
f.savefig('OTE_throughput_plot.pdf')
plt.close(f)
