#! /usr/bin/env/ python

'''plot the components of the optics file...for the TR'''

from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

def replace_nan(items):
    for index, item in enumerate(items):
        if (item == '---'):
            items[index] = float('nan')
    return items


ofile = 'NIRCam_optics_transmission_29Oct2015.csv'

opttab = ascii.read(ofile,header_start=1,data_start=2,format='csv')


wave = opttab['Wavelength'].data.data
nvr_thru = opttab['NVR_Transmission'].data.data
nvr_wave = opttab['NVR_Wavelength'].data.data
collimator = opttab['Collimator'].data.data
sw_triplet = replace_nan(opttab['SW_Triplet'].data.data).astype('float')
sw_mirrors = replace_nan(opttab['SW_Mirrors'].data.data).astype('float')
lw_triplet = replace_nan(opttab['LW_triplet'].data.data).astype('float')
lw_mirrors = replace_nan(opttab['LW_Mirrors'].data.data).astype('float')
sw_particulates = replace_nan(opttab['SW_Particulates'].data.data).astype('float')
lw_particulates = replace_nan(opttab['LW_Particulates'].data.data).astype('float')

#remove extra entries in NVR columns
good = np.where(nvr_wave != 0.)[0]
nvr_thru = nvr_thru[good]
nvr_wave = nvr_wave[good]


#interpolate NVR to the same wavelength scale as the other columns
nvr_interp = np.interp(wave,nvr_wave,nvr_thru)

#combine the elements to produce a SW optics curve and a LW optics curve
#The 0.98 factor is a 'contingency factor' John Stansberry included in a 
#previous version. He said we can keep it out for this version
sw_optics = collimator * sw_triplet * sw_mirrors * sw_particulates #* 0.98
lw_optics = collimator * lw_triplet * lw_mirrors * lw_particulates #* 0.98


f,a = plt.subplots()
a.plot(wave,collimator,color='red',label='Collimator')

a.plot(wave,sw_triplet,color='blue',label='SW Triplet')
a.plot(wave,sw_mirrors,color='black',label='SW Mirrors')
a.plot(wave,sw_particulates,color='green',label='SW Particulates')

a.plot(wave,lw_triplet,color='blue',linestyle='--',label='LW Triplet')
a.plot(wave,lw_mirrors,color='black',linestyle='--',label='LW Mirrors')
a.plot(wave,lw_particulates,color='green',linestyle='--',label='LW Particulates')

a.plot(nvr_wave,nvr_thru,color='orange',label='NVR')

a.plot(wave,sw_optics,color='magenta',label='Total SW')
a.plot(wave,lw_optics,color='magenta',linestyle='--',label='Total LW')

a.set_xlim(0.5,5.5)
a.set_ylabel('Throughput')
a.set_xlabel('Wavelength (microns)')
a.legend(loc='lower right')

#f.savefig('Optics_components_plot.pdf')
plt.close(f)


#2 panel version
#f,(a,b) = plt.subplots(2,sharex=True)
#f = plt.figure()
#gs = gridspec.GridSpec(1,2,height_ratios=[3,1])
#a = plt.subplot(gs[0])
#b = plt.subplot(gs[1])
f, (a, b) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[1,3]},sharex=True)
#dummy plots so that we can get the curves from the upper panel in the
#legend of the lower panel
b.plot(wave,sw_mirrors+10,color='black',label='SW Mirrors')
b.plot(wave,lw_mirrors+10,color='black',linestyle='--',label='LW Mirrors')
b.plot(wave,sw_particulates+10,color='lightgreen',label='SW Particulates')
b.plot(wave,lw_particulates+10,color='lightgreen',linestyle='--',label='LW Particulates')


b.plot(nvr_wave,nvr_thru,color='orange',label='NVR')
b.plot(wave,collimator,color='orchid',linestyle='-.',linewidth=2,label='Collimator')

b.plot(wave,sw_triplet,color='darkcyan',label='SW Triplet')
a.plot(wave,sw_mirrors,color='black',label='SW Mirrors')
a.plot(wave,sw_particulates,color='lightgreen',label='SW Particulates')

b.plot(wave,lw_triplet,color='darkcyan',linestyle='--',label='LW Triplet')
a.plot(wave,lw_mirrors,color='black',linestyle='--',label='LW Mirrors')
a.plot(wave,lw_particulates,color='lightgreen',linestyle='--',label='LW Particulates')

b.plot(wave,sw_optics,color='red',label='Total SW')
b.plot(wave,lw_optics,color='red',linestyle='--',label='Total LW')



a.set_xlim(0.5,5.25)
b.set_ylim(0,1)
a.set_ylim(0.9,1.0)
a.set_ylabel('Throughput')
b.set_ylabel('Throughput')
b.set_xlabel('Wavelength (microns)')
b.legend(loc='lower right',fontsize=11)

f.savefig('Optics_components_plot_2panel.pdf')
plt.close(f)
