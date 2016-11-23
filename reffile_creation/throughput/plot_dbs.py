#! /usr/bin/env python

'''plot DBS curves for TR'''

lwafile = 'DBS_LW_ModA_highres.txt'
lwbfile = 'DBS_LW_ModB_highres.txt'
swabadfile = 'DBS_SW_ModA_highres.txt'
swagoodfile = 'DBS_SW_modA_highres_1_minus_trans_plus_absorb.txt'
swbfile = 'DBS_SW_ModB_highres.txt'


import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

lwa = Table.read(lwafile,header_start=0,data_start=1,format='ascii')
lwb = Table.read(lwbfile,header_start=0,data_start=1,format='ascii')

swabad = Table.read(swabadfile,header_start=0,data_start=1,format='ascii')
swagood = Table.read(swagoodfile,header_start=0,data_start=1,format='ascii')
swb = Table.read(swbfile,header_start=0,data_start=1,format='ascii')

f,a = plt.subplots()
a.plot(lwa['microns'],lwa['transmission'],color='red',label='LWA')

good = swabad['microns'] < 2.71
a.plot(swabad['microns'][good],swabad['Reflection'][good],color='orange',label='SWA Bad')
#a.plot(swabad['microns'],swabad['Reflection'],color='orange',linestyle='--',label='SWA Bad')

good2 = swagood['microns'] < 2.71
a.plot(swagood['microns'][good2],swagood['reflection'][good2],color='blue',label='SWA Good')

a.legend(loc='lower right')
a.set_xlim(0.5,5.5)
a.set_ylabel('Transmission')
a.set_xlabel('Wavelength (microns)')
a.set_title('Module A DBS Transmission')
#f.savefig('DBS_ModA_curves.pdf')
plt.close(f)

absb = 1.0 - (swb['reflection'] + np.interp(swb['microns'],lwb['microns'],lwb['transmission']))

#zoomed in view of shoulder area in mod A data, including absorption
f,a = plt.subplots()
a.plot(lwa['microns'],lwa['transmission'],color='red',label='LWA')
good = swabad['microns'] < 2.71
good2 = swagood['microns'] < 2.71
a.plot(swagood['microns'][good2],swagood['reflection'][good2],color='blue',linewidth=1.5,label='SWA Good')
a.plot(swabad['microns'][good],swabad['Reflection'][good],color='orange',label='SWA Bad')
good3 = swb['microns'] < 2.71
a.plot(swb['microns'][good3],swb['reflection'][good3],color='magenta',linestyle='--',label='SWB')
a.plot(swb['microns'],absb,color='green',label='Absorption')
a.legend(loc='center right')
a.set_xlim(2.0,3.0)
a.set_ylabel('Transmission')
a.set_xlabel('Wavelength (microns)')
f.savefig('DBS_ModA_curves_zoom.pdf')
plt.close(f)



f,a = plt.subplots()
a.plot(lwb['microns'],lwb['transmission'],color='red',label='LWB')
a.plot(swb['microns'],swb['reflection'],color='blue',label='SWB')

a.legend(loc='lower right')
a.set_xlim(0.5,5.5)
a.set_ylabel('Transmission')
a.set_xlabel('Wavelength (microns)')
a.set_title('Module B DBS Transmission')
f.savefig('DBS_modB_curves.pdf')
plt.close(f)

#make a separate plot of absorption
#absb = 1.0 - (swb['reflection'] + np.interp(swb['microns'],lwb['microns'],lwb['transmission']))
#f,a = plt.subplots()
#a.plot(swb['microns'],absb,color='red',label='Absorption')
#a.set_xlabel('Wavelength (microns)')
#a.set_ylabel('Absorption')
#a.set_xlim(0.5,5.5)
#f.savefig("DBS_absorption_curve.pdf")
#plt.close(f)

                                        
                                    
