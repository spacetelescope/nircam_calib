#!/usr/bin/env python

'''Create the DBS curves to use in the NIRCam throughput tables.
This means extracting the high-resolution tables for mod A and
mod B for SW and LW. The high resolution SW curve that does not
have the anomalous knee in transmission between 2.3-2.6um is for 
mod B. This is safe to use. For mod A, use the high-resolution 
curve, but in the 2.3-2.6um range, replace the high-res curve with
values interpolated from the low-resolution copy. '''

#interpolation region:
#2.418um onward
#no - now 2.30 onward.
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from scipy.interpolate import interp1d
from astropy.table import Table

#dbsdir = '/grp/jwst/wit/nircam/reference_files/SpectralResponse_2015-Final/DBS_QE_Optics/'
dbsdir = '/grp/jwst/wit/nircam/hilbert/filter_throughputs/nircam_system_throughputs/'
swa_hres_file = 'DBS_SW_ModA_highres.txt'
swa_lres_file = 'DBS_SW_ModA_LowResSmoothed.txt'
swb_hres_file = 'DBS_SW_ModB_highres.txt'
lwb_hres_file = 'DBS_LW_ModB_highres.txt'
lwa_hres_file = 'DBS_LW_ModA_highres.txt'

swa_hres = ascii.read(dbsdir+swa_hres_file,header_start=0,data_start=1)
swa_hres_w = swa_hres.columns[0].data
swa_hres_t = swa_hres.columns[1].data

#lres = ascii.read(dbsdir+swa_lres_file,header_start=0,data_start=1)
#lres_w = lres.columns[0].data
#lres_t = lres.columns[1].data

swb_hres = ascii.read(dbsdir+swb_hres_file,header_start=0,data_start=1)
swb_hres_w = swb_hres.columns[0].data
swb_hres_t = swb_hres.columns[1].data

lwa_hres = ascii.read(dbsdir+lwa_hres_file,header_start=0,data_start=1)
lwa_hres_w = lwa_hres.columns[0].data
lwa_hres_t = lwa_hres.columns[1].data

lwb_hres = ascii.read(dbsdir+lwb_hres_file,header_start=0,data_start=1)
lwb_hres_w = lwb_hres.columns[0].data
lwb_hres_t = lwb_hres.columns[1].data

#good = hres_w >= 2.418
good = lwa_hres_w >= 2.275
short = swa_hres_w < 2.275
#interp_w = hres_w[good]

sw_to_lw_b = np.interp(lwb_hres_w,swb_hres_w,swb_hres_t)
absorption = 1.0 - sw_to_lw_b - lwb_hres_t
abs_interp = np.interp(lwa_hres_w,lwb_hres_w,absorption)
solution = 1 - lwa_hres_t - abs_interp

new_swa_w = np.concatenate((swa_hres_w[short],lwa_hres_w[good]))
new_swa = np.concatenate((swa_hres_t[short],solution[good]))

#interp_function = interp1d(lres_w,lres_t,kind='linear')
#interp_t = interp_function(interp_w)

#final_t = np.concatenate((hres_t[short],interp_t))

oldgood = swa_hres_w < 2.67

f,a = plt.subplots()
a.plot(lwa_hres_w,lwa_hres_t,color='red',label='LWA',linewidth=1)
#a.plot(swb_hres_w,swb_hres_t,color='black',label='SWB',linewidth=1)
#a.plot(lwb_hres_w,lwb_hres_t,color='black',label='LWB',linewidth=1)
a.plot(new_swa_w[oldgood],new_swa[oldgood],color='blue',label='SWA Good',linewidth=2)
a.plot(swa_hres_w,swa_hres_t,color='orange',label='SWA Bad',linewidth=1)
a.plot(lwa_hres_w,abs_interp,color='green',label='Absorption')
#a.scatter(lres_w,lres_t,color='blue',marker='o',label='Low Res')
#a.plot(interp_w,interp_t,color='green',linewidth=2,linestyle='--',label='Interpolated')
#a.plot(hres_w,final_t,color='red',label='Final Curve')
a.set_xlim(2,3)
#a.set_xlim(.5,5.5)
a.set_ylim(-0.05,1.0)
a.set_xlabel('Wavelength (microns)')
a.set_ylabel("DBS Throughput")
a.legend(loc='center right')
f.savefig('new_DBS_SWA_highres_1_minus_trans_plus_absorb_zoom.png')
f.savefig('new_DBS_SWA_highres_1_minus_trans_plus_absorb_zoom.pdf')
plt.close(f)

#save the new Mod A high-res DBS curve
#outtab = Table([new_swa_w,new_swa],names=('microns','reflection'))
#ascii.write(outtab,'DBS_SW_modA_highres_1_minus_trans_plus_absorb.txt')


