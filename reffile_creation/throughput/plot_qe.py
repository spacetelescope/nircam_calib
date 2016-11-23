#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

sw_wavelength = np.arange(0.5,2.8,0.001)
lw_wavelength = np.arange(2.25,5.5,0.001)

sw_coeffs = np.array([0.65830,-0.05668,0.25580,-0.08350])
sw_exponential = 100.
sw_wavecut = 2.38
lw_coeffs_a = np.array([0.934871,0.051541,-0.281664,0.243867,-0.086009,0.014509,-0.001])
lw_factor_a = 0.88
lw_coeffs_b = np.array([2.9104951,-2.182822,0.7075635,-0.071767])

sw_qe = sw_coeffs[0] + sw_coeffs[1]*sw_wavelength + sw_coeffs[2]*sw_wavelength**2 + sw_coeffs[3]*sw_wavelength**3
red = sw_wavelength > sw_wavecut
sw_qe[red] = sw_qe[red] * np.exp((sw_wavecut-sw_wavelength[red])*sw_exponential)


lw_qe_a = lw_factor_a * (lw_coeffs_a[0] + lw_coeffs_a[1]*lw_wavelength + lw_coeffs_a[2]*lw_wavelength**2 + lw_coeffs_a[3]*lw_wavelength**3 + lw_coeffs_a[4]*lw_wavelength**4 + lw_coeffs_a[5]*lw_wavelength**5 + lw_coeffs_a[6]*lw_wavelength**6)
lw_qe_b = lw_coeffs_b[0] + lw_coeffs_b[1]*lw_wavelength + lw_coeffs_b[2]*lw_wavelength**2 + lw_coeffs_b[3]*lw_wavelength**3

f,a = plt.subplots()
a.plot(sw_wavelength,sw_qe,color='green',label='SW')
a.plot(lw_wavelength,lw_qe_a,color='red',label='LWA')
a.plot(lw_wavelength,lw_qe_b,color='blue',label='LWB')
a.set_xlim(0.5,5.5)
a.set_ylim(0,1)
a.set_ylabel('QE')
a.set_xlabel('Wavelength (microns)')
a.legend(loc='lower left')
f.savefig('QE_curves.pdf')
plt.close(f)
