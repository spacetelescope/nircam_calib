#! /usr/bin/env python

"""Create imaging TSO catalog, modify point source catalog, create lightcurve file, etc
"""
import os
import glob
import shutil
import yaml
import numpy as np
import matplotlib.pyplot as plt
import batman
from mirage.catalogs.hdf5_catalog import save_tso
from mirage.catalogs.catalog_generator import ImagingTSOCatalog

output_dir = '/ifs/jwst/wit/nircam/simulationWG/Imaging/CAR-19/Input_Catalogs/'

# CREATE LIGHTCURVE------------------------------------------------

# In this case, we don't want the lightcurve to change at all during the observation
params = batman.TransitParams()       # object to store transit parameters
params.t0 = 2500.                        # time of inferior conjunction
params.per = 31620.24                       # orbital period
params.rp = 0.723                       # planet radius (in units of stellar radii)
params.a = 9.37                        # semi-major axis (in units of stellar radii)
params.inc = 83.3                      # orbital inclination (in degrees)
params.ecc = 0.                       # eccentricity
params.w = 90.                        # longitude of periastron (in degrees)
params.limb_dark = "nonlinear"        # limb darkening model
params.u = [0.5, 0.1, 0.1, -0.1]      # limb darkening coefficients [u1, u2, u3, u4]

times = np.linspace(0, 780, 1000)  # times at which to calculate light curve

# Create lightcurve file for the TSO object with batman. In this case, lightcurve is flat
m = batman.TransitModel(params, times)
flux = m.light_curve(params)

f, a = plt.subplots()
a.scatter(times, flux, color='red', marker='v')
a.set_xlabel('Time (sec)')
a.set_ylabel('Normalized Signal')
plt.show()

lightcurve_file = os.path.join(output_dir, 'example_lightcurve.hdf5')

contents = {}
contents['1'] = {'times': times,
                 'fluxes': flux}

save_tso(contents, lightcurve_file, time_unit='second')
#-------------------------------------------------------------


# CREATE MIRAGE CATALOG---------------------------------------------------
imaging_tso_catalog = os.path.join(output_dir, 'tso_imaging_source.cat')

object_ra = 80.48259735107422
object_dec = -69.49395751953125
tsimg_cat = ImagingTSOCatalog(ra=[object_ra], dec=[object_dec], lightcurve_file=[lightcurve_file])

tsimg_cat.add_magnitude_column([13.910049438476562], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f200w')
tsimg_cat.add_magnitude_column([13.898650169372559], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f356w')
tsimg_cat.add_magnitude_column([13.903349876403809], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f164n')
tsimg_cat.add_magnitude_column([13.893850326538086], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f480m')
tsimg_cat.add_magnitude_column([13.89954948425293], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f322w2')
tsimg_cat.add_magnitude_column([13.913049697875977], magnitude_system='vegamag',
                               instrument='nircam', filter_name='wlp4')
tsimg_cat.save(imaging_tso_catalog)


# Insert tso catalog, and modified point source catalog (car19_ptsrc_minus_tso.cat),
# into the yaml files for the extended source subarrays as well as the stripe subarray
files = sorted(glob.glob('yaml_files/jw0106800[567]*yaml'))

for filename in files:
    # Read in
    with open(filename,'r') as f:
        data = yaml.safe_load(f)

    # Update catalog names
    data['simSignals']['pointsource'] = os.path.join(output_dir, 'car19_ptsrc_minus_tso.cat')
    data['simSignals']['tso_imaging_catalog'] = imaging_tso_catalog

    # Move the original yaml to a new name
    d, f = os.path.split(filename)
    yaml_copy = os.path.join(d, 'orig_{}'.format(f))
    shutil.copy2(filename, yaml_copy)
    os.remove(filename)

    # Save
    with open(filename, 'w') as output:
        yaml.dump(data, output, default_flow_style=False)


