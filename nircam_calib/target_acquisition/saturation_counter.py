#! /usr/bin/env python

'''
Count the number of pixels at hard saturation (65535)
in each group of the given files'''


from glob import glob
from astropy.io import fits
from astropy.table import Table
import numpy as np

unsorted_files = glob('*/nrca5_TA_timeseries_NRCFLATA5GRTS*uncal.fits')
print("Number of uncal files found: {}".format(len(unsorted_files)))

center = []
for file in unsorted_files:
    split = file.split('_')
    c = split.index('center')
    center.append(split[c+1] + '_' + split[c+2])
center = np.array(center)
centers = np.unique(center)
centers = sorted(centers)

files = []
for c in centers:
    good = [f for f in unsorted_files if c in f]
    good = sorted(good)
    files += good

satall = []
sat23 = []
sat3 = []
for file in files:
    with fits.open(file) as h:
        data = h[1].data
    #tot = []
    for grp in range(data.shape[1]):
        if grp == 0:
            satall.append(np.sum(data[0,grp,:,:] == 65535))
        elif grp > 0:
            sat = ((data[0,grp,:,:] == 65535) & (data[0,grp-1,:,:] != 65535))
            if grp == 1:
                sat23.append(np.sum(sat))
            elif grp == 2:
                sat3.append(np.sum(sat))
        #tot.append(sum(sat))
    #print("{}, {}".format(file,tot))

sattab = Table()
sattab['File'] = files
sattab['Sat_in_all_groups'] = satall
sattab['Sat_in_grps_2_3'] = sat23
sattab['Sat_in_grp_3_only'] = sat3
#print(sattab)
sattab.write('saturation_table.tab', format='ascii', overwrite=True)
