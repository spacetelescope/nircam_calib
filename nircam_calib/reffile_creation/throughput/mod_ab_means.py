#!/usr/bin/env python

import mean_across_modules
from itertools import izip

filta = 'filteronly_modA.list'
filtb = 'filteronly_modB.list'
nrca = 'nrc_only_moda.list'
nrcb = 'nrc_only_modb.list'
nrcotea = 'nrc_plus_ote_moda.list'
nrcoteb = 'nrc_plus_ote_modb.list'

def get_names(file):
    names = []
    with open(file) as f:
        for line in f:
            if len(line) > 2:
                names.append(line.strip())
    return names

#create means for filter only tables
afiles = get_names(filta)
bfiles = get_names(filtb)
for a,b in izip(afiles,bfiles):
    print('Working on files {} and {}.'.format(a,b))
    m = mean_across_modules.ModMean()
    m.file1 = a
    m.file2 = b
    under = a.find('_')
    filt = a[0:under]
    m.outfile = filt + '_filteronly_ModAB_mean.txt'
    m.run()

#create means for nircam-only files
afiles = get_names(nrca)
bfiles = get_names(nrcb)
for a,b in izip(afiles,bfiles):
    m = mean_across_modules.ModMean()
    m.file1 = a
    m.file2 = b
    under = a.find('_')
    filt = a[0:under]
    m.outfile = filt + '_NRConly_ModAB_mean.txt'
    m.run()

#create means for nircam+ote files
afiles = get_names(nrcotea)
bfiles = get_names(nrcoteb)
for a,b in izip(afiles,bfiles):
    m = mean_across_modules.ModMean()
    m.file1 = a
    m.file2 = b
    under = a.find('_')
    filt = a[0:under]
    m.outfile = filt + '_NRC_and_OTE_ModAB_mean.txt'
    m.run()



