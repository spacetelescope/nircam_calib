#! /usr/bin/env python

'''
Create copies of source catalogs for TA data. Assume we have all
the catalogs for a given peak position, and make copies for all 
other peak positions
'''

import os
#from glob import glob
from copy import copy
from astropy.io import ascii
import numpy as np

#instring = '15.5_15.5'
#incats = glob('*_15.5_15.5*cat')
#incats = 'A5_TA_timeseries_ptsrc_mag04.0_center_15.5_15.5.cat'
#xpos = np.arange(1525,1601,25) / 100.
#ypos = np.arange(1525,1601,25) / 100.

xpos = [15.0,15.25,15.5,15.75,16.0]
ypos = [15.0,15.25,15.5,15.75,16.0]
maglist = np.arange(57,101,1)/10   # basic mag coverage

#still need to determine this range below...
#maglist2 = np.arange(65,80,1)/10 # higher res. near area of interest

#maglist = np.append(maglist,maglist2)
maglist = sorted(np.unique(maglist))

inlines = ["#position_pixels",
           "#",
           "#",
           "#",
           "x_or_RA y_or_Dec magnitude",
           "15.5 15.5 4.0"]

prefix = 'A5_TA_timeseries_ptsrc_mag'

for xp in xpos:
    for yp in ypos:
        #posstr = str(xp)+'_'+str(yp)
        for mag in maglist:
            magstr = str(np.around(mag,decimals=1))
            xstr = str(np.around(xp,decimals=2))
            ystr = str(np.around(yp,decimals=2))
            lines = copy(inlines)            
            lines[-1] = xstr + ' ' + ystr + ' ' + magstr
            if mag < 10:
                magpre = '0'
            else:
                magpre = ''
            newfilename = prefix + magpre + magstr + '_center_'+ xstr +'_' + ystr+'.cat'
            magint = str(np.int(mag))
            newfilename = os.path.join('mag'+magint,newfilename)
            with open(newfilename,'w') as h:
                for ll in lines:
                    h.write(ll+'\n')

        
