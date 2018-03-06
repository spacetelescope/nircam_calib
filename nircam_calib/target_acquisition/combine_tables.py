#! /usr/bin/env python

'''
Combine the centroid results table and the saturation count tables
'''

import numpy as np
from astropy.table import Table,join
import magnitude_conversions as conv

ctabfile = 'centroid_table_full_dataset.tab'
stabfile = 'saturation_table.tab'

ctab = Table.read(ctabfile, format='ascii')
stab = Table.read(stabfile, format='ascii')
fulltab = join(ctab,stab)

# Now convert the magnitudes, which are in units of F335M
# ABMAGS, to K-band VEGAMAGS
fmags = fulltab['Source_Magnitude_in_Simulator'].data
converter = conv.MagConvert()
kmags = [converter.f335abmag_to_kvega(f).value for f in fmags]
fulltab['K_G2V_Mag'] = kmags

outtab = 'Final_results_table.tab'
fulltab.write(outtab, format='ascii', overwrite=True)
