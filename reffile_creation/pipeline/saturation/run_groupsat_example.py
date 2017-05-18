#! /usr/bin/env python

'''
Example runs of group_saturation.py, both with and
without the inverse linearity reference file.
'''

import group_saturation as g

#First the faster case, where a reference file of
#coefficients is provided to convert linear to
#non-linear data
coeffcase = g.GroupSat()
coeffcase.satfile = 'my_satfile.fits'
coeffcase.linfile = 'my_linfile.fits'
coeffcase.superbiasfile = 'my_superbias.fits'
coeffcase.detector = 'NRCA1'
coeffcase.readpatt = 'SHALLOW4'
coeffcase.nframe = 4
coeffcase.nskip = 1
coeffcase.inverse_lin_file = 'inverse_linearity_coeffs.fits'
coeffcase.outfile = 'groupsat_inverse_lin_from_file.fits'
coeffcase.make_file()
print('Done with case #1')

#Now the slower case, where the inverse linearity coefficients
#are calculated on the fly
coeffcase.inverse_lin_file = None
coeffcase.outfile = 'groupsat_inverse_lin_onthefly.fits'
coeffcase.make_file()
print('Done with case #2')
