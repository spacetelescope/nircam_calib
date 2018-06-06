#! /usr/bin/env python

'''
Modify a superbias reference file for pipeline testing. Currently
all we should need to do is make sure that the reference file
contains at least one pixel flagged with UNRELIABLE_BIAS, so that
in testing we can make sure that the superbias is still subtracted
from this pixel
'''
from jwst.datamodels import SuperBiasModel

infile = '../reffiles/A1_superbias_from_list_of_biasfiles.list.fits'
outfile = 'modified_superbias_reffile.fits'


#read in reffile
sb = SuperBiasModel(infile)

#set the unreliable_bias flag for one pixel
sb.dq[50,50] = 2

#save
sb.save(outfile)
