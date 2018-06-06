#! /usr/bin/env python

'''
Test the persistence step of the pipeline. Written during 
testing of build 7.1

Validation Part 1:
Check that trapsfilled file is generated correctly
Check that it’s (the step?) used correctly


SSB documentation:
Based on a model, this step computes the number of traps 
that are expected to have captured or released a charge 
during an exposure. The released charge is proportional 
to the persistence signal, and this will be subtracted 
(group by group) from the science data. An image of the 
number of filled traps at the end of the exposure will 
be written as an output file, in order to be used as input 
for correcting the persistence of a subsequent exposure.

Input
The input science file is a RampModel.

A trapsfilled file (TrapsFilledModel) may optionally be 
passed as input as well. This normally would be specified 
unless the previous exposure with the current detector was 
taken more than several hours previously, that is, so long 
ago that persistence from that exposure could be ignored.

Output
The output science file is a RampModel, a persistence-corrected 
copy of the input data.

A second output file will be written, with suffix “_trapsfilled”. 
This is a TrapsFilledModel, the number of filled traps at each 
pixel at the end of the exposure. This takes into account the 
capture of charge by traps due to the current science exposure, 
as well as the release of charge from traps shown in the input 
trapsfilled file, if one was specified.

If the user specified save_persistence=True, a third output file 
will be written, with suffix “_output_pers”. This is a RampModel 
matching the output science file, but this gives the persistence 
that was subtracted from each group in each integration.

input file -> run persistence step -> output hdu and file -> 
run tests against...what truth?

'''

import pytest
import os
import numpy as np
from astropy.io import fits
from jwst import datamodels
from jwst.persistence import PersistenceStep
from jwst.datamodels import TrapsFilledModel
from jwst.datamodels import dqflags
from persistence_utils import subtracted_persist
from persistence_utils import dq_flagged_pix

@pytest.fixture(scope='function')
def in_hdul(input_file):
    # open the input_file defined above once for each module
    yield fits.open(input_file)


@pytest.fixture(scope='function')
def out_hdul(in_hdul):
    outname = '_persist.'.join(in_hdul[0].header['filename'].split('.'))
    yield fits.open(outname)
    # Delete output fits file after module is finished
    #os.remove(outname)

    
@pytest.fixture(scope='function')
def traps_hdul(in_hdul):
    tfile = '_trapsfilled.'.join(in_hdul[0].header['filename'].split('.'))
    yield fits.open(tfile)
    #os.remove(tfile)

    
@pytest.fixture(scope='function')
def pers_hdul(in_hdul):
    pfile = '_output_pers.'.join(in_hdul[0].header['filename'].split('.'))
    yield fits.open(pfile)
    #os.remove(pfile)
    

def test_run_persist_step(in_hdul,trapsfilled):
    outfile = in_hdul[0].header['FILENAME'].replace('.fits','_persist.fits')
    PersistenceStep.call(in_hdul,save_persistence=True,\
                         output_file=outfile,save_results=True,\
                         input_trapsfilled=trapsfilled)


def test_trapsfilled_shape(in_hdul,traps_hdul):
    '''Check to see that the trapsfilled
    file was created. 
    input filename is the name of the 
    file input into the persistence step'''
    x,y = in_hdul['SCI'].data.shape[-2:]
    print("Science data shape (x,y) = ({},{})".format(x,y))
    assert traps_hdul['SCI'].data.shape == (3,y,x)

    
def test_output_pers_shape(in_hdul,pers_hdul):
    '''Check that the optional output file
    "_output_pers.fits" was created if
    the save_persistence option in the persistence
    step was set to True. (Assume this test will
    only be called in instances when save_persistence
    is True'''
    opshape = pers_hdul['SCI'].data.shape
    print("Output_pers data shape: {}".format(opshape))
    assert opshape == in_hdul['SCI'].data.shape


def test_subtracted_persist(in_hdul, out_hdul, pers_hdul):
    '''Check that the signal values contained in the 
    output_pers file are indeed subtracted from the original
    input file.'''
    assert subtracted_persist(in_hdul, out_hdul, pers_hdul)


def test_dq_flagged_pix(out_hdul,pers_hdul,flagthresh):
    '''Pixels that have more persistence signal than flag_pers_cutoff
    should be flagged in the DQ array of the output file. The default
    value of flag_pers_cutoff is 40 DN'''
    assert dq_flagged_pix(out_hdul,pers_hdul,flagthresh)
    
