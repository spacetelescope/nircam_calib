"""
This file contains the functions which will be used to test the persistence step
of the JWST Calibration Pipeline
"""

import numpy as np
from jwst.datamodels import dqflags


def subtracted_persist(in_hdul, out_hdul, pers_hdul):
    '''
    Check that the signal values contained in the 
    output_pers file are indeed subtracted from the original
    input file.
    '''
    result = np.allclose(out_hdul[1].data,in_hdul[1].data - pers_hdul[1].data)
    return result


def dq_flagged_pix(out_hdul,pers_hdul,flagthresh):
    '''Pixels that have more persistence signal than flag_pers_cutoff
    should be flagged in the DQ array of the output file. The default
    value of flag_pers_cutoff is 40 DN'''
    # Check only integration #1
    pdata = pers_hdul['SCI'].data[0,:,:,:]
    # Keep only the maximum persistence value
    # for each pixel
    if ((flagthresh is not None) and (flagthresh > 0)):
        collapsed = np.max(pdata,axis=0)
        flagged = collapsed > flagthresh
        dq_data = out_hdul['PIXELDQ'].data
        print(("{} pixels have persistence values above the threshold "
               "of {}.".format(np.sum(flagged),flagthresh)))
        result = np.all(dq_data[flagged] & dqflags.pixel['DO_NOT_USE'] > 0)
        return result
    else:
        print("Flagthresh is {}".format(flagthresh))
        return True
