#! /usr/bin/env python

"""This module contains functions concerning file i/o
"""

def integration_times(hdulist):
    """Retrieve information on integration start, mid, and end times from
    the INT_TIMES extension of the input file

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        HDU List from an opened FITS file

    Returns
    -------
    starting : list
        List of integration start times, in MJD UTC

    mid : list
        List of integration mid times, in MJD UTC

    ending : list
        List of integration ending times, in MJD UTC
    """
    int_times = hdulist['INT_TIMES'].data
    starting = int_times['int_start_MJD_UTC']
    mid = int_times['int_mid_MJD_UTC']
    ending = int_times['int_end_MJD_UTC']
    return starting, mid, ending