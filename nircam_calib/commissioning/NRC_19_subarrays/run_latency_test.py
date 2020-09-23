#! /usr/bin/env python

"""Run the latency test code in nircam_calib
"""

from glob import glob
from nircam_calib.commissioning.NRC_19_subarrays import latency_investigation as li
import os

files = sorted(glob('Pipeline_Level1/jw01068001001*nrcb5_rateints.fits'))
li.check(files)

