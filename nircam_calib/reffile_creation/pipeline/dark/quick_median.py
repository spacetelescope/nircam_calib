#! /usr/bin/env python

'''
From a list of input ramps, generate a quick median ramp
'''

import numpy as np
from astropy.io import fits

class Median():
    def __init__(self):
        self.files = []
        self.ngroupopen = 15
        
    def calculate(self):
        with fits.open(self.files[0]) as h:
            nint,ngrp,ny,nx = h[1].data.shape

        for file in self.files:
            with fits.open(file) as h:
                data = h[1].data[:,,:,:]
