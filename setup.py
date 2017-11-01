#!/usr/bin/env python

from setuptools import setup
from setuptools import find_packages

setup(
    name='nircam_calib',
    version = '0.1',
    description = 'NIRCam calibration tools',
    long_description = ('A collection of tools for creating NIRCam'
                        'reference files and performing basic manipulation'
                        'of NIRCam data.'),
    author = 'STScI NIRCam Team',
    author_email = 'hilbert@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python'],
    packages = find_packages(exclude=["examples"]),
    install_requires = [],
    include_package_data = True
    )
