#!/usr/bin/env python

from distutils.core import setup

setup(name='nircam_calib',
      version='0.1',
      description='NIRCam calibration tools',
      author='Bryan Hilbert',
      author_email='hilbert@stsci.edu',
      packages=['nircam_calib', 'nircam_calib.pipeline_testing','nircam_calib.pipeline_testing.assign_wcs','nircam_calib.pipeline_testing.dq_init','nircam_calib.pipeline_testing.refpix_corr','nircam_calib.reffile_creation','nircam_calib.reffile_creation','nircam_calib.reffile_creation.ETC','nircam_calib.reffile_creation.pipeline','nircam_calib.reffile_creation.throughput','nircam_calib.reffile_creation.pipeline.badpix_map','nircam_calib.reffile_creation.pipeline.crosstalk','nircam_calib.reffile_creation.pipeline.dark','nircam_calib.reffile_creation.pipeline.distortion','nircam_calib.reffile_creation.flats','nircam_calib.reffile_creation.pipeline.gain','nircam_calib.reffile_creation.pipeline.ipc','nircam_calib.reffile_creation.pipeline.linearity','nircam_calib.reffile_creation.pipeline.readnoise','nircam_calib.reffile_creation.pipeline.saturation','nircam_calib.reffile_creation.pipeline.superbias','nircam_calib.tools','nircam_calib.tools.WCS','nircam_calib.tools.image_orientation','nircam_calib.tools.math','nircam_calib.tools.misc','nircam_calib.tools.ramp_plots','nircam_calib.tools.reorganize_integrations'])
