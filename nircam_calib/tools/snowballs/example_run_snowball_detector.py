#! /usr/bin/env python

import os
from glob import glob
import snowball_detector as sd

# Collect files to examine
sw_files = glob('/path/to/my/nircam/sw/data/*nrc??_jump*fits')
lw_files = glob('/path/to/my/nircam/lw/data/*nrc?long_jump*fits')

# Run on the SW files
sd.SnowballDetector(sw_files, min_snowball_area=25, histogram_filebase='sw_snowball_data',
                    gallery_filename='sw_snowball_data_gallery_minarea25.png',
                    summary_file='sw_snowball_data_summary_minarea25.ecsv',
                    snowball_table='sw_snowball_data_table_minarea25.ecsv',
                    max_ellipticity=1.0, min_total_signal=100000.,
                    save_segment_map=True, make_gallery_view=True,
                    output_dir='/save/results/here/sw/')


# Run on the LW files
sd.SnowballDetector(lw_files, min_snowball_area=25, histogram_filebase='lw_snowball_data',
                    gallery_filename='lw_snowball_data_gallery_minarea25.png',
                    summary_file='lw_snowball_data_summary_minarea25.ecsv',
                    snowball_table='lw_snowball_data_table_minarea25.ecsv',
                    max_ellipticity=1.0, min_total_signal=100000.,
                    save_segment_map=True, make_gallery_view=True,
                    output_dir='/save/results/here/lw/')
