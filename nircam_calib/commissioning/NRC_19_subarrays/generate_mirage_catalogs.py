#! /usr/bin/env python

"""
Generate Mirage-formatted catalogs for CAR-19. Include 2MASS sources as well as a representative
sample from the Besancon model
"""

from mirage.catalogs import create_catalog

# First we need to query the Besancon model from mirage.catalogs import create_catalog
ra = 80.4875  # degrees
dec = -69.4975  # degrees
box_width = 360  # arcseconds
#create_catalog.besancon(ra, dec, box_width, username='hilbert', kmag_limits=(17, 30))


# Once the above is finished and the results saved into a file, the command below can be used
# to create the final catalog
filter_list = ['F200W', 'F356W', 'F164N', 'F480M', 'F322W2', 'WLP4']
cat, mag_column_names = create_catalog.get_all_catalogs(ra, dec, box_width, besancon_catalog_file='besancon_results.cat',
                                                        instrument='NIRCAM', filters=filter_list
                                                        )