#! /usr/bin/env python

"""Try generating a 2MASS catalog for the region to be viewed in car-19"""

import copy

from astroquery.vizier import Vizier
import astropy.units as u
from astropy.io import ascii
import astropy.coordinates as coord
from astropy.table import Table, Column

central_ra = 80.482577916
central_dec = -69.493958333
cat = Vizier.query_region(coord.SkyCoord(ra=central_ra, dec=central_dec,
                                         unit=(u.deg, u.deg), frame='icrs'),
                          width="2m", catalog=["2MASS"])


ra = cat[0]['RAJ2000']
dec = cat[0]['DEJ2000']
jmag = cat[0]['Jmag']
hmag = cat[0]['Hmag']
kmag = cat[0]['Kmag']

# Create Mirage catalog. Put all filters in one catalog
mirage_cat = Table()
mirage_cat.add_column(ra)
mirage_cat['RAJ2000'].name = 'x_or_RA'
mirage_cat.add_column(dec)
mirage_cat['DEJ2000'].name = 'y_or_Dec'

# F200W catalog uses K magnitudes with no changes
mirage_cat.add_column(kmag)
mirage_cat['Kmag'].name = 'nircam_f200w_magnitude'

# F480M catalog uses K magnitudes also
mirage_cat.add_column(kmag)
mirage_cat['Kmag'].name = 'nircam_f480m_magnitude'

# F356W magnitudes - subtract 0.75 from Kmags
mags = copy.deepcopy(kmag.data.data)
mags -= 0.75
cat[0]['Kmag'] = mags
mirage_cat.add_column(cat[0]['Kmag'])
mirage_cat['Kmag'].name = 'nircam_f356w_magnitude'

# F162M magnitudes - add 0.75 to H band mags
mags = copy.deepcopy(hmag.data.data)
mags += 0.75
cat[0]['Hmag'] = mags
mirage_cat.add_column(cat[0]['Hmag'])
mirage_cat['Hmag'].name = 'nircam_f162m_magnitude'

ascii.write(mirage_cat, "../catalogs/2MASS_v2_pointsource_all_filters.cat", overwrite=True)
