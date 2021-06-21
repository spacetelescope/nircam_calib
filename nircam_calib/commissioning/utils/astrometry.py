#! /usr/bin/env python

"""This module contains functions related to astrometry that may be generally
useful to analysis of multiple CARs

NOTE: In order to create an attitude matrix, see mirage.siaf_interface.get_siaf_information()
"""
import pysiaf


def RADec_To_XY(ra, dec, array_name, attitude_matrix):
        """Translate backwards, RA, Dec to V2, V3. If a distortion reference file is
        provided, use that. Otherwise fall back to pysiaf.
        Parameters:
        -----------
        ra : float
            Right ascention value, in degrees, to be translated.
        dec : float
            Declination value, in degrees, to be translated.
        Returns:
        --------
        pixelx : float
            X coordinate value in the aperture corresponding to the input location
        pixely : float
            Y coordinate value in the aperture corresponding to the input location
        """
        siaf = pysiaf.Siaf('nircam')[array_name]
        loc_v2, loc_v3 = pysiaf.utils.rotations.getv2v3(attitude_matrix, ra, dec)

        pixelx, pixely = siaf.tel_to_sci(loc_v2, loc_v3)

        # Subtract 1 from SAIF-derived results since SIAF works in a 1-indexed coord system
        pixelx -= 1
        pixely -= 1
        return pixelx, pixely


def XY_To_RADec(pixelx, pixely, array_name, attitude_matrix):
        """Translate a given x, y location on the detector to RA, Dec. If a
        distortion reference file is provided, use that. Otherwise fall back to
        using pysiaf.
        Parameters:
        -----------
        pixelx : float
            X coordinate value in the aperture
        pixely : float
            Y coordinate value in the aperture
        Returns:
        --------
        ra : float
            Right ascention value in degrees
        dec : float
            Declination value in degrees
        ra_str : str
            Right ascention value in HH:MM:SS
        dec_str : str
            Declination value in DD:MM:SS
        """
        siaf = pysiaf.Siaf('nircam')[array_name]

        # Use SIAF to do the calculations
        #In this case, add 1 to the input pixel values
        # since SIAF works in a 1-indexed coordinate system.
        loc_v2, loc_v3 = siaf.sci_to_tel(pixelx + 1, pixely + 1)

        ra, dec = pysiaf.utils.rotations.pointing(attitude_matrix, loc_v2, loc_v3)

        return ra, dec

