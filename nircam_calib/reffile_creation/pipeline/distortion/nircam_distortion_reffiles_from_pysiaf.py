"""
This module contains functions to create NIRCAM reference files, using distortion
information in pysiaf (rather then the previous version, which used excel spreadsheets)

NIRCAM model:
im.meta.instrument.name :'NIRCAM'
im.meta.instrument.channel : 'SHORT'
im.meta.instrument.module : 'B'
im.meta.instrument.detector : 'NRCB1'
im.meta.exposure.type : 'NRC_IMAGE'

???
im.meta.instrument.pupil : 'FLAT' (would be GRISMR or GRISMC for slitless)

Transform Paths for Imaging mode:

science --> ideal --> V2V3
V2V3 --> ideal --> science

Where the "science" frame has units of distorted pixels. The "ideal" frame
is distortion-free and is the distance in arcseconds from 0,0. V2/V3 is also
in units of arcseconds from the refrence pixel.


"""
from asdf import AsdfFile
from astropy.io import ascii
from astropy.modeling.models import Polynomial2D, Mapping, Shift
import astropy.units as u
from jwst.datamodels import DistortionModel
from stdatamodels import util
from mirage.utils.siaf_interface import sci_subarray_corners
import numpy as np
import pysiaf

#import read_siaf_table


def create_nircam_distortion(detector, aperture, outname, sci_pupil,
                             sci_subarr, sci_exptype, history_entry,
                             author=None, descrip=None, pedigree=None,
                             useafter=None, dist_coeffs_file=None,
                             siaf_xml_file=None):
    """
    Create an asdf reference file with all distortion components for the NIRCam imager.

    NOTE: The IDT has not provided any distortion information. The files are constructed
    using ISIM transformations provided/(computed?) by the TEL team which they use to
    create the SIAF file.
    These reference files should be replaced when/if the IDT provides us with distortion.

    Parameters
    ----------
    detector : str
        NRCB1, NRCB2, NRCB3, NRCB4, NRCB5, NRCA1, NRCA2, NRCA3, NRCA4, NRCA5

    aperture : str
        Name of the aperture/subarray. (e.g. FULL, SUB160, SUB320, SUB640, GRISM_F322W2)

    outname : str
        Name of output file.
    siaf_xml_file : str
        Name of SIAF xml file to use in place of the default SIAF version from pysiaf.
        If None, the default version in pysiaf will be used.

    sci_pupil : list
        Pupil wheel values for which this distortion solution applies

    sci_subarr : list
        List of subarray/aperture names to which this distortion solution applies

    sci_exptype : list
        List of exposure types to which this distortion solution applies

    history_entry : str
        Text to be added as a HISTORY entry in the output reference file

    author : str
        Value to place in the output file's Author metadata entry

    descrip : str
        Text to place in the output file's DECRIP header keyword

    pedgree : str
        Value to place in the output file's PEDIGREE header keyword

    useafter : str
        Value to place in the output file's USEAFTER header keyword (e.g. "2014-10-01T00:00:01")
    dist_coeffs_file : str
        Name of ascii file (nominally output by jwst_fpa package) containing distortion
        coefficients. If this is provided, the coefficients in this file are used, rather
        than those in pysiaf.

    Examples
    --------

    """
    degree = 5  # distotion in pysiaf is a 5th order polynomial
    numdet = detector[-1]
    module = detector[-2]
    channel = 'SHORT'
    if numdet == '5':
        channel = 'LONG'

    full_aperture = detector + '_' + aperture

    # Get Siaf instance for detector/aperture
    if siaf_xml_file is None:
        print('Using default SIAF version in pysiaf.')
        inst_siaf = pysiaf.Siaf('nircam')
    else:
        print(f'SIAF to be loaded from {siaf_xml_file}...')
        inst_siaf = pysiaf.Siaf(filename=siaf_xml_file, instrument='nircam')

    siaf = inst_siaf[full_aperture]


    # Find the distance between (0,0) and the reference location
    xshift, yshift = get_refpix(inst_siaf, full_aperture)

    # *****************************************************
    # If the user provides files containing distortion coefficients
    # (as output by the jwst_fpa package), use those rather than
    # retrieving coefficients from siaf.
    if dist_coeffs_file is not None:
        coeff_tab = read_distortion_coeffs_file(dist_coeffs_file)
        xcoeffs = convert_distortion_coeffs_table(coeff_tab, 'Sci2IdlX')
        ycoeffs = convert_distortion_coeffs_table(coeff_tab, 'Sci2IdlY')
        inv_xcoeffs = convert_distortion_coeffs_table(coeff_tab, 'Idl2SciX')
        inv_ycoeffs = convert_distortion_coeffs_table(coeff_tab, 'Idl2SciY')
    elif dist_coeffs_file is None:
        xcoeffs, ycoeffs = get_distortion_coeffs('Sci2Idl', siaf)
        inv_xcoeffs, inv_ycoeffs = get_distortion_coeffs('Idl2Sci', siaf)

    # V3IdlYAngle and V2Ref, V3Ref should always be taken from the latest version
    # of SIAF, rather than the output of jwst_fpa. Separate FGS/NIRISS analyses must
    # be done in order to modify these values.
    v3_ideal_y_angle = siaf.V3IdlYAngle * np.pi / 180.

    # *****************************************************
    # "Forward' transformations. science --> ideal --> V2V3
    #label = 'Sci2Idl'
    ##from_units = 'distorted pixels'
    ##to_units = 'arcsec'

    #xcoeffs, ycoeffs = get_distortion_coeffs(label, siaf)

    sci2idlx = Polynomial2D(degree, **xcoeffs)
    sci2idly = Polynomial2D(degree, **ycoeffs)

    # Get info for ideal -> v2v3 or v2v3 -> ideal model
    parity = siaf.VIdlParity
    #v3_ideal_y_angle = siaf.V3IdlYAngle * np.pi / 180.
    idl2v2v3x, idl2v2v3y = v2v3_model('ideal', 'v2v3', parity, v3_ideal_y_angle)

    # Finally, we need to shift by the v2,v3 value of the reference
    # location in order to get to absolute v2,v3 coordinates
    v2shift, v3shift = get_v2v3ref(siaf)

    # *****************************************************
    # 'Reverse' transformations. V2V3 --> ideal --> science
    #label = 'Idl2Sci'
    ##from_units = 'arcsec'
    ##to_units = 'distorted pixels'

    #xcoeffs, ycoeffs = get_distortion_coeffs(label, siaf)

    idl2scix = Polynomial2D(degree, **inv_xcoeffs)
    idl2sciy = Polynomial2D(degree, **inv_ycoeffs)

    # Get info for ideal -> v2v3 or v2v3 -> ideal model
    #parity = siaf.VIdlParity
    #v3_ideal_y_angle = siaf.V3IdlYAngle * np.pi / 180.
    v2v32idlx, v2v32idly = v2v3_model('v2v3', 'ideal', parity, v3_ideal_y_angle)

    ##"Forward' transformations. science --> ideal --> V2V3
    #sci2idlx, sci2idly, sciunit, idlunit = read_siaf_table.get_siaf_transform(coefffile,full_aperture,'science','ideal', 5)
    #idl2v2v3x, idl2v2v3y = read_siaf_table.get_siaf_v2v3_transform(coefffile,full_aperture,from_system='ideal')

    ##'Reverse' transformations. V2V3 --> ideal --> science
    #v2v32idlx, v2v32idly = read_siaf_table.get_siaf_v2v3_transform(coefffile,full_aperture,to_system='ideal')
    #idl2scix, idl2sciy, idlunit, sciunit = read_siaf_table.get_siaf_transform(coefffile,full_aperture,'ideal','science', 5)

    # Now create a compound model for each with the appropriate inverse
    sci2idl = Mapping([0, 1, 0, 1]) | sci2idlx & sci2idly
    sci2idl.inverse = Mapping([0, 1, 0, 1]) | idl2scix & idl2sciy

    idl2v2v3 = Mapping([0, 1, 0, 1]) | idl2v2v3x & idl2v2v3y
    idl2v2v3.inverse = Mapping([0, 1, 0, 1]) | v2v32idlx & v2v32idly

    # Now string the models together to make a single transformation

    # We also need
    # to account for the difference of 1 between the SIAF
    # coordinate values (indexed to 1) and python (indexed to 0).
    # Nadia said that this shift should be present in the
    # distortion reference file.

    core_model = sci2idl | idl2v2v3

    # Now add in the shifts to create the full model
    # including the shift to go from 0-indexed python coords to
    # 1-indexed

    # SIAF coords
    index_shift = Shift(1)
    model = index_shift & index_shift | xshift & yshift | core_model | v2shift & v3shift

    # Since the inverse of all model components are now defined,
    # the total model inverse is also defined automatically

    # In the reference file headers, we need to switch NRCA5 to
    # NRCALONG, and same for module B.
    if detector[-1] == '5':
        detector = detector[0:4] + 'LONG'

    # Save using the DistortionModel datamodel
    d = DistortionModel(model=model, input_units=u.pix,
                        output_units=u.arcsec)

    #Populate metadata

    # Keyword values in science data to which this file should
    # be applied
    p_pupil = ''
    for p in sci_pupil:
        p_pupil = p_pupil + p + '|'

    p_subarr = ''
    for p in sci_subarr:
        p_subarr = p_subarr + p + '|'

    p_exptype = ''
    for p in sci_exptype:
        p_exptype = p_exptype + p + '|'

    d.meta.instrument.p_pupil = p_pupil
    d.meta.subarray.p_subarray = p_subarr
    d.meta.exposure.p_exptype = p_exptype

    #d.meta.instrument.p_pupil = "CLEAR|F162M|F164N|F323N|F405N|F470N|"
    #d.meta.p_subarray = "FULL|SUB64P|SUB160|SUB160P|SUB320|SUB400P|SUB640|SUB32TATS|SUB32TATSGRISM|SUB8FP1A|SUB8FP1B|SUB96DHSPILA|SUB96DHSPILB|SUB64FP1A|SUB64FP1B|"
    #d.meta.exposure.p_exptype = "NRC_IMAGE|NRC_TSIMAGE|NRC_FLAT|NRC_LED|NRC_WFSC|"

    # metadata describing the reference file itself
    d.meta.title = "NIRCam Distortion"
    d.meta.instrument.name = "NIRCAM"
    d.meta.instrument.module = module
    d.meta.instrument.channel = channel
    d.meta.instrument.detector = detector
    d.meta.telescope = 'JWST'
    d.meta.subarray.name = 'FULL'

    if pedigree is None:
        d.meta.pedigree = 'GROUND'
    else:
        #if pedigree.upper() not in ['DUMMY', 'GROUND', 'INFLIGHT']:
        #    raise ValueError("Bad PEDIGREE value.")
        d.meta.pedigree = pedigree.upper()

    d.meta.reftype = 'DISTORTION'

    if author is None:
        author = "B. Hilbert"
    d.meta.author = author

    d.meta.litref = "https://github.com/spacetelescope/nircam_calib/nircam_calib/reffile_creation/pipeline/distortion/nircam_distortion_reffiles_from_pysiaf.py"

    if descrip is None:
        d.meta.description = "TEST OF UPDATED CODE"
    else:
        d.meta.description = descrip

    #d.meta.exp_type = exp_type
    if useafter is None:
        d.meta.useafter = "2014-10-01T00:00:01"
    else:
        d.meta.useafter = useafter

    # To be ready for the future where we will have filter-dependent solutions
    d.meta.instrument.filter = 'N/A'

    # Create initial HISTORY ENTRY
    sdict = {'name': 'nircam_distortion_reffiles_from_pysiaf.py',
             'author': author,
             'homepage': 'https://github.com/spacetelescope/nircam_calib',
             'version': '0.0'}

    entry = util.create_history_entry(history_entry, software=sdict)
    d.history = [entry]

    #Create additional HISTORY entries
    #entry2 = util.create_history_entry(history_2)
    #d.history.append(entry2)

    d.save(outname)
    print("Output saved to {}".format(outname))


def convert_distortion_coeffs_table(tab, label):
    """Convert the one set of coefficients output by read_distortion_coeffs_file into
    the proper dictionary format to be then saved in the reference file.

    Parameters
    ----------
    tab : astropy.table.Table
        Table containing all coefficients

    label : str
        e.g. 'Sci2IdlX'

    Returns
    -------
    coeffs : dict
        Dictionary of one set of coefficients (e.g. Sci2IdlX) from the input table.
        Keys list the polynomial order numbers (e.g. c3_4)
    """
    coeffs = {}
    for row in tab:
        i = int(row["siaf_index"][0])
        j = int(row["siaf_index"][1])
        key = f'c{i-j}_{j}'
        coeffs[key] = row[label]

    return coeffs


def get_distortion_coeffs(direction_label, siaf):
    """Retrieve the requested set of distortion coefficients from Siaf
    and package into a dictionary

    Paramters
    ---------

    direction_label ; str
        Either 'Sci2Idl' or 'Idl2Sci'

    Returns
    -------
    x_coeffs : dict
        Dictionary containing x coefficients

    y_coeffs : dict
        Dictionary containing y coefficients
    """
    # Create dictionaries of distortion coefficients
    x_coeffs = {}
    y_coeffs = {}

    for i in range(6):
        for j in range(0, i+1):
            xcolname = '{}X{}{}'.format(direction_label, i, j)
            ycolname = xcolname.replace('X', 'Y')
            key = 'c{}_{}'.format(i-j, j)
            x_coeffs[key] = eval('siaf.{}'.format(xcolname))
            y_coeffs[key] = eval('siaf.{}'.format(ycolname))
    return x_coeffs, y_coeffs


def get_refpix(siaf_instance, apername):
    """Return the reference location within the given aperture

    Parameters
    ----------
    siaf_instance : pysiaf.Siaf('nircam')
    """
    siaf_aperture = siaf_instance[apername]
    xref = siaf_aperture.XSciRef
    yref = siaf_aperture.YSciRef
    #return Shift(-xref) & Shift(-yref)

    # Check to see if we can use coeffs from a subarray aperture
    # and have them apply to all apertures. Need to update the shift
    # in that case by adding the distance from detector (0, 0) to the
    # lower left corner of the aperture
    #siaf = pysiaf.Siaf('nircam')
    xc, yc = sci_subarray_corners('nircam', apername, siaf=siaf_instance, verbose=False)
    llx, urx = xc
    lly, ury = yc
    print('Lower left corner x and y:', llx, lly)

    return Shift(-xref-llx) & Shift(-yref-lly)


def get_v2v3ref(siaf_instance):
    """Return v2 and v3 at the reference location
    These are arcsec in the SIAF file. Convert to degrees

    Parameters
    ----------
    siaf_instance : pysiaf.Siaf[aperture]

    Returns
    -------
    v2 : astropy.modeling.models.Shift
        Shift between V2 at reference location and V2=0

    v3 : astropy.modeling.models.Shift
        Shift between V3 at reference location and V3=0
    """
    v2ref = siaf_instance.V2Ref
    v3ref = siaf_instance.V3Ref
    return Shift(v2ref) & Shift(v3ref)


def read_distortion_coeffs_file(filename):
    """Read the file containing the table of distortion coefficients

    Example:

    # NIRCAM distortion coefficient file

    # Source file: jw01144001001_01101_00001_nrcb4_cal.fits
    # Aperture: NRCB4_FULL
    # Filter/Pupil: F200W/CLEAR
    # Generated 2022-01-25T16:16:04.533 utc
    # by verap
    #
      AperName , siaf_index , exponent_x , exponent_y ,                Sci2IdlX ,                Sci2IdlY ,                Idl2SciX ,                Idl2SciY
    NRCB4_FULL ,         00 ,          0 ,          0 ,                     0.0 ,                     0.0 ,                     0.0 ,                     0.0
    NRCB4_FULL ,         10 ,          1 ,          0 ,    0.031281790934487304 ,  0.00014142457551002174 ,      31.967141087158005 ,    -0.14404661727118445
    NRCB4_FULL ,         11 ,          0 ,          1 ,                     0.0 ,    0.031447520345431045 ,   3.469446951953614e-18 ,       31.79851632221204
    NRCB4_FULL ,         20 ,          2 ,          0 ,  -6.709581883542899e-08 ,    6.38422037163669e-08 ,     0.00215373180436267 ,  -0.0020935324940174927
    NRCB4_FULL ,         21 ,          1 ,          1 , -2.1509448922459775e-07 ,  -9.112311025594254e-08 ,     0.00702920876108879 ,   0.0028750871441249734

    Parameters
    ----------
    filename : str
        Name of text file containing the data.

    Returns
    -------
    tab : astropy.table.Table
        Table containing distortion coefficients.
    """
    converters = {'siaf_index': [ascii.convert_numpy(str)]}
    tab = ascii.read(filename, format='csv', header_start=7, data_start=8, converters=converters)

    # Catch if the file format changes
    if 'Sci2IdlX' not in tab.colnames:
        raise ValueError("distortion_coeffs_file was not read correctly. You may need to adjust header and data starting lines.")

    return tab


def read_siaf_params_file(filename):
    """Read the file containing V2/V3 Ref and Angle data.

    Currently, this file is not used when populating a new distortion
    reference file.

    Example file:

    V2Ref       = 0.028307615938974295
    V3Ref       = 0.03990523039978199
    V3SciXAngle = -90.2219128530577
    V3SciYAngle = 0.03711886628887312

    Parameters
    ----------
    filename : str
        Name of text file containing the data.

    Returns
    -------
    results : dict
        Dictionary of values with parameter names for keys
    """
    results = {}

    tab = ascii.read(filename)
    for row in tab:
        results[row['col1']] = row['col3']

    return results


def v2v3_model(from_sys, to_sys, par, angle):
    """
    Creates an astropy.modeling.Model object
    for the undistorted ("ideal") to V2V3 coordinate translation
    """
    if from_sys != 'v2v3' and to_sys != 'v2v3':
        raise ValueError("This function is designed to generate the transformation either to or from V2V3.")

    # Cast the transform functions as 1st order polynomials
    xc = {}
    yc = {}
    if to_sys == 'v2v3':
        xc['c1_0'] = par * np.cos(angle)
        xc['c0_1'] = np.sin(angle)
        yc['c1_0'] = (0.-par) * np.sin(angle)
        yc['c0_1'] = np.cos(angle)

    if from_sys == 'v2v3':
        xc['c1_0'] = par * np.cos(angle)
        xc['c0_1'] = par * (0. - np.sin(angle))
        yc['c1_0'] = np.sin(angle)
        yc['c0_1'] = np.cos(angle)

    #0,0 coeff should never be used.
    xc['c0_0'] = 0
    yc['c0_0'] = 0

    xmodel = Polynomial2D(1, **xc)
    ymodel = Polynomial2D(1, **yc)

    return xmodel, ymodel
