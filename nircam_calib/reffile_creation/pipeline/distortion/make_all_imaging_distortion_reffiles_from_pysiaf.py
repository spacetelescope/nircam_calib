#import nircam_reftools as ref
import nircam_distortion_reffiles_from_pysiaf as ref
import numpy as np

nrc_a_apertures = ['NRCA{}_FULL' for i in range(5)]
nrc_b_apertures = ['NRCB{}_FULL' for i in range(5)]
nrc_apertures = nrc_a_apertures + nrc_b_apertures

hist = ("This reference file was created from the distortion coefficients contained in pysiaf "
        "as of 21 Oct 2019. This includes the fix for the round trip (V2,V3 -> x,y -> V2, V3) error "
        "present in previous versions of the coefficients. This update is described in JIRA issues: "
        "JWSTSIAF-161, JWSTSIAF-123.")

for aperture in nrc_apertures:
    detector, apname = aperture.split('_')
    outname = os.path.join('reffiles_Oct2019/', '{}_distortion.asdf'.format(aperture))
    ref.create_nircam_distortion(detector, apname, apname, outname)

#for entry in siaf:
#    apername = entry['col4']
#    opgs = str(entry['col18'])
#    if apername != 'na' and opgs != 'na' and 'GRISMR' not in apername and 'GRISMC' not in apername:
#        aper_minus_det = apername[6:]
#        det = apername[0:5]
#        try:
#            detnum = np.int(det[-1])
#            outname = 'reffiles_27Oct2016/' + det + '_' + aper_minus_det + '_distortion.asdf'
#            print('running {},{}'.format(apername,opgs))
#            ref.create_nircam_distortion('NIRCam_SIAF_2016-09-29.csv',det,aper_minus_det,opgs,outname)
#        except:
#            #print('skipping {}.'.format(entry))
#            pass
