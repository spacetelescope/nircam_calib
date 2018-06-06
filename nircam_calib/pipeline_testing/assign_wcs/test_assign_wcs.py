import pytest
import sys, os
import numpy as np
from jwst.assign_wcs import AssignWcsStep
from jwst import datamodels
from asdf import AsdfFile
from astropy.coordinates import SkyCoord
from astropy.io import fits


@pytest.fixture(scope='module')
def in_datamodel(fits_input):
    '''open input file as a datamodel'''
    yield datamodels.open(fits_input[0].header['FILENAME'])
    
    
@pytest.fixture(scope='module')
def out_datamodel(in_datamodel):
    '''open output file'''
    outname = '_assign_wcs.'.join(in_datamodel.meta.filename.split('.'))
    yield datamodels.open(outname)
    os.remove(outname)

    
@pytest.fixture(scope='module')
def dist_reffile(out_datamodel):
    '''determine the distortion reference file used'''
    origname = out_datamodel.meta.ref_file.distortion.name
    if 'crds://' in origname:
        origname = origname.replace('crds://','/grp/crds/cache/references/jwst/')
    yield origname
    
    
def test_assign_wcs_step(in_datamodel):
    '''Run the assign_wcs pipeline step'''
    outfile = in_datamodel.meta.filename.replace('.fits','_assign_wcs.fits')
    AssignWcsStep.call(in_datamodel,output_file = outfile,save_results=True)


def test_asdf_extension(fits_input):
    '''Make sure there is a new ASDF extension in the output file'''
    print("Test to be sure the ASDF extension was added to the file")
    outname = '_assign_wcs.'.join(fits_input[0].header['FILENAME'].split('.'))
    h = fits.open(outname)
    exts = []
    for i in range(1,len(h)):
        exts.append(h[i].header['EXTNAME'])
    assert 'ASDF' in exts
    
        
def test_inserted_wcs_model(out_datamodel):
    '''Check RA, Dec values from output file WCS model 
    are correct for refrence location. Check pixel scale implied
    by the WCS model'''

    expectedra = out_datamodel.meta.wcsinfo.ra_ref
    expecteddec = out_datamodel.meta.wcsinfo.dec_ref

    refloc_x = out_datamodel.meta.wcsinfo.crpix1 - 1
    refloc_y = out_datamodel.meta.wcsinfo.crpix2 - 1

    # pixel scale in arcsec per pixel
    xpixscale = out_datamodel.meta.wcsinfo.cdelt1 * 3600.
    ypixscale = out_datamodel.meta.wcsinfo.cdelt2 * 3600.
    #distscale = np.sqrt(xpixscale**2 + ypixscale**2)
    
    # Tolerance to use when checking if values are close enough
    atol_pix = 0.01
    xatol_arcsec = atol_pix * xpixscale
    yatol_arcsec = atol_pix * ypixscale
    
    #check pixel scale by reporting RA,Dec of adjacent pixels
    exp_type = out_datamodel.meta.exposure.type
    if 'GRISM' not in exp_type: # imaging data
        refra,refdec = out_datamodel.meta.wcs(refloc_x,refloc_y)
        adra,addec   = out_datamodel.meta.wcs(refloc_x+1,refloc_y)
        ad2ra,ad2dec = out_datamodel.meta.wcs(refloc_x,refloc_y+1)

        pos1 = SkyCoord(ra=adra, dec=addec, unit='deg')
        refpos_skycoords = SkyCoord(ra=refra, dec=refdec, unit='deg')
        dist = refpos_skycoords.separation(pos1).value * 3600

        pos2 = SkyCoord(ra=ad2ra, dec=ad2dec, unit='deg')
        dist2 = refpos_skycoords.separation(pos2).value * 3600

        print('Expected RA, Dec at reference location (deg):',expectedra,expecteddec)
        print('Ref loc. RA, Dec (deg):',refra,refdec)
        print('RA, Dec of adjacent pixel (deg):',adra,addec)
        print('Delta Distance, horiz. adjacent pix: (arcsec)',dist)
        print('Delta Distance, vertical adjacent pix: (arcsec)',dist2)
        
        assert np.allclose(expectedra,refra,atol=xatol_arcsec,rtol=0.)
        assert np.allclose(expecteddec,refdec,atol=yatol_arcsec,rtol=0)
        assert np.allclose(dist,xpixscale,atol=yatol_arcsec,rtol=0)
        assert np.allclose(dist2,ypixscale,atol=yatol_arcsec,rtol=0)


    else: # WFSS data
        pupil = out_datamodel.meta.instrument.pupil
        filter = out_datamodel.meta.instrument.filter #to support NIRISS grism
        if pupil[0] == 'G':
            pass
        elif filter[0] == 'G':
            pupil = filter
        grisms_r = ['GRISMR','G150R']
        grisms_c = ['GRISMC','G150C']
        adjpix = (None,None)
        adjwave = (None,None)
        if pupil in grisms_r:
            adjpix = (0,1) # delta x, delta y
            adjwave = (1,0)
        elif pupil in grism_c:
            adjpix = (1,0)
            adjwave = (0,1)
        else:
            print("Grism value of {} is not recognized. Skipping testing.".format(pupil))

        if adjwave[0] is not None:
            # RA, Dec, Wavelength, Order for the reference location
            # pixel in the dispersed image
            refra,refdec,refwave,reforder = out_datamodel.meta.wcs(refloc_x
                                                                   ,refloc_y
                                                                   ,refloc_x
                                                                   ,refloc_y,1)

            # Move one pixel in the dispersion direction in the dispersed
            # image (refloc in direct image). Should have same RA, Dec as
            # reference location in direct image.
            adwavera,adwavedec,adwavewave,adwaveord = out_datamodel.meta.wcs(refloc_x+adjwave[0]
                                                                             ,refloc_y+adjwave[1]
                                                                             ,refloc_x,refloc_y,1)

            # Move one pixel perp to dispersion direction in dispersed
            # image (refloc in direct image). Should have same RA, Dec as
            # reference location in direct image. Should also have same
            # wavelength as reference location.
            adpixra,adpixdec,adpixwave,adpixord = out_datamodel.meta.wcs(refloc_x+adjpix[0],
                                                                         refloc_y+adjpix[1],
                                                                         refloc_x,refloc_y,1)
        
            # Move one pixel in the dispersion direction in both the dispersed
            # AND direct image. The resulting wavelength should be the same as
            # refwave case above.
            dra,ddec,dwave,dord = out_datamodel.meta.wcs(refloc_x+adjwave[0],refloc_y+adjwave[1],
                                                         refloc_x+adjwave[0],refloc_y+adjwave[1],1)

            # Calculate distances to compare with the stated pixel scale
            pos1 = SkyCoord(ra=dra, dec=ddec, unit='deg')
            refpos_skycoords = SkyCoord(ra=refra, dec=refdec, unit='deg')
            ddist = refpos_skycoords.separation(pos1).value * 3600

            print('Expected RA, Dec at reference location (deg):',expectedra,expecteddec)
            print('Ref loc. RA, Dec (deg):',refra,refdec)
            print('RA, Dec of adjacent pixel (deg):',dra,ddec)
            print('Delta Distance, horiz. adjacent pix: (arcsec)',ddist)

            assert np.allclose(expectedra,refra,atol=xatol_arcsec,rtol=0.)
            assert np.allclose(expecteddec,refdec,atol=yatol_arcsec,rtol=0)
            assert np.allclose(refra,adwavera,atol=1e-8,rtol=0.)
            assert np.allclose(refdec,adwavedec,atol=1e-8,rtol=0.)
            assert np.allclose(refra,adpixra,atol=1e-8,rtol=0.)
            assert np.allclose(refdec,adpixdec,atol=1e-8,rtol=0.)
            assert np.allclose(refwave,adpixwave,atol=1e-10,rtol=0.)
            assert np.allclose(refwave,dwave,atol=1e-10,rtol=0.)
            assert np.allclose(ddist,xpixscale,atol=0.0005,rtol=0.)
        
        
def test_wcs_vs_reffile(out_datamodel,dist_reffile):
    '''Test WCS model in the distortion reference file
    matches that in the output file'''

    print("Distortion reference file used: {}".format(dist_reffile))
    distortion = AsdfFile.open(dist_reffile).tree['model']
    reverse = distortion.inverse

    refx = out_datamodel.meta.wcsinfo.crpix1 - 1
    refy = out_datamodel.meta.wcsinfo.crpix2 - 1

    shape = out_datamodel.data.shape
    if len(shape) == 2:
        ylen,xlen = shape
    elif len(shape) == 3:
        nint,ylen,xlen = shape

    inx = [refx]
    iny = [refy]
    fractions = [0.1,0.3,0.5,0.7,0.9]
    for f in fractions:
        xup = int(xlen * f)
        yup = int(ylen * f)
        ydown = ylen - yup
        inx.append(xup)
        iny.append(yup)
        inx.append(xup)
        iny.append(ydown)

    # Convert x,y positions to RA, Dec and back, to see if
    # you recover the same x,y
    reffilex = np.array([])
    reffiley = np.array([])
    modx = np.array([])
    mody = np.array([])
    exp_type = out_datamodel.meta.exposure.type
    if 'GRISM' not in exp_type: # imaging data
        for x,y in zip(inx,iny):
            refra,refdec = distortion(x,y)
            refnewx,refnewy = reverse(refra,refdec)
            reffilex = np.append(reffilex,refnewx)
            reffiley = np.append(reffiley,refnewy)
            fra, fdec = out_datamodel.meta.wcs(x,y)
            fnewx, fnewy = out_datamodel.meta.wcs.backward_transform(fra,fdec)
            modx = np.append(modx,fnewx)
            mody = np.append(mody,fnewy)

        print("Input x coords:")
        print(inx)
        print("X coords from reference file WCS model:")
        print(reffilex)
        print("X coords from output file WCS model:")
        print(modx)

        print("Input y coords:")
        print(iny)
        print("Y coords from reference file WCS model:")
        print(reffiley)
        print("Y coords from output file WCS model:")
        print(mody)
        assert np.allclose(reffilex,modx,atol=1e-8,rtol=0.)
        assert np.allclose(reffiley,mody,atol=1e-8,rtol=0.)

    else:
        print("This test not implemented for GRISM data.")
        
        # The GRISM comparison below would need more work before
        # it could be used. The creation of the WCS model to go into
        # the output file is more complicated than in the imaging case
        # and it's not clear how to create this WCS without simply
        # copying the code in the JWST pipeline
    #else:  # GRISM data
    #    refdirectx = np.array([])
    #    refdirecty = np.array([])
    #    reffilewave = np.array([])
    #    reforder = np.array([])
    #    moddirectx = np.array([])
    #    moddirecty = np.array([])
    #    modwave = np.array([])
    #    modorder = np.array([])
    #
    #    pupil = out_datamodel.meta.instrument.pupil
    #    filter = out_datamodel.meta.instrument.filter #to support NIRISS grism
    #    if pupil[0] == 'G':
    #        pass
    #    elif filter[0] == 'G':
    #        pupil = filter
    #    grisms_r = ['GRISMR','G150R']
    #    grisms_c = ['GRISMC','G150C']
    #            
    #    for x,y in zip(inx,iny):
    #        #refra,refdec,refwave,reford = distortion(x,y,refx,refy,1)
    #        #refnewx,refnewy,refdirx,refdiry,refneworder = reverse(refra,refdec,refwave,reford)
    #        #reffilex = np.append(reffilex,refnewx)
    #        #reffiley = np.append(reffiley,refnewy)
    #        #refdirectx = np.append(refdirectx,refdirx)
    #        #refdirecty = np.append(refdirecty,refdiry)
    #        #reffilewave = np.append(reffilewave,refwave)
    #        #reforder = np.append(reforder,refneworder)
    #        refra,refdec = distortion(x,y)
    #        refnewx,refnewy = reverse(refra,refdec)
    #        if pupil in grisms_r:
    #            reffilex = np.append(reffilex,refnewx)
    #            reffiley = np.append(reffiley,refy)
    #        elif pupil in grisms_c:
    #            reffilex = np.append(reffilex,refx)
    #            reffiley = np.append(reffiley,refnewy)
    #
    #        fra, fdec, fwave, forder = out_datamodel.meta.wcs(x,y,refx,refy,1)
    #        fdispx, fdispy, fdirectx, fdirecty, ford = out_datamodel.meta.wcs.backward_transform(fra,fdec,fwave,forder)
    #
    #        print(x,fdispx,y,fdispy)
    #
    #        
    #        if pupil in grisms_r:
    #            modx = np.append(modx,fdispx)
    #            mody = np.append(mody,refy)
    #            moddirectx = np.append(moddirectx,fdirectx)
    #            moddirecty = np.append(moddirecty,fdirecty)
    #            modwave = np.append(modwave,fwave)
    #            modorder = np.append(modorder,ford)
    #        elif pupil in grisms_c:
    #            modx = np.append(modx,refx)
    #            mody = np.append(mody,fdispy)
    #            moddirectx = np.append(moddirectx,fdirectx)
    #            moddirecty = np.append(moddirecty,fdirecty)
    #            modwave = np.append(modwave,fwave)
    #            modorder = np.append(modorder,ford)
    #            
    #            
    #    print('reffilex, modx, reffiley, mody')
    #    print(reffilex)
    #    print(modx)
    #    print(reffiley)
    #    print(mody)
    #        
    #    assert np.allclose(reffilex,modx,atol=1e-8,rtol=0.)
    #    assert np.allclose(reffiley,mody,atol=1e-8,rtol=0.)
    #    #assert np.allclose(refdirectx,moddirectx,atol=1e-8,rtol=0.)
    #    #assert np.allclose(refdirecty,moddirecty,atol=1e-8,rtol=0.)
    #    #assert np.allclose(reffilewave,modwave,atol=1e-8,rtol=0.)
    #    #assert np.allclose(reforder,modorder,atol=1e-8,rtol=0.)

