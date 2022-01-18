#! /usr/bin/env python

'''
Create NIRCam photom (flux calibration) and
pixel area map (PAM) reference files using refactored synphot
'''

from synphot import SpectralElement, Observation, SourceSpectrum, units
import matplotlib.pyplot as plt
import glob,argparse
import copy
import numpy as np
from math import hypot, atan2, sin
from scipy.stats import sigmaclip
from astropy.io import fits,ascii
from astropy.table import Table,vstack
from astropy import units as u
import pysiaf
from jwst.datamodels import NrcImgPhotomModel, NrcWfssPhotomModel
from scipy.integrate import simps
import sys

class Photom:
    def __init__(self):
        """
        Set up NIRCam-specific information
        """
        self.verbose = True
        self.imaging_throughput_files = None #total sys. throughput: include JWST+NIRCam+filter
        self.grism_throughput_files = None #assume these include JWST+NIRCam+grism+crossing filter
        self.detector = None #'A1-A5, B1-B5'
        self.telescope_area = 25.326 * 10000. #JWST primary area in cm^2
        self.detector_width = 2040 #pixels
        self.detector_length = 2040 #pixels
        self.xcoeffs = None
        self.ycoeffs = None
        self.pixel_area_a2 = None
        self.pixel_area_sr = None
        self.vega = SourceSpectrum.from_vega()
        self.maxlen = 3000 #wavelength and relresponse arrays have to be 3000 elements long
        self.resolving_power = 1525 #for NIRCam at 4.0 microns
        #dictionary of pupil wheel filter:filter wheel filter
        self.filter_dict = {'F164N':'150W2','F323N':'F322W2','F405N':'F444W',
                            'F466N':'F444W','F470N':'F444W'}
        self.crossing_filters = ['F250M','F277W','F300M','F322W2','F335W',
                                 'F356W','F360M','F410M','F430M','F444W',
                                 'F460M','F480M']
        self.disp = {'1':10.04,'2':5.02} #angstroms/pixel - from Nor's code
        self.zeropoint_base = 'NIRCam_zeropoints'
        self.photom_base = 'NIRCam_photom'
        self.grism_thruput = '/grp/jwst/wit/nircam/reference_files/SpectralResponse_2015-Final/Grisms/NIRCam_LW_grism_efficiencies_from_Tom_Greene.csv'
        self.pam_outfile = None
        self.photom_outfile = None
        self.author = 'Nircam Team'
        self.useafter = '2017-01-01'
        self.pam_descrip = 'NIRCam Pixel Area Map'
        self.photom_descrip = 'NIRCam flux conversion reference file'
        self.pedigree = 'GROUND'
        self.pam_history = None
        self.photom_history = None


    def get_pixel_area(self):
        """Use SIAF to get nominal pixel area at reference location"""
        xscale = hypot(self.siaf.Sci2IdlX10, self.siaf.Sci2IdlY10)
        yscale = hypot(self.siaf.Sci2IdlX11, self.siaf.Sci2IdlY11)
        bx = atan2(self.siaf.Sci2IdlX10, self.siaf.Sci2IdlY10)
        return xscale * yscale * sin(bx) * u.arcsecond * u.arcsecond


    def find_gain_pre_launch(self):
        """Use gain values determined by the IDT from bootstrapping
        CV3 data
        """
        if (('5' in self.detector) or ('LONG' in self.detector.upper())):
            self.gain = 1.82
        else:
            self.gain = 2.05


    def find_gain(self):
        """Read in gain data from a fits file"""
        with fits.open(self.gain_file) as h:
            gaindet = h[0].header['DETECTOR'][3:]
            gain2d = h[1].data[4:2044,4:2044]

        # Make sure the gain file detector matches that
        # for the throughput files
        if 'LONG' in gaindet:
            gaindet = gaindet.replace('LONG','5')

        if gaindet != self.detector:
            print("WARNING!! Detector for the gain file ({}) does".format(gaindet))
            print("not match the input detector for the throughput")
            print("files ({}). Quitting.".format(self.detector))
            sys.exit()

        # Calculate the sigma-clipped mean of the gain data
        good = np.isfinite(gain2d)
        clipped,lo,hi = sigmaclip(gain2d[good],low=3,high=3)
        self.gain = clipped.mean()


    def make_photom(self):
        """MAIN FUNCTION"""
        # Check inputs
        if self.imaging_throughput_files is None:
            raise ValueError("No imaging throughput files specified! Quitting.")
        if self.detector is None:
            raise ValueError("No detector specified! Quitting.")

        self.siaf = pysiaf.Siaf('nircam')['NRC{}_FULL'.format(self.detector)]

        # Lists of polynomial coefficients from SIAF
        self.xcoeffs, self.ycoeffs = self.get_coefficients()

        # Find the nominal pixel area from the SIAF
        self.pixel_area_a2 = self.get_pixel_area()

        # Required header keyword outputs
        sterrad_per_arcsec2 = (1. / 3600. * np.pi / 180.)**2
        self.str_per_detector = (self.detector_width * self.detector_length *
                                 self.pixel_area_a2) * sterrad_per_arcsec2 * u.sr
        self.pixel_area_sr = self.pixel_area_a2.to(u.sr)

        # Calculate the appropriate gain value to use
        #self.find_gain()
        self.find_gain_pre_launch()

        # Calculations for imaging filters
        img_tab = self.imaging_calibrations(self.imaging_throughput_files)

        # Now do the calculations for the grisms
        if self.grism_throughput_files is not None:
            gfiles = self.read_listfile(self.grism_throughput_files)
            grism_table = self.grism_cal(gfiles)

            # Combine the imaging and grism photom tables
            print(grism_table.shape,grism_table[0].shape)
            print(img_tab.shape,img_tab[0].shape)

            photom_table = np.append(img_tab,grism_table)
            print(photom_table.shape,photom_table[0].shape)
        else:
            photom_table = img_tab

        # Get module name from the detector name input
        self.module = self.detector[0].upper()

        # Now you should be ready to write out the reference file
        if self.photom_outfile is None:
            self.photom_outfile = 'NIRCam_{}_photom.fits'.format(self.detector)
        self.save_photom_model(photom_table,self.photom_outfile)


    def read_listfile(self,file):
        """
        Read in a file containing a list of files

        Arguments:
        ----------
        file -- name of a text file containing a list of filenames

        Returns:
        --------
        List of filenames
        """
        flist = []
        with open(file) as f:
            for line in f:
                if line.strip() != '':
                    flist.append(line.strip())
        return flist


    def imaging_calibrations(self,listfile):
        """
        Calculate flux cal information for imaging mode

        Arguments:
        ----------
        listfile -- A text file that lists all of the filter
                    throughput curves

        Returns:
        --------
        Astropy table containing flux calibration information in
        the photom reference file format
        """
        self.model = NrcImgPhotomModel()

        results = Table() #Table of filter vegamag zeropoints
        dets = []
        filts = np.array([])
        realfilt = np.array([])
        pupil = np.array([])
        mods = []
        megajy_sterradian = np.array([]) * u.megajansky / u.sr
        unc = np.array([]) * u.megajansky / u.sr
        vegazp = []
        stzp = []
        abzp = []
        flam = []
        fnu = []
        mjy_str = []
        pivot = []
        wave = []
        relresp = []
        photom_results = []

        # Read in imaging throughput file list
        files = []
        with open(listfile) as f:
            for line in f:
                if line.strip() != '':
                    files.append(line.strip())

        nelem = np.zeros(len(files))
        order = np.zeros(len(files)).astype(np.int)
        zeros = np.zeros(self.maxlen)
        waves = np.append([zeros],(len(files)-1)*[zeros],axis=0)
        relresp = waves

        for file in files:
            # Get filter, module name from bandpass filename
            compname = fits.getval(file,'COMPNAME')
            under = compname.find('_')
            filter = compname[0:under]
            module = compname[-1:].upper()
            filts = np.append(filts,filter)
            mods.append(module)

            # Inputs for photom reference file
            if filter in self.filter_dict:
                f = self.filter_dict[filter]
                p = filter
            else:
                f = filter
                p = 'CLEAR'
            realfilt = np.append(realfilt,f)
            pupil = np.append(pupil,p)

            # Calculate zeropoints in various photometric systems
            # as well as pivot wavelength.
            # Uncertainty blindly set to 20%
            vzp, azp, szp, flambda, fn, mjy, pivotwave = self.synphot_calcs(file)
            mjy = mjy  * u.megajansky / self.pixel_area_sr
            megajy_sterradian = np.append(megajy_sterradian,mjy)
            unc = np.append(unc,mjy * 0.2)

            vegazp.append(vzp)
            abzp.append(azp)
            stzp.append(szp)
            flam.append(flambda)
            fnu.append(fn)
            pivot.append(pivotwave)

        # print('megajansky per sterradian: {}'.format(megajy_sterradian))
        # print(realfilt)
        # print(pupil)
        # print(order)
        # print(megajy_sterradian)
        # print(unc)
        # print(nelem)
        # print(waves.shape)
        # print(relresp.shape)

        #photom_response = np.array(zip(realfilt,pupil,megajy_sterradian,unc,order,nelem,waves,relresp),dtype=self.model.phot_table.dtype)
        photom_response = np.array(list(zip(realfilt,pupil,megajy_sterradian.value,unc.value)),dtype=self.model.phot_table.dtype)


        print('Data placed into zipped array:')
        for i in range(4):
            print(photom_response[0][i])

        print("unc is not being brought into photom_response!! Why not????")

        print(unc)
        print(photom_response[0][2])
        print(photom_response[0][3])

        # Results below will be placed into a table and written out in ascii
        # It is not used in the construction of the reference file, and
        # contains additional info that the reference file doesn't need
        results['Filter'] = np.array(filts)
        results['Module'] = np.array(mods)
        results['Detector'] = np.repeat(self.detector,len(mods))
        results['VEGAMAG'] = np.array(vegazp)
        results['ABMAG'] = np.array(abzp)
        results['STMAG'] = np.array(stzp)
        results['FLAMBDA'] = np.array(flam)
        results['FNU'] = np.array(fnu)
        results['Pivot_Wavelength'] = pivot

        # Save table of filter vegamag, stmag, abmag, photflambda, photfnu
        outzpa = self.zeropoint_base + '_mod'+self.detector[0].upper()+'.txt'
        results.write(outzpa,names=['Filter','Module','Detector','VEGAMAG','ABMAG',
                                    'STMAG','FLAM','FNU','Pivot_wavelength'],
                      format='ascii',overwrite=True)
        #outzpb = self.zeropoint_base + '_modB' + '.txt'
        #bmod.write(outzpb,names=['Filter','Module','VEGAMAG','ABMAG','STMAG',
        #                         'FLAM','FNU'],format='ascii')

        #photom_results['filter'] = np.array(realfilt)
        #photom_results['pupil'] = np.array(pupil)
        #photom_results['photmjsr'] = np.array(mjy_str)
        #photom_results['uncertainty'] = np.array(mjy_str) * 0.2
        #photom_results['nelem'] = np.zeros(len(filts),dtype=np.int)
        #photom_results['order'] = np.zeros(len(filts),dtype=np.int)

        #zeroarr = np.zeros(self.maxlen,dtype=np.float32)
        #zerolist = np.append([zeroarr],(len(pupil)-1)*[zeroarr],axis=0)

        #photom_results['wavelength'] = zerolist
        #photom_results['relresponse'] = zerolist
        return photom_response


    def package_grism_data(self):
        #take as input all the information that needs to go into the
        #grism photom file, and put it into CRDS format
        pass


    def grism_cal(self,throughput_files):
        """
        Calculate flux cal outputs for grism mode. Input files should
        contain the total system throughput for the grism plus crossing
        filter. The input file list should contain throughputs for all
        orders of just a single grism+crossing filter.

        Arguments:
        ----------
        throughput_files -- list of text files containing filter throughputs

        Returns:
        --------
        Astropy table containing flux calibration information for grism mode
        """
        print('need to update this for synphot')
        self.model = NrcWfssPhotomModel()

        allrows = []
        filters = np.array([])
        pupils = np.array([])
        orders = np.array([])
        fluxes = np.array([])
        uncs = np.array([])
        nelems = np.array([])
        waves = np.array([])
        resps = np.array([])
        #waves = np.expand_dims(waves,axis=0)
        #resps = np.expand_dims(resps,axis=0)

        for file in throughput_files:
            # Read in necessary file info
            with fits.open(file) as h:
                cname = h[0].header['COMPNAME']

            junk,filter,mod = cname.split('_')
            mod = mod[-1].upper()
            pupil = 'GRISMR'
            print('Eventually need to be smarter about pupil value!!!!!')

            print("opening {}".format(file))
            bp = SpectralElement.from_file(file)

            # Now reduce the PCE curve by a factor equal to
            # the gain, in order to get the output to translate
            # from ADU/sec rather than from e/sec
            bp.model.lookup_table = bp.model.lookup_table / self.gain

            # Pivot wavelength in microns
            pwave = bp.pivot()
            #pwave = pwave.to(u.micron)

            # Need to find the order here, so that the proper
            # dispersion value is used.
            ord = file.find('order')
            order = file[ord+5]
            dispersion = self.disp[order]

            #find pivot? mean? center? effective? wavelength
            #denom = self.h * self.c / eff_lambda
            #countratedensity = self.telescope_area * tp['Throughput'] * vegaflux / denom

            #countratedensityflux,countratedensitywave,pwave = self.getcountsflux(bp)
            #totalrate = self.integrate(countratedensity)

            #Is the variable below used at all?
            #grism_total_rate += totalrate

            obs = self.getcountsflux(bp)

            # From observation, get the flux in counts
            countratedensityflux = obs(obs.binset, flux_unit='count', area=self.telescope_area)

            print('countratedensityflux',countratedensityflux)

            # Multiply by dispersion
            countratedensityflux *= dispersion

            # Countrate density at the pivot wavelength
            print('pivot',pwave.value)
            print('countratedensityflux*dispersion',countratedensityflux)
            print('obs.binset',obs.binset)
            #cnorm = np.interp(pwave.value,obs.binset,countratedensityflux)
            cnorm = obs(pwave,flux_unit='count',area=self.telescope_area) * dispersion
            print('cnorm',cnorm)

            # Vega flux value at pivot wavelength, convert to Jansky
            vega_pivot = self.vega(pwave)
            j0 = units.convert_flux(pwave,vega_pivot,'MJy')

            # Now we need to divide by the area of a pixel in
            # sterradians, so we can eventually get MJy/str per ADU/sec
            j0 /= self.pixel_area_sr

            #vega_pivotorig = np.interp(pwave.value,self.vega.waveset,self.vega(self.vega.waveset))
            #print("units of vega flux are: {}".format(self.vega(self.vega.waveset).unit))
            #j0 = self.toJansky(vega_pivot.value,pwave.value) / 1.e6

            # Ratio of Vega flux to the countrate density at pivot wavelength
            ratio = j0 / cnorm

            print('')
            print('PHOTMJSR',ratio)
            print('NIRISS values are 0.01 to 2.0. I would expect ours to be similar!')
            print('')

            # Define a set of wavelengths to evaluate relative fluxcal
            goodpts = bp(bp.waveset) > 0.0001
            minwave = np.min(bp.waveset[goodpts])
            maxwave = np.max(bp.waveset[goodpts])
            w = minwave.value
            allwaves = np.array([w])
            while (w < maxwave.value):
                delt = w / (np.absolute(np.int(order))*self.resolving_power)
                allwaves = np.append(allwaves,w+delt)
                w += delt

            # Calculate Vega flux at each wavelength
            nelem = len(allwaves)
            allfluxes = self.vega(allwaves)
            alljansky = units.convert_flux(allwaves,allfluxes,'MJy')
            allcounts = np.interp(allwaves,obs.binset,countratedensityflux)
            #allfluxes = np.interp(allwaves,self.vega.waveset,self.vega(self.vega.waveset))
            #alljansky = self.toJansky(allfluxes,allwaves) / 1.e6
            # Normalize the Vega counts at all wavelengths by the value at the pivot
            # wavelength
            allratio = alljansky.value / allcounts / ratio / self.pixel_area_sr

            #print(np.min(allwaves),np.max(allwaves),allwaves[0],allwaves[-1])
            #f,a = plt.subplots()
            #a.plot(allwaves,allratio,color='red')
            #a.set_xlabel('Wavelength ({})'.format(bp.waveset.unit))
            #a.set_ylabel('Normalized Ratio (MJy/str per count/sec)')
            #a.set_title("{}, Order {}".format(filter.upper(),order))
            #f.savefig(os.path.split(file)[-1]+'_allratios.pdf')

            #f,a = plt.subplots()
            #a.plot(allwaves,alljansky,color='red')
            #a.set_xlabel('Wavelength ({})'.format(bp.waveset.unit))
            #a.set_ylabel('Vega Flux (MJy)')
            #a.set_title("{}, Order {}".format(filter.upper(),order))
            #f.savefig(os.path.split(file)[-1]+'_alljansky.pdf')

            #f,a = plt.subplots()
            #a.plot(allwaves,allcounts,color='red')
            #a.set_xlabel('Wavelength ({})'.format(bp.waveset.unit))
            #a.set_ylabel('Counts (count/sec)')
            #a.set_title("{}, Order {}".format(filter.upper(),order))
            #f.savefig(os.path.split(file)[-1]+'_allcounts.pdf')

            if np.min(allcounts) < 0:
                print('WARNING: counts<0 for {},{}'.format(filter.upper(),order))
                stop

            if '444' in filter:
                print("")
                print('444!!!!!!')
                print(allratio.value)
                print(len(allratio.value))
                print(type(allratio.value))
                #bad = allratio.value > 1e5
                bad = np.where(allratio.value > 1e5)
                print(len(bad[0]))
                print(type(bad[0]))
                print(alljansky.value[bad[0][0:10]])
                print(allcounts[bad[0][0:10]])
                print(ratio)
                print(self.str_per_detector)
                print(allratio.value[bad[0][0:10]])
                print(allwaves)
                stop


            # Pad allwaves and allratio to be the length needed by the
            # photom reference file. Also convert wavelengths to microns.
            allwaves = np.pad(allwaves/1.e4,(0,self.maxlen-nelem),'constant')
            allratio = np.pad(allratio,(0,self.maxlen-nelem),'constant')

            print(allwaves[100:105])
            print(allcounts[100:105])
            print(alljansky[100:105])
            print(alljansky[100:105]/allcounts[100:105])
            print(ratio)
            print(allratio[100:105])

            print("******************")
            print("******************")
            print("******************")
            print("need real conversion factor and uncertainty!!!")
            conversion_factor = 1000.
            uncertainty = ratio*.1
            #row = [filter,pupil,np.int(order),ratio*conversion_factor,uncertainty,nelem,allwaves,allratio]

            # Populate lists that will be used to create the final table
            filters = np.append(filters,filter)
            pupils = np.append(pupils,pupil)
            orders = np.append(orders,np.int(order))
            fluxes = np.append(fluxes,ratio*conversion_factor)
            uncs = np.append(uncs,uncertainty)
            nelems = np.append(nelems,nelem)

            print(allwaves.shape)

            if len(waves) == 0:
                waves = allwaves
                waves = np.expand_dims(waves,axis=0)
                resps = allratio
                resps = np.expand_dims(resps,axis=0)
            else:
                waves = np.append(waves,np.expand_dims(allwaves,axis=0),axis=0)
                resps = np.append(resps,np.expand_dims(allratio,axis=0),axis=0)


        print('waves.shape',waves.shape)
        print('resps.shape',resps.shape)

        print(filters)
        print(pupils)
        print(orders)
        print(fluxes)
        print(uncs)
        print(nelems)
        print(waves[0,40:45])
        print(resps[0,40:45])

        # Zip all data together
        alldata = np.array(zip(np.array(filters),np.array(pupils),np.array(fluxes),np.array(uncs),np.array(orders),np.array(nelems),np.array(waves),np.array(resps)),dtype=self.model.phot_table.dtype)

        return alldata


    def save_photom(self,tab,outfile):
        """
        Save reference file using astropy since the nircam phot model isn't finalized yet?
        and because we need astroconda for pysynphot but jwst is not in astroconda yet.
        OBSOLETE!!!!! Use save_photom_model below.
        """
        print('We need to feed in a header!!!!')
        h0 = fits.PrimaryHDU()
        h1 = fits.table_to_hdu(tab)
        h1.name = 'PHOTOM'
        hdulist = fits.HDUList([h0,h1])
        hdulist.writeto(outfile,clobber=True)
        hdulist.close()


    def save_photom_model(self,tab,outfile):
        """
        Save flux cal information using the NIRCam photom data model

        Arguments:
        ----------
        tab -- Astropy table contiaining flux cal information to be saved
        outfile -- Name of file in which to save flux cal info

        Returns:
        --------
        Fits file in photom data model format containing the JWST-
        reference file formatted flux calibration information
        """
        self.model.phot_table = tab

        print("Within model at time of saving:")
        print(self.model.phot_table[0][0])
        print(self.model.phot_table[0][1])
        print(self.model.phot_table[0][2])
        print(self.model.phot_table[0][3])
        print('')
        print(self.model.phot_table['filter'])
        print(self.model.phot_table['pupil'])
        print(self.model.phot_table['photmjsr'])
        print(self.model.phot_table['uncertainty'])

        # Set metadata
        self.model.meta.author = self.author
        self.model.meta.telescope = 'JWST'
        self.model.meta.instrument.name = 'NIRCAM'
        if '5' in self.detector:
            self.detector = self.detector[0] + 'LONG'
        self.model.meta.instrument.detector = 'NRC' + self.detector
        self.model.meta.instrument.module = self.detector[0].upper()

        self.model.meta.description = self.photom_descrip
        self.model.meta.filetype = 'PHOTOM'
        self.model.meta.reftype = 'PHOTOM'
        self.model.meta.pedigree = self.pedigree
        self.model.meta.useafter = self.useafter

        self.model.meta.photometry.pixelarea_arcsecsq = self.pixel_area_a2.value
        self.model.meta.photometry.pixelarea_steradians = self.pixel_area_sr.value

        if isinstance(self.model, jwst.datamodels.photom.NrcWfssPhotomModel):
            self.model.meta.exposure.type = 'NRC_WFSS'
            self.model.meta.exposure.p_exptype = 'NRC_WFSS|NRC_TSGRISM|NRC_GRISM|'
        elif isinstance(self.model, jwst.datamodels.photom.NrcImgPhotomModel):
            self.model.meta.exposure.type = 'NRC_IMAGE'
            self.model.meta.exposure.p_exptype = 'NRC_IMAGE|NRC_TSIMAGE|NRC_FLAT|NRC_CORON|NRC_TACONFIRM|NRC_TACQ|NRC_FOCUS|'

        if outfile is None:
            outfile = 'NIRCam_{}_photom.fits'.format(self.detector)
        self.model.save(outfile)


    def spec_phot(self,file):
        """OBSOLETE. Replaced by grism_cal"""
        #read in necessary file info
        with fits.open(file) as h:
            cname = h[0].header['COMPNAME']

        junk,filter,mod = cname.split('_')
        mod = mod[-1].upper()
        pupil = 'GRISMR'
        print('Eventually need to be smarter about pupil value!!!!!')

        print("opening {}".format(file))
        bp = S.FileBandpass(file)

        #convert to angstroms, to match the source spectrum
        bp.convert('Angstrom')

        #need to find the order here, so that the proper
        #dispersion value is used.
        #order = h[0].header['SOMETHING']
        ord = file.find('order')
        order = file[ord+5]
        dispersion = self.disp[order]

        #find pivot? mean? center? effective? wavelength
        #denom = self.h * self.c / eff_lambda
        #countratedensity = self.telescope_area * tp['Throughput'] * vegaflux / denom

        countratedensity,pwave = self.getcountsflux(file)
        #totalrate = self.integrate(countratedensity)

        #Is the variable below used at all?
        #grism_total_rate += totalrate

        countratedensityflux = countratedensity.flux
        countratedensityflux *= dispersion

        #countrate density at the pivot wavelength
        cnorm = np.interp(pwave,countratedensity.wave,countratedensityflux)

        #Vega flux value at pivot wavelength, convert to Jansky
        vega_pivot = np.interp(pwave,self.src_spectrum.wave,self.src_spectrum.flux)
        j0 = self.toJansky(vega_pivot,pwave)

        ratio = j0 / cnorm

        #Define a set of wavelengths to evaluate relative fluxcal
        goodpts = bp.throughput > 0.0001
        minwave = np.min(bp.wave[goodpts])
        maxwave = np.max(bp.wave[goodpts])
        deltawave = minwave / (np.absolute(np.int(order))*self.resolving_power)
        allwaves = np.arange(minwave+deltawave/2,maxwave+deltawave,deltawave)
        nelem = len(allwaves)

        allcounts = np.interp(allwaves,countratedensity.wave,countratedensityflux)
        allfluxes = np.interp(allwaves,self.src_spectrum.wave,self.src_spectrum.flux)
        alljansky = self.toJansky(allfluxes,allwaves)
        allratio = alljansky / allcounts / ratio

        #order will be needed eventually. At the moment there is no place for it in the NircamPhotomModel, which is why we'll use the NirissPhotomModel
        print("******************")
        print("******************")
        print("******************")
        print("need real conversion factor and uncertainty!!!")
        conversion_factor = 1000.
        uncertainty = ratio*.1
        newrow = [filter,pupil,np.int(order),ratio*conversion_factor,uncertainty,nelem,allwaves,allratio]
        return newrow


    def toJansky(self, fluxdensity, wavelength):
        """
        Convert flux density in flam units (erg s-1 cm-2 A-1)
        to Jansky (10^-23 erg s-1 cm-2 Hz-1)

        Arguments:
        ----------
        fluxdensity -- flux density value
        wavelength -- wavelength corresponding to the flux density

        Returns:
        --------
        Flux density converted to Jansky
        """
        c = 299792458.
        fnu = fluxdensity * wavelength**2 / c / 1e6
        return fnu * 1.e26


    def integrate(self, obs):
        """
        Integrate a synphot Observation object
        to get total countrate density

        Arguments:
        ----------
        obs -- synphot observation object

        Returns:
        --------
        Integrated countrate density
        """
        return simps(obs.flux, obs.wave)


    def getcountsflux(self, bandpass):
        """
        Calculate the total flux through a bandpass

        Arguments:
        ----------
        bandpass -- synphot SpectralElement object that contains the
                    bandpass to use

        Returns:
        A synphot Observation object that combines the bandpass
        with a Vega spectrum
        """
        obs = Observation(self.vega, bandpass)
        return obs


    def getcountsflux_kevin(self,something):
        """NOT USED"""
        #given an input flux calibrated spectrum and throughput curve, calculate flux in a bandpass
        self.src_spectrum

        #filter throughput curve is probably smoother than vega curve, so interpolate the filter
        #passband to match the vega wavelength points

        #multiply filter passband by the Vega spectrum
        aval = filter * spectrum

        #translate to photon flux
        photon= h * speedoflight * 1.e+06 / spectrumwl
        value = aperture * aval * spectrumfl / photon
        return value


    def save_photom_reffile_astropy(self, photomtab, descrip):
        """NOT USED"""
        #eliminates the need for jwst library modules
        #Currently we need astroconda for pysynphot, but jwst
        #is not in astroconda
        h = fits.PrimaryHDU()

        h[0].header['TELESCOP'] = 'JWST'
        h[0].header['INSTRUME'] = 'NIRCam'
        h[0].header['FILETYPE'] = 'PHOTOM'
        h[0].header['AUTHOR'] = 'Bryan Hilbert'
        h[0].header['DESCRIP'] = descrip
        h[0].header['PEDIGREE'] = pedigree
        h[0].header['USEAFTER'] = useafter
        h[0].header['EXPTYPE'] = 'ANY'
        h[0].header['SUBARRAY'] = 'GENERIC'

        h[0].header['PIXAR_SR'] = self.pixar_sr
        h[0].header['PIXAR_A2'] = self.pixar_a2

        #insert data
        #model.phot_table['filter'] = photomtab['FILTER']
        #model.phot_table['pupil'] = photomtab['PUPIL']
        #model.phot_table['photmjsr'] = photomtab['PHOTMJSR']
        #model.phot_table['uncertainty'] = 0.
        #model.phot_table['nelem'] = photomtab['NELEM']
        #model.phot_table['wavelength'] = 0.
        #model.phot_table['relresponse'] = 0.


    def save_photom_reffile(self,photomtab,descrip):
        """NOT USED"""
        #save an input table of photom results to a file
        #using the appropriate datamodel
        model = NircamPhotomModel()

        #necessary metadata
        model.meta.telescope = 'JWST'
        model.meta.instrument.name = 'NIRCam'
        model.meta.filetype = 'PHOTOM'
        model.meta.reffile.author = 'Bryan Hilbert'
        model.meta.reffile.description = descrip
        model.meta.reffile.pedigree = pedigree
        model.meta.reffile.useafter = useafter
        model.meta.exposure.type = 'ANY'
        model.meta.subarray.name = 'GENERIC'

        model.meta.photometry.pixelarea_steradians = self.pixar_sr
        model.meta.photometry.pixelarea_arcsecsq = self.pixar_a2

        #insert data
        model.phot_table['filter'] = photomtab['FILTER']
        model.phot_table['pupil'] = photomtab['PUPIL']
        model.phot_table['photmjsr'] = photomtab['PHOTMJSR']
        model.phot_table['uncertainty'] = 0.
        model.phot_table['nelem'] = photomtab['NELEM']
        model.phot_table['wavelength'] = 0.
        model.phot_table['relresponse'] = 0.
        model.history.append('String to add')
        model.history.append('another string to add')


    def synphot_calcs(self, bpfile):
        """
        Calculate zeropoint for a given bandpass in several
        photometric systems

        Arguments:
        ----------
        bpfile -- Text file containing the throuput table for
                  a single bandpass

        Returns:
        --------
        Bandpass zeropoint value in Vega mags, AB mags, ST mags,
        photflam, photfnu, and megajansky. Also returns pivot wavelength
        """
        # Define wavelength list to use
        wave_bins = np.arange(0.5, 5, 0.01) * u.micron

        # Use refactored synphot to calculate zeropoints
        orig_bp = SpectralElement.from_file(bpfile)

        # Now reduce the PCE curve by a factor equal to
        # the gain, in order to get the output to translate
        # from ADU/sec rather than e/sec
        bp = orig_bp / self.gain
        #bp.model.lookup_table = bp.model.lookup_table / self.gain

        photflam = bp.unit_response(self.telescope_area)
        photplam = bp.pivot()
        st_zpt = -2.5 * np.log10(photflam.value) - 21.1
        ab_zpt = (-2.5 * np.log10(photflam.value) - 21.1 - 5 * np.log10(photplam.value) + 18.6921)
        #mjy_zpt = units.convert_flux(photplam,photflam, 'MJy')
        mjy_zpt = photflam.to(u.MJy, u.spectral_density(photplam))
        obs = Observation(self.vega, bp, binset=wave_bins)
        vega_zpt = -obs.effstim(flux_unit='obmag', area=self.telescope_area)
        photfnu = units.convert_flux(photplam, photflam, units.FNU)
        return (vega_zpt.value, ab_zpt, st_zpt, photflam.value, photfnu.value,
                mjy_zpt.value, photplam.value)


    def pysynphot_calcs(self,bpfile):
        """NOT USED"""
        #Use pysynphot to calculate zeropoints

        #read in bandpass file
        bp = S.FileBandpass(bpfile)

        #convert to angstroms, to match the source spectrum
        bp.convert('Angstrom')

        #calculate the pivot wavlength, return to microns
        pivot = round(bp.pivot())/10000.

        #renormalize the source spectrum to 1 count/sec in the bandpass
        sp_zeropoint = self.src_spectrum.renorm(1,'counts',bp)

        #create observation
        obs_zeropoint = S.Observation(sp_zeropoint, bp)

        #verify that the countrate is 1 count/sec
        countrateFullAperture_zeropoint = obs_zeropoint.countrate()
        cr = np.absolute(1.-countrateFullAperture_zeropoint)
        if cr > 0.00015:
            print('WARNING: error in normalization of source spectrum')
            print('for file {}'.format(bpfile))
            print('expecting countrate of 1.0, found countrate of {}.'.
                  format(countrateFullAperture_zeropoint))
            sys.exit()

        #zeropoint = magnitude that gives 1 count per second
        vzp = obs_zeropoint.effstim('vegamag')
        azp = obs_zeropoint.effstim('abmag')
        szp = obs_zeropoint.effstim('stmag')
        flambda = obs_zeropoint.effstim('flam')
        fn = obs_zeropoint.effstim('fnu')
        jy = obs_zeropoint.effstim('Jy')

        return (vzp,azp,szp,flambda,fn,jy)


    def make_pam(self):
        """Create pixel area map reference file"""
        import polynomial #from Colin Cox

        # Lists of polynomial coefficients
        if ((self.xcoeffs is None) | (self.ycoeffs is None)):
            print("Getting coefficients from SIAF")
            self.xcoeffs, self.ycoeffs = self.get_coefficients()

        # Get the pixel area if it was not
        # previously defined
        if self.pixel_area_a2 is None:
            print('Pixel area not defined')
            print('calculating from SIAF.')
            self.pixel_area_a2 = self.get_pixel_area()

        if self.pixel_area_sr is None:
            self.pixel_area_sr = self.pixel_area_a2.to(u.steradian)

        # Get order
        order = self.get_order(self.xcoeffs)

        # Run Jacobian to get area
        coords = np.mgrid[-1024:1024,-1024:1024]
        ycoords = coords[0,:,:]
        xcoords = coords[1,:,:]
        map = polynomial.jacob(self.xcoeffs, self.ycoeffs, xcoords, ycoords, order)
        yd,xd = map.shape
        print('Pixel area from Jacob at {},{}: {}'.format(yd/2, xd/2, map[yd/2, xd/2]))
        print('Manually calculated reference location pixel area: {}'.format(self.pixel_area_a2))

        test_ratio = map[yd/2, xd/2] / self.pixel_area_a2.value #should be 1.0
        print('Ratio of areas at 1024,1024 (should be 1.0): {}'.format(test_ratio))

        map /= self.pixel_area_a2.value
        self.save_pam(map)


    def save_pam(self,pixmap):
        """
        Save pixel area map in correct reffile format

        Arguments:
        ----------
        pixmap -- 2D array containing the pixel area map

        Returns:
        --------
        Nothing (but saves pixmap in PixelAreaModel file)
        """
        from jwst.datamodels import PixelAreaModel

        detector = self.detector
        if '5' in self.detector:
            detector = self.detector.replace('5','LONG')

        model = PixelAreaModel(pixmap)

        model.meta.author = self.author
        model.meta.filetype = 'AREA'
        model.meta.instrument.name = 'NIRCAM'
        model.meta.instrument.detector = 'NRC'+detector
        model.meta.instrument.filter = 'ANY'
        model.meta.instrument.pupil = 'ANY'
        model.meta.description = self.pam_descrip
        model.meta.pedigree = self.pedigree
        model.meta.useafter = self.useafter
        model.meta.exposure.type = 'NRC_IMAGE'
        model.meta.reftype = 'AREA'
        model.meta.telescope = 'JWST'
        model.meta.photometry.pixelarea_arcsecsq = self.pixel_area_a2.value
        model.meta.photometry.pixelarea_steradians = self.pixel_area_sr.value
        #model.meta.subarray.fastaxis = 1
        #model.meta.subarray.slowaxis = -2
        #model.meta.subarray.name = 'FULL'
        #model.meta.subarray.xsize = 2048
        #model.meta.subarray.ysize = 2048
        #model.meta.subarray.xstart = 1
        #model.meta.subarray.ystart = 1
        # Now split up the pam_history entry into a list of 60-character
        # strings
        if self.pam_history is not None:
            stringlist = self.chop_string(self.pam_history, 60)
        else:
            print("WARNING: No HISTORY info for PAM file. This will")
            print("not be accepted when it comes time to deliver files.")

        for s in stringlist:
            model.history.append(s)

        if self.pam_outfile is None:
            self.pam_outfile = 'NIRCam_{}_PAM.fits'.format(detector)
        model.save(self.pam_outfile)


    def chop_string(self, full, length):
        """
        Chop a long string into a list of shorter strings
        of size length.

        Arguments:
        ----------
        full -- String to be chopped
        length -- Length of substrings into which the full string
                  is chopped

        Returns:
        --------
        List of substrings
        """
        fulllen = len(full)
        s = []
        i = 0
        while i < fulllen:
            s.append(full[i:i+length])
            i += length
        return s


    def get_order(self, clist):
        """
        Determine the order of the polynomial coefficients
        retrieved from SIAF

        Arguments:
        ----------
        clist -- List of coefficients

        Returns:
        --------
        The (integer) polynomial order of the coefficients
        """
        ord_dict = {}
        ord_dict['21'] = 5
        ord_dict['15'] = 4
        ord_dict['10'] = 3
        ord_dict['6'] = 2
        ord_dict['3'] = 1
        return ord_dict[str(len(clist))]


    def get_coefficients(self):
        """Retrieve polynomial coefficeints from SIAF"""
        xc = self.siaf.get_polynomial_coefficients()['Sci2IdlX']
        yc = self.siaf.get_polynomial_coefficients()['Sci2IdlY']
        return xc, yc


    def get_siaf_row(self):
        """Extract the row in SIAF corresponding to the
           aperture being used"""
        # Read in SIAF
        siaf = self.read_siaf()

        # Extract the appropriate row
        aperture = 'NRC{}_FULL'.format(self.detector)
        match = siaf['AperName'] == aperture
        if np.sum(match) < 1:
            print("No matching row found in SIAF for {}".format(aperture))
            sys.exit()
        if np.sum(match) > 1:
            print("WARNING: multiple rows in SIAF match {}".format(aperture))
            sys.exit()
        self.siaf_row = siaf[match]


    def read_siaf(self):
        """Read in SIAF"""
        return ascii.read(self.siaf,header_start=1,data_start=2)


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='make photom reffile')
        parser.add_argument("detector",help='Name of detector the throughput files correspond to')
        parser.add_argument("imaging_throughput_files",help='File containing list of imaging throughput files')
        parser.add_argument("--grism_throughput_files",help='File containing list of throughput files',default=None)
        parser.add_argument("--siaf",help="SIAF from which to get pixel area values")
        parser.add_argument("--gain_file",help="Gain reference file for the input detector")
        parser.add_argument("--author",help="Name of author for reffiles",default="Nircam Team")
        parser.add_argument('--pedigree',default='GROUND')
        parser.add_argument('--useafter',default='2017-01-01T00:00:00')
        parser.add_argument('--photom_descrip',help="Descrip for photom reffile",default='Flux calibration reffile')
        parser.add_argument('--photom_outfile',help="Name of output photom reference file",default=None)
        parser.add_argument('--photom_history',help="Text to be added to the HISTORY field of the photom file",default=None)
        parser.add_argument('--pam_descrip',help="Descrip for PAM reffile",default='Pix area map for NIRCam')
        parser.add_argument('--pam_outfile',help="Name of output pixel area map file",default=None)
        parser.add_argument('--pam_history',help="Text to be added to the HISTORY field",default=None)
        parser.add_argument('--pam_only',help='If True, only the PAM is calculated and the photom file is skipped',action='store_true')
        parser.add_argument('--photom_only',help='If True, only the photom file is calculated and the PAM is skipped.',action='store_true')
        return parser

if __name__ == '__main__':
    usagestring = 'USAGE: photom.py listfile'

    photom = Photom()
    parser = photom.add_options(usage = usagestring)
    args = parser.parse_args(namespace=photom)
    #photom.get_siaf_row()

    if args.pam_only:
        photom.make_pam()

    elif args.photom_only:
        photom.make_photom()

    else:
        photom.make_photom()
        photom.make_pam()

