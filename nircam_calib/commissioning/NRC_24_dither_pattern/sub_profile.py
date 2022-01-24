name = "sub_profile"
import inspect
import re
#
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse
import numpy as np
#
from astropy.io import fits
import astropy.io.fits as fits
from astropy.table import Table,Column, MaskedColumn
from astropy.stats import sigma_clipped_stats
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.nddata import NDData
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.visualization import simple_norm

from astropy.stats import sigma_clipped_stats
    #
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
#
import photutils
from photutils.detection import DAOStarFinder
from photutils.detection import IRAFStarFinder
from photutils.detection import find_peaks
#
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.psf import BasicPSFPhotometry
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
#
import sep
#
#------------------------------------------------------------------------
#
def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno
#
#
#------------------------------------------------------------------------
#
def scale_from_filter(filter):
    if(filter == 'F070W' or filter == 'F090W'  or filter == 'F115W' or
       filter == 'F140M' or filter == 'F150W2'  or filter == 'F150W' or
       filter == 'F162M' or filter == 'F164N'  or filter == 'F182M' or
       filter == 'F187N' or filter == 'F200W'  or filter == 'F210M' or
       filter == 'F212N'):
        pixel_scale = 0.0310
    else:
        pixel_scale = 0.0630
    return pixel_scale
#
#------------------------------------------------------------------------
#
def scale_from_filename(file):
    pixel_scale = 0.0310
    filter = 'null'
    # Mirage output 
    if(re.search('b5',file) != None or re.search(file,'a5') ) :
        pixel_scale = 0.0630
    
    # Webbpsf templates used by guitarra
    if(re.search('F277W', file) != None):
        pixel_scale = 0.0630
        filter = 'F277W'

    if(re.search('F250M', file) != None):
        pixel_scale = 0.0630
        filter = 'F250M'
    
    if(re.search('F300M', file) != None):
        pixel_scale = 0.0630
        filter = 'F300M'
    
    if(re.search('F322W2', file) != None):
        pixel_scale = 0.0630
        filter = 'F322W2'
    
    if(re.search('F323N', file) != None):
        pixel_scale = 0.0630
        filter = 'F323N'
    
    if(re.search('F335M', file) != None):
        pixel_scale = 0.0630
        filter = 'F335M'
    
    if(re.search('F356W', file) != None):
        pixel_scale = 0.0630
        filter = 'F356W'
    
    if(re.search('F360M', file) != None):
        pixel_scale = 0.0630
        filter = 'F360M'
    
    if(re.search('F405N', file) != None):
        pixel_scale = 0.0630
        filter = 'F405N'
    
    if(re.search('F410M', file) != None):
        pixel_scale = 0.0630
        filter = 'F410M'
    
    if(re.search('F430M', file) != None):
        pixel_scale = 0.0630
        filter = 'F430M'
    
    if(re.search('F444W', file) != None):
        pixel_scale = 0.0630
        filter = 'F444W'
    
    if(re.search('F460M', file) != None):
        pixel_scale = 0.0630
        filter = 'F460M'
    
    if(re.search('F466N', file) != None):
        pixel_scale = 0.0630
        filter = 'F466N'
    
    if(re.search('F470N', file) != None):
        pixel_scale = 0.0630
        filter = 'F470N'
    
    if(re.search('F480M', file) != None):
        pixel_scale = 0.0630
        filter = 'F480M'

    return pixel_scale, filter
#
#------------------------------------------------------------------------
#   push(@x, $radius);
#    push(@y, $flux);
#    push(@sb, $flux/$sarea);
#    push(@area,$sarea);
#    if($i == 0) {
#        push(@z, $flux);
#        $df = $flux;
#        $xx = $df/$sarea;
#    } else {
#        $df = $flux - $y[$i-1];
#        push(@z, $df);
#        $xx = $df/($sarea-$area[$i-1]);
#    }
#    push(@intensity, log10($xx));

def differential(apmag,area, norm):
    diff         = np.zeros(len(apmag))
    encircled    = np.zeros(len(apmag))
    diff[0]      = apmag[0]
    encircled[0] = apmag[0]/area[0]
#    ii = 0
#    print(apmag[ii], diff[ii], area[ii], encircled[ii])
    for ii in range(1,len(apmag)):
        diff[ii]= apmag[ii]-apmag[ii-1]
        encircled[ii] = diff[ii]/(area[ii]-area[ii-1])
#        print(apmag[ii], diff[ii], area[ii], encircled[ii])
    if(norm == True or norm == 1) :
        diff = diff/diff[0]
        encircled = encircled/encircled[0]
    return diff, encircled
#
#------------------------------------------------------------------------
#
#
def read_ascii(file, debug):
    with open(file) as f:
        length = len(f.readlines())
    f.close()
    radius = np.zeros(shape=[length])
    flux   = np.zeros(shape=[length])
    area   = np.zeros(shape=[length])
    
    data = open(file,"r")
    ii = -1
    for line in data:
        line = line.strip()
        columns = line.split()
        if (columns[0] != '#') :
            ii= ii +1
            radius[ii] = float(columns[0])
            flux[ii]   = float(columns[2])
            area[ii] = float(columns[4])
    data.close()
    radius = radius[0:ii]
    flux   = flux[0:ii]
    area = area[0:ii]
    return radius, flux, area
        
##
#------------------------------------------------------------------------
#
def read_iraf_isim_cv3(file, debug):
    with open(file) as f:
        length = len(f.readlines())
    f.close()
    radius = np.zeros(shape=[length])
    flux   = np.zeros(shape=[length])
    area   = np.zeros(shape=[length])
    data = open(file,"r")
    ii = -1
    for line in data:
        line = line.strip()
        columns = line.split()
        if (columns[0] != '#') :
            ii= ii +1
            radius[ii] = float(columns[0])
            flux[ii]   = float(columns[3])
            area[ii]   = float(columns[4])
    data.close()
    radius = radius[0:ii]
    flux   = flux[0:ii]
    return radius, flux, area
        
#
#------------------------------------------------------------------------
#
def read_image(file, debug):
    hdulist = fits.open(file)
    nhdu  = len(hdulist)
    if(debug == 1):
        hdulist.info()
    samp_rate = 1
    filter    = 'null'
    pixel_scale = 0.0
    for ii in range(0, nhdu):
#        print("at line ", lineno()," ii ", ii, " nhdu", nhdu)
        header = hdulist[ii].header
        if('FILTER' in header):
            filter = header['FILTER']
#            print("at line ", lineno()," filter ", filter)
        if('OVERSAMP' in header):
            samp_rate = int(header['OVERSAMP'])
        if('PIXELSCL' in header):
            pixel_scale = float(header['PIXELSCL'])
        #        if(re.search('.slp.fits',file) != None ):
        naxis  = int(header['NAXIS'])
        if(naxis > 0 ) :
            if(naxis == 2):
                image  = hdulist[ii].data
            if(naxis == 3):
                image  = hdulist[ii].data[0,:,:]
            break

        if('EXTNAME' in header):
            if(header['EXTNAME'] =='SCI' ):
                naxis  = int(header['NAXIS'])
                print("sub_profile at line ", lineno()," naxis ", naxis)
                if(naxis == 2) :
                    image  = hdulist[ii].data
                else:
                    image  = hdulist[ii].data[0,:,:]
            break
            
        if('XTENSION' in header):
            if(header['XTENSION'] == 'IMAGE'): 
                naxis  = int(header['NAXIS'])
                if(naxis == 2) :
                    image  = hdulist[ii].data
                else:
                    image  = hdulist[ii].data[0,:,:]
                break
#        print("at line ", lineno(), " ii ", ii)

    print("sub_profile line ",lineno()," samp_rate, filter, pixel_scale",samp_rate, filter, pixel_scale)
    return image, header, samp_rate, filter, pixel_scale
#
#------------------------------------------------------------------------
#
def read_webbpsf_template(file, ext, debug):
    hdulist = fits.open(file)
    header = hdulist[ext].header
    over_sampling_rate = int(header['OVERSAMP'])
    naxis  = int(header['naxis'])
    naxis1 = int(header['naxis1'])
    naxis2 = int(header['naxis2'])
    print('EXT', ext, ' naxis ', naxis, ' naxis1 ', naxis1,' naxis2 ', naxis2,
          'over_sampling_rate',over_sampling_rate,' file ',file)
    if(debug == 1) :
        hdulist.info()
    xc     = naxis1/2.0 - 0.5
    yc     = naxis2/2.0 - 0.5
    print("line ", lineno()," xc, yc ", xc, yc)
    (xc, yc) = run_sep(file,ext, 200.0, debug,plot_results=False)
    print("line ", lineno()," SEP xc, yc ", xc, yc)
    if(naxis == 2):
        image  = hdulist[ext].data
        error  = np.sqrt(np.abs(image))
    
    return image, xc, yc, over_sampling_rate, naxis1, naxis2
#
#-----------------------------------------------------------------------------

def read_webbpsf_ext(tempfile, ext, template_r, norm_profile, debug):
    template, xct, yct, temp_samp, naxis1, naxi2 = read_webbpsf_template(tempfile, ext, debug)
    tempmask = np.zeros(template.shape, dtype=bool)
    
    rskyt= [xct-11., xct-1.]
    #
    (apmagt, areat, difft) = aper_phot(template, tempmask, xct, yct, template_r, rskyt, debug)
    difft, encircledt = differential(apmagt, areat, norm_profile)
    # print(lineno()," radii ", radii)
    rt = []
    for index in range(0,len(template_r)):
        rt.append((template_r[index]/temp_samp))
        if(debug == 1):
            print(lineno(), index, rt[index], apmagt[index], areat[index], encircledt[index])

    return rt, apmagt, areat, difft, encircledt, temp_samp
#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------
#
def find_objects(image, fwhm, nsigma, debug, mask=None):
    mean, median, std_dev = sigma_clipped_stats(image, mask=mask, sigma=nsigma)
    threshold  = nsigma * std_dev
    if(debug == 1) :
        print("mean, median, std_dev ",mean, median, std_dev)
        print("nsigma, detection threshold ",nsigma, threshold)

    iraf_find = IRAFStarFinder(threshold, fwhm, minsep_fwhm=25.0,
                               sharplo=0.5, sharphi = 2.0,
                               roundlo=0.0, roundhi=0.05,
                               peakmax=1000,exclude_border=True,
                               brightest=100)
    sources = iraf_find(image-median, mask= mask)
    return sources
#
#------------------------------------------------------------------------
#
def find_centroid(image, fwhm, nsigma, debug):

    sources = find_objects(image, fwhm, nsigma, debug)
    if(sources == None):
        print("find_centroid line #", lineno(), " sources", sources)
        sources = first_run(image, fwhm, nsigma)
    print("find_centroid # sources: ",len(sources))
    xc    = sources['xcentroid']
    yc    = sources['ycentroid']
    flux  = sources['flux']
    brightest = -1
    flux_max  = -1.0
    for ii in range(0,len(flux)):
        if(flux[ii] > flux_max) :
            brightest = ii
            flux_max = flux[ii]
    print("find_centroid : xc, yc, flux", xc[brightest], yc[brightest], flux_max)
    return xc[brightest], yc[brightest]
#
#------------------------------------------------------------------------
#
def first_run(image, fwhm, nsigma, debug=0, mask=None):
    mean, median, std_dev = sigma_clipped_stats(image, mask=mask, sigma=nsigma)

#    daofind = DAOStarFinder(fwhm=fwhm, threshold=5.*std_dev) 
#    sources = daofind(image - median)
    threshold  = nsigma * std_dev
    print("first_run: mean, median, std_dev, threshold ",mean, median, std_dev, threshold)
    iraf_find = IRAFStarFinder(threshold, fwhm)
    sources = iraf_find(image-median)

    if(debug == 1):
        print("mean, median, std_dev ",mean, median, std_dev)
        print("detection threshold ",threshold)
        for col in sources.colnames:
            sources[col].info.format = "%.8g"
        print(sources)

        np.fwhm      = sources['fwhm']
        np.area      = np.log10(sources['npix'])
        np.flux      = sources['flux']
        np.mag      = sources['mag']
        np.roundness = sources['roundness']
        np.sharpness = sources['sharpness']

        fig, ax = plt.subplots()
        plt.plot(np.area, np.mag,'bo',markersize=1)
        ax.set_xlabel('Log(area)')
        ax.set_ylabel('mag')
        plt.show()
        plt.plot(np.roundness, np.sharpness,'bo',markersize=1)
        ax.set_xlabel('roundness')
        ax.set_ylabel('sharpness')
        plt.show()
    return sources
#
#------------------------------------------------------------------------
#
def run_sep(file, ext, nsigma, debug, plot_results=False): 
    hdulist = fits.open(file)
    if(debug == 1 or debug == True):
        print("run_sep : hdulist.info() ")
        hdulist.info()
    header  = hdulist[ext].header
    naxis   = int(header['NAXIS'])
    if(naxis == 2):
        data    = hdulist[ext].data
    if(naxis == 3):
        data    = hdulist[ext].data[0,:,:]
    data    = data.byteswap(inplace=True).newbyteorder()
# measure a spatially varying background on the image
    bkg = sep.Background(data)
    # bkg = sep.Background(data, mask=mask, bw=64, bh=64, fw=3, fh=3)

    # get a "global" mean and noise of the image background:
    #    if(debug == 1 or debug == True):'
    print("run_sep: background : mean ", bkg.globalback, " rms ", bkg.globalrms)

    # evaluate background as 2-d array, same size as original image
    bkg_image = bkg.back()
    # bkg_image = np.array(bkg) # equivalent to above
# evaluate the background noise as 2-d array, same size as original image
    bkg_rms = bkg.rms()

# subtract the background
    data_sub = data  - bkg

    objects = sep.extract(data_sub, nsigma, err=bkg.globalrms)
# how many objects were detected
    len(objects)

    radius = 2.5 
    flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
                                    radius, err=bkg.globalrms, gain=1.0)

    if(debug == 1) :
        print("object     x        y     flux      +/-")
        for i in range(len(objects)):
            print("{:d}     {:8.2f} {:8.2f}  {:f}  {:f}".format(i, objects['x'][i], objects['y'][i], flux[i], fluxerr[i]))

    brightest = -1
    flux_max  = -1.0
    for ii in range(0,len(flux)):
        if(flux[ii] > flux_max) :
            brightest = ii
            flux_max = flux[ii]
            xc       = objects['x'][ii]
            yc       = objects['y'][ii]

    print("run_sep: xc, yc, flux_max ", xc, yc, flux_max)
    if(plot_results == True):
        rcParams['figure.figsize'] = [10., 8.]
        fig, ax = plt.subplots()
        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       vmin=m-s, vmax=m+s, origin='lower')
    
        for i in range(len(objects)):
            e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                        width=6*objects['a'][i],
                        height=6*objects['b'][i],
                        angle=objects['theta'][i] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        plt.xlabel("X (pixels)")
        plt.ylabel("Y (pixels)")
        plt.title(file)
        plt.show()

    return xc, yc
#
#------------------------------------------------------------------------
#
def aper_phot(image, mask, xc, yc, radii, rsky, debug):

    positions = [(xc, yc)]

    # Define apertures 
    apertures = [CircularAperture(positions, r=r) for r in radii ]
    if(debug == 1):
#        print("line ", lineno()," apertures  : ", apertures)
        print("line ", lineno()," aper_phot: positions: " , positions)
        print("line ", lineno()," aper_phot: sky aperture ", rsky)
#        for rr in range(0,len(radii)):
#            print("line ", lineno(), " apertures[rr].r, apertures[rr].area :", apertures[rr].r, apertures[rr].area )

# Background, masking bad pixels 
    annulus_aperture = CircularAnnulus(positions, r_in=rsky[0], r_out=rsky[1])
    annulus_masks = annulus_aperture.to_mask(method='center')
    bkg_median    = []
    for anm in annulus_masks:
        annulus_data = anm.multiply(image)
        annulus_data_1d = annulus_data[anm.data > 0]
        if(debug == 1) :
            print("line ", lineno()," aper_phot: annulus_data_1d.shape ",annulus_data_1d.shape)

        # Remove NaNs, Infs
        annulus_data_1d = annulus_data_1d[np.isfinite(annulus_data_1d)]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip) 
        if(debug == 1) :
            print("line ", lineno(), " aper_phot: annulus_data_1d.shape ",annulus_data_1d.shape)
            print("line ", lineno(), " aper_phot: median sigclip" , median_sigclip)
    if(debug == 1) :
        print("line ", lineno(), " aper_phot: bkg_median ", bkg_median)
        
    phot_table = aperture_photometry(image, apertures, mask=mask)
#
    junk = []
    area_list = []
    n   = -1
    for index  in phot_table.colnames:
        if('aperture_sum' in index):
            n = n +1
            array = phot_table[index].data[0]
            flux  = array.tolist()
            bkg = apertures[n].area *  bkg_median[0]
            junk.append(flux-bkg)
            area_list.append(apertures[n].area)
    apflux = np.array(junk)
    area  = np.array(area_list)
    diff  = differential(apflux,area,True)
    return apflux, area, diff
#
#------------------------------------------------------------------------
#
