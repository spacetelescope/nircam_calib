#! /usr/bin/env python

import numpy as np
from astropy.io import fits
import os,sys,argparse
from jwst_lib.models import SuperBiasModel
from sigmacut import calcaverageclass


def sigmacut_meandev(pixels,mask=None):
    sigmacut = calcaverageclass()
    pixels = np.ravel(pixels)
    initmask = np.isnan(pixels)
    if mask != None:
        initmask = np.logical_or(initmask,np.ravel(mask))
    sigmacut.calcaverage_sigmacutloop(pixels,mask=initmask,Nsigma=3.0,verbose=0)
    return sigmacut.mean,sigmacut.stdev

def run(listfile,outfile,origfilelist,overwrite=True):
    #check for the existence of outfile and quit if it exists and overwrite is false
    if os.path.isfile(outfile):
        if overwrite == True:
            os.remove(outfile)
        else:
            print(("WARNING: {} exists and overwrite is set to False. Quitting.".format(outfile)))
            sys.exit()

    #if origfilelist is not provided, use listfile
    if origfilelist == None:
        origfilelist = listfile

    #read in file with list of files
    files = []
    with open(listfile) as f:
        for line in f:
            files.append(line.strip())

    print(("Files to use:",files))

    #read in individual superbiases
    superbiasdqdef = []
    for i,file in enumerate(files):
        #ind_sb = SuperBiasModel(file)
        #bias = ind_sb.data
        #print("bias shape: {}".format(bias.shape))
        
        #err = ind_sb.err
        #dq = ind_sb.dq
        #dqdef = ind_sb.dq_def

        with fits.open(file) as h:
            bias = h[1].data
            err = h[2].data
            dq = h[3].data
            dqdef = []

        if i == 0:
            yd,xd = bias.shape
            superbiases = np.zeros((len(files),yd,xd))
            supererrors = np.zeros((len(files),yd,xd))
            superbiasdq = np.zeros((len(files),yd,xd))
            

        superbiases[i,:,:] = bias
        supererrors[i,:,:] = err
        superbiasdq[i,:,:] = dq
        superbiasdqdef.append(dqdef)
        

    #combine into a single superbias
    meanbias, meanerrs, meandq = calc_avg_superbias(superbiases,supererrors)

    #for a full frame bias, explicitly set the reference pixel values to zero
    #if yd == 2048:
    #    superbias[0:4,:] = 0.
    #    superbias[2044:,:] = 0.
    #    superbias[:,0:4] = 0.
    #    superbias[:,2044:] = 0
    #    errors[0:4,:] = 0.
    #    errors[2044:,:] = 0.
    #    errors[:,0:4] = 0.
    #    errors[:,2044:] = 0

    #dq definition
    dq_def = []
    dqssb_bpm={'DO_NOT_USE':np.uint8(1),'UNRELIABLE_BIAS':np.uint8(2)}
    for bitname in dqssb_bpm:
        bitvalue = dqssb_bpm[bitname]
        bitnumber = int(np.log(bitvalue)/np.log(2))
        newrow = (bitnumber,bitvalue,bitname,'')
        dq_def.append(newrow)

    #save as a properly formatted CRDS file
    save_lin_file(meanbias,meanerrs,np.uint8(meandq),dq_def,outfile,origfilelist)

    #save the averaged superbias and uncertainties
    #hdu0 = fits.PrimaryHDU()
    #hdu1 = fits.ImageHDU(meanbias)
    #hdu2 = fits.ImageHDU(meanerrs)
    #hdulist = fits.HDUList([hdu0,hdu1,hdu2])
    #h1 = hdulist[1].header
    #h1['DATATYPE'] = 'superbias'
    #h1['HISTORY'] = 'Files used to create this superbais:'
    #for file in files:
    #    h1['HISTORY'] = file
    #hdulist[1].header = h1
    #hdulist[2].header['DATATYPE'] = 'uncertainties'
    #hdulist.writeto(outfile)


def save_lin_file(bias,err,dq,dq_def,outfile,inlist):
    '''write out a properly formatted CRDS reference file'''

    #create empty instance of superbias model
    f = SuperBiasModel()
    
    #populate data
    f.data = bias
    f.err = err
    f.dq = dq
    f.dq_def = dq_def

    #header values
    f = make_proper_header(f,inlist)

    #save the reference file
    f.save(outfile)


#def calc_avg_superbias(biases,errs):
#    #combine into a single superbias using pixel-by-pixel weighted average
#    superbias = np.average(biases,axis=0,weights=1./errs)
#    errors = np.std(biases,axis=0)/np.sqrt(biases.shape[0])
#    return superbias,errors

def make_proper_header(model,filelist):
    '''set the proper header keywords for this reffile type'''
    files = []
    with open(filelist) as f:
        for line in f:
            if len(line) > 2:
                files.append(line.strip())

    #read in headers from one of the input files
    with fits.open(files[0]) as h:
        header0 = h[0].header
        header1 = h[1].header

    model.meta.reffile.type = 'SUPERBIAS'

    model.meta.subarray.xsize = header0['SUBSTRT1']
    model.meta.subarray.ysize = header0['SUBSTRT2']
    model.meta.subarray.xstart = header0['SUBSIZE1']
    model.meta.subarray.ystart = header0['SUBSIZE2']
    model.meta.subarray.name = header0['SUBARRAY']
    model.meta.instrument.name = 'NIRCAM'
    detector = header0['DETECTOR']
    if detector == 'NRCA5':
        detector = 'NRCALONG'
    if detector == 'NRCB5':
        detector = 'NRCBLONG'
    model.meta.instrument.detector = detector


    try:
        model.meta.subarray.fastaxis = header0['FASTAXIS']
        model.meta.subarray.slowaxis = header0['SLOWAXIS']
    except KeyError:
        print('===============================================')
        print("FASTAXIS and SLOWAXIS header keywords not found in the input data.")
        print("Assuming they are in native (fitswriter) orientation, and adding the")
        print("native orientation values for those keywords to the static pixel mask.")
        print('===============================================')
        model.meta.subarray.fastaxis = 1
        model.meta.subarray.slowaxis = 2
            
    model.meta.reffile.author = 'Hilbert'
    model.meta.reffile.description = 'Superbias reffile from CV3 data'
    model.meta.reffile.pedigree = 'GROUND'
    model.meta.reffile.useafter = '2015-10-01'

    #HISTORY keyword
    model.history.append('Description of Reference File Creation')

    model.history.append('DOCUMENT:')
    model.history.append('JWST-STScI-TR-XXXX')

    model.history.append('SOFTWARE:')
    model.history.append('/ifs/jwst/wit/witserv/data4/nrc/')
    model.history.append('hilbert/superbias/cv3/B1/superbias_create_ngroups_parallel')
    model.history.append('ized.py and average_superbiases.py')
        
    #put the list of input files into the HISTORY keyword
    model.history.append('DATA USED:')
    for file in files:
        totlen = len(file)
        div = np.arange(0,totlen,60)
        for val in div:
            if totlen > (val+60):
                model.history.append(file[val:val+60])
            else:
                model.history.append(file[val:])
        
            model.history.append('DIFFERENCES:')
            model.history.append('N/A. No previous version.')
    return model


def calc_avg_superbias(biases,errs):
    #combine into a single superbias using pixel-by-pixel sigma clipped average
    nframe,ny,nx = biases.shape

    superbias = np.zeros((ny,nx))
    supererrs = np.zeros((ny,nx))
    superdq = np.zeros((ny,nx))
    for y in range(ny):
        for x in range(nx):
    #for y in xrange(500,505):
    #    for x in xrange(500,505):
            if ((x==0) & (y % 100 == 0)):
                print(('Working on row {}'.format(y)))
            mm,dd = sigmacut_meandev(biases[:,y,x])
            superbias[y,x] = mm

            #calculate errors
            #use rms of input superbiases
            supererrs[y,x] = dd

            #use errors from superbiases divided by the number of frames
            #mme,dde = sigmacut_meandev(errs[:,y,x])
            #supererrs[y,x] = mme / np.sqrt(nframe)

    #look for pixels which have a much larger uncertainty
    #than average. Flag as UNRELIABLE_BIAS
    meanerr,deverr = sigmacut_meandev(supererrs)

    #save for evaluation
    #h=fits.PrimaryHDU()
    #h0=fits.ImageHDU(superbias)
    #h1=fits.ImageHDU(supererrs)
    #hl = fits.HDUList([h,h0,h1])
    #hl.writeto('TEST_superbias_supererrs.fits',clobber=True)

    #for code testing
    ##sigs = [3,5,7,9,11]
    ##for sig in sigs:
    ##    bad = len(np.where(supererrs > (meanerr+sig*deverr))[0])
    ##    print('{} sigma = {} unreliable bias pixels.'.format(sig,bad))

    #bad = supererrs > (meanerr + 5.*deverr)
    #superdq[bad] = np.uint8(2)
    #print("Found {} UNRELIABLE BIAS pixels.".format(len(np.where(superdq == 2)[0])))

    return superbias,supererrs,superdq

#def combine_dqarrays(dqs):
#    #Inputs to this will probably be blank, so don't waste much time on this.
    


def create_dq_array(biases):
    #Looking at the input data, flag pixels that have UNRELIABLE_BIAS
    zd,yd,xd = biases.shape
    superbias = np.zeros((yd,xd))
    uncertainty = np.zeros((yd,xd))
    badmap = np.zeros((yd,xd))
    for y in range(yd):
        for x in range(xd):
            mn,dev = sigmacut_meandev(biases[:,y,x])
            superbias[y,x] = mn
            uncertainty[y,x] = dev
            maxdiff = np.max(biases[:,y,x]) - np.min(biases[:,y,x])
            if maxdiff > threshold:
                badmap[y,x] = np.uint8(1)

    #How to define what an unreliable bias is? A sigma-clipped stdev that is larger 
    #than the average sigma-clipped stdev by some amount? Or larger than some number
    #of ADU? A pixel with one of the inputs from which the superbias is created with a 
    #value more than N-sigma from the mean?
    
    




def add_options(parser=None,usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage,description='Average together a list of biases into a superbias')
    parser.add_argument("listfile",help="file containing a list of bias files")
    parser.add_argument("outfile",help="name of file to output superbias into")
    parser.add_argument("--origfilelist",help='name of file containing the list of original files used.',default=None)
    parser.add_argument('-o','--overwrite',help='overwrite existing output file',action="store_true")
    return parser



if __name__ == '__main__':
    usagestring = 'average_biases.py files.list output.fits'
    parser = add_options(usage=usagestring)
    args = parser.parse_args()

    run(args.listfile,args.outfile,args.origfilelist,overwrite=args.overwrite)

