#! /usr/bin/env python

'''Create a CRDS-formatted dark current reference file from a list of
individual dark current integrations

HOT PIXEL SEARCH HAS NOT BEEN TESTED. IT WAS TAKEN FROM THE BAD PIXEL 
GENERATOR SCRIPT WHERE IT WAS TESTED AND WORKED.
HAS NOT BEEN TESTED ON FILES WHICH CONTAIN >1 INTEGRATION. SPECIFICALLY THE 
INTERACTION BETWEEN THE CR MAPS AND THE HOT PIXEL SEARCH.

'''


import sys, os
from astropy.io import fits
from astropy.table import Table,Column
import numpy as np
import argparse
import datetime
from itertools import izip
import subprocess
import matplotlib.pyplot as plt

# put the tools directory into the path
# add pythonmodules to PATH
#if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
#    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
#else:
#    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
#    sys.exit(0)

from sigmacut import calcaverageclass
from jwst.datamodels import DarkModel,MaskModel,dqflags


class make_nrc_dark_reffile():
    def __init__(self):
        self.verbose = False


    def read_listfile(self,file):
        '''Read in a listfile and return list'''
        files = []
        with open(file) as f:
            for line in f:
                if len(line) > 2:
                    files.append(line.strip())
        return files

    def input_consistency_check(self,files):
        '''check basic input file attributes to make sure everything
        is consistent'''
        readpatt = set()
        ngroups = set()
        nx = set()
        ny = set()
        pupil = set()
        detector = set()
        instrument = set()

        for file in files:
            with fits.open(file) as h:
                head0 = h[0].header

            pupil.add(head0['PUPIL'].strip())
            instrument.add(head0['INSTRUME'].strip())
            detector.add(head0['DETECTOR'].strip())
            readpatt.add(head0['READPATT'].strip())
            ngroups.add(head0['NGROUPS'])
            nx.add(head0['SUBSIZE1'])
            ny.add(head0['SUBSIZE2'])

        if len(pupil) > 1:
            print("WARNING!! MORE THAN ONE PUPIL SETTING USED FOR INPUT DATA!!")
            sys.exit(0)

        if len(instrument) > 1:
            print("WARNING!! MORE THAN ONE INSTRUMENT USED FOR INPUT DATA!!")
            sys.exit(0)

        if len(detector) > 1:
            print("WARNING!! MORE THAN ONE DETECTOR USED FOR INPUT DATA!!")
            sys.exit(0)

        if len(readpatt) > 1:
            print("WARNING!! MORE THAN ONE READ PATTERN USED FOR INPUT DATA!!")
            sys.exit(0)

        #dark current reference files must be RAPID. When the reference file is 
        #applied to data with a different readpattern, the groups in the reference
        #file are averaged and dropped in order to "create" a reference file that
        #matches the readpattern of the data being corrected
        if list(readpatt)[0] != 'RAPID':
            print("WARNING!! Input data do not use the RAPID read pattern! Quitting")
            sys.exit(0)

        if len(ngroups) > 1:
            print("WARNING!! MORE THAN ONE NUMBER OF GROUPS IN INPUT DATA!!")
            sys.exit(0)

        if len(nx) > 1:
            print("WARNING!! MORE THAN ONE ARRAY X-SIZE IN INPUT DATA!!")
            sys.exit(0)

        if len(ny) > 1:
            print("WARNING!! MORE THAN ONE ARRAY Y-SIZE IN INPUT DATA!!")
            sys.exit(0)
            
        return list(detector)[0],list(readpatt)[0],list(ngroups)[0],list(nx)[0],list(ny)[0]

    def hot_search(self,file,hot_sigfrac=10,in_mask=0):
        #mean and noise values used to determine the hot pixel threshold.
        #Better to use a consistent threshold ramp-to-ramp, rather than
        #relying on calculated mean and stdev for each.
        mnvals = {'SW':0.,'LW':20}
        devvals = {'SW':8.,'LW':8}

        slash = file.rfind('/')
        if slash != -1:
            filepath = file[0:slash+1]
            filefile = file[slash+1:]
        else:
            filepath = './'
            filefile = file

        #search for hot pixels in the input file.
        h = fits.open(file)
        data = h[1].data
        tint = h[0].header['NGROUPS'] * h[0].header['TGROUP']
        channel = h[0].header['CHANNEL'].strip()

        if channel == 'SHORT':
            mn = mnvals['SW'] / tint
            dev = devvals['SW'] / tint
        if channel == 'LONG':
            mn = mnvals['LW'] / tint
            dev = devvals['LW'] / tint
        if (channel != 'SHORT') and (channel != 'LONG'):
            print('Unable to use the channel header keyword to determine which channel the data are from. Quitting.')
            sys.exit()

        #increase dev to account for the use of the difference between two frames
        dev = dev * np.sqrt(2.)

        #assume 4 rows and cols of refpix which we will ignore
        #qstart = [4,512,1024,1536,2044]
        yd = data.shape[2]
        xd = data.shape[3]

        #use the mask to set CR-hit pixels to NaN
        if len(in_mask.shape) == 3:
            for integ in xrange(data.shape[0]):
                intmask = in_mask[integ,:,:]
                for y in xrange(yd):
                    for x in xrange(xd):
                        if intmask[y,x] != -1:
                            data[integ,intmask[y,x]:,y,x] = np.nan


        ngroup = data.shape[1]
        #readtimes = np.arange(ngroup) * tgroup #seconds, assuming full-frame
        hotmap = np.zeros((data.shape[0],yd,xd))
        for integration in xrange(data.shape[0]):
            #we want a good measure of the signal accumulated up the ramp. Under the assumption
            #that darks are being used in the search, we can use last-first read.
            diff = (data[integration,-1,:,:] - data[integration,0,:,:]) / tint
            
            #histograms
            allfig = plt.figure()
            hotfig = plt.figure()
            allax = allfig.add_subplot(1,1,1)
            #allaxright = allfig.add_subplot(111, sharex=allax, frameon=False)
            hotax = hotfig.add_subplot(1,1,1)

            #Use a single threshold across the detector. Doesn't make sense to have
            #a separate threshold for each quadrant.

            #if a mask is given, apply it. Don't check CR-hit pixels to see if they are hot.
            #if np.max(in_mask) > -1:
            #    bad = in_mask != -1
            #    diff[bad] = np.nan

            #calculate the mean of the science pixels
            #mn,dev = self.sigma_clipped_mean(diff[4:yd-4,4:xd-4])
            #RELY ON PAST MEASUREMENTS. THIS WAY WE ARE CONSISTENT FROM RAMP TO RAMP

            #determine threshold and locate hot pixels
            hotpix_threshold = mn + hot_sigfrac*dev
            #print('Hot pix search for {}, mean and stdev of data are {} and {}.'.format(file,mn,dev)) 
            #print('Hot pixel threshold is {}'.format(hotpix_threshold))
            hot = diff > hotpix_threshold 
            numhot = len(np.where(hot == True)[0])
            print('Number of hot pixels found: {}'.format(numhot))

            #histograms - put in units of DN/sec
            if numhot > 2:
                allbins = np.arange(mn-3*dev,hotpix_threshold+dev,dev*0.9)
                alln, allbs, allpatches = allax.hist(np.ravel(diff),bins=allbins,facecolor='blue',label='Good Pix') #,range=(-100,1000))
                #allrightbins = np.arange(mn-3*dev,mn+20*dev,dev*0.9)
                #allrightn,allrightbs,allrightpatches = allaxright.hist(np.ravel(diff),bins=allrightbins,normed=1,cumulative=True,histtype='step',color='green')
                allax.plot([mn,mn],[0,np.max(alln)],color='black',linestyle='--')
                allax.plot([hotpix_threshold,hotpix_threshold],[0,np.max(alln)],color='black',linestyle='-.')
                threshstr =  "{:.3f}".format(hotpix_threshold)
                allax.text(hotpix_threshold*1.05,0.1*np.max(alln),'Threshold = '+threshstr+' DN/sec',color='black',fontsize=10)
                hotbins = np.arange(hotpix_threshold,mn+20*dev,dev*0.9)
                hotn, hotbs, hotpatches = allax.hist(np.ravel(diff[hot]),bins=hotbins,facecolor='red',label='Hot Pix') #,range=(-100,1000))    

                hotzbins = np.arange(hotpix_threshold,np.nanmax(diff),(np.nanmax(diff)-hotpix_threshold)/20)
                hotzn, hotzbins, hotzpatches = hotax.hist(np.ravel(diff[hot]),bins=hotzbins,facecolor='red',label='Hot Pix') #,range=(-100,1000))


            #print('all pix: bins min, max, and delta: ',mn-3*dev,hotpix_threshold+dev,dev*0.9)
            #print('hot pix: bins min, max, and delta: ',hotpix_threshold,mn+20*dev,dev*0.9)
            #print('x and y of text: ',hotpix_threshold+3,0.1*np.max(alln))

            #flag the hot pixels in the hot pixel map
            hotmap[integration,:,:] = hot

            
            ##Calculate the mean signal in each quadrant
            #quadcolors = ['red','blue','green','magenta']
            #for i in xrange(len(qstart)-1):
            #    qdata = diff[4:yd-4,qstart[i]:qstart[i+1]]
            #
            #    if np.max(in_mask) > -1:
            #        qmask = in_mask[4:yd-4,qstart[i]:qstart[i+1]]
            #        bad = qmask != -1
            #        qdata[bad] = np.nan
            #    #mn,dev = self.sigma_clipped_mean(qdata)
            #    mn,dev = self.multi_step_sigma_clipped_mean(qdata)
            #
            #    #hot = qdata > np.absolute(mn*hot_sigfrac)
            #    hotpix_threshold = hot_sigfrac*24 #ADU. 24 is from adding ktc noise (23 ADU) and readnoise (7ADU)
            #    hotpix_threshold = mn + hot_sigfrac*dev
            #    print('hot pix search, mean and stdev of data ',mn,dev) 
            #    hot = qdata > hotpix_threshold 
            #    print("or should we define hotpix level based on stdev of input data?")
            #
            #    alln, allbins, allpatches = allax.hist(np.ravel(qdata),bins = np.arange(-24,24,0.1),facecolor=quadcolors[i],alpha=0.25,label='Amp '+str(i+1))
            #    allax.plot([mn,mn],[0,np.max(alln)],color='black',linestyle='--')
            #    allax.plot([hotpix_threshold,hotpix_threshold],[0,np.max(alln)],color='black',linestyle='-.')
            #    hotn, hotbins, hotpatches = hotax.hist(np.ravel(qdata[hot]),facecolor=quadcolors[i],alpha=0.25,label='Amp '+str(i+1))
            #
            #    #hot = qdata > (mn+hot_nsig*dev)
            #    hotmap[integration,4:yd-4,qstart[i]:qstart[i+1]] = hot
                
            #Plot a histogram of the pixel values from the first file, along with the limits of what
            #defines a hot pixel. 
            #Plot a histogram of only hot pixel values? There will be some odd entries because of pixels
            #that aren't hot in every integration.
            allax.set_yscale('log')
            hotax.set_yscale('log')
            allax.set_ylim(bottom=1)
            hotax.set_ylim(bottom=1)
            allax.set_ylabel('Number of Pixels')
            allax.set_xlabel('Signal Rate from (Last Group - First Group) (DN/sec)')
            allax.legend(loc='upper right')

            #allaxright.yaxis.tick_right()
            #allaxright.yaxis.set_label_position("right")
            #allaxright.set_ylabel("Cumulative Percentage of Pixels")

            hotax.set_ylabel('Number of Pixels')
            hotax.set_xlabel('Signal Rate from (Last Group - Frist Group) (DN/sec)')
            hotax.legend(loc='upper right')
            slash = file.rfind('/')
            #dir = file[0:slash+1]
            allhname = filepath+'Allpix_hist_'+file[slash+1:]+'_Integ'+str(integration)+'.pdf'
            hothname = filepath+'Hotpix_hist_'+file[slash+1:]+'_Integ'+str(integration)+'.pdf'
            allfig.savefig(allhname)
            hotfig.savefig(hothname)
            plt.close(allfig)
            plt.close(hotfig)
            print("All pixel and hot pixel histograms for {} saved to {} and {}.".format(file,allhname,hothname))


        #save the hotpixel map
        dot = file.rfind('.')
        slash = file.rfind('/')
        hotmapfile = filepath + file[slash+1:dot] + '_HOTMAP.fits'
        h[1].data = hotmap
        h.writeto(hotmapfile,clobber=True)
        print("Hotmap for {} written to {}.".format(file,hotmapfile))

        return hotmap


    def update_cr_map(self,cr_map,groupdq,groupstart):
        '''update the Cosmic Ray map for a given file'''
        nint = groupdq.shape[0]

        #loop over integrations
        for integration in xrange(nint):
            #extract the appropriate cr map for the integration 
            cr_int_map = cr_map[integration,:,:]

            #loop over the groups in the groupdq array, which contains
            #dq values for some set of groups within the integration
            for group in xrange(groupdq.shape[1]):
                
                #get appropriate groupdq group
                gdq = groupdq[integration,group,:,:]

                #find CR hits
                hits = (gdq & dqflags.pixel['JUMP_DET'] > 0)

                #for pixels with a CR hit in this group, we need to see if they
                #suffered a CR hit in any previous group. If so, keep the lower value.
                #If not, insert the current group number
                oldvals = cr_int_map[hits]
                newvals = np.zeros(oldvals.shape) + (group+groupstart)

                #keep the minimum value in each pixel
                minvals = np.minimum(oldvals,newvals)

                #set the CR map values to the minimum values
                #and place these back in the main mask
                cr_int_map[hits] = minvals
                cr_map[integration,group,:,:] = cr_int_map

        return cr_map

    def convert_dq_arr_to_cr_map(self,file):
        #read in dq array from a given file, extract cr hits from groupdq
        #and return a 2D array listing group number of initial CR hit for each pixel
        with fits.open(file) as h:
            groupdq = h[3].data

        #variable to hold CR map. Be ready for files with more than one integration
        map = np.zeros((groupdq.shape[0],groupdq.shape[2],groupdq.shape[3]))
        
        #variable that will hold the number of CR hits in each pixel
        numhits = np.zeros_like(map)

        #extract only CR bits at the moment
        crmask = (groupdq & dqflags.pixel['JUMP_DET'] > 0).astype(int)

        #loop over integrations
        for integration in xrange(groupdq.shape[0]):
            #loop over pixels
            for y in xrange(groupdq.shape[2]):
                for x in xrange(groupdq.shape[3]):
                    crvals = crmask[integration,:,y,x]
                    if max(crvals) > 0:
                        hits = np.where(crvals > 0)[0]
                        map[integration,y,x] = hits[0]

                    numhits[integration,y,x] = np.sum(crvals)

        return map,numhits


    def combine_ind_maps(self,maps,minfrac=0.4):
        #combine the hotpixel maps for all of the input files into a final hot pixel map
        #If a pixel is flagged as hot in more than minfrac fraction of the inputs, then the pixel 
        #is considered hot in the final map.
        nmap,yd,xd = maps.shape
        mapfrac = np.sum(maps,axis=0) / nmap
        
        finalmap = np.zeros((yd,xd))
        finalmap = mapfrac >= minfrac
        finalmap = finalmap.astype(int)
        return finalmap

    def table_of_badpix(self,array):
        #given a hot pixel table/map of some sort, create a table that lists the number
        #of hot pixels in each quadrant. Return a table.
        bad_stats = Table()
        #for j in xrange(4):
        #quad = array[:,qstart[j]:qstart[j+1]]
        #header = 'Quad'+str(j+1)
        header = 'Entire_Detector'
        bad_stats[header] = [len(np.where(array != 0)[0])]
        return bad_stats


    def run(self,infile,outfile,maxopen):
        print("USING np mean rather than sigma clipping!!!!")
        print("see below in code")                        

        #if the input filename includes a path, separate that out
        slash = infile.rfind('/')
        if slash != -1:
            infilepath = infile[0:slash+1]
            infilefile = infile[slash:]
        else:
            infilepath = './'
            infilefile = infile

        #check for the DQ array to use in the final file
        #if self.read_dq == None:
        #    print("FIXME! Need to find DQ values manually (import badpix script)")
        #    sys.exit(0)
            
        #read in list of dark current integrations
        files = self.read_listfile(infile)

        #check to be sure all input files have the same readpattern, array
        #size, etc
        detector,readpatt,ngroups,nx,ny = self.input_consistency_check(files)
        #print(detector,readpatt,ngroups,nx,ny)

        #get CR hit maps for all files, so we know which pixels to avoid
        #Use the DQ arrays from the input files
        if ((self.load_cr == False) & (self.empty_cr == False)):
            print("Locating CR hits using the input files' DQ arrays")
            full_cr_maps = np.zeros((1,ny,nx))
            full_cr_nhits = np.zeros((1,ny,nx))

            for file in files:
                cr_maps, cr_num = self.convert_dq_arr_to_cr_map(file)
                full_cr_maps = np.vstack((full_cr_maps,cr_maps))
                full_cr_nhits = np.vstack((full_cr_nhits,cr_num))

            #remove empty initial frame of each
            full_cr_maps = full_cr_maps[1:,:,:]
            full_cr_nhits = full_cr_nhits[1:,:,:]

            print("Done creating cosmic ray impact maps.")
            print("Saved to CR_hit_maps.fits")
            hh=fits.PrimaryHDU(full_cr_maps)
            hh1=fits.ImageHDU(full_cr_nhits)
            hl = fits.HDUList([hh,hh1])
            hl.writeto('CR_hit_maps.fits',clobber=True)

        elif ((self.load_cr == True) & (self.empty_cr == False)):
            print("Reading in CR maps from a previous run.")
            print("Hardwired to CR_hit_maps.fits")
            print("--------------------------------------")
            print("Still hardwired!! FIXME!!")
            with fits.open('CR_hit_maps.fits') as h:
                full_cr_maps = h[0].data
                full_cr_nhits = h[1].data
        elif self.empty_cr == True:
            print("Ignoring CR hits when calcualting means.")
            full_cr_maps = np.zeros((1,ny,nx)) + 999


        #find hot pixels for the DQ array if requested
        #For this, we read in one entire file at a time
        if ((self.read_dq == None) & (self.empty_dq == False)):
            print("read_dq was not set, so manually finding hot pixels for the DQ array")

            hotmaps = np.zeros((1,ny,nx))
            ctr = 0
            for i,darkfile in enumerate(files):
                numints = fits.getval(darkfile,'NINTS')
                hmap = self.hot_search(darkfile,hot_sigfrac=self.hot_sigfrac,in_mask = full_cr_maps[ctr:ctr+numints,:,:])
                hotmaps = np.concatenate((hotmaps,hmap),axis=0)
                ctr += numints
            hotmaps = hotmaps[1:,:,:]
        
            #combine the hot maps from all the files into one final hot map
            hot_map = self.combine_ind_maps(hotmaps,minfrac = self.min_hot_frac)

            #set the pixel values according to expected DQ values
            hot_map[hot_map != 0] = dqflags.pixel['HOT']
            
            #create a table of how many pixels in each quadrant are hot
            stats_hot = self.table_of_badpix(hot_map)
            source = Column(['Hot'],name='Source')
            stats_hot.add_column(source,index=0)

            #save table
            hot_table_file = infilepath+'Hot_stats_table_from_'+infilefile+'.txt'
            stats_hot.write(hot_table_file,format='ascii')
            print("Hot pixel statistics saved to {}".format(hot_table_file))

            #save the final hot map
            h = fits.PrimaryHDU()
            hi = fits.ImageHDU(hot_map)
            hdu = fits.HDUList([h,hi])
            if self.hotfile == None:
                self.hotfile = infilepath+'Hot_Pixel_Mask_from_'+infilefile+'.fits'
            hdu.writeto(self.hotfile,clobber=True)
            print('Hot pixel mask saved as {}'.format(self.hotfile))


        #output variables 
        meandark = np.zeros((ngroups,ny,nx)) - 999.
        meandarkerr = np.zeros((ngroups,ny,nx)) - 999.

        #index numbers of groups to be opened
        deltagroup = self.maxopen / len(files)
        groupstart = np.arange(0,ngroups,deltagroup)
        groupend = np.arange(deltagroup,ngroups,deltagroup)
        groupend = np.append(groupend,ngroups)

        #maps to hold cosmic ray information for all integrations
        #cr_maps = np.zeros((len(files),

        #instance of sigma clipping mean
        sigmacut = calcaverageclass()

        #total cosmic ray map. One frame per integration.
        #add frames as we loop over files
        #full_cr_map = np.zeros((1,ny,nx))
        filexts = {}

        #loop over groups of groups
        for gs,ge in izip(groupstart,groupend):
            
            #set up variables to hold data from groups n to n+deltagroup
            gnum = ge - gs
            chunk = np.zeros((1,gnum,ny,nx))
            chunk_err = np.zeros((1,gnum,ny,nx))

            #loop over files
            for file in files:
            
                #open with astropy so that we can read in only part of the file
                with fits.open(file) as h:
                    data = h[1].data[:,gs:ge,:,:]
                    err = h[4].data[:,gs:ge,:,:]
                    #pixeldq = h[2].data
                    #groupdq = h[3].data[:,gs:ge,:,:]
                    nint = h[0].header['NINTS']
                    ny = h[0].header['SUBSIZE2']
                    nx = h[0].header['SUBSIZE1']

                #throw out any data with less than 4 dimensions (i.e. we want to use
                #SSB-formatted data)
                if len(data.shape) < 4:
                    print("WARNING. INPUT DATA SHAPE HAS LESS THAN 4 DIMENSIONS.")
                    sys.exit(0)

                #find groups that have been impacted by a cosmic ray and set
                #them so that they are exluded from the averaging
                #The CR flag only appears in the group where the CR hit,
                #but we also want to exlude all groups after the CR hit,
                #so we need to keep a running CR map
                chunk = np.vstack((chunk,data))
                chunk_err = np.vstack((chunk_err,err))


            #remove the 0th extension of the chunk variables
            chunk = chunk[1:,:,:,:]
            chunk_err = chunk_err[1:,:,:,:]
            
            #calculate the mean for each pixel in each group
            for group in xrange(gnum):
            #for group in xrange(1):
                for y in xrange(ny):
                    for x in xrange(nx):
                        #for y in xrange(300,350):
                        #for x in xrange(300,350):
                        if ((y==x) & (y % 100 == 0) & (self.verbose == True)):
                            print("Group {}, Pixel ({},{})".format(group+gs,x,y))
                        pix = chunk[:,group,y,x]
                        pixerr = chunk_err[:,group,y,x]

                        crgroups = full_cr_maps[:,y,x]

                        bad = crgroups < group
                        good = crgroups >= group
                        numbad = len(np.where(bad == True)[0])
                        percbad = numbad / len(full_cr_maps[:,y,x])

                        #if the pixel shows a CR hit in > 50% of the input ramps,
                        #then the pixel must be weird, as opposed to there being 
                        #real CR hits. In this case, keep all of the input values when
                        #calculating the mean. If you don't do this, then there will be pixels 
                        #flagged as CR hit in every integration which will lead to no data
                        #being put in the mean calculator, and a 0 or NaN for that pixel in the 
                        #resulting mean dark current ramp.
                        if percbad > 0.5:
                            bad = bad * 0
                            good = good * 1

                        #if ((y == 100) & (x == 100)):
                        #    print(pix)
                        #    print(crgroups)
                        #    print(bad)

                        #calculate mean dark
                        #sigmacut.calcaverage_sigmacutloop(pix,mask=bad,Nsigma=3.0,verbose=0)
                        #if sigmacut.converged:
                        #    meandark[gs+group,y,x] = sigmacut.mean
                        #    meandarkerr[gs+group,y,x] = sigmacut.stdev

                        #do we trust the numpy mean to do it? 
                        #assuming we have done a good job filtering out CR hits, maybe 
                        #we can just use numpy, which is much faster than the sigma clipping
                        #above.
                        meandark[gs+group,y,x] = np.mean(pix[good])
                        meandarkerr[gs+group,y,x] = np.std(pix[good])


        #save as a reference file
        #if no hot pixel map is to be read in and it was not calculated, pass in 0
        if ((self.read_dq != None) | (self.empty_dq == True)):
            hot_map = 0
            dq_def = None
        self.save_reffile(meandark,meandarkerr,hot_map,files)

        #file format checks
        subprocess.call(['fitsverify',self.outfile])


    def save_reffile(self,dark,err,passedmap,files):
        '''save the mean dark current in CDBS format'''
        finaldark = DarkModel()
        finaldark.data = dark
        finaldark.err = err

        ngroup,ny,nx = dark.shape

        #if an empty dq array is requested, do that 
        if self.empty_dq == True:
            print("Using empty DQ array in output reference file!!")

        #read in the appropriate DQ mask if read_dq is set
        if self.read_dq != None:
            #looks like the MaskModel bug is still present. Read in with astropy
            #dq = MaskModel(self.read_dq)
            #dqmap = dq.dq
            #dqdef = dq.dq_def
            with fits.open(self.read_dq) as h:
                dqmap = h[1].data
                dqdef = h[2].data
        else:
            dqdef = []
            dqssb_dark={'DO_NOT_USE':np.uint8(1),'HOT':np.uint8(2),'WARM':np.uint8(4),'UNRELIABLE_DARK':np.uint8(8),'UNRELIABLE_SLOPE':np.uint8(16)}
            for bitname in dqssb_dark:
                bitvalue = dqssb_dark[bitname]
                bitnumber = int(np.log(bitvalue)/np.log(2))
                newrow = (bitnumber,bitvalue,bitname,'')
                dqdef.append(newrow)
            if self.empty_dq == True:
                dqmap = np.zeros((ny,nx))
            else:
                dqmap = passedmap

        #insert dq array into dark reffile
        finaldark.dq = dqmap
        finaldark.dq_def = dqdef

        #create the proper header
        finaldark = self.make_proper_header(finaldark,files)

        #save the reference file
        finaldark.save(self.outfile)

    def make_proper_header(self,model,files):
        '''Insert CRDS-required keywords'''
        
        #read in headers from one of the input files
        with fits.open(files[0]) as h:
            header0 = h[0].header
            header1 = h[1].header

        model.meta.reffile.type = 'DARK'
        model.meta.exposure.groupgap = 1
        model.meta.exposure.nframes = 1

        model.meta.subarray.xsize = header0['SUBSIZE1']
        model.meta.subarray.ysize = header0['SUBSIZE2']
        model.meta.subarray.xstart = header0['SUBSTRT1']
        model.meta.subarray.ystart = header0['SUBSTRT2']
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
        model.meta.reffile.description = 'Dark Current reffile from CV3 data'
        model.meta.reffile.pedigree = 'GROUND'
        model.meta.reffile.useafter = '2015-10-01'

        #HISTORY keyword
        model.history.append('Description of Reference File Creation')

        model.history.append('DOCUMENT:')
        model.history.append('JWST-STScI-TR-XXXX')

        model.history.append('SOFTWARE:')
        model.history.append('/ifs/jwst/wit/witserv/data4/nrc/')
        model.history.append('hilbert/darks/A1/make_NIRCAM_dark_reffile.py')
        
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



    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='Make NIRCam dark current reference file.')
        parser.add_argument("infile",help="List file of input dark current ramps.")
        parser.add_argument("--outfile",help="Name of output dark current file",default=None)
        parser.add_argument("--maxopen",help="Maximum number of groups to have open at one time.",default=200,type=int)
        parser.add_argument("--read_dq",help="File to read in that contains the DQ map to use.", default=None)
        parser.add_argument("--empty_dq",help="If set, the reference file will have a DQ array of all zeros.",action='store_true',default=False)
        parser.add_argument("--load_cr",help="Read in CR hit map from a previous run of the script.", default=False, action='store_true')
        parser.add_argument("--empty_cr",help="If True, CR hits are ignored when calculating means.",action='store_true',default=False)
        parser.add_argument("--hot_sigfrac",help="Multiple of local mean signal above which a pixel is flagged as hot",default=2,type=int)
        parser.add_argument("--min_hot_frac",help="Pixels found to be hot in more than min_hot_frac fraction of the input integrations are defined as hot",default=0.4)
        parser.add_argument("--hotfile",help="Name of file to save hot pixel mask to.",default=None)
        parser.add_argument("-v","--verbose",default=False,action='store_true')
        return parser


if __name__ == '__main__':
    usagestring = 'make_NIRCAM_dark_reffile.py inputs.list'

    starting = datetime.datetime.now()

    dark = make_nrc_dark_reffile()
    parser = dark.add_options(usage=usagestring)
    args = parser.parse_args(namespace=dark)

    dark.run(args.infile,args.outfile,args.maxopen)

    ending = datetime.datetime.now()
    difftime = ending - starting
    print("DONE. Elapsed time: {}".format(difftime))
