#!/usr/bin/env python



#script to compare: WFC3 cycle 19 bad pixel generator:
#/grp/hst/wfc3d/hilbert/cycle19/badpixtab/cdbs_ir_badpixtab_generator_cy19.pro


'''
notes from confluence page on bad pix values from build 5
hot, unreliable_slope go in dark reference files as 2,16       unreliable_slope == unstable?? set unstable as do_not_use?
dead, lowQE go in the BPM as 4,8 RC in BPM as ?? set dead,lowQE as do_not_use
bad in 0th ???
weird???

karl has
hot
dead
rc
noisy
weird == nonlinear in lin reffile?? not really. they are only nonlinear above a certain threshold
'''         




import numpy as np
import sys
from astropy.io import fits,ascii
sys.path.append('/user/hilbert/detector_test/python_procs/nircam/utils/')
import sigmacut
import argparse
import datetime,os
from itertools import izip
import even_odd_correct
import matplotlib.pyplot as plt
from astropy.table import Table,Column,vstack
from jwst_lib.models import MaskModel,dqflags
import subprocess

#qstart = [0,508,1020,1532,2040]
qstart = [4,512,1024,1536,2044]

lindir = '/grp/jwst/wit/nircam/CV2_reffile_delivery_v1/Linearity_delivery/'

class Badpix_Map:
    def __init__(self):
        self.verbose = False
        self.ron = 7.  #readnoise approximation 

    def get_filelist(self,file):
        files = []
        with open(file) as f:
            for line in f:
                if len(line) > 3:
                    files.append(line.strip())

        return files


    def sigma_clipped_mean(self,pix):
        #calculate and return the sigma-clipped mean and stdev of the input data
        pix = np.ravel(pix)
        nanpts = np.isnan(pix)
        pixclass = sigmacut.calcaverageclass()
        pixclass.calcaverage_sigmacutloop(pix,Nsigma=3,saveused=True,mask=nanpts)
        finalmask = np.logical_or(nanpts,pixclass.clipped)
        pixclass.calcaverage_sigmacutloop(pix,Nsigma=0,mask=finalmask)

        infpoints = ~np.isfinite(pix)
        if len(np.where(infpoints == True)[0]) > 0:
            #print('At least some infinite data point values present for sigma_clipped_mean. Ignoring and repeating calculation.')
            nanpts = np.logical_or(np.isnan(pix),infpoints)
            pixclass = sigmacut.calcaverageclass()
            pixclass.calcaverage_sigmacutloop(pix,Nsigma=3,saveused=True,mask=nanpts)
            finalmask = np.logical_or(nanpts,pixclass.clipped)
            pixclass.calcaverage_sigmacutloop(pix,Nsigma=0,mask=finalmask)

        return pixclass.mean,pixclass.stdev


    def cr_map_from_dq(self,filename,crval=4):
        #generate a CR map using the data quality arrays of the input file
        h = fits.open(filename)
        dq = h[3].data
        h.close()

        dims = dq.shape
        
        #I don't think it's possible to have only 3 dimensions with SSB
        #data, but just in case...
        if dims == 3:
            dq = np.expand_dims(dq,axis=0)

        crmap = np.zeros((dq.shape[0],dq.shape[2],dq.shape[3]))
        yd = dq.shape[2]
        xd = dq.shape[3]
        
        for integration in xrange(dq.shape[0]):
            for y in xrange(yd):
                for x in xrange(xd):
                    integ = dq[integration,:,y,x]
                    deltadq = integ - np.roll(integ,1)
                    deltadq = deltadq[1:]
                    change = np.where(deltadq != 0)[0]
                    if len(change) > 0:
                        for chg in change:
                            
                            #easy case, dq array changes by 4
                            if deltadq[chg] == crval:
                                crmap[integration,y,x] = chg+1
                                break #once we find the first CR, stop looking
                            
                            #harder case, dq array changes by crval+another bit
                            #at the same time
                            else:
                                power = np.log2((deltadq[chg]-crval))
                                intcheck = power - int(power)
                                if np.abs(intcheck) < 1.**-10:
                                    crmap[integration,y,x] = chg+1
                                    break #once we find the first CR, stop looking
        return crmap


    def multi_step_sigma_clipped_mean(self,pix):
        #calculate and return the sigma-clipped mean and stdev of the input data
        #use Armin's more careful method of determining the mean
        pix = np.ravel(pix)
        nanpts = np.isnan(pix)
        good = ~np.isnan(pix)
        goodpix = pix[good]

        med = np.median(goodpix)
        dev = np.std(goodpix)
        #cut = ((goodpix < (med-3*dev)) or (goodpix > (med+3*dev))) 
        cut = np.logical_or(goodpix < (med-3*dev),goodpix > (med+3*dev))

        #pixclass = sigmacut.calcaverageclass()
        #pixclass.calcaverage_sigmacutloop(goodpix,Nsigma=3,saveused=True,mask=cut)
        #cut2 = ((goodpix < (pixclass.mean-3*pixclass.stdev)) or (goodpix > (pixclass.mean+3*pixclass.stdev))) 
        #
        #pixclass2 = sigmacut.calcaverageclass()
        #pixclass2.calcaverage_sigmacutloop(goodpix,Nsigma=3,saveused=True,mask=cut2)
        #cut3 = ((goodpix < (pixclass2.mean-3*pixclass2.stdev)) or (goodpix > (pixclass2.mean+3*pixclass2.stdev))) 
        
        #print('median and dumb stdev: ',med,dev)

        counter = 0
        while True:
            #print('sigma mean counter = ',counter)
            pixclass = sigmacut.calcaverageclass()
            pixclass.calcaverage_sigmacutloop(goodpix,Nsigma=3,saveused=True,mask=cut)
            #print('pixclass mean and stdev: ',pixclass.mean,pixclass.stdev)
            #cut2 = ((goodpix < (pixclass.mean-3*pixclass.stdev)) | (goodpix > (pixclass.mean+3*pixclass.stdev)))
            cut2 = np.logical_or(goodpix < (pixclass.mean-3*pixclass.stdev),goodpix > (pixclass.mean+3*pixclass.stdev))
            if (cut == cut2).all():
                break
            else:
                cut = cut2
                counter+=1
            if counter == 100:
                print('multi_step_sigma_clipped_mean not converging!')
                sys.exit()


        #pixclass.calcaverage_sigmacutloop(pix,Nsigma=0,mask=finalmask)
        return pixclass.mean,pixclass.stdev
        

    def cr_search(self,file,cr_sigma=5,flag_neighbors=False):
        #look for cosmic ray hits in the input file
        #return a map where the value in each pixel is the
        #group number where a cosmic ray was found. Return 
        #-1 for pixels with no cosmic ray
        
        #if the flag_neighbors flag is set, then for all found
        #CR hits, flag the 4 nearest neighbors in addition to the CR-imacted pixel

        #get data
        h = fits.open(file)
        data = h[1].data


        #work-around: crappy linearity correction is screwing the cr-search at high 
        #signal levels. so if we are looking at a flat from our test dataset, limit it
        #to the first 120 reads
        #if data.shape[1] > 121:
        #    print('TRUNCATING FLAT RAMP FOR CODE TESTING DUE TO POOR LINEARITY CORRECTIONS')
        #    data = data[:,0:120,:,:]

        nint = data.shape[0]
        ngroup = data.shape[1]
        ydim = data.shape[2]
        xdim = data.shape[3]

        #loop over integrations if present
        cr_map = np.zeros((nint,ydim,xdim)) - 1
        booger = 0
        if booger != 1:

            for integration in xrange(nint):

                #calculate differences between adjacent groups
                diffs = np.zeros((ngroup-1,ydim,xdim))
                for group in xrange(ngroup-1):
                    diffs[group,:,:] = data[integration,group+1,:,:] - data[integration,group,:,:]

                #calculate the mean and stdev of the rates for each pixel
                #for y in xrange(ydim):
                #    for x in xrange(xdim):
                for y in xrange(4,2044):
                    for x in xrange(4,2044):
                        ramp = diffs[:,y,x]
                        nancheck = ~np.isnan(ramp)
                        if any(nancheck):
                            mn,dev = self.sigma_clipped_mean(ramp)
                    
                            #now look for groups where the rate is more than N sigma from the mean.
                            #bad = (ramp > (mn+cr_sigma*dev)) | (ramp < (mn-cr_sigma*dev))
                            #print(x,y,mn,dev,cr_sigma,ramp)
                            bad = np.where((ramp > (mn+cr_sigma*dev)) | (ramp < (mn-cr_sigma*dev)))[0]
                            #bad = np.where((ramp > (mn+cr_sigma*dev)) | (ramp < (mn-cr_sigma*dev)),1,0)
                            if len(bad) > 0:
                                cr_map[integration,y,x] = bad[0]
                                #if flag_neighbors is set, then flag the 4 neighboring pixels along with
                                #the CR-impacted signal. But, if a neighboring pixel is already flagged
                                #because of an earlier CR hit, then keep the original flag.
                                if flag_neighbors == True:
                                    square = cr_map[integration,y-1:y+2,x-1:x+2]
                                    m1 = square == -1
                                    square[m1] = 10000
                                    square[0,1] = np.min([square[0,1],bad[0]])
                                    square[1,0] = np.min([square[1,0],bad[0]])
                                    square[1,2] = np.min([square[1,2],bad[0]])
                                    square[2,1] = np.min([square[2,1],bad[0]])
                                    square[m1] = -1
                                    cr_map[integration,y-1:y+2,x-1:x+2] = square
                        else:
                            pass
        else:
            print('skipping actual cr search for code testing')

        #create a table of how many pixels in each quadrant suffered a CR hit
        stats_cr = self.table_of_badpix(cr_map)
        source = Column(['CR_Hit'],name='Source')
        stats_cr.add_column(source,index=0)

        #save table
        slash = self.darklist.rfind('/')
        #dir = self.darklist[0:slash+1]
        darklistfileonly = self.darklist[slash+1:]
        cr_table_file = self.outdir+'CR_stats_table_from_'+darklistfileonly+'.txt'
        stats_cr.write(cr_table_file,format='ascii')

        #create a histogram of the groups of the CR hits, as a sanity check
        fcr,acr = plt.subplots()
        n,bins,patches = acr.hist(np.ravel(cr_map),facecolor='blue',range=(0,ngroup+1))
        crhistfile = self.outdir+'CR_hits_histogram_from_'+darklistfileonly+'.eps'
        fcr.savefig(crhistfile,clobber=True)

        #save the resulting CR maps
        dot = file.rfind('.')
        slash = file.rfind('/')
        croutfile = self.outdir+file[slash+1:dot] + '_CRMAPS.fits'
        h[1].data = cr_map
        h.writeto(croutfile,clobber=True)

        print("Cosmic ray hit map for {} written to {}".format(file,croutfile))

        return cr_map

    def badin0th(self,darklist,sigma=5):
        #identify pixels that are bad in the 0th read.
        #The original strategy from Don Hall was to sum the 0th
        #reads from two ramps, and then sigma-clip the results. 
        #In this case, it seems appropriate to make a mean 0th read
        #from all the input ramps, and then clip.

        h = fits.open(darklist[0])
        xd = h[1].header['NAXIS1']
        yd = h[1].header['NAXIS2']
        mean0 = np.zeros((yd,xd))

        #create mean 0th read
        for darkfile in darklist:
            h = fits.open(darkfile)
            data = h[1].data
            h.close()

            #even/odd correct the 0th read
            print('Applying even/odd correction to 0th read of '+darkfile+' in preparation for bad in 0th read search.')
            zero = data[0,0,:,:]
            eo_corr = even_odd_correct.Oddeven()
            evenmn,oddmn = eo_corr.calc_averages(zero)
            corrzero = eo_corr.submeans(zero,evenmn,oddmn)

            mean0 = mean0 + corrzero
        mean0 = mean0 / len(darklist)

        #set up the bad pixel map that will be returned
        badmap = np.zeros(mean0.shape)

        #Calculate means and rms values separately for each quadrant, and separately for 
        #the science pixels vs reference pixels.
        for i in xrange(4):
            sector = mean0[4:2044,qstart[i]:qstart[i+1]]
            topref = mean0[2044:,qstart[i]:qstart[i+1]]
            bottomref = mean0[0:4,qstart[i]:qstart[i+1]]
            sideref = []

            if i == 0:
                sideref = mean0[:,0:4]
            if i == 3:
                sideref = mean0[:,2044:]
            refpix = np.concatenate([np.ravel(topref),np.ravel(bottomref),np.ravel(sideref)])


            #Now find the bad pixels in the science pixels
            scimn,scidev = self.sigma_clipped_mean(sector)
            badsector = ((sector < (scimn-sigma*scidev)) | (sector > (scimn+sigma*scidev)))
            badmap[4:2044,qstart[i]:qstart[i+1]] = badsector

            #Now find the bad pixels in the reference pixels
            refmn,refdev = self.sigma_clipped_mean(refpix)
            badtop = ((topref < (refmn-sigma*refdev)) | (topref > (refmn+sigma*refdev)))
            badbottom = ((bottomref < (refmn-sigma*refdev)) | (bottomref > (refmn+sigma*refdev)))
            badmap[2044:,qstart[i]:qstart[i+1]] = badtop
            badmap[0:4,qstart[i]:qstart[i+1]] = badbottom
            if i == 0:
                badside = ((sideref < (refmn-sigma*refdev)) | (sideref > (refmn+sigma*refdev)))
                badmap[:,0:4] = badside
            if i == 3:
                badside = ((sideref < (refmn-sigma*refdev)) | (sideref > (refmn+sigma*refdev)))
                badmap[:,2044:] = badside
                
        #create a table of how many pixels in each quadrant are bad in the 0th read
        stats_badin0 = self.table_of_badpix(badmap)
        source = Column(['Bad_in_0_Group'],name='Source')
        stats_badin0.add_column(source,index=0)

        #save table
        slash = self.darklist.rfind('/')
        #dir = self.darklist[0:slash+1]
        darkfilenameonly = self.darklist[slash+1:]
        bad0_table_file = self.outdir+'Badin0_stats_table_from_'+darkfilenameonly+'.txt'
        stats_badin0.write(bad0_table_file,format='ascii')
        print("Bad in 0th read stats file for {} written to {}".format(self.darklist,bad0_table_file))

        return badmap.astype(int)


    def dead_search(self,file,dead_sigfrac=0.5,in_mask=None):
        #search for dead pixels in the input file.
        h = fits.open(file)
        data = h[1].data
        tgroup = h[0].header['TGROUP']
        header0 = h[0].header
        h.close()

        #assume 4 rows and cols of refpix which we will ignore
        #qstart = [4,512,1024,1536,2044]
        yd = data.shape[2]
        xd = data.shape[3]

        #use the mask to set CR-hit pixels to NaN
        for integ in xrange(data.shape[0]):
            for y in xrange(yd):
                for x in xrange(xd):
                    if in_mask[y,x] != -1:
                        data[integ,in_mask[y,x]:,y,x] = np.nan

        ngroup = data.shape[1]
        readtimes = np.arange(ngroup) * tgroup #seconds, assuming full-frame
        dead1 = np.zeros((data.shape[0],yd,xd))
        deadmap = np.zeros((data.shape[0],yd,xd))
        for integration in xrange(data.shape[0]):
            #we want a good measure of the signal accumulated up the ramp, but it might
            #be wise to avoid saturation due to pixels saturating at different levels.
            #Let's try to keep the mean signal under 30000 DN.

            #So first get a quick estimate of the signal rate
            rate = np.nanmean((data[integration,3,:,:] - data[integration,0,:,:])) / (tgroup*3)

            #At that rate, how long to get to 30000DN
            timetohighsig = 30000. / rate
            tdiff = timetohighsig - readtimes

            #The group number to use. This is the last group taken before timeto25k 
            if np.max(tdiff) > 0:
                final_group = np.where(tdiff > 0)[0][-1]
            else:
                final_group = data.shape[1]

            #print('rate:',rate)
            #print('final group is {}.'.format(final_group))
            #Now take difference
            diff = (data[integration,final_group,:,:] - data[integration,1,:,:]) / (tgroup*(final_group-1))
            
            #print(len(np.where(np.isnan(diff))[0]))
            #print(len(np.where(~np.isnan(diff))[0]))

            #dead search #1, pixels with signal <= 0
            #exclude reference pixels
            tmp = diff[4:yd-4,4:xd-4]
            dead0 = tmp <= 0.
            dead1[integration,4:yd-4,4:xd-4] = dead0

            #Work on one quad at a time
            for i in xrange(len(qstart)-1):
                #print("Dead search starting on quad {}".format(i+1))
                #print("Extracting pixels (x1-x2):",qstart[i],qstart[i+1])
                #print("(y1-y2) ",4,yd-4)
                qdata = diff[4:yd-4,qstart[i]:qstart[i+1]]
                #print('qdata shape is ',qdata.shape)
                
                #if np.max(in_mask) > -1:
                #    qmask = in_mask[4:yd-4,qstart[i]:qstart[i+1]]
                #    bad = qmask != -1
                #    #print(i,len(np.where(bad ==True)[0])*1.)
                #    qdata[bad] = np.nan

                #loop over quad and make 53x53-pixel boxes
                #53 in rather than 50 so that we don't have a tiny box
                #along the edges of the quad
                totx = qdata.shape[1] / 53
                toty = qdata.shape[0] / 53
                quaddead = np.zeros(qdata.shape)
                for xbox in xrange(totx+1):
                    xstrt = 53 * xbox
                    xend = 53 * (xbox+1)
                    if xend > qdata.shape[1]:
                        xend = qdata.shape[1]+1

                    for ybox in xrange(toty+1):
                        ystrt = 53 * ybox
                        yend = 53 * (ybox+1)
                        if yend > qdata.shape[0]:
                            yend = qdata.shape[0]+1

                        box = qdata[ystrt:yend,xstrt:xend]

                        mnbox,devbox = self.sigma_clipped_mean(box)
                        #print('Dead search on box {},{} to {},{}. Mean of box is {}'.format(xstrt+qstart[i],xend+qstart[i],ystrt,yend,mnbox))
                        #within each box, declare any pixels with signal less than
                        #dead_sigfrac times the mean to be dead.

                        #nans = np.isnan(box)
                        #numnan = len(np.where(nans == True)[0])*1.
                        #totalbox = box.shape[0]*box.shape[1]*1.
                        #print(ystrt,yend,xstrt,xend,mnbox,dead_sigfrac,numnan/totalbox)
                        deadbox = (box < (mnbox*dead_sigfrac)) & (box > 0)
                        #nlow = len(np.where(deadbox == True)[0])
                        quaddead[ystrt:yend,xstrt:xend] = deadbox

                #insert search results fo the quad into the full dead map
                deadmap[integration,4:yd-4,qstart[i]:qstart[i+1]] = quaddead
                
            
            #save the two types of deadmaps.
            #g = fits.PrimaryHDU()
            #g1 = fits.ImageHDU(dead1.astype(int))
            #g2 = fits.ImageHDU(deadmap.astype(int))
            #ghdu = fits.HDUList([g,g1,g2])
            #ghdu.writeto('deadmaps.fits',clobber=True)
            #booger

            #combine the two dead maps
            #deadmap[integration,:,:] = np.any([deadmap[integration,:,:],dead1],axis=0)


        ndead = len(np.where(dead1 == 1)[0])
        nlow = len(np.where(deadmap == 1)[0])
        print('number of bad pix in dead map: ',ndead)
        print('number of low QE: ',nlow)


        #save the deadpixel map
        dot = file.rfind('.')
        slash = file.rfind('/')
        deadmapfile = self.outdir + file[slash+1:dot] + '_DEADMAP.fits'
        h = fits.PrimaryHDU()
        header0['MAPTYPE'] = 'DEAD'
        h1 = fits.ImageHDU(dead1,header0)
        header0['MAPTYPE'] = 'LOWQE'
        h2 = fits.ImageHDU(deadmap,header0)
        hdu = fits.HDUList([h,h1,h2])
        hdu.writeto(deadmapfile,clobber=True)
        print("Dead map for {} written to {}".format(file,deadmapfile))
        return dead1,deadmap
        

    def table_of_badpix(self,array):
        #given a bad pixel table/map of some sort, create a table that lists the number
        #of bad pixels in each quadrant. Return a table.
        bad_stats = Table()
        for j in xrange(4):
            quad = array[4:2044,qstart[j]:qstart[j+1]]
            header = 'Quad'+str(j+1)
            bad_stats[header] = [len(np.where(quad != 0)[0])]
        return bad_stats


    def combine_ind_maps(self,maps,minfrac=0.4):
        #combine the dead pixel maps for all of the input files into a final dead pixel map
        #If a pixel is flagged as dead in more than minfrac fraction of the inputs, then the pixel 
        #is considered dead in the final map.
        nmap = maps.shape[0]
        yd = maps.shape[1]
        xd = maps.shape[2]
        mapfrac = np.sum(maps,axis=0) / nmap
        
        finalmap = np.zeros((yd,xd))
        finalmap = mapfrac >= minfrac
        finalmap = finalmap.astype(int)
        return finalmap

    def hot_search(self,file,hot_sigfrac=10,in_mask=0):
        #mean and noise values used to determine the hot pixel threshold.
        #Better to use a consistent threshold ramp-to-ramp, rather than
        #relying on calculated mean and stdev for each.
        mnvals = {'SW':0.,'LW':20}
        devvals = {'SW':8./np.sqrt(2),'LW':8./np.sqrt(2)}

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
        for integ in xrange(data.shape[0]):
            for y in xrange(yd):
                for x in xrange(xd):
                    if in_mask[y,x] != -1:
                        data[integ,in_mask[y,x]:,y,x] = np.nan


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

            #histograms - put in units of DN/sec
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
            allhname = self.outdir+'Allpix_hist_'+file[slash+1:]+'_Integ'+str(integration)+'.pdf'
            hothname = self.outdir+'Hotpix_hist_'+file[slash+1:]+'_Integ'+str(integration)+'.pdf'
            allfig.savefig(allhname)
            hotfig.savefig(hothname)
            plt.close(allfig)
            plt.close(hotfig)
            print("All pixel and hot pixel histograms for {} saved to {} and {}.".format(file,allhname,hothname))


        #save the hotpixel map
        dot = file.rfind('.')
        slash = file.rfind('/')
        hotmapfile = self.outdir + file[slash+1:dot] + '_HOTMAP.fits'
        h[1].data = hotmap
        h.writeto(hotmapfile,clobber=True)
        print("Hotmap for {} written to {}.".format(file,hotmapfile))

        return hotmap
        
    def make_mean_img(self,fileimgs): #,crimgs):
        #assume for the moment that inputs are arrays. fileimgs will have to be an array made from images 
        #(i.e. ramps have been converted to slopes or last-first, etc)
        
        #set any pixels hit by cosmic rays to NaN.
        #filimgs[crimgs == 1] = np.nan

        #calculate the sigma-clipped mean and stdev through the stack
        yd = fileimgs.shape[1]
        xd = fileimgs.shape[2]
        mean_im = np.zeros((yd,xd))
        dev_im = np.zeros((yd,xd))
        for x in xrange(xd):
            for y in xrange(yd):
                mn,dev = self.sigma_clipped_mean(fileimgs[:,y,x])
                mean_im[y,x] = mn
                dev_im[y,x] = dev
        return mean_im,dev_im

    def ramp_to_im(self,filename,crmask):
        #take the input file name of a ramp, open, and translate into an image.
        #For high signal ramps (flats), use read m - read n.
        #For SW darks, just calculate the mean. (??? or should we line-fit???)
        #For LW darks...(line-fit???)

        #read in ramp
        h = fits.open(filename)
        data = h[1].data
        time_bet_groups = h[0].header['TGROUP']
        h.close()

        ngroup = data.shape[1]
        yd = data.shape[2]
        xd = data.shape[3]
        images = np.zeros((data.shape[0],yd,xd))

        #loop over integrations
        for integration in xrange(data.shape[0]):
            cr = crmask[integration,:,:]

            #apply CR masks
            for y in xrange(yd):
                for x in xrange(xd):
                    if cr[y,x] != -1:
                        data[cr[y,x]:,y,x] = np.nan

            #find signal level and decide which strategy to use.
            #quickrate = np.nanmean(data[integration,4,:,:]-data[integration,0,:,:]) / (time_bet_groups*4)
            #time_to_10k = 10000. / quickrate
            

            #just linefit regardless. Even with the high signal ramps, we would have to check for
            #NaNs and then adjust which difference we take to avoid the CR hits.
            #if time_to_10k < 0:  #SW darks may fall in here
            im,imunc = self.linefit(data,time_bet_groups)
            #if time_to_10k < (time_bet_groups*ngroup): #if ramps hits saturation, just do image difference/time
            #    dt = time_to_10k - (time_bet_groups * np.arange(ngroup))
            #    m = np.where(dt > 0)[0]
            #    m = m[-1]
            #    im = (data[integration,m,:,:] - data[integration,0,:,:]) / (dt*m)
            #    imunc = np.sqrt(im)
            #if time_to_10k >= (time_bet_groups*ngroup):
            #    im,imunc = self.linefit(data,time_bet_groups)
        return im,unc


    def slope_from_last_minus_first(self,filename,crmask):
        '''
        Attempt to speed things up a bit. Get an approximation of slopes
        from last-first, taking into account satruated groups. Probably not
        a good idea to use this for darks.
        '''
        #read in ramp
        h = fits.open(filename)
        data = h[1].data
        time_bet_groups = h[0].header['TGROUP']
        groupdq = h[3].data
        h.close()

        #pull out the saturation bit from the DQ masks. That's all we care about here.
        satval = dqflags.pixel['SATURATED']
        #satval = 2.
        #print("WARNING: assuming a value of 2 for saturation in the groupdq array.")
        groupdqb = np.expand_dims(groupdq,axis=0)
        groupdqb = np.unpackbits(groupdqb,axis=0)
        groupdqb = groupdqb[7-np.log2(satval)]
        groupdqb = groupdqb[0,:,:,:,:]

        ngroup = data.shape[1]
        yd = data.shape[2]
        xd = data.shape[3]
        images = np.zeros((data.shape[0],yd,xd))
        uncs = np.zeros((data.shape[0],yd,xd))

        #loop over integrations
        for integration in xrange(data.shape[0]):
            cr = crmask[integration,:,:]
            #dq = groupdqb[integration,:,:,:]
            
            #loop over pix
            for y in xrange(yd):
                for x in xrange(xd):

                    #apply CR masks
                    if cr[y,x] != -1:
                        data[cr[y,x]:,y,x] = np.nan

                    #apply saturation flags
                    dqpix = dq[integration,:,y,x]
                    data[dqpix == 1] = np.nan
                    
                    #take last - first difference
                    good = ~np.isnan(data[:,y,x]) 
                    goodindex = np.where(good == True)[0]
                    if len(goodindex) > 0:
                        diff = (data[goodindex[-1],y,x] - data[goodindex[0],y,x]) / ((goodindex[-1]-goodindex[0])*time_bet_groups)
                        err = np.sqrt( (data[goodindex[-1],y,x]+self.ron**2) + (data[goodindex[0],y,x]+ self.ron**2) )
                    else:
                        diff = np.nan
                        err = np.nan

                    images[integration,y,x] = diff
                    uncs[integration,y,x] = err

        return images,uncs
        

    def linefit(self,data,dt):
        #perform line fitting to transform a ramp into a slope image.
        yd = data.shape[1]
        xd = data.shape[2]
        slopeim = np.zeros((yd,xd))
        
        time = np.arange(0,dt*data.shape[0],dt)

        for y in xrange(4,2044):
            for x in xrange(4,2044):
                signals = data[:,y,x]

                #if CR hits have introduced nan's into the data, then
                #line-fit only up to the first nan.
                firstnan = np.nonzero(np.isnan(signals))[0]
                if len(firstnan) != 0:
                    signals = signals[0:firstnan[0]]
                    itime = time[0:firstnan[0]]

                coeffs = np.polyfit(itime,signals,1,cov=True)
                slopeim[y,x] = coeffs[0][1]
                slopeunc[y,x] = coeffs[1][1,1]
        
        return slopeim,slopeunc



    def noise_in_meanimg(self,meanimg,thresh):
        #search a mean noise image (noise/signal) for pixels which exhibit anomalously
        #high noise levels. This is designed to catch pixels that are pretty wildly 
        #inconsistent, with signals that vary quite a bit from ramp-to-ramp.

        yd = meanimg.shape[0]
        xd = meanimg.shape[1]
        bad  = np.zeros((yd,xd))

        #Work one amp at a time in case the mean noise characteristics vary by quadrant
        badmap = np.zeros(meanimg.shape)
        #qstart = [4,512,1024,1536,2044]
        for i in xrange(4):
            amp = meanimg[4:2044,qstart[i]:qstart[i+1]]
            mn,dev = self.sigma_clipped_mean(amp)

            #do we care about pixels with anomalously LOW noise? Probably not....
            ampnoise = (amp > (mn+thresh*dev)) #| (amp < (mn-thresh*dev))
            badmap[4:2044,qstart[i]:qstart[i+1]] = ampnoise
            badmap = badmap.astype(int)
        return badmap


    def inconsistent_signal_search(self,indslope,indslopeunc,meanslope,devslope,thresh):
        #Compare individual files to the mean file, and mark as bad pixels with signal that 
        #is inconsistent with that in the mean file.
        #This is designed to catch pixels that are self-consistent most of the time, but then
        #show an occasionnal bad ramp. (i.e. the sigma-clipped stdev will be small, so they won't
        #be caught by noise_in_meanimg above)
        yd = indslope.shape[1]
        xd = indslope.shape[2]
        totalbad = np.zeros((yd,xd))

        for i in xrange(indslope.shape[0]):
            ##for frame,frameunc in izip(indslope,indslopeunc):
            deltamean = thresh * devslope
            deltaind = thresh * indslopeunc[i,:,:]
            deltahigh = (indslope[i,:,:]+deltaind) - (meanslope-deltamean)
            bad1 = np.where(deltahigh < 0,1,0)
            deltalow = (indslope[i,:,:]-deltaind) - (meanslope+deltamean)
            bad2 = np.where(deltalow > 0,1,0)
            totalbad = totalbad + bad1 + bad2
            
        totalbad[totalbad > 0] = 1
        return totalbad


    def unstable_search(self,darkfiles,flatfiles,dark_cr_map,flat_cr_map):
        #multiple parts
        #create mean flat field
        #STEP 1 - PIXEL WHOSE RELATIVE NOISE LEVEL IN THE OVERALL MEAN FLAT
        #FIELD (NOISE/SIGNAL), IS MORE THAN X-SIGMA FROM THE MEAN
        #STEP 2 - PIXEL WITH FLT VALUES THAT VARY IN TIME BY MORE THAN Y-SIGMA
        #;DARK CURRENT
        #;need to do for the darks what we did for the flats in the dead pixel
        #;search. need to create an overall mean dark image and uncertainty, 
        #;and compare each individual dark to that. no need to normalize though
        #;since there's no filter dependence here.
        #
        #;STEP 1 - PIXEL WHOSE RELATIVE NOISE LEVEL IN THE OVERALL MEAN DARK
        #;         FIELD (NOISE/SIGNAL), IS MORE THAN X-SIGMA FROM THE MEAN
        #;STEP 2 - PIXEL WITH FLT VALUES THAT VARY IN TIME BY MORE THAN Y-SIGMA'


        #create a mean flat field integration. If the flats span multiple filters, then create a mean
        #for each filter. Then normalize each and combine into a single overall mean flat.

        #if the user provides a list of slope image files for the flats, use those. If not, we need
        #to create our own slope images
        flatslope = np.zeros((len(flatfiles),2048,2048))
        flatslopeunc = np.zeros((len(flatfiles),2048,2048))
        
        if self.flatslopes == None:
            i = 0
            for ffile,thisflatcr in izip(flatfiles,flat_cr_map):
                flatslope[i,:,:],flatslopeunc[i,:,:] = self.ramp_to_im(ffile,thisflatcr)
                i = i + 1

        else:
            flatslopefiles = []
            with open(self.flatslopes) as fsl:
                for line in fsl:
                    if len(line) > 3:
                        fslopefile = line.strip()
                        flatslopefiles.append(fslopefile)
            
            i=0
            for ffile in flatslopefiles:
                h = fits.open(ffile)
                f = h[1].data
                ferr = h[2].data
                h.close()
            
                flatslope[i,:,:] = f
                flatslopeunc[i,:,:] = ferr
                i = i + 1
        #print('XXXXXXXXXXXXXskipping the calculation of the mean flatXXXXXXXX')
        meanflat,devflat = self.make_mean_img(flatslope) #,flat_cr_map)

        #print('do we want to deal with making separate means for separate files')
        #print('then normalizing them and combining? or should we assume a single filter?')

        darkslope = np.zeros((len(darkfiles),2048,2048))
        darkslopeunc = np.zeros((len(darkfiles),2048,2048))

        if self.darkslopes == None:
            i = 0
            for dfile,thisdarkcr in izip(darkfiles,dark_cr_map):
                darkslope[i,:,:],darkslopeunc[i,:,:] = self.ramp_to_im(dfile,thisdarkcr)
                i = i + 1

        else:
            darkslopefiles = []
            with open(self.darkslopes) as dsl:
                for line in dsl:
                    if len(line) > 3:
                        dslopefile = line.strip()
                        darkslopefiles.append(dslopefile)
            
            i=0
            for dfile in darkslopefiles:
                h = fits.open(dfile)
                d = h[1].data
                derr = h[2].data
                h.close()
            
                darkslope[i,:,:] = d
                darkslopeunc[i,:,:] = derr
                i = i + 1

        #print('XXXXXXXXXXXXXskipping the calculation of the mean darkXXXXXXXX')

        meandark,devdark = self.make_mean_img(darkslope) #dark_cr_map) - don't need CR maps here.


        #save the meandark and meanflat
        #save the individual unstable maps for later analysis
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(meandark)
        h2 = fits.ImageHDU(devdark)
        hdu = fits.HDUList([h0,h1,h2])
        slash = self.darklist.rfind('/')
        meandarkfile = 'Meandark_from_'+self.darklist[slash+1:]+'.fits'
        hdu.writeto(meandarkfile,clobber=True)
        print("Mean dark slope image for {} saved to {}.".format(self.darklist,meandarkfile))
        
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(meanflat)
        h2 = fits.ImageHDU(devflat)
        hdu = fits.HDUList([h0,h1,h2])
        slash = self.flatlist.rfind('/')
        meanflatfile = 'Meanflat_from_'+self.flatlist[slash+1:]+'.fits'
        hdu.writeto(meanflatfile,clobber=True)
        print("Mean flat slope image for {} saved to {}.".format(self.flatlist,meanflatfile))



        #-------------NOTES ON TESTS---------------------------
        #darknoise peak of hist is about 2. stdev is 2.5. values go out to almost 100.
        #5-sigma cutoff results in ~4.7% of all pixels flagged as unstable
        #7-sigma cutoff results in ~3.4%
        #9-sigma cutoff results in ~2.6%
        #
        #stdev is 2.5. readnoise is 8 per group, so 8/sqrt(108)=0.77 per ramp. this is much smaller than even 1-sigma...
        #meanwhile, 667,000 pixels have a negative slope in the mean dark...
        #pixels with very negative values of unc/dark should be unstable also, right?
        #maybe better to look only at devdark:
        #devdark - hist peak at 0.0025, stdev 0.00037
        #in this case, only 0.8% of all pixels would be unstable using a 5-sigma limit
        #              about 0.4% of all pixels                          7-sigma limit
        #              about 0.25% of all pixels                         9-sigma limit
        #
        #flatnoise peak of hist is about 0.004. values go out to almost 0.02 for flatdev/flat
        #Looking at flatdev only, the histogram peak is at about 0.16, with a stdev of 0.04
        #3-sigma cut then finds 0.6% of pixels as unstable (3-sigma, expect to catch 0.27% for gaussian dist.)
        #5-sigma cut finds 0.1% of pixels as unstable (5-sigma, expect to catch 0.00006% (aka 2 pixels) of pixels for gaussian dist.)
        #
        #This suggests we should just look at the noise in the mean dark and flat, and not worry about making a noise/signal image.
        #--------------------------------------------------------



        #print('XXXXXXXXXXXXXxReading in mean dark and flat to save time in code developmentXXXXXXXXXXX')

        #hh = fits.open('Meandark_from_dark_ramps.list_DMSorient.fits')
        #meandark = hh[1].data
        #devdark = hh[2].data
        #hh.close

        #with fits.open('Meanflat_from_flat_ramps.list_DMSorient.fits') as hhh:
        #    meanflat = hhh[1].data
        #    devflat = hhh[2].data




        #UNSTABLE, TYPE 1: PIXELS WHOSE RELATIVE NOISE LEVEL IN THE OVERALL MEAN FLAT
        #FIELD (NOISE/SIGNAL), IS MORE THAN X-SIGMA FROM THE MEAN -- DO EACH AMP SEPARATELY
        flatnoise = devflat #/ meanflat
        bad_noise_in_meanflat = self.noise_in_meanimg(flatnoise,self.flatmeannoisethresh)
        
        darknoise = devdark #/ meandark
        bad_noise_in_meandark = self.noise_in_meanimg(darknoise,self.darkmeannoisethresh)

        #histogram of the noise images to be searched for unstable pixels
        dnoisefig = plt.figure()
        dnoiseax = dnoisefig.add_subplot(1,1,1)
        #dbins = np.arange(np.min(darknoise),np.max(darknoise),(np.max(darknoise)-np.min(darknoise))/25.)
        dbins = np.arange(0,0.006,0.00002)
        dnhist, dnbins,dnlpatch = dnoiseax.hist(np.ravel(darknoise),bins=dbins,facecolor='blue')
        dnoiseax.set_xlabel('Noise of Mean Dark Slopes')
        dnoiseax.set_ylabel('Number of Pixels')
        roughmean,roughdev = self.sigma_clipped_mean(darknoise)
        dnoiseax.plot([roughmean+roughdev*self.darkmeannoisethresh,roughmean+roughdev*self.darkmeannoisethresh],[0,np.max(dnhist)],color='black',linestyle='--',label='Threshold')
        dnoiseax.legend(loc='upper right')
        dnoisehistfile = 'Darknoise_histogram.png'
        dnoisefig.savefig(dnoisehistfile)
        plt.close(dnoisefig)

        #histogram of the noise images to be searched for unstable pixels
        fnoisefig = plt.figure()
        fnoiseax = fnoisefig.add_subplot(1,1,1)
        #fbins = np.arange(np.min(flatnoise),np.max(flatnoise),(np.max(flatnoise)-np.min(flatnoise))/25.)
        fbins = np.arange(0,0.4,0.001)
        fnhist, fnbins,fnlpatch = fnoiseax.hist(np.ravel(flatnoise),bins=fbins,facecolor='blue')
        fnoiseax.set_xlabel('Noise of Mean Flat Slopes')
        fnoiseax.set_ylabel('Number of Pixels')
        roughmean,roughdev = self.sigma_clipped_mean(flatnoise)
        fnoiseax.plot([roughmean+roughdev*self.flatmeannoisethresh,roughmean+roughdev*self.flatmeannoisethresh],[0,np.max(fnhist)],color='black',linestyle='--',label='Threshold')
        fnoiseax.legend(loc='upper right')
        fnoisehistfile = 'Flatnoise_histogram.png'
        fnoisefig.savefig(fnoisehistfile)
        plt.close(fnoisefig)
        

        #UNSTABLE, TYPE 2: PIXELS IN SLOPE IMAGES WITH SIGNAL INCONSISTENT WITH THAT IN THE MEAN IMAGE
        bad_ind_flat_signal = self.inconsistent_signal_search(flatslope,flatslopeunc,meanflat,devflat,self.flat_signal_thresh)
        bad_ind_dark_signal = self.inconsistent_signal_search(darkslope,darkslopeunc,meandark,devdark,self.dark_signal_thresh)

        #combine the various flavors of unstable pixel masks into one final unstable bad pixel map for darks and one for flats
        unstable_darks = bad_noise_in_meandark + bad_ind_dark_signal
        unstable_flats = bad_noise_in_meanflat + bad_ind_flat_signal
        unstable_darks[unstable_darks > 1] = 1
        unstable_flats[unstable_flats > 1] = 1
        #unstable = bad_noise_in_meanflat + bad_noise_in_meandark + bad_ind_flat_signal + bad_ind_dark_signal
        #unstable[unstable > 1] = 1

        #stats on the number of bad pixels in each quadrant
        stats_noise_in_meanflat = self.table_of_badpix(bad_noise_in_meanflat)
        source = Column(['Unstab_from_mean_flat'],name='Source')
        stats_noise_in_meanflat.add_column(source,index=0)

        stats_noise_in_meandark = self.table_of_badpix(bad_noise_in_meandark)
        source = Column(['Unstab_from_mean_dark'],name='Source')
        stats_noise_in_meandark.add_column(source,index=0)

        stats_noise_in_indflat = self.table_of_badpix(bad_ind_flat_signal)
        source = Column(['Unstab_from_ind_flat'],name='Source')
        stats_noise_in_indflat.add_column(source,index=0)

        stats_noise_in_inddark = self.table_of_badpix(bad_ind_dark_signal)
        source = Column(['Unstab_from_ind_dark'],name='Source')
        stats_noise_in_inddark.add_column(source,index=0)
            

        #save the individual unstable maps for later analysis
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(bad_noise_in_meanflat)
        h2 = fits.ImageHDU(bad_noise_in_meandark)
        h3 = fits.ImageHDU(bad_ind_flat_signal)
        h4 = fits.ImageHDU(bad_ind_dark_signal)
        hdu = fits.HDUList([h0,h1,h2,h3,h4])
        indunstablefile = 'Individual_types_of_unstable_pix.fits'
        hdu.writeto(indunstablefile,clobber=True)

        print("Unstable maps for {} and {} saved to {}.".format(self.darklist,self.flatlist,indunstablefile))

        unstab_table = vstack([stats_noise_in_meanflat,stats_noise_in_meandark,stats_noise_in_indflat,stats_noise_in_inddark])
        
        print("Number of unstable pixels found, broken down by type of unstable pixel and quadrant:")
        print(unstab_table)
        unstab_table_file = 'Table_num_of_unstab_types_per_quad_from_'+self.darklist+'_'+self.flatlist+'.txt'
        unstab_table.write(unstab_table_file,format='ascii')
        print("Unstable pixel statistics saved to {}".format(unstab_table_file))

        return unstable_darks, unstable_flats


    def convert_dq_to_cr(self,file,crval=4):
        #grab the DQ array of an input ramp, extract the cosmic ray flags, and convert into a map image of CR hits,
        #where the value in each pixel is equal to the first group number where a cosmic ray was found. Groups later than
        #this will be ignored in the various bad pixel searches.
        print("Starting CR conversion step: ",datetime.datetime.now())

        h = fits.open(file)
        dq = h[3].data
        h.close()

        #need to find the definition of the DQ flags. Need to confirm the bit value corresponding to a CR hit
        crval = np.uint8(crval)
        crval = np.unpackbits([crval])
        badbit = np.where(crval == 1)[0][0]

        #convert the DQ integer flags to binary
        dq = np.expand_dims(dq,axis=4)
        dqb = np.unpackbits(dq,axis=4)
        
        #keep only the dq bits for the 'jumps'
        dqb = dqb[:,:,:,:,badbit]

        #create cr_map
        cr_map = np.zeros((dq.shape[0],2048,2048))
        
        #search for cr-hits 
        for integration in xrange(dq.shape[0]):
            for y in xrange(dq.shape[2]):
                for x in xrange(dq.shape[3]):
                    ramp = dqb[integration,:,y,x]
                    onebad = np.where(ramp != 0)[0]
                    if len(onebad) > 0:
                        cr_map[integration,y,x] = onebad[0]
                    else:
                        cr_map[integration,y,x] = -1
        return cr_map




    def weird_search_OLDXXXX(self,ramp,ngroup=5,nsigma=9):
        #find the "weird" pixels first identified by Karl
        #Try finding the slope for every N groups going up the ramp, and look
        #for a decrease followed by an increase just before hard saturation

        #if it doesn't exist, create a directory to hold the example plots of the 
        #weird pixels
        if not os.path.exists("weird_pix"):
            os.makedirs("weird_pix")

        weird = np.zeros((2048,2048))

        #assume that the input ramp is really that...a single ramp (i.e. integration)
        if len(ramp.shape) == 4:
            print("Weird pixel search works on one integration at a time. You sent > 1.")
            stophere


        #This is designed to work best with high SNR flats...
        
        zd = ramp.shape[0]
        yd = ramp.shape[1]
        xd = ramp.shape[2]
        
        #make a list of endpoints for the blocks of ramp to determine the slope of
        boundaries = np.arange(0,zd-1,ngroup)

        iii = 0
        #for y in xrange(4,yd-4):
        #    for x in xrange(4,xd-4):
        for y in xrange(1022,1223):
            for x in xrange(1240,1441):
                pix = ramp[:,y,x]
                
                #if the pixel has limited signal, skip it
                if (np.nanmax(pix)-np.nanmin(pix)) > 5000:
                    slopes = np.zeros(len(boundaries)-1)
                    for i in xrange(len(slopes)):
                        #use a simple last-first group within the block as our slope measurement
                        slopes[i] = pix[boundaries[i+1]] - pix[boundaries[i]]

                    #find where the slope is maximized
                    maxslopeloc = np.where(slopes == np.max(slopes))[0][-1]
                    maxslope = slopes[maxslopeloc]
                    #find the minimum slope prior to the location of the max slope
                    #minslopeloc = np.where(slopes[0:maxslopeloc+1] == np.min(slopes[0:maxslopeloc+1]))[0][-1]
                    #minslope = slopes[minslopeloc]

                    nominalslope,nominalslopeerr = self.sigma_clipped_mean(slopes[0:zd/ngroup/2])

                    #print(x,y,minslope,maxslope,nominalslope,nominalslopeerr,minslopeloc)
                    #print(slopes)

                    #find where slope goes to zero at saturation
                    zeroslopeloc = np.where(slopes < (0.+3.*nominalslopeerr))[0]
                    if len(zeroslopeloc) != 0:
                        zeroslopeloc = zeroslopeloc[0] - 4
                    else:
                        zeroslopeloc = zd/ngroup
                    #print(x,y,nominalslope,nominalslopeerr,maxslope,(nominalslope+nsigma*nominalslopeerr),(maxslope-nominalslope)/nominalslopeerr,maxslopeloc,zeroslopeloc,30./ngroup)
                    #print(maxslope,(nominalslope+nsigma*nominalslopeerr),zeroslopeloc-maxslopeloc,30./ngroup)


                    #conditions that must be true to have a "weird" pixel
                    #if ((maxslope > (nominalslope+nominalslopeerr)) and (minslope < (nominalslope-nominalslopeerr)) and (minslopeloc > (zd/2.))):
                    #if minslopeloc > (zd/ngroup/2.):
                    if ((maxslope > (nominalslope+nsigma*nominalslopeerr)) and (zeroslopeloc-maxslopeloc < 30./ngroup)):
                        #if maxslope > 0:   #get all pix, to see what we should be flagging
                        weird[y,x] = 1
                        #print("found weird pixel at ({},{})".format(x,y))
                        #print(x,y,maxslope,nominalslope,nominalslopeerr,maxslopeloc,zd/ngroup)#,minslopeloc)
                        #print(slopes)
                        #print('==========================')
                        if iii % 100 == 0:
                            print(x,y)
                            f,a = plt.subplots()
                            a.plot(np.arange(zd),pix,'ro')
                            a.set_xlabel('Group')
                            a.set_ylabel('Signal (DN)')
                            a.set_title('Pixel ('+str(x)+','+str(y)+')')
                            f.savefig('weird_pix/weirdpix_'+str(x)+'_'+str(y)+'.pdf')
                            plt.close(f)
                        iii = iii + 1
                else:
                    pass
                
        #save for code testing
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(weird)
        h = fits.HDUList([h0,h1])
        h.writeto('weirdmap.fits',clobber=True)
        return weird


    def alternate_findmax(self,slopes,zeroslopeloc,meanval,meanunc,ngroup):
        #find location of local maximum for weird pixel search
        slopes = slopes[zeroslopeloc/2:zeroslopeloc]

        #f,a = plt.subplots()
        #a.plot(np.arange(len(slopes)),slopes,'ro')

        for step in np.arange(-.25,5,.25):
            #find all pixels with slopes higher than given value
            above = np.where(slopes > (meanval-step*meanunc))[0]
            #a.plot([0,len(slopes)],[meanval-step*meanunc,meanval-step*meanunc])
            #print(np.roll(above,-1))
            #print(above)
            #shift and add indices, in order to look for gaps
            diff = np.roll(above,-1) - above
            diff = diff[0:-1]
            #print(diff)
            #Look for gaps. Any gap of more than ?? implies the slope has decreased, then increased?
            gaps = np.where(diff > 25/ngroup)[0]
            #print(gaps)
            if len(gaps) > 0:
                break

        #f.savefig('testnewmethod.pdf')
        if len(gaps) == 0:
            #case where we have a goodpix
            gaploc = -1
        else:
            gaploc = above[gaps[-1]+1] + zeroslopeloc/2
            #print(above[gaps[-1]+1],zeroslopeloc/2)
        return gaploc

            
#pixels to check with weird pixel search
#1306,1286 - looks normal. nominalslopeerr/nominalslope cutoff of 0.1 was too aggressive
#1291,1275 - early saturation, looks odd. saturates, then later signal increase, then saturated again. 
     #perhaps this is a symptom of whatever is creating the weird pixels? a change in well depth somehow??
#1293,1258 - high early slope, relaxes, then odd saturation - not weird - PROBABLY GOOD to not call this weird?????? - check karls map
#1296, 1257 - changing slope early in ramp- not called weird. Shows same well depth change as 1291,1275
#1294, 1257 - saturation only in last couple reads - no significant decrease in slope, just a sudden increase in the final ~20 reads - caught as weird.
#1293, 1257 - clearly a weird pixel but early changes in slope in first 10 reads. caught as weird. GOOD
#1258,1218 - correctly found as good pixel. Exhibits the change in well-depth seen in other pixels.
#1291,1275 - large changing slope early on. correctly flagged as high noise. Shows change in well-depth seen in other pixels.
#1240,1022 - weird.
#1324,1022 - weird
#1250,1022 - good
#1410,1028 - weird, but doesn't saturate at all
#check the pixel below in detail!
#1275,1087 - flagged as weird, but doesn't seem like it really is.
#1409,1139 - initialzero is way too early because the weird jump has such a large slope. gotta deal with this


    def weird_search(self,ramp,ngroup=2,nsigma=9):
        #find the "weird" pixels first identified by Karl
        #Try finding the slope for every N groups going up the ramp, and look
        #for a decrease followed by an increase just before hard saturation

        #if it doesn't exist, create a directory to hold the example plots of the 
        #weird pixels
        if not os.path.exists("weird_pix"):
            os.makedirs("weird_pix")

        weird = np.zeros((2048,2048))
        high_noise = np.zeros((2048,2048))

        #assume that the input ramp is really that...a single ramp (i.e. integration)
        if len(ramp.shape) == 4:
            print("Weird pixel search works on one integration at a time. You sent > 1.")
            stophere


        #This is designed to work best with high SNR flats...
        
        zd = ramp.shape[0]
        yd = ramp.shape[1]
        xd = ramp.shape[2]
        
        #make a list of endpoints for the blocks of ramp to determine the slope of
        boundaries = np.arange(0,zd-1,ngroup)

        iii = 0
        for y in xrange(4,yd-4):
            for x in xrange(4,xd-4):
        #for y in xrange(1644,1822):
        #    for x in xrange(4,xd-4):
                if x % 100 == 0:
                    if y % 100 == 0:
                        print(x,y)
                pix = ramp[:,y,x]
                #if the pixel has limited signal, skip it
                if (np.nanmax(pix)-np.nanmin(pix)) > 5000:

                    slopes = np.zeros(len(boundaries)-1)
                    for i in xrange(len(slopes)):
                        #use a simple last-first group within the block as our slope measurement
                        slopes[i] = (pix[boundaries[i+1]] - pix[boundaries[i]]) / (10.7*ngroup)

                    #initial guess at saturation location
                    #initialzero = np.where(slopes < 20.)[0]
                    initialzero = np.where(slopes/np.max(slopes) < 0.1)[0]
                    if len(initialzero) != 0:
                        initialzero = initialzero[0]
                        if initialzero < len(slopes)*0.1:
                            #large jump could be causing trouble
                            minslope = np.min(slopes)
                            minloc = np.where(slopes < (minslope+2.5))[0]
                            initialzero = minloc[0]
                    else:
                        #no saturation
                        initialzero = zd/ngroup
                    
                    #if initialzero == 0:
                    #    print('initialzero is 0!!',x,y)
                    #    stophere

                    
                    #find the nominal mean slope
                    #nominalslope,nominalslopeerr = self.sigma_clipped_mean(slopes[0:zd/ngroup/2])
                    #nominalslope,nominalslopeerr = self.sigma_clipped_mean(slopes[0:initialzero/2])
                    if initialzero/2 >= 10:
                        nominalslope,nominalslopeerr = self.sigma_clipped_mean(slopes[initialzero/2-10:initialzero/2+10])
                    if initialzero/2 < 10:
                        nominalslope,nominalslopeerr = self.sigma_clipped_mean(slopes[0:10]) #initialzero/2])

                    #print(slopes)
                    #print('nominal slope, err, initialzero, total points',nominalslope,nominalslopeerr,initialzero,zd/ngroup)
                    

                    #print('nominal slope and error, and location/2',nominalslope,nominalslopeerr,initialzero/2)
                    if nominalslopeerr/nominalslope > 0.2:
                        print("high noise pixel {},{}".format(x,y))
                        #print(x,y,nominalslope,nominalslopeerr,nominalslopeerr/nominalslope)
                        high_noise[y,x] = 1
                    else:                        
                        #find where slope goes to zero at saturation
                        zeroslopeloc = np.where(slopes < (0.+7.*nominalslopeerr))[0]
                        #zeroslopeloc = initialzero
                        if len(zeroslopeloc) != 0:
                            zeroslopeloc = zeroslopeloc[0] - 2/ngroup
                            pass
                        else:
                            zeroslopeloc = zd/ngroup

                        if zeroslopeloc < zd/ngroup/6:
                            print("moving zeroslopeloc later in ramp for pixel ({},{}). Moving to {}".format(x,y,initialzero))
                            zeroslopeloc = initialzero

                            #if zeroslopeloc is still at the beginning of the ramp after all this, the pixel is very strange and
                            #will be flagged in other ways, and we don't need to worry about it from the standpoint of it being weird or not
                            if zeroslopeloc == 0:
                                print("Pixel ({},{}) is still giving very odd results, with zeroslopeloc at zero. Setting to the end of the ramp and proceeding.".format(x,y))
                                zeroslopeloc = zd/ngroup

                                #f,a = plt.subplots()
                                #a.plot(np.arange(len(pix)),pix)
                                #f.savefig('test.pdf')
                                #stophere


                        #now back up from the point where the slope went to zero, and clip all
                        #points back to where the slope reaches the mean level...
                        #backup = np.where(slopes[zeroslopeloc/2:zeroslopeloc] < (nominalslope + 1.*nominalslopeerr))[0]
                        #hmmmmm, not sure about this strategy

                        #print(slopes)
                        #print('nominalslope,nominalslopeerr',nominalslope,nominalslopeerr)
                        #print('x,y,zeroslopeloc,initialzero,totalpoints',x,y,zeroslopeloc,initialzero,zd/ngroup)
                        ##find the minimum slope in the second half of the ramp
                        ##print(x,y,zeroslopeloc)
                        ##if ((y == 1257) and (x == 1296)) or zeroslopeloc < 22:
                        ##    print(zeroslopeloc,nominalslope,nominalslopeerr,slopes)
                        ##    stophere

                        #find the maximum slope in the part of the ramp between the minimum and saturation
                        #maxslopeloc = self.alternate_findmax(slopes,zeroslopeloc,nominalslope,nominalslopeerr,ngroup)
                        #print('new maxslopeloc ',maxslopeloc)

                        maxslopeloc = np.where(slopes[zeroslopeloc/2:zeroslopeloc] == np.max(slopes[zeroslopeloc/2:zeroslopeloc]))[0][-1]
                        maxslopeloc = maxslopeloc + zeroslopeloc/2
                        maxslope = slopes[maxslopeloc]
                        #print('old maxslopeloc ',maxslopeloc)
                        #stop

                        #print(zeroslopeloc/2,maxslopeloc)
                        #if maxslopeloc == zeroslopeloc/2:
                        if maxslopeloc < (zeroslopeloc/2 + 0.33*(zeroslopeloc - zeroslopeloc/2)):
                            #maxslopeloc = np.where(slopes[zeroslopeloc-10:zeroslopeloc] == np.max(slopes[zeroslopeloc-10:zeroslopeloc]))[0][-1]
                            #maxslopeloc = maxslopeloc + zeroslopeloc-10

                            maxslopeloc = self.alternate_findmax(slopes,zeroslopeloc,nominalslope,nominalslopeerr,ngroup)
                            if maxslopeloc == -1:
                                maxslope = -999


                            #print('maxslope at same place as zeroslope/2 ',x,y,maxslopeloc,zeroslopeloc/2)
                            #print('second try for maxslopeloc gave {}'.format(maxslopeloc))
                            #f,a = plt.subplots()
                            #xs = np.arange(len(slopes))
                            #a.plot(xs,slopes,'ro')
                            #a.plot(xs[zeroslopeloc:],slopes[zeroslopeloc:],'bo')
                            #a.set_xlabel('Group')
                            #a.set_ylabel('Signal (DN)')
                            #a.set_title('Pixel ('+str(x)+','+str(y)+')')
                            #f.savefig('test_{}_{}.pdf'.format(str(x),str(y)))
                            #plt.close(f)
                            #stop
                            
                                

                        minslopeloc = np.where(slopes[zeroslopeloc/2:maxslopeloc] == np.min(slopes[zeroslopeloc/2:maxslopeloc]))[0][-1]
                        minslopeloc = minslopeloc + zeroslopeloc/2
                        minslope = slopes[minslopeloc]


                        #print('x,y,zd/ngroup,nominalslope,nominalslopeerr,minslope,minslopeloc,maxslope,maxslopeloc,zeroslope,zeroslopeloc,nominalslopeerr/nominalslope',x,y,zd/ngroup,nominalslope,nominalslopeerr,minslope,minslopeloc,maxslope,maxslopeloc,slopes[zeroslopeloc],zeroslopeloc,nominalslopeerr/nominalslope)
                        
                        #print(maxslope,(nominalslope+nsigma*nominalslopeerr),zeroslopeloc-maxslopeloc,30./ngroup)

                        #PLOT OF SLOPE VS TIME SHOWING ZEROSLOPELOC
                        #f,a = plt.subplots()
                        #xs = np.arange(len(slopes))
                        #a.plot(xs,slopes,'ro')
                        #a.plot(xs[zeroslopeloc:],slopes[zeroslopeloc:],'bo')
                        #a.set_xlabel('Group')
                        #a.set_ylabel('Signal (DN)')
                        #a.set_title('Pixel ('+str(x)+','+str(y)+')')
                        #f.savefig('test_{}_{}.pdf'.format(str(x),str(y)))
                        #plt.close(f)
                        #stop


                        #conditions that must be true to have a "weird" pixel
                        #if ((maxslope > (nominalslope+nominalslopeerr)) and (minslope < (nominalslope-nominalslopeerr)) and (minslopeloc > (zd/2.))):
                        #if minslopeloc > (zd/ngroup/2.):
                        if ((maxslope > (minslope+nsigma*nominalslopeerr)) and (nominalslopeerr/nominalslope < 0.2)):  #and (zeroslopeloc-maxslopeloc < 30./ngroup)):
                            #if maxslope > 0:   #get all pix, to see what we should be flagging
                            weird[y,x] = 1
                            #print('WERID!',x,y)
                            #print("found weird pixel at ({},{})".format(x,y))
                            #print(x,y,maxslope,nominalslope,nominalslopeerr,maxslopeloc,zd/ngroup)#,minslopeloc)
                            #print(slopes)
                            #print('==========================')
                            if iii % 100 == 0:
                                #print(x,y)
                                f,a = plt.subplots()
                                a.plot(np.arange(zd),pix,'ro')
                                a.set_xlabel('Group')
                                a.set_ylabel('Signal (DN)')
                                a.set_title('Pixel ('+str(x)+','+str(y)+')')
                                f.savefig('weird_pix/weirdpix_'+str(x)+'_'+str(y)+'.pdf')
                                plt.close(f)
                                #stop
                            iii = iii + 1
                        else:
                            #print('GOODPIX')
                            pass
                #stop
        #stop
        #save for code testing
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(weird)
        h = fits.HDUList([h0,h1])
        h.writeto('my_weirdmap.fits',clobber=True)
        #stop
        return weird


    def get_weird_from_karl(self,file):
        with fits.open(file) as h:
            data = h[1].data

        weird = data >= 64
        map = np.zeros((2048,2048))
        map[weird] = 1
        return map

    def import_weird(self,file):
        with fits.open(file) as h:
            data = h[1].data

        weird = data >= 1
        map = np.zeros((2048,2048))
        map[weird] = 1
        return map


    def run(self):
        
        #get the lists of input files
        darkfiles = self.get_filelist(self.darklist)
        flatfiles = self.get_filelist(self.flatlist)

        #for output filenames later
        slash = self.darklist.rfind('/')
        darklistdir = self.darklist[0:slash+1]
        darklistfileonly = self.darklist[slash+1:]
        slash = self.flatlist.rfind('/')
        flatlistdir = self.flatlist[0:slash+1]
        flatlistfileonly = self.flatlist[slash+1:]

        #run cosmic ray search on each input file, (darks and flats)------------------------
        dark_cr_maps = np.zeros((1,2048,2048))
        flat_cr_maps = np.zeros((1,2048,2048))

        #dark_cr_maps = np.zeros((len(darkfiles),2048,2048))  #--FOR CODE TESTING ONLY!!!!!!!!!!
        #print('CHANGE ME BEFORE RUNNING FOR REAL!!!!!!!!!!!!!!!!!!!!!!')
        #USE THE ORIGINAL LINE WITH (1,2048,2048) WHEN RUNNING FOR REAL

        if self.find_cr == True:
            print('Beginning CR search block...',datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

            if self.darkcrmapfiles == None:
                print('Beginning cosmic ray search on the darks...',datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

                for darkfile in darkfiles:
                    cmap = self.cr_search(darkfile,cr_sigma=self.cr_sigma,flag_neighbors=True)
                    dark_cr_maps = np.concatenate((dark_cr_maps,cmap),axis=0)
                dark_cr_maps = dark_cr_maps[1:,:,:]

                #now use all the maps together to filter out repeat CR hits, which are likely pixels
                #with bad linearity corrections that are being collected here as false positives. Just like
                #with the hot pixels, let's declare any pixel that is flagged as CR hit in more than 10% of the
                #input files to be CR-free.
                cr_tot = np.sum(dark_cr_maps,axis=0) / dark_cr_maps.shape[0]
                false_positives = cr_tot > 0.10
                h = fits.PrimaryHDU()
                hi = fits.ImageHDU(false_positives.astype(int))
                hdu = fits.HDUList([h,hi])
                darkfalse = 'CR_search_false_positives_darks.fits'
                hdu.writeto(darkfalse,clobber='True')
                print("False positives from running the CR search on the darks written to {}".darkfalse)


                for i in xrange(dark_cr_maps.shape[0]):
                    frame = dark_cr_maps[i,:,:]
                    frame[false_positives] = -1
                    dark_cr_maps[i,:,:] = frame

                if self.darkcrfile == None:
                    self.darkcrfile = 'CosmicRay_Map_from_'+self.darklist+'.fits'
                h = fits.PrimaryHDU()
                hi = fits.ImageHDU(dark_cr_maps)
                hdu = fits.HDUList([h,hi])
                hdu.writeto(self.darkcrfile,clobber='True')
                print("Cosmic Ray Maps from the dark current files written to {}".format(self.darkcrfile))

                #user-provided CR maps
            else:
                if 'fits' in self.darkcrmapfiles:
                    print("Reading in dark current CR maps from {}".format(self.darkcrmapfiles))
                    with fits.open(self.darkcrmapfiles) as h:
                        dark_cr_maps = h[1].data

                else:
                    print("Using dark current CR maps listed in {}".format(self.darkcrmapfiles))
                    darkcrmaps = []
                    with open(self.darkcrmapfiles) as f:
                        for line in f:
                            if len(line) > 3:
                                darkcrmaps.append(line.strip())
            
                    dark_cr_maps = np.zeros((1,2048,2048))
                    for dmap in darkcrmaps:
                        h = fits.open(dmap)
                        dark_cr_maps = np.concatenate((dark_cr_maps,h[1].data),axis=0)
                        h.close()
                    dark_cr_maps = dark_cr_maps[1:,:,:]



            if self.flatcrmapfiles == None:
                print('beginning cosmic ray search on the flats...',datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
                #flat_cr_maps = np.zeros((1,2048,2048))
                for flatfile in flatfiles:
                    cmap = self.cr_search(flatfile,cr_sigma=self.cr_sigma,flag_neighbors=True)
                    flat_cr_maps = np.concatenate((flat_cr_maps,cmap),axis=0)
                flat_cr_maps = flat_cr_maps[1:,:,:]

                #now use all the maps together to filter out repeat CR hits, which are likely pixels
                #with bad linearity corrections that are being collected here as false positives. Just like
                #with the hot pixels, let's declare any pixel that is flagged as CR hit in more than 10% of the
                #input files to be CR-free.
                cr_tot = np.sum(flat_cr_maps,axis=0) / flat_cr_maps.shape[0]
                false_positives = cr_tot > 0.10
                h = fits.PrimaryHDU()
                hi = fits.ImageHDU(false_positives.astype(int))
                hdu = fits.HDUList([h,hi])
                flatfalse = 'CR_search_false_positives_flats.fits'
                hdu.writeto(flatfalse,clobber='True')
                print("False positives from running the CR search on the flats written to {}".flatfalse)

                for i in xrange(flat_cr_maps.shape[0]):
                    frame = flat_cr_maps[i,:,:]
                    frame[false_positives] = -1
                    flat_cr_maps[i,:,:] = frame

                if self.flatcrfile == None:
                    self.flatcrfile = 'CosmicRay_Map_from_'+self.flatlist+'.fits'
                h = fits.PrimaryHDU()
                hi = fits.ImageHDU(flat_cr_maps)
                hdu = fits.HDUList([h,hi])
                hdu.writeto(self.flatcrfile,clobber='True')
                print("Cosmic Ray Maps from the flat field files written to {}".format(self.flatcrfile))

                #user-provided CR maps
            else:
                if 'fits' in self.flatcrmapfiles:
                    print("Reading in flat field CR maps from {}".format(self.flatcrmapfiles))
                    with fits.open(self.flatcrmapfiles) as h:
                        flat_cr_maps = h[1].data
                else:
                    print("Using flat field CR maps listed in {}".format(self.flatcrmapfiles))
                    flatcrmaps = []
                    with open(self.flatcrmapfiles) as f:
                        for line in f:
                            if len(line) > 3:
                                flatcrmaps.append(line.strip())
            
                    #flat_cr_maps = np.zeros((1,2048,2048))
                    for fmap in flatcrmaps:
                        h = fits.open(fmap)
                        h1 = h[1].data
                        if len(h1.shape) == 2:
                            h1 = np.expand_dims(h1,axis=0)
                        flat_cr_maps = np.concatenate((flat_cr_maps,h1),axis=0)
                        h.close()
                    flat_cr_maps = flat_cr_maps[1:,:,:]


        else:
            #use the DQ extensions of the input files to locate CR hits
            print("Creating CR maps from the DQ arrays of the input data.")
            if self.crval == None:
                self.crval = dqflags.pixel['JUMP_DET']
            for darkfile in darkfiles:
                crmap = self.convert_dq_to_cr(darkfile,self.crval)
                dark_cr_maps = np.concatenate((dark_cr_maps,crmap),axis=0)
            dark_cr_maps = dark_cr_maps[1:,:,:]

            #save the CR maps from the darks
            if self.darkcrfile == None:
                self.darkcrfile = 'CosmicRay_Maps_from_'+self.darklist+'.fits'
            h = fits.PrimaryHDU()
            hi = fits.ImageHDU(dark_cr_maps)
            hdu = fits.HDUList([h,hi])
            hdu.writeto(self.darkcrfile,clobber='True')
            print("Converted CR maps for the darks saved to: "+self.darkcrfile)

            for flatfile in flatfiles:
                crmap = self.convert_dq_to_cr(flatfile,self.crval)
                flat_cr_maps = np.concatenate((flat_cr_maps,crmap),axis=0)
            flat_cr_maps = flat_cr_maps[1:,:,:]

            #save the CR maps from the flats
            if self.flatcrfile == None:
                self.flatcrfile = 'CosmicRay_Maps_from_'+self.flatlist+'.fits'
            h = fits.PrimaryHDU()
            hi = fits.ImageHDU(flat_cr_maps)
            hdu = fits.HDUList([h,hi])
            hdu.writeto(self.flatcrfile,clobber='True')
            print("Converted CR maps for the flats saved to: "+self.flatcrfile)

        #end cosmic ray search----------------------------------------------------------


        #find pixels that are bad in the 0th read, a la Don Hall------------------------
        #My CR search is not senstive to CRs that hit between reset and the 0th read, so
        #no need to read in and use cosmic ray maps here. For this step we need raw files,
        #or at least, files without the superbias subtraction done. For the moment, let's
        #just generate filenames on the fly from the names of the calibrated files.
        if self.skip_bad0 != True:
            print('beginning bad-in-0th-read pixel search...',datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            darkuncalfiles = []
            for darkfile in darkfiles:
                u = darkfile.rfind('ramp')
                darkuncalfiles.append(darkfile[0:u]+'uncal.fits')

            badin0_map = self.badin0th(darkuncalfiles,sigma=self.bad0_sigma)

            #save the bad in 0 map:
            #save the final hot map
            h = fits.PrimaryHDU()
            hi = fits.ImageHDU(badin0_map)
            hdu = fits.HDUList([h,hi])
            if self.zerofile == None:
                self.zerofile = 'Bad_in_0Read_Mask_from_'+self.flatlist+'.fits'
            hdu.writeto(self.zerofile,clobber=True)
            print('Mask of pixels which are bad in the 0th group saved as {}'.format(self.zerofile))
        else:
            print(self.skip_bad0)
            print("Skipping bad in 0th search.")
            badin0_map = np.zeros((2048,2048))
        #end bad-in-0th-read search-----------------------------------------------------
        

        #find hot pixels. Feed in CR maps so that CR-impacted pixels are ignored--------
        #print('skipping hot pix search during code testing')
        if self.skip_hot != True:
            print('beginning hot pixel search....',datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            hotmaps = np.zeros((1,2048,2048))
            #print(dark_cr_maps.shape)
            #print('hot search turned off for code testing')
            for i,darkfile in enumerate(darkfiles):
                hmap = self.hot_search(darkfile,hot_sigfrac=self.hot_sigfrac,in_mask = dark_cr_maps[i,:,:])
                hotmaps = np.concatenate((hotmaps,hmap),axis=0)
            hotmaps = hotmaps[1:,:,:]
        
            #combine the hot maps from all the files into one final hot map
            hot_map = self.combine_ind_maps(hotmaps,minfrac = self.min_hot_frac)
        

            #create a table of how many pixels in each quadrant are hot
            stats_hot = self.table_of_badpix(hot_map)
            source = Column(['Hot'],name='Source')
            stats_hot.add_column(source,index=0)

            #save table
            hot_table_file = darklistdir+'Hot_stats_table_from_'+darklistfileonly+'.txt'
            stats_hot.write(hot_table_file,format='ascii')
            print("Hot pixel statistics saved to {}".format(hot_table_file))

            #save the final hot map
            h = fits.PrimaryHDU()
            hi = fits.ImageHDU(hot_map)
            hdu = fits.HDUList([h,hi])
            if self.hotfile == None:
                self.hotfile = darklistdir+'Hot_Pixel_Mask_from_'+darklistfileonly+'.fits'
            hdu.writeto(self.hotfile,clobber=True)
            print('Hot pixel mask saved as {}'.format(self.hotfile))
        else:
            print(self.skip_hot)
            print("Skipping hot pixel search.")
            hot_map = np.zeros((2048,2048))
        #end hot pixel search---------------------------------------------------------

        
        #find dead pixels in the flats, feed in CR maps again-------------------------
        if self.skip_dead != True:
            print('Beginning dead pixel search...',datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            flat_dead_maps = np.zeros((len(flatfiles),2048,2048))
            flat_low_qe_maps = np.zeros((len(flatfiles),2048,2048))
            #print('dead search turned off for code testing.')
            for i,flatfile in enumerate(flatfiles):
                uncal = flatfile.rfind('ramp')
                uncalfile = flatfile[0:uncal] + 'uncal.fits'
                print('dead search: ',uncalfile)
                dmap,lowqemap = self.dead_search(uncalfile,dead_sigfrac=self.dead_sigfrac,in_mask = flat_cr_maps[i,:,:])
                #flat_dead_maps = np.concatenate((flat_dead_maps,dmap),axis=0)
                #flat_low_qe_maps = np.concatenate((flat_low_qe_maps,lowqemap),axis=0)
                flat_dead_maps[i,:,:] = dmap
                flat_low_qe_maps[i,:,:] = lowqemap
            #flat_dead_maps = flat_dead_maps[1:,:,:]
            #flat_low_qe_maps = flat_low_qe_maps[1,:,:,:]
        

            for i in xrange(flat_dead_maps.shape[0]):
                nnnbad = len(np.where(flat_low_qe_maps[i,:,:] != 0)[0])
                nnngood = len(np.where(flat_low_qe_maps[i,:,:] == 0)[0])
                print("number of non-zero pix in low QE map {}: {}".format(i,nnnbad))
                print("number of good pix: {}".format(nnngood))

            #combine the dead maps from all the files into one final dead map
            dead_map = self.combine_ind_maps(flat_dead_maps,minfrac = self.min_dead_frac)
            lowqe_map = self.combine_ind_maps(flat_low_qe_maps,minfrac = self.min_dead_frac)
            nlow = len(np.where(lowqe_map != 0)[0])
            print("Number of low qe pix is {}".format(nlow))
            nlow = len(np.where(dead_map != 0)[0])
            print("Number of dead pix is {}".format(nlow))

            #create a table of how many pixels in each quadrant have low QE
            stats_lowqe = self.table_of_badpix(lowqe_map)
            source = Column(['LowQE'],name='Source')
            stats_lowqe.add_column(source,index=0)

            #save table
            slash = self.flatlist.rfind('/')
            flatlistdir = self.flatlist[0:slash+1]
            flatlistfileonly = self.flatlist[slash+1:]
            lowqe_table_file = flatlistdir+'LowQE_stats_table_from_'+flatlistfileonly+'.txt'
            stats_lowqe.write(lowqe_table_file,format='ascii')        
            print("Low QE pixel statistics saved to {}".format(lowqe_table_file))

            #create a table of how many pixels in each quadrant are dead
            stats_dead = self.table_of_badpix(dead_map)
            source = Column(['Dead'],name='Source')
            stats_dead.add_column(source,index=0)
            
            #save table
            dead_table_file = flatlistdir+'Dead_stats_table_from_'+flatlistfileonly+'.txt'
            stats_dead.write(dead_table_file,format='ascii')
            print("Dead pixel statistics saved to {}".format(dead_table_file))

            #save the final dead map
            h = fits.PrimaryHDU()
            hi = fits.ImageHDU(dead_map)
            hdu = fits.HDUList([h,hi])
            if self.deadfile == None:
                self.deadfile = 'Dead_Pixel_Mask_from_'+self.flatlist+'.fits'
            hdu.writeto(self.deadfile,clobber=True)
            print('Dead pixel mask saved as {}'.format(self.deadfile))

            #save the low QE map
            h = fits.PrimaryHDU()
            hi = fits.ImageHDU(lowqe_map)
            hdu = fits.HDUList([h,hi])
            if self.lowqefile == None:
                self.lowqefile = 'LowQE_Pixel_Mask_from_'+self.flatlist+'.fits'
            hdu.writeto(self.lowqefile,clobber=True)
            print('Low QE pixel mask saved as {}'.format(self.lowqefile))
        else:
            print(self.skip_dead)
            print("Skipping dead pixel search.")
            dead_map = np.zeros((2048,2048))
            lowqe_map = np.zeros((2048,2048))
        #end dead pixel search----------------------------------------------------------

        #search for Karl's weird pixels, which have a sudden slope increase just before saturation------
        flatuncalfiles = []
        for ff in flatfiles:
            u = ff.rfind('ramp')
            flatuncalfiles.append(ff[0:u]+'uncal.fits')


        if self.skip_weird == False:
            if self.karlweird == None:
                print("Beginning search for 'weird' pixels...",datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

                with fits.open(flatuncalfiles[0]) as h:
                    wdata = h[1].data
                wdata = wdata[0,:,:,:]

                weird_map = self.weird_search(wdata,nsigma=self.weirdsigma)

            else:  #if Karl's BPM name is given, use those weird pix instead of searching for them
                print('Skipping the search for weird pixels. Taking the weird pixels from '+self.karlweird+' instead.')
                weird_map = self.get_weird_from_karl(self.karlweird)
        else:
            if self.karlweird == None:
                weird_map = np.zeros((2048,2048))
            else:
                print('Skipping the search for weird pixels. Taking the weird pixels from '+self.karlweird+' instead.')
                weird_map = self.get_weird_from_karl(self.karlweird)

        #create a table of how many pixels in each quadrant are weird
        stats_weird = self.table_of_badpix(weird_map)
        source = Column(['Weird'],name='Source')
        stats_weird.add_column(source,index=0)

        #save table
        base = flatuncalfiles[0]
        slash = base.rfind('/')
        basedir = base[0:slash+1]
        basefileonly = base[slash+1:]
        weird_table_file = basedir+'Weird_stats_table_from_'+basefileonly+'.txt'
        stats_weird.write(weird_table_file,format='ascii')        
        print("Weird pixel statistics saved to {}".format(weird_table_file))

        #save the final weird map
        h = fits.PrimaryHDU()
        hi = fits.ImageHDU(weird_map)
        hdu = fits.HDUList([h,hi])
        if self.weirdfile == None:
            self.weirdfile = 'Weird_Pixel_Mask_from_'+basefileonly+'.fits'
        hdu.writeto(self.weirdfile,clobber=True)
        print('Weird pixel mask saved as {}'.format(self.weirdfile))
        #end weird pixel search-------------------------------------------------------------------------



        #find unstable pixels-----------------------------------------------------------
        if self.skip_unstable == False:
            print('beginning unstable pixel search...',datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            unstable_darkmap, unstable_flatmap = self.unstable_search(darkfiles,flatfiles,dark_cr_maps,flat_cr_maps)
            all_unstable = unstable_darkmap + unstable_flatmap
            all_unstable[all_unstable > 1] = 1

            #get a basic header so you can add to it for the unstable file extensions
            with fits.open(darkfiles[0]) as ll:
                head0 = ll[0].header
                head1 = ll[1].header

        

            #save the unstable pixel map
            h = fits.PrimaryHDU()
            head1['EXTNAME'] = 'UNS_DARK'
            hi = fits.ImageHDU(unstable_darkmap,head1)
            head1['EXTNAME'] = 'UNS_FLAT'
            hii = fits.ImageHDU(unstable_flatmap,head1)
            head1['EXTNAME'] = 'UNS_BOTH'
            hiii = fits.ImageHDU(all_unstable,head1)
            hdu = fits.HDUList([h,hi,hii,hiii])
            if self.unstablefile == None:
                self.unstablefile = darklistdir+'Unstable_Pixel_Mask_from_'+darklistfileonly+'_and_'+flatlistfileonly+'.fits'
            hdu.writeto(self.unstablefile,clobber=True)
            print('Unstable pixel mask saved as {}'.format(self.unstablefile))
        else:
            print("Skipping unstable pixel search.")
            unstable_darkmap = np.zeros((2048,2048))
            unstable_flatmap = np.zeros((2048,2048))
            all_unstable = np.zeros((2048,2048))
        #end unstable pixel search------------------------------------------------------

        #combine the various maps into a single bad pixel mask--------------------------

        #SET UP DQ DEFINITIONS, FOLLOWING SSB CONVENTIONS
        #WE NEED TO CREATE *TWO* OUTPUT FILES HERE. ONE IS THE STATIC BAD PIXEL MASK,
        #AND THE OTHER IS THE DQ ARRAY WHICH WILL GO INTO THE DARK CURRENT REFERENCE PIXELS.
        #THIS IS BECAUSE SSB SAYS THAT EACH REFERENCE FILE CONTAINS ONLY THE TYPES OF BAD PIXELS
        #RELEVANT TO IT.

        #So, unstable pixels, which will be separated by unstable in the flat and unstable in the dark, will be 
        #renamed unreliable_slope and unreliable_dark, respectively, and be put into the dark current dq array?
        #(I don't really agree with this.....)

        #Weird pixels will be recast as NONLINEAR in the non-lin reference file?? This doesn't seem especially appropriate either

 
        #create dictionaries that will be used to create dq_def extension definitions: static bad pixel mask and dark current dq array
        #unstable pix are re-cast as unreliable dark and unreliable slope depending on how they were detected. weird pix, and bad_0group
        #pixels are currently kept with those mnemonics. That means the flags will be thrown out if the bpm is read in via MaskModel
        #(as is done in the pipeline), but kept if the file is read in via astropy.
        dqdef_bpm = []
        dqssb_bpm={'DO_NOT_USE':np.uint8(1),'NON_SCIENCE':np.uint8(2),'DEAD':np.uint8(4),'LOW_QE':np.uint8(8),'NO_GAIN_VALUE':np.uint8(16),'BAD_0GROUP':np.uint8(64),'WEIRD':np.uint8(128)}
        dqdef_dark = []
        dqssb_dark={'DO_NOT_USE':np.uint8(1),'HOT':np.uint8(2),'WARM':np.uint8(4),'UNRELIABLE_DARK':np.uint8(8),'UNRELIABLE_SLOPE':np.uint8(16)}

        #now create the dq_def data
        for bitname in dqssb_bpm:
            bitvalue = dqssb_bpm[bitname]
            bitnumber = int(np.log(bitvalue)/np.log(2))
            newrow = (bitnumber,bitvalue,bitname,'')
            dqdef_bpm.append(newrow)

        for bitname in dqssb_dark:
            bitvalue = dqssb_dark[bitname]
            bitnumber = int(np.log(bitvalue)/np.log(2))
            newrow = (bitnumber,bitvalue,bitname,'')
            dqdef_dark.append(newrow)

        #now combine the individual bad pixel maps into the final maps and save
        #print(badin0_map.shape,dead_map.shape,hot_map.shape,lowqe_map.shape,unstable_map.shape)
        #print(bad0val,deadval,lowqeval,hotval,unstabval)

       
        xxb,yyb = 1650,160
        xxd,yyd = 1296,959
        print(badin0_map[yyb,xxb],dead_map[yyb,xxb],lowqe_map[yyb,xxb],weird_map[yyb,xxb])
        print(hot_map[yyd,xxd],unstable_darkmap[yyd,xxd],unstable_flatmap[yyd,xxd])

        badmapb = np.zeros((8,2048,2048))
        badmapb[7-np.log2(dqssb_bpm['BAD_0GROUP']),:,:] = np.uint8(badin0_map)
        badmapb[7-np.log2(dqssb_bpm['DEAD']),:,:] = np.uint8(dead_map)
        badmapb[7-np.log2(dqssb_bpm['LOW_QE']),:,:] = np.uint8(lowqe_map)
        badmapb[7-np.log2(dqssb_bpm['WEIRD']),:,:] = np.uint8(weird_map)

        darkmapb = np.zeros((8,2048,2048))
        darkmapb[7-np.log2(dqssb_dark['HOT']),:,:] = np.uint8(hot_map)
        darkmapb[7-np.log2(dqssb_dark['UNRELIABLE_DARK']),:,:] = np.uint8(unstable_darkmap)
        darkmapb[7-np.log2(dqssb_dark['UNRELIABLE_SLOPE']),:,:] = np.uint8(unstable_flatmap)


        #for code development------
        checky,checkx = 598,314
        tempbadmap = np.packbits(np.int8(badmapb),axis=0)
        tempdarkmap = np.packbits(np.int8(darkmapb),axis=0)
        print("After adding my values, before importing anything from Karl: ",np.max(tempbadmap))
        print(tempbadmap.shape)
        print(np.max(tempdarkmap))
        print(tempbadmap[0,checky,checkx])
        #for code development------


        #badmap = np.uint8(badin0_map)*dqssb_bpm['BAD_0GROUP'] + np.uint8(dead_map)*(dqssb_bpm['DEAD']+dqssb_bpm['DO_NOT_USE']) + np.uint8(lowqe_map)*(dqssb_bpm['LOW_QE']+dqssb_bpm['DO_NOT_USE']) + np.uint8(weird_map)*dqssb_bpm['WEIRD']
        #darkmap = np.uint(hot_map)*dqssb_dark['HOT'] + np.uint8(unstable_darkmap)*dqssb_dark['UNRELIABLE_DARK'] + np.uint8(unstable_flatmap)*(dqssb_dark['UNRELIABLE_SLOPE']+dqssb_dark['DO_NOT_USE'])

        
        #if requested, pull over desired bad pixel types from one of Karl's bad pixel masks
        #I'm thinking abot RC pixels, which this script doesn't check for.
        if self.Karl_file != None:
            print("Importing bad pixel mask bit values {} from {}.".format(self.Karl_vals_to_keep,self.Karl_file))

            vals_to_keep = np.array(self.Karl_vals_to_keep)
            karl_mask, karl_dict = self.import_from_Karl(self.Karl_file,vals_to_keep)

            print(karl_mask[0,yyb,xxb],karl_mask[0,yyd,xxd])


            #unpack karl's extracted bpm
            #karl_maskb = np.expand_dims(karl_mask,axis=0)
            karl_maskb = np.unpackbits(np.uint8(karl_mask),axis=0)

            #unpack the SSB bpm and dark dq arrays
            #badmapb = np.expand_dims(np.uint8(badmap),axis=0)
            #badmapb = np.unpackbits(badmapb,axis=0)
            #darkmapb = np.expand_dims(np.uint8(darkmap),axis=0)
            #darkmapb = np.unpackbits(darkmapb,axis=0)
 
            #now go through Karl's dictionary and equate his types of bad pixels
            #to those in SSB, and change the bit numbers accordingly
            #Keep the scope of this limited for now...
            for key in karl_dict:
                #get the bit number from the bit value
                badbit = np.log2(np.uint8(key))
                #extract the appropriate plane
                kplane = karl_maskb[7-badbit,:,:]
                #if the plane has bad pix, continue
                if np.max(kplane == 1):
                    #get the type of bad pixel
                    badtype = karl_dict[key]
                    #print("working on "+badtype+" with badbit equal to {}".format(badbit))
                    #print("before any work: ",dqdef_bpm)

                    #check to see if badtype is already in dqdef_bpm:
                    present_bpm = False
                    for s in dqdef_bpm:
                        if badtype in s:
                            present_bpm = True
                            break

                    present_dark = False
                    for s in dqdef_dark:
                        if badtype in s:
                            present_dark = True
                            break


                    #if the type matches a type already defined in SSB
                    #if badtype in dqdef_bpm:
                    if present_bpm == True:

                        print("{} is already in the bad pixel mask definition. Adding these pixels in.".format(badtype))
                        
                        #get the ssb bit number
                        newbit = np.log2(np.uint8(dqssb_bpm[badtype]))
                    
                        #add karl's bitplane into the appropriate ssb bit plane
                        tmpplane = badmapb[7-newbit,:,:] 
                        tmpplane = tmpplane + kplane
                        tmpplane[tmpplane > 1] = 1
                        badmapb[7-newbit,:,:] = tmpplane


                    #if the type matches a type in the dark dq_def, then add the plane in there
                    #(with the exception of DO_NOT_USE, which we will put only in the bad pix mask)
                    if ((badtype != 'DO_NOT_USE') and (present_dark)): #(badtype in dqdef_dark)):
                    
                        print("{} is already in the dark reffile DQ definition. Adding these pixels in.".format(badtype))

                        #get the ssb bit number
                        newbit = np.log2(np.uint8(dqssb_dark[badtype]))
                    
                        #add karl's bitplane into the appropriate ssb bit plane
                        tmpplane = darkmapb[7-newbit,:,:] 
                        tmpplane = tmpplane + kplane
                        tmpplane[tmpplane > 1] = 1
                        darkmapb[7-newbit,:,:] = tmpplane
                       
                    #if the type is not present in the SSB bad pixel mask nor the dark dq array
                    #then add it to the bad pixel mask with the next unoccupied bit value
                    #if ((badtype not in dqdef_bpm) and (badtype not in dqdef_dark)):
                    if ((present_bpm == False) & (present_dark == False)):
   
                        #find the bit value to use
                        assigned = np.sort(np.array([s[0] for s in dqdef_bpm]))
                        comp = np.arange(8)
                        for s in comp:
                            if s not in assigned:
                                newbit = s
                                break
                        print("{} is not present in the bad pixel mask nor dark current DQ definition.".format(badtype))
                        print("Adding to the bad pixel mask with a bit value of {}".format(2**newbit))

                        if newbit > 7:
                            print("WARNING: CANNOT IMPORT {} FLAGS FROM {}, AS THE FINAL BAD PIXEL MASK".format(badtype,self.Karl_file))
                            print("MUST CONTAIN ONLY 8-BIT INTEGERS, AND ALL 8 BITS ARE ALREADY DEFINED.")
                            print("Continuing, but skipping the import.")

                        else:
                            #add karl's bitplane into the appropriate ssb bit plane
                            tmpplane = badmapb[7-newbit,:,:] 
                            tmpplane = tmpplane + kplane
                            tmpplane[tmpplane > 1] = 1


                            #check
                            #diff = tmpplane-badmapb[7-newbit,:,:]
                            #print(len(np.where(diff != 0)[0]))
                            #check = np.where(diff != 0)
                            #checky,checkx = check[0][2],check[1][2]
                            #print(checkx,checky)

                            badmapb[7-newbit,:,:] = tmpplane
                            #add the new definition to dqdef_bpm
                            bval = np.uint8(2**newbit)

                            #print('bval is {}'.format(bval))

                            newrow = (newbit,bval,badtype,'')
                            dqdef_bpm.append(newrow)
                            #also add to the definition dictionary
                            dqssb_bpm[badtype] = bval

                else: #plane has no bad pix to convert
                    pass


        #for code development------
        tempbadmap = np.packbits(np.int8(badmapb),axis=0)
        tempdarkmap = np.packbits(np.int8(darkmapb),axis=0)
        print("After importing anything from Karl: ",np.max(tempbadmap))
        print(np.max(tempdarkmap))
        print(tempbadmap[0,checky,checkx])
        #for code development------




        #add the DO_NOT_USE flag where requested
        badmapb = self.add_donotuse(badmapb,dqdef_bpm)
        darkmapb = self.add_donotuse(darkmapb,dqdef_dark)

        #pack the bits back up
        darkmap = np.packbits(np.uint8(darkmapb),axis=0)
        darkmap = darkmap[0,:,:]
        badmap = np.packbits(np.uint8(badmapb),axis=0)
        badmap = badmap[0,:,:]

        print('imports from karl done,bits packed back up, do not use added',badmap[yyb,xxb],darkmap[yyd,xxd])

        #for code development------
        tempbadmap = np.packbits(np.int8(badmapb),axis=0)
        tempdarkmap = np.packbits(np.int8(darkmapb),axis=0)
        print("After adding donotuse values: ",np.max(tempbadmap))
        print(np.max(tempdarkmap))
        print(tempbadmap[0,checky,checkx])
        #for code development------



        #flag the reference pixels in the static pixel mask
        badmap[0:4,:] = badmap[0:4,:] + dqssb_bpm['NON_SCIENCE']
        badmap[2044:,:] = badmap[2044:,:] + dqssb_bpm['NON_SCIENCE']
        badmap[4:2044,0:4] = badmap[4:2044,0:4] + dqssb_bpm['NON_SCIENCE']
        badmap[4:2044,2044:] = badmap[4:2044,2044:] + dqssb_bpm['NON_SCIENCE']


        #for code development------
        tempbadmap = np.uint8(badmap)
        tempdarkmap = darkmap
        print("After flagging refpix: ",np.max(tempbadmap))
        print(np.max(tempdarkmap))
        print(tempbadmap[checky,checkx])
        #save manually as a fits file to see if odd large values
        #are still present
        h=fits.PrimaryHDU()
        h1=fits.ImageHDU(badmap)
        hdulist = fits.HDUList([h,h1])
        hdulist.writeto('tempbadmap.fits',clobber=True)
        #for code development------



        #print('after importing RC from karl')
        #print(badmapb[:,yyb,xxb],darkmapb[:,yyd,xxd])
        #print(badmap[yyb,xxb],darkmap[yyd,xxd])
        
        #check
        #diff = badmap2 - badmap
        #print(np.max(diff),np.min(diff))
        #print(badmap[checky,checkx])
        #print(karl_mask[0,checky,checkx])
        #print(badmap2[checky,checkx])
        #stop


        if self.outfile == None:
            self.outfile = darklistdir+'BadPixMask_'+darklistfileonly+'_and_'+flatlistfileonly+'.fits'

        #save files using JWST MaskModel
        #create empty instance of MaskModel and then populate via attributes
        #to get around current bug in MaskModel that pops up when creating
        #an instance of MaskModel and populating via arguments.
        finalmask = MaskModel()
        darkmask = MaskModel()
        finalmask.dq = np.uint8(badmap)
        finalmask.dq_def = dqdef_bpm
        darkmask.dq = np.uint8(darkmap)
        darkmask.dq_def = dqdef_dark

        #grab info for header
        detector = fits.getval(darkfiles[0],'DETECTOR')

        #set appropriate header values
        finalmask = self.make_proper_header(finalmask,detector,darkfiles,flatfiles)
        darkmask = self.make_proper_header(darkmask,detector,darkfiles,flatfiles)

        #save masks
        finalmask.save(self.outfile)
        darkmask_outfile = self.outfile[0:-5] + '_DARKREFDQ.fits'
        darkmask.save(darkmask_outfile)


        #fits file format checks
        print("Quick check of fits file format. Printing file info, and reading in DQ_DEF extension")
        check_ssb = fits.open(self.outfile)
        print(check_ssb.info())
        print(check_ssb['DQ_DEF'].data)
        baddq = check_ssb[1].data
        print('BPM pixel check',baddq[yyb,xxb])

        check_ssb = fits.open(darkmask_outfile)
        print(check_ssb.info())
        print(check_ssb['DQ_DEF'].data)
        baddq = check_ssb[1].data
        print('DARK DQ PIXEL CHECK',baddq[yyb,xxb])


        #redcat team checks
        print("Running fitsverify on bad pixel mask and dark reffile DQ array file")
        subprocess.call(['fitsverify',self.outfile])
        subprocess.call(['fitsverify',darkmask_outfile])


        #insert badpix definitions into header
        #with fits.open(flatuncalfiles[0]) as hh:
        #    head = hh[0].header
        #
        #head['bad0val'] = bad0val
        #head['deadval'] = deadval
        #head['lowqeval'] = lowqeval
        #head['hotval'] = hotval
        #head['unstable'] = unstabval
        #head['donotuse'] = 1
        #head['refpxval'] = refpixval
        #head['weirdval'] = weirdval
        #
        #h = fits.PrimaryHDU()
        #hi = fits.ImageHDU(badmap,head)
        #hdu = fits.HDUList([h,hi])
        #hdu.writeto(self.outfile,clobber=True)
        print('Final bad pixel mask saved as {}'.format(self.outfile))
        print('Accompanying DQ values to be put into dark current reference file saved to {}'.format(darkmask_outfile))
        
        
    def add_donotuse(self,map,dqdef):
        '''
        Insert the DO_NOT_USE flag for pixels which are flagged as bad in the required ways.
        At the moment the DO_NOT_USE flag is set for pix which are flagged as DEAD, LOW_QE, and UNRELIABLE_SLOPE
        '''
        #newrow = (bitnumber,bitvalue,bitname,'')
        #example row in dqdef: (1, 2, 'DEAD','')
        dnu_list = ['DEAD','LOW_QE','UNRELIABLE_SLOPE','RC']
        dnu = map[7,:,:]

        for row in dqdef:
            badtype = row[2]
            if badtype in dnu_list:
                source = map[7-row[0],:,:]
                dnu[source == 1] = np.uint8(1)

        map[7,:,:] = dnu
        return map

        

        #badmap = np.uint8(badin0_map)*dqssb_bpm['BAD_0GROUP'] + np.uint8(dead_map)*(dqssb_bpm['DEAD']+dqssb_bpm['DO_NOT_USE']) + np.uint8(lowqe_map)*(dqssb_bpm['LOW_QE']+dqssb_bpm['DO_NOT_USE']) + np.uint8(weird_map)*dqssb_bpm['WEIRD']
        #darkmap = np.uint(hot_map)*dqssb_dark['HOT'] + np.uint8(unstable_darkmap)*dqssb_dark['UNRELIABLE_DARK'] + np.uint8(unstable_flatmap)*(dqssb_dark['UNRELIABLE_SLOPE']+dqssb_dark['DO_NOT_USE'])
        pass

    
    def make_proper_header(self,model,detector,darkfiles,flatfiles):
        #set header keywords that are specific to the MaskModel
        model.meta.reffile.type = 'MASK'
        model.meta.subarray.name = 'FULL'
        model.meta.subarray.xstart = 1
        model.meta.subarray.xsize = 2048
        model.meta.subarray.ystart = 1
        model.meta.subarray.ysize = 2048
        model.meta.instrument.name = 'NIRCAM'
        if detector == 'NRCA5':
            detector = 'NRCALONG'
        if detector == 'NRCB5':
            detector = 'NRCBLONG'
        model.meta.instrument.detector = str(detector)


        #look for the fastaxis and slowaxis keywords in the input data.
        #if they are present propogate those values into the bad pixel
        #mask. If they are not present, then you must be working with 
        #native orientation data, so use the appropriate values
        inhdu = fits.open(darkfiles[0])
        try:
            model.meta.subarray.fastaxis = inhdu[0].header['FASTAXIS']
            model.meta.subarray.slowaxis = inhdu[0].header['SLOWAXIS']
        except KeyError:
            print('===============================================')
            print("FASTAXIS and SLOWAXIS header keywords not found in the input data.")
            print("Assuming they are in native (fitswriter) orientation, and adding the")
            print("native orientation values for those keywords to the static pixel mask.")
            print('===============================================')
            model.meta.subarray.fastaxis = 1
            model.meta.subarray.slowaxis = 2
            
        model.meta.reffile.author = 'Hilbert'
        model.meta.reffile.description = 'Bad Pixel Mask reffile from CV3 data'
        model.meta.reffile.pedigree = 'GROUND'
        model.meta.reffile.useafter = '2015-10-01'

        #HISTORY keyword
        model.history.append('Description of Reference File Creation')

        model.history.append('DOCUMENT:')
        model.history.append('JWST-STScI-TR-XXXX')

        model.history.append('SOFTWARE:')
        model.history.append('/ifs/jwst/wit/witserv/data4/nrc/')
        model.history.append('hilbert/bad_pixel_map/cv3/badpix_maps/B1/nircam_find_badpix_cv3.py')
        
        #put the list of input files into the HISTORY keyword
        model.history.append('DATA USED:')
        for file in darkfiles+flatfiles:
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


    def import_from_Karl(self,file,dqval):
        #import bad pix of the specified flavor from one of Karl's bad pixel masks
        #designed so that you specify the dq value you want to bring over

        #get data
        with fits.open(file) as h:
            dq = h[1].data
            header = h[0].header
        #dq = np.array(dq,dtype=np.uint8)
        dqb = np.expand_dims(dq,axis=0)

        #translate dq map into binary
        dqb = np.unpackbits(np.uint8(dqb),axis=0)

        #translate the requested dq value to bit number
        dqvalb = np.array([np.unpackbits(np.uint8(s)) for s in dqval])
        badbitplanes = [np.where(s==0)[0] for s in dqvalb]

        if len(badbitplanes) == 1:
            intersect = badbitplanes[0]
        else:
            intersect = list(set(badbitplanes[0]) & set(badbitplanes[1]))
        if len(badbitplanes) > 2:
            for i in xrange(2,len(badbitplanes)):
                intersect = list(set(intersect) & set(badbitplanes[i]))
            

        #zero out all pixels in the bit planes that you are not 
        #interested in
        for plane in intersect:
            dqb[plane,:,:] = 0
            
        #re-pack the bits to form your extracted mask
        newdq = np.packbits(dqb,axis=0)
        
        #get Karl's header values to go with the mask
        karl_dq_def = self.get_karl_dq_defs(header)

        return newdq,karl_dq_def


    def get_karl_dq_defs(self,header):
        #extract the header keywords from Karl's bad pixel file that 
        #definte which bits correspond to which types of bad pixels.
        dict = {}

        for i in xrange(7):
            planestr = str(2**i)
            hdrkwd = 'DQ_'+planestr
            dict[2**i] = header[hdrkwd]

        return dict




    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description="Generate a NIRCam bad pixel mask from a set of input files")
        parser.add_argument("darklist",help="File containing a list of dark current integrations to use.")
        parser.add_argument("flatlist",help="File containing a list of flat field integrations to use.")
        parser.add_argument("-v","--verbose",help="detailed outputs displayed",action="store_true",default="False")
        parser.add_argument("-d","--outdir",help="output directory",default='./')
        parser.add_argument("-o","--outfile",help="file to save bias frame and uncertainty to.",default=None)
        parser.add_argument("-s","--saveplots",help="save output plots",action="store_true",default=False)
        parser.add_argument("--cr_sigma",help="Sigma level to use when searching for CR hits.",default=5)
        parser.add_argument("--hot_sigfrac",help="Multiple of local surrounding mean signal above which a pixel is flagged as hot",default=10)
        parser.add_argument("--dead_sigfrac",help="Fraction of local surrounding mean signal below which a pixel is flagged as dead",default=0.5)
        parser.add_argument("--min_dead_frac",help="Pixels found to be dead in more than min_dead_frac fraction of the input integrations are defined as dead",default=0.4)
        parser.add_argument("--min_hot_frac",help="Pixels found to be hot in more than min_hot_frac fraction of the input integrations are defined as hot",default=0.4)
        parser.add_argument("--bad0_sigma",help="Sigma value to use for clipping when looking for pixels that are bad in the 0th read.",default=5)
        parser.add_argument("--hotfile",help="Name of file to hold final hot pixel map.",default=None)
        parser.add_argument("--deadfile",help="Name of file to hold final dead pixel map.",default=None)
        parser.add_argument("--lowqefile",help="Name of file to hold final low QE pixel map.",default=None)
        parser.add_argument("--unstablefile",help="Name of file to hold final unstable pixel map.",default=None)
        parser.add_argument("--zerofile",help="Name of file to hold final bad-in-0th-group pixel map.",default=None)
        parser.add_argument("--flatmeannoisethresh",help="Number of sigma to use as cutoff when looking for unstable pix in mean flat",default=5)  #5 arrived at through examination of data
        parser.add_argument("--darkmeannoisethresh",help="Number of sigma to use as cutoff when looking for unstable pix in mean dark",default=9)  #9 arrived at through examination of data
        parser.add_argument("--flat_signal_thresh",help="Number of sigma to use as cutoff when looking for unstable pix in individual flats",default=3)
        parser.add_argument("--dark_signal_thresh",help="Number of sigma to use as cutoff when looking for unstable pix in individual darks",default=9)
        #parser.add_argument("--badpix_defs",help="Text file containing the definition for bad pixel types. (i.e. dead pix = 4, hot pix = 16, etc)",default=None)
        parser.add_argument("--flatcrfile",help="Name of file to save the CR maps from the flat field inputs within.",default=None)
        parser.add_argument("--darkcrfile",help="Name of file to save the CR maps from the dark current inputs within.",default=None)
        parser.add_argument("--darkslopes",help="Name of file that contains a list of slopefiles for the dark current ramps.",default=None)
        parser.add_argument("--flatslopes",help="Name of file that contains a list of slopefiles for the flat field ramps.",default=None)
        parser.add_argument("--flatcrmapfiles",help="Name of a list file that contains cosmic ray maps for all of the flat field input files.",default=None)
        parser.add_argument("--darkcrmapfiles",help="Name of a list file that contains cosmic ray maps for all of the dark current input files.",default=None)
        parser.add_argument("--skip_hot",help="Skip hot pixel search.",action='store_true',default=False)
        parser.add_argument("--skip_dead",help="Skip dead pixel search",action='store_true',default=False)
        parser.add_argument("--skip_bad0",help="Skip bad in 0th search",action='store_true',default=False)
        parser.add_argument("--skip_cr",help="Skip CR search",action='store_true',default=False)
        parser.add_argument("--skip_weird",help="Skip weird pixel search",action='store_true',default=False)
        parser.add_argument("--skip_unstable",help="Skip unstable pixel search",action='store_true',default=False)
        parser.add_argument("--crval",help="Bit value corresponding to a CR hit. Needed while the DQ definition extension is not in the data.",default=None)
        parser.add_argument("--find_cr",help="If set, the scripts attempts its own search for CRs. Implemented before the SSB pipeline was doing this.",action='store_true',default=False)
        parser.add_argument("--weirdsigma",help="Sigma value above the nominal mean slope that an elevated slope at the end of the ramp will cause that pixel to be flagged as weird.",default=5)
        parser.add_argument("--weirdfile",help="Filename to save the map of weird pixels into.",default=None)
        parser.add_argument("--karlweird",help="Filename of Karl's bad pixel map, from which the weird pixels will be taken and used.")
        parser.add_argument("--Karl_file",help="Filename of Karl's bad pixel map from which to import bad pixels.",default=None)
        parser.add_argument("--Karl_vals_to_keep",help="DQ flag numbers to import from Karl's bad pixel file",default=[],type=int,nargs='+')
        return parser

if __name__ == '__main__':
    usagestring = 'python nircam_find_badpix.py input_darks.list input_flats.list'

    badpix = Badpix_Map()
    parser = badpix.add_options(usage=usagestring)
    args = parser.parse_args()

    args = parser.parse_args()

    badpix.darklist = args.darklist
    badpix.flatlist = args.flatlist
    badpix.verbose = args.verbose
    badpix.outdir = args.outdir
    badpix.outfile = args.outfile
    badpix.saveplots = args.saveplots
    badpix.cr_sigma = args.cr_sigma
    badpix.hot_sigfrac = int(args.hot_sigfrac)
    badpix.dead_sigfrac = float(args.dead_sigfrac)
    badpix.min_hot_frac = args.min_hot_frac
    badpix.min_dead_frac = args.min_dead_frac
    badpix.bad0_sigma = args.bad0_sigma
    badpix.hotfile = args.hotfile
    badpix.deadfile = args.deadfile
    badpix.lowqefile = args.lowqefile
    badpix.unstablefile = args.unstablefile
    badpix.zerofile = args.zerofile
    badpix.flatmeannoisethresh = args.flatmeannoisethresh
    badpix.darkmeannoisethresh = args.darkmeannoisethresh
    badpix.flat_signal_thresh = args.flat_signal_thresh
    badpix.dark_signal_thresh = args.dark_signal_thresh
    ##badpix.badpix_defs = args.badpix_defs
    badpix.flatcrfile = args.flatcrfile
    badpix.darkcrfile = args.darkcrfile
    badpix.darkslopes = args.darkslopes
    badpix.flatslopes = args.flatslopes
    badpix.flatcrmapfiles = args.flatcrmapfiles
    badpix.darkcrmapfiles = args.darkcrmapfiles
    badpix.skip_hot = args.skip_hot
    badpix.skip_dead = args.skip_dead
    badpix.skip_bad0 = args.skip_bad0
    badpix.skip_unstable = args.skip_unstable
    badpix.crval = args.crval
    badpix.find_cr = args.find_cr
    badpix.weirdsigma = args.weirdsigma
    badpix.skip_weird = args.skip_weird
    badpix.weirdfile = args.weirdfile
    badpix.karlweird = args.karlweird
    badpix.Karl_file = args.Karl_file
    badpix.Karl_vals_to_keep = args.Karl_vals_to_keep
    badpix.run()


    
