#! /usr/bin/env python

'''verify that the refpix correction works as intended by comparing
pipeline-output ramps to manually calculated ramps.'''

import numpy as np
from scipy.stats import sigmaclip
from astropy.io import fits
import copy,os,sys
from astropy.convolution import convolve, Box1DKernel, Box2DKernel
import matplotlib.pyplot as plt
if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
else:
    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
    sys.exit(0)
import sigmacut

#frame numbers to extract and use in the testing
frames = np.arange(10)

#files to compare
origfile = 'NRCNRCB1-DARK-53500942331_1_486_SE_2015-12-16T11h02m18_uncal_dq_init_superbias.fits'
refqfile = 'NRCNRCB1-DARK-53500942331_dq_init_superbias_refpixquadsonly.fits'
refeofile = 'NRCNRCB1-DARK-53500942331_dq_init_superbias_refpixquadsevenodd.fits'
ref1ffile = 'NRCNRCB1-DARK-53500942331_dq_init_superbias_refpixinclude1overf.fits'

def sigma_clipped_mean(pix):
    pix = np.ravel(pix)
    goodpts = ~np.isnan(pix)
    nanpts = np.isnan(pix)
    pixclass = sigmacut.calcaverageclass()
    pixclass.calcaverage_sigmacutloop(pix,Nsigma=3,saveused=True,mask=nanpts)
    finalmask = np.logical_or(nanpts,pixclass.clipped)
    pixclass.calcaverage_sigmacutloop(pix,Nsigma=0,mask=finalmask)
    return pixclass.mean,pixclass.stdev


#read in ramps and extract a subset of frames to test
with fits.open(origfile) as h:
    orig = h[1].data
orig = orig[0,frames,:,:].astype('float')

with fits.open(refqfile) as h:
    refq = h[1].data
refq = refq[0,frames,:,:].astype('float')

with fits.open(refeofile) as h:
    refeo = h[1].data
refeo = refeo[0,frames,:,:].astype('float')

with fits.open(ref1ffile) as h:
    ref1f = h[1].data
ref1f = ref1f[0,frames,:,:].astype('float')

#starting coordinates of each quadrant
qstart = [0,512,1024,1536,2044]

#-------------TEST THE QUAD MEAN ONLY SUBTRACTION--------------------------
#manually calculate the mean values in each quad
#and subtract from the data, to simulate the quad-average subtraction
qmeans = np.zeros((len(frames),4))
refq_manual = copy.deepcopy(orig)
for frameno in xrange(len(frames)):
    for i in xrange(4):
        rpix = np.append(orig[frameno,0:4,qstart[i]:qstart[i+1]],orig[frameno,2044:,qstart[i]:qstart[i+1]])
        #qmeans[frameno,i] = np.mean(rpix)
        #qmeans[frameno,i],junk = sigma_clipped_mean(rpix)
        goodrpix,lower,higher = sigmaclip(rpix,4.,4.)
        qmeans[frameno,i] = np.mean(goodrpix)
        #qmeans[frameno,i] = np.median(rpix)
        
        #manually subtact these mean values
        #refq_manual[frameno,4:2044,qstart[i]:qstart[i+1]] -= qmeans[frameno,i]
        refq_manual[frameno,4:2044,qstart[i]:qstart[i+1]] = refq_manual[frameno,4:2044,qstart[i]:qstart[i+1]] - qmeans[frameno,i]


print('Reference pixel means per quad in the non refpix subtracted file, as calculated manually:',qmeans)

#save the manually subtracted version
h0 = fits.PrimaryHDU()
h1 = fits.ImageHDU(refq_manual)
hdulist = fits.HDUList([h0,h1])
hdulist.writeto(refqfile[0:-5] + '_MANUAL.fits',clobber=True)

#save the extracted groups from the pipeline version
h0 = fits.PrimaryHDU()
h1 = fits.ImageHDU(refq)
hdulist = fits.HDUList([h0,h1])
hdulist.writeto(refqfile[0:-5] + '_PIPELINE.fits',clobber=True)


#compare to the values in refq, which come from the pipeline
diff1 = (refq - refq_manual) / refq

#this histogram is decieving. The variation shown is just the variation already present
#in the signal, rather than variations in the refpix subtraction differences
#ff, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#nn,bb,pp = ax1.hist(np.ravel(diff1[0,4:2044,0:512]),bins=np.arange(-0.1,0.005,.002))
#nn,bb,pp = ax2.hist(np.ravel(diff1[0,4:2044,512:1024]),bins=np.arange(-0.1,0.005,.002))
#nn,bb,pp = ax3.hist(np.ravel(diff1[0,4:2044,1024:1536]),bins=np.arange(-0.1,0.005,.002))
#nn,bb,pp = ax4.hist(np.ravel(diff1[0,4:2044,1536:2044]),bins=np.arange(-0.1,0.005,.002))
#ax1.set_title('Sec 1')
#ax2.set_title('Sec 2')
#ax3.set_title('Sec 3')
#ax4.set_title('Sec 4')
#ax3.set_xlabel('(Pipeline-Manual) Fractional Difference')
#ax3.set_ylabel('# Pixels')
#ax3.set_yscale('log')
#ax1.set_yscale('log')
#ax4.set_xlim(-0.1,0.0)
#ax3.set_xlim(-0.1,0.0)
#ax2.set_xlim(-0.1,0.0)
#ax1.set_xlim(-0.1,0.0)
#ff.savefig('diff_quadonly_hist_fractional.pdf')

#compare the manual subtraction to the pipeline subtraction, and check for agreement
diff = (refq - refq_manual)
print("Index, Group, Diffs:Q1 through Q4")
for i in xrange(diff.shape[0]):
    print(i,frames[i],diff[i,12,12],diff[i,12,600],diff[i,12,1030],diff[i,12,1600])

#print("the 4 quadrants are {}, {}, {}, {} DN.".format(q1diff[1],q2diff[1],q3diff[1],q4diff[1]))
#print("These values represent {}, {}, {}, and {}%".format(q1percmx,q2percmx,q3percmx,q4percmx))
#print("Doing quad-only calcs for the moment, skipping others.")
#sys.exit()


#-----------------TEST QUAD MEAN + EVEN ODD SUBTRACTION DIFFERENCE---------------
#same as above, but rather than calculating a single mean value for each quadrant, 
#find the mean for the even columns and the odd columns separately
qmeans = np.zeros((len(frames),4,2))
refeo_manual = copy.deepcopy(orig)
for frameno in xrange(len(frames)):
    for i in xrange(4):
        evenrpix = np.append(orig[frameno,0:4,qstart[i]:qstart[i+1]:2],orig[frameno,2044:,qstart[i]:qstart[i+1]:2])
        oddrpix = np.append(orig[frameno,0:4,qstart[i]+1:qstart[i+1]:2],orig[frameno,2044:,qstart[i]+1:qstart[i+1]:2])
        
        #qmeans[frameno,i,0],junk = sigma_clipped_mean(evenrpix)
        #qmeans[frameno,i,1],junk = sigma_clipped_mean(oddrpix)
        goodeven,lower,higher = sigmaclip(evenrpix,4.,4.)
        qmeans[frameno,i,0] = np.mean(goodeven)
        goododd,lower,higher = sigmaclip(oddrpix,4.,4.)
        qmeans[frameno,i,1] = np.mean(goododd)

        #qmeans[frameno,i,0] = np.mean(evenrpix)
        #qmeans[frameno,i,1] = np.mean(oddrpix)

        #manually subtact these mean values
        refeo_manual[frameno,4:2044,qstart[i]:qstart[i+1]:2] -= qmeans[frameno,i,0]
        refeo_manual[frameno,4:2044,qstart[i]+1:qstart[i+1]:2] -= qmeans[frameno,i,1]

#save the manually subtracted version
h0 = fits.PrimaryHDU()
h1 = fits.ImageHDU(refeo_manual)
hdulist = fits.HDUList([h0,h1])
hdulist.writeto(refeofile[0:-5] + '_MANUAL.fits',clobber=True)

#compare to the values from the pipeline version of the even/odd correction
diff = (refeo - refeo_manual)
print("Quad+Even/Odd: Index, Group, Diffs:Q1 through Q4")
for i in xrange(diff.shape[0]):
    print("Even:",i,frames[i],diff[i,12,12],diff[i,12,600],diff[i,12,1030],diff[i,12,1600])
    print("Odd:",i,frames[i],diff[i,12,13],diff[i,12,601],diff[i,12,1031],diff[i,12,1601])


#q1ediff = np.min(diff[0,4:2044,4:512:2]),np.max(diff[0,4:2044,4:512:2])
#q1odiff = np.min(diff[0,4:2044,5:512:2]),np.max(diff[0,4:2044,5:512:2])
#q2ediff = np.min(diff[0,4:2044,512:1024:2]),np.max(diff[0,4:2044,512:1024:2])
#q2odiff = np.min(diff[0,4:2044,513:1024:2]),np.max(diff[0,4:2044,513:1024:2])
#q3ediff = np.min(diff[0,4:2044,1024:1536:2]),np.max(diff[0,4:2044,1024:1536:2])
#q3odiff = np.min(diff[0,4:2044,1025:1536:2]),np.max(diff[0,4:2044,1025:1536:2])
#q4ediff = np.min(diff[0,4:2044,1536:2044:2]),np.max(diff[0,4:2044,1536:2044:2])
#q4odiff = np.min(diff[0,4:2044,1537:2044:2]),np.max(diff[0,4:2044,1537:2044:2])
#q1eperc = q1ediff[0]/qmeans[0,0,0]*100.
#q2eperc = q2ediff[0]/qmeans[0,1,0]*100.
#q3eperc = q3ediff[0]/qmeans[0,2,0]*100.
#q4eperc = q4ediff[0]/qmeans[0,3,0]*100.
#q1operc = q1odiff[0]/qmeans[0,0,1]*100.
#q2operc = q2odiff[0]/qmeans[0,1,1]*100.
#q3operc = q3odiff[0]/qmeans[0,2,1]*100.
#q4operc = q4odiff[0]/qmeans[0,3,1]*100.
#print("Mean-quadrant plus even/odd refpix subtraction, differences between pipeline and manual subtractions for")
#print("the 4 quadrants are (even refpix) {}, {}, {}, {} DN.".format(q1ediff[0],q2ediff[0],q3ediff[0],q4ediff[0]))
#print("These values represent {}, {}, {}, and {}%".format(q1eperc,q2eperc,q3eperc,q4eperc))
#print("the 4 quadrants are (odd refpix) {}, {}, {}, {} DN.".format(q1odiff[0],q2odiff[0],q3odiff[0],q4odiff[0]))
#print("These values represent {}, {}, {}, and {}%".format(q1operc,q2operc,q3operc,q4operc))


#----------------------TEST Quadmean, even/odd, and 1/f-----------------------
#Now test the case where the even/odd correction is performed, and followed by the 1/f correction using the
#side reference pixels. For the manual case, begin with the pipeline version of the even/odd corrected file, 
#then do a manual 1/f correction
refeo1overf_manual = copy.deepcopy(refeo)
mncollalll = np.zeros((2040,4))
mncollallr = np.zeros((2040,4))
mncoll = np.zeros(2040)
mncolr = np.zeros(2040)
mncol = np.zeros(2040)
rpix = np.append(refeo[0,:,0:4],refeo[0,:,2044:],axis=1)
smoothed = np.zeros((2040,8))


#Robert says the SSB pipeline takes the median of the 44 (11rows x 4 columns) pixels 
#surrounding the pixel in question for the left and right sides, and then takes the average
#of the left and right sides and uses that as the smoothed 1/f value to subtract.
#Ignore the rows near the top and bottom for now, so we don't need to worry about reflection...
for frameno in xrange(len(frames)):
    for y in xrange(6,2030):
        mncoll[y] = np.median(refeo[frameno,y-5:y+6,0:4])
        mncolr[y] = np.median(refeo[frameno,y-5:y+6,2044:])
        mncol[y] = np.mean([mncoll[y],mncolr[y]])
        #subtract from the science pixels in that row
        refeo1overf_manual[frameno,y,4:2044] = refeo1overf_manual[frameno,y,4:2044] - mncol[y]

#save the manually subtracted version
h0 = fits.PrimaryHDU()
h1 = fits.ImageHDU(refeo1overf_manual)
hdulist = fits.HDUList([h0,h1])
hdulist.writeto(ref1ffile[0:-5] + '_MANUAL.fits',clobber=True)

diff = (ref1f - refeo1overf_manual)
print("1/f Removal: Index, Group, Diffs:Q1 through Q4")
for i in xrange(diff.shape[0]):
    print("Row 1:",i,frames[i],diff[i,12,12],diff[i,12,600],diff[i,12,1030],diff[i,12,1600])
    print("Row 2:",i,frames[i],diff[i,500,12],diff[i,500,600],diff[i,500,1030],diff[i,500,1600])

#yay!! my 1/f correction matches Robert's almost perfectly. differences are 1e-7DN

q1diff = np.min(diff[0,6:2030,4:512]),np.max(diff[0,6:2030,4:512])
q2diff = np.min(diff[0,6:2030,512:1024]),np.max(diff[0,6:2030,512:1024])
q3diff = np.min(diff[0,6:2030,1024:1536]),np.max(diff[0,6:2030,1024:1536])
q4diff = np.min(diff[0,6:2030,1536:2044]),np.max(diff[0,6:2030,1536:2044])
q1perc = q1diff[0]/qmeans[0,0]*100.
q2perc = q2diff[0]/qmeans[0,1]*100.
q3perc = q3diff[0]/qmeans[0,2]*100.
q4perc = q4diff[0]/qmeans[0,3]*100.
print("Mean-quadrant, even/odd, 1/f, differences between pipeline and manual subtractions for")
print("the 4 quadrants are {}, {}, {}, {} DN.".format(q1diff[0],q2diff[0],q3diff[0],q4diff[0]))
print("These values represent {}, {}, {}, and {}%".format(q1perc,q2perc,q3perc,q4perc))
print("Skipping frame reset calcs for the moment.")
sys.exit()

#---------------------------------------------------------------------
#Next, quantify how well the quadrant-to-quadrant offsets were removed
#by looking at the quad-to-quad offsets in the subtracted data
with fits.open(origfile) as h:
    orig = h[1].data
with fits.open(refeofile) as h:
    refeo = h[1].data

ints,groups,yd,xd = orig.shape
origmeans = np.zeros((ints,groups,4))
origmeanerrs = np.zeros((ints,groups,4))
refeomeans = np.zeros((ints,groups,4))
refeomeanerrs = np.zeros((ints,groups,4))
qstart = [0,512,1024,1536,2040]
sig = sigmacut.calcaverageclass()
for integration in xrange(ints):
    for group in xrange(groups):
        for quad in xrange(4):
            #use sigma-clipped mean
            sig.calcaverage_sigmacutloop(orig[integration,group,4:2040,qstart[quad]:qstart[quad+1]],Nsigma=3.0,verbose=0)
            origmeans[integration,group,quad] = sig.mean
            origmeanerrs[integration,group,quad] = sig.mean_err
            sig.calcaverage_sigmacutloop(refeo[integration,group,4:2040,qstart[quad]:qstart[quad+1]],Nsigma=3.0,verbose=0)
            refeomeans[integration,group,quad] = sig.mean
            refeomeanerrs[integration,group,quad] = sig.mean_err

origmeanmean = np.zeros(4)
origmeanmeanerr = np.zeros(4)
refeomeanmean = np.zeros(4)
refeomeanmeanerr = np.zeros(4)
for quad in xrange(4):
    sig.calcaverage_sigmacutloop(origmeans[0,:,quad],Nsigma=3.0,verbose=0)
    origmeanmean[quad] = sig.mean
    origmeanmeanerr[quad] = sig.mean_err
    sig.calcaverage_sigmacutloop(refeomeans[0,:,quad],Nsigma=3.0,verbose=0)
    refeomeanmean[quad] = sig.mean
    refeomeanmeanerr[quad] = sig.mean_err
print('Looking at quad to quad offsets in mean signal:')
print('Uncorrected data, mean scipix values, 4 quads: {}'.format(origmeanmean))
print('Refpix corrected data, mean scipix values, 4 quads: {}'.format(refeomeanmean))

ff,aa = plt.subplots()#2,sharex=True)
xs = np.arange(groups)
#aa.plot(xs,origmeans[0,:,0],color='black',marker='8',label='Uncorrected')
aa.plot(xs,origmeans[0,:,1]-origmeans[0,:,0],color='blue',marker='8',label='Uncorr, Q2')
aa.plot(xs,origmeans[0,:,2]-origmeans[0,:,0],color='red',marker='8',label='Uncorr, Q3')
aa.plot(xs,origmeans[0,:,3]-origmeans[0,:,0],color='green',marker='8',label='Uncorr, Q4')
#aa.legend(loc='best')
#aa.plot(xs,refeomeans[0,:,0],color='lightgrey',marker='8',label='Corrected')
aa.plot(xs,refeomeans[0,:,1]-refeomeans[0,:,0],color='lightskyblue',marker='8',label='Corr, Q2')
aa.plot(xs,refeomeans[0,:,2]-refeomeans[0,:,0],color='salmon',marker='8',label='Corr, Q3')
aa.plot(xs,refeomeans[0,:,3]-refeomeans[0,:,0],color='lightgreen',marker='8',label='Corr, Q4')
aa.plot([0,108],[0,0],color='black',linestyle='--')
aa.legend(loc='best',fontsize=9)
aa.set_xlim(0,108)
aa.set_xlabel('Group Number')
aa.set_ylabel('Mean Rel to Sector 0 (DN)')
ff.savefig('quad_to_quad_levels.pdf')
plt.close(ff)

ff,aa = plt.subplots()
xs = np.arange(4)+1
aa.plot(xs,origmeanmean,color='black',marker='8',label='Uncorrected')
aa.plot(xs,refeomeanmean,color='red',marker='8',label='Corrected')
labels=['1','2','3','4']
plt.xticks(xs, labels)
aa.set_xlabel('Sector Number')
aa.set_ylabel('Mean Scipix Signal (DN)')
aa.legend(loc='best')
ff.savefig('quad_to_quad_levels_meanof108groups.pdf')
plt.close(ff)
