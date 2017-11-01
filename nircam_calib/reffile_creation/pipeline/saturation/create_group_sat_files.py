#! /usr/bin/env python

'''
Create data files for grouped data saturation

Basic idea: 

Begin with standard saturation reference file.
For each NIRCam readout pattern/subarray/detector, 
create a grouped saturation reference file.

WARNING: This adds up to over 1300 files given
all of NIRCam's readout modes.

'''


from jwst.datamodels import SaturationModel,LinearityModel,SuperBiasModel
from glob import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


class GroupSat():
    def __init__(self):
        self.verbose = False
        satdir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/welldepth/'
        self.satfiles = glob(satdir + '*wfact_DMSorient.fits')

        lindir = '/ifs/jwst/wit/witserv/data7/nrc/reference_files/SSB/CV3/cv3_reffile_conversion/linearity/'
        self.linfiles = glob(lindir + '*v2_DMSorient.fits')
        self.maxiter = 10
        self.accuracy = 0.000001

        sbdir = '/ifs/jwst/wit/witserv/data4/nrc/hilbert/superbias/cv3/deliver_to_CRDS/'
        self.sbfiles = glob(sbdir + '*superbias*fits')
                            
        #self.readpatts = {'BRIGHT2':(2,0),'SHALLOW2':(2,3),'SHALLOW4':(4,1),'MEDIUM2':(2,8),'MEDIUM8':(8,2),'DEEP2':(2,18),'DEEP8':(8,12)}
        self.readpatts = {'SHALLOW4':(4,1)} # for testing
        self.maxgroups = 10
        self.frametime = 10.73676
        #self.detectors = ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5']
        self.detectors = ['A1'] #for testing
        self.plot = True
        
        
    def linearize(self,ramparr,coeffarr):
        #remove non-linearity from the input data
        #taken from jwst pipeline
        
        ncoeffs = coeffarr.shape[0]
        
        # Accumulate the polynomial terms into the corrected counts
        scorr = coeffarr[ncoeffs - 1] * ramparr
        for j in range(ncoeffs - 2, 0, -1):
            scorr = (scorr + coeffarr[j]) * ramparr
        scorr = coeffarr[0] + scorr
        return scorr


    def unlinearize(self,image,coeffs,sat):
        #insert non-linearity into the linear synthetic sources

        #find pixels with "good" signals, to have the nonlin applied. Negative pix or pix with
        #signals above the requested max value will not be changed.
        x = np.copy(image)
        i1 = np.where(image > 0.) 
        dev = np.zeros_like(image,dtype=float)
        dev[i1] = 1.
        i2 = np.where(image <= 0.) 
        print("In unlinearize, {} pixels in the input have signal values < 0".format(len(i2)))
        print("and will be set to 1.0")
        
        i = 0

        #initial run of the nonlin function - when calling the non-lin function, give the
        #original satmap for the non-linear signal values
        val = self.linearize(image,coeffs)#,sat)
        val[i2]=1.

        #set any pix where val = 0 to 1.
        #z = val == 0.
        #val[z] = 1.
        neg = val <= 0
        if np.sum(neg) > 0:
            print("Linearized signal has {} pix with zero or negative signals. Setting to 1.".format(np.sum(neg)))
            val[neg] = 1.

        
        x[i1] = (image[i1]+image[i1]/val[i1]) / 2. 

        while i < self.maxiter:
            i=i+1
            val = self.linearize(x,coeffs)#,sat)
            val[i2]=1.

            #if val <= 0, set to 1
            z = val <= 0
            val[z] = 1.

            
            dev[i1] = abs(image[i1]/val[i1]-1.)
            inds = np.where(dev[i1] > self.accuracy)
            if inds[0].size < 1:
                break
            val1 = self.nonLinDeriv(x,coeffs,sat)
            val1[i2] = 1.
            zv1 = val1 == 0
            if np.sum(zv1) > 0:
                print("{} pixels with zero val1. Setting to 1.".format(np.sum(zv1)))
                val1[zv1] = 1.
            x[i1] = x[i1] + (image[i1]-val[i1])/val1[i1]

        #if we max out the number of iterations, save the array of accuracy values
        #Spot checks reveal the pix that don't meet the accuracy reqs are randomly
        #located on the detector, and don't seem to be correlated with point source
        #locations.
        if i == self.maxiter:
            ofile = 'doNonLin_accuracy.fits'
            devcheck = np.copy(dev)
            devcheck[i2] = -1.
            h0 = fits.PrimaryHDU()
            h1 = fits.ImageHDU(devcheck)
            hl = fits.HDUList([h0,h1])
            hl.writeto(ofile,overwrite=True)
                
        return x


    def nonLinDeriv(self,image,coeffs,limits):
        values = np.copy(image)

        bady = 0
        badx = 1
        if len(image.shape) == 3:
            bady = 1
            badx = 2
        
        bad = np.where(values > limits)
        values[bad] = limits[bad[bady],bad[badx]]
        ncoeff = coeffs.shape[0]
        t = (ncoeff-1) * np.copy(coeffs[-1,:,:])
        for i in range(ncoeff-3,-1,-1):
            t = (i+1) * coeffs[i+1,:,:] + values*t
        #t = coeffs[0,:,:] + values*t
        return t


    
    def savefile(self,ramp,detector,readpatt):
        outfile = 'GroupedData_Satfile_Det{}_Readpatt{}.fits'.format(detector,readpatt)
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(ramp)
        h1.header['DETECTOR'] = detector
        h1.header['READPATT'] = readpatt
        hlist = fits.HDUList([h0,h1])
        hlist.writeto(outfile,clobber=True)
    
    
    def make_files(self):

        #loop over detector
        for det in self.detectors:

            satfile = [s for s in self.satfiles if det in s][0]
            sat_model = SaturationModel(satfile)
            sat = sat_model.data

            #pixels with saturation level of 0.
            #set sat level to 65535
            bad = sat == 0
            sat[bad] = 65535

            linfile = [s for s in self.linfiles if det in s][0]
            lin_model = LinearityModel(linfile)
            lin = lin_model.coeffs

            #pixels with bad linearity coeffs
            #set so no correction is applied
            nans = np.isnan(lin[1,:,:])
            tmp = lin[1,:,:]
            tmp[nans] = 1.
            lin[1,:,:] = tmp
            for i in range(2,7):
                tmp = lin[i,:,:]
                nans = np.isnan(tmp)
                tmp[nans] = 0.
                lin[i,:,:] = tmp
            
            #superbias file
            sbfile = [s for s in self.sbfiles if det in s][0]
            sb_model = SuperBiasModel(sbfile)
            superbias = sb_model.data

            #linearize the saturation values
            sat_lin = self.linearize(sat-superbias,lin)
            
            #loop over readout patterns
            for readpat in self.readpatts:

                nframe,nskip = self.readpatts[readpat]

                #optional output plot
                if self.plot:
                    xx = 400
                    yy = 400
                    fsize=12
                    f = plt.figure()
                    a = f.add_subplot(111)
                    a2 = a.twiny()
                    f.subplots_adjust(bottom=0.2)
                    xs = np.arange(self.maxgroups*(nframe+nskip))
                    a.plot(xs,np.repeat(sat[yy,xx],len(xs)),linestyle=':',color='blue',linewidth=2,label='Original Saturation')
                    a.plot(xs,np.repeat(sat_lin[yy,xx]+superbias[yy,xx],len(xs)),linestyle=':',color='red',linewidth=2,label='Linearized Saturation')
                    a.plot(xs,np.repeat(superbias[yy,xx],len(xs)),linestyle=':',color='black',linewidth=2,label='Superbias Level')
                
                #loop over groups
                satramp = np.zeros((self.maxgroups,2048,2048))
                linsatramp = np.zeros((self.maxgroups,2048,2048))
                xfmeans = []
                for i in range(self.maxgroups):

                    #exposure time to the final frame in the group
                    exptime = self.frametime * (nframe+nskip) * (i+1)
                    satslope = sat_lin / exptime

                    #now calculate signals for each frame within the
                    #group, by reducing exposure time by one frametime
                    #for each
                    fsigs = np.zeros((nframe,2048,2048))
                    fsigs[0,:,:] = sat_lin 
                    for frame in range(1,nframe):
                        fsigs[frame,:,:]  = satslope * (exptime-(self.frametime*frame))
                        linsatramp[i,:,:] = np.mean(fsigs,axis=0) + superbias
                        
                    #non-linearize fsigs
                    fsigs_nl = self.unlinearize(fsigs,lin,sat-superbias)
                    
                    #add superbias back in
                    fsigs_nl += superbias
                    
                    satramp[i,:,:] = np.mean(fsigs_nl,axis=0)

                    #print("Group: {}".format(i))
                    #print("Exptime: {}".format(exptime))
                    #print("Sat: {}, Lin Sat: {}".format(sat[yy,xx],sat_lin[yy,xx]))
                    #print("Satslope: {}".format(satslope[yy,xx]))
                    #print("Frame exps: {}".format(exptime-(self.frametime*np.arange(1,nframe))))
                    #print("signals: {}".format(fsigs[:,yy,xx]))
                    #print("mean signal: {}".format(satramp[i,yy,xx]))
                    #if i == 1:
                    #    stop
                    
                    if self.plot:
                        xf = np.array([-1])
                        xf = np.append(xf,np.arange((i+1)*(nframe+nskip)-nframe,(i+1)*(nframe+nskip)))
                        yf = np.array([superbias[yy,xx]])
                        yf = np.append(yf,fsigs[:,yy,xx][::-1]+superbias[yy,xx])
                        yfnl = np.array([superbias[yy,xx]])
                        yfnl = np.append(yfnl,fsigs_nl[:,yy,xx][::-1])
                        if i == 9:
                            a.plot(xf,yf,mfc='red',mec='red',color='black',marker='o',linestyle='-',label='Linearized Frames',alpha=0.5)
                            #a.plot(np.mean(xf[1:]),linsatramp[i,yy,xx],color='red',marker='8',mfc='none',mec='red',markersize=10,label='Linearized Frames Mean')
                            a.plot(xf,yfnl,mfc='black',mec='black',color='black',marker='o',linestyle='-',label='Non-Linearized Frames',alpha=0.5)
                            #a.plot(np.mean(xf[1:]),linsatramp[i,yy,xx],color='red',marker='8',mfc='none',mec='red',markersize=10,label='Linearized Frames Mean')                  
                        else:
                            a.plot(xf,yf,mfc='red',mec='red',color='black',marker='o',linestyle='-',alpha=0.5)
                            #a.plot(np.mean(xf[1:]),linsatramp[i,yy,xx],color='red',mfc='none',mec='red',marker='8',markersize=10)
                            a.plot(xf,yfnl,mfc='black',mec='black',color='black',marker='o',linestyle='-',alpha=0.5)
                        xfmeans.append(np.mean(xf[1:]))


                #optional plot
                if self.plot:
                    new_tick_locations = np.array(xfmeans)
                    a.plot(xfmeans,satramp[:,yy,xx],color='blue',marker='8',markersize=10,label='Non-Linear Frames Mean')
                    
                    a.legend(loc='lower right',numpoints=1,fontsize=fsize)
                    a.set_xlabel('Frame Number')
                    a.set_ylabel('Signal (DN)')
                    a.set_title('NRC{}, Readpattern: {}, Pixel ({},{})'.format(det,readpat,xx,yy))
                    a.set_xlim(-1,np.max(xf)+1)


                    a.set_ylim(0,60000)
                    
                    # Move twinned axis ticks and label from top to bottom
                    a2.xaxis.set_ticks_position("bottom")
                    a2.xaxis.set_label_position("bottom")

                    # Offset the twin axis below the host
                    a2.spines["bottom"].set_position(("axes", -0.15))
                    
                    # Turn on the frame for the twin axis, but then hide all 
                    # but the bottom spine
                    a2.set_frame_on(True)
                    a2.patch.set_visible(False)
                    for sp in a2.spines.itervalues():
                        sp.set_visible(False)
                    a2.spines["bottom"].set_visible(True)

                    a2.set_xticks(new_tick_locations)
                    a2.set_xticklabels(np.arange(len(new_tick_locations)))
                    a2.set_xlabel("Group Number")
                    
                    f.savefig('GroupSatPlot_forTR_{}_{}_C_filledcircs.png'.format(det,readpat))

                #save output file
                self.savefile(satramp,det,readpat)


    def tick_function(self,X):
        #V = 1/(1+X)
        #return ["%.3f" % z for z in V]
        return np.arange(len(X))
    
if __name__ == '__main__':

    gs = GroupSat()
    gs.make_files()
