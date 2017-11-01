#! /usr/bin/env python

'''
Create data files for grouped data saturation

Basic idea: 

Begin with standard saturation reference file
For each readout pattern....


'''


from jwst.datamodels import SaturationModel,LinearityModel,SuperBiasModel
from glob import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
#import datetime as d

class GroupSat():
    """
    Create satruation values for grouped data

    Given an ungrouped saturation reference file
    (i.e. a standard JWST saturation reference file)
    and a description of the grouped data characteristics,
    calculate the grouped saturation values, which are 
    specific to the specified readout pattern

    Parameters
    ----------
    satfile
       Name of ungrouped saturation reference file
    linfile 
       Name of linearity correction reference file
    superbiasfile 
       Name of superbias reference file
    maxiter
       Maximum number of iterations of Newton's method
       when calculating inverse linearity coefficients
    accuracy
       Threshold accuracy value at which to stop Newton's
       method and declare success.
    maxgroups
       The maximum number of groups per integration for the 
       specified readout pattern. For NIRCam, this is 10 for
       all readout patterns except the DEEP patterns, where it
       is 20.
    frametime 
       Readout time for a single frame. For full frame NIRCam
       observations, this is 10.73676 seconds.
    detector 
       Name of detector. (e.g. 'NRCA1')
    readpatt 
       Readout pattern of output. (e.g. 'SHALLOW4')
    nframe
       Number of frames per group for the specified readout
       pattern.
    nskip
       Number of skipped frames per group for the specified
       readout pattern.
    xsize
       Number of pixels in the x-dimension of the output. Use
       'full' for full-frame, otherwise use a number.
    ysize
       Number of pixels in the y-dimension of the output. Use
       'full' for full-frame, otherwise use a number.
    xstart 
       Column number (in full-frame coords) to treat as left-
       most column in the case that the output is a subarray.
       Works in concert with xsize to extract subarray from
       the input reference files.
    ystart 
       Row number (in full-frame coords) to treat as bottom-
       most row in the case that the output is a subarray.
       Works in concert with ysize to extract subarray from
       the input reference files.       
    inverse_lin_file
       Name of reference file containing inverse linearity 
       coefficients. These are coefficients for translateing
       from linear to non-linear signal (opposite of the 
       standard non-linearity reference file). If left as 
       None, these values will be calculated on the fly.
       This significantly slows the speed of the grouped 
       saturation values.
    outfile
       Name of fits file to contain the grouped saturation 
       values.

    Returns
    -------
    satramp
       3D array of grouped saturation values for the specified
       detector/readout pattern/subarray
    """

    
    def __init__(self):
        self.satfile = None
        self.linfile = None
        self.superbiasfile = None
        self.maxiter = 10
        self.accuracy = 0.000001
        self.maxgroups = 10
        self.frametime = 10.73676
        self.detector = None
        self.readpatt = None
        self.nframe = 4
        self.nskip = 1
        self.xsize = 'full'
        self.ysize = 'full'
        self.ystart = 0
        self.xstart = 0
        self.inverse_lin_file = None
        self.outfile = None
        
        
    def linearize(self,ramparr,coeffarr):
        """Remove non-linearity from the input data
        code taken from jwst pipeline 
        """
        ncoeffs = coeffarr.shape[0]
        
        #Calculate corrected counts
        scorr = coeffarr[ncoeffs - 1] * ramparr
        for j in range(ncoeffs - 2, 0, -1):
            scorr = (scorr + coeffarr[j]) * ramparr
        scorr = coeffarr[0] + scorr
        return scorr


    def unlinearize(self,image,coeffs,sat):
        """Insert non-linearity into input data using
        the nonlin correction coefficients and
        Newton's method"""
        #This function is the bottleneck in terms of speed

        #find pixels with "good" signals. Negative pix or pix with
        #signals above the requested max value will not be changed.
        x = np.copy(image)
        i1 = np.where(image > 0.) 
        dev = np.zeros_like(image,dtype=float)
        dev[i1] = 1.
        i2 = np.where(image <= 0.) 
        
        i = 0
        #initial run of the nonlin function 
        val = self.linearize(image,coeffs)
        val[i2]=1.

        #set any pix where val <= 0 to 1.
        neg = val <= 0
        if np.sum(neg) > 0:
            val[neg] = 1.
        
        x[i1] = (image[i1]+image[i1]/val[i1]) / 2. 

        while i < self.maxiter:
            i=i+1
            val = self.linearize(x,coeffs)
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
        return x


    def unlinearize_polynomial(self,s,c):
        """Apply coeffs to convert linear->non-linear signal"""
        return np.array(c[0] + c[1]*s + c[2]*s**2 + c[3]*s**3
                        + c[4]*s**4 + c[5]*s**5 + c[6]*s**6)
    

    def nonLinDeriv(self,image,coeffs,limits):
        """1st derivative of non-lin function"""
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
        return t

    
    def savefile(self,ramp,detector,readpatt):
        """save output as a generic fits file"""
        if self.outfile is None:
            self.outfile = 'GroupedDataSaturation_' + self.infile
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(ramp)
        h1.header['DETECTOR'] = detector
        h1.header['READPATT'] = readpatt
        hlist = fits.HDUList([h0,h1])
        hlist.writeto(self.outfile,clobber=True)
    

    def get_lin_coeffs(self,file):
        """Read in non-linearity coefficients, and check
        for bad pixels
        """
        lin_model = LinearityModel(file)
        data = lin_model.coeffs

        #set pixels with bad linearity coeffs
        #so that no correction is applied
        numcof,ys,xs = data.shape
        bad = np.where(np.isnan(data))
        for y,x in zip(bad[1],bad[2]):
            data[0,y,x] = 0.
            data[1,y,x] = 1.
            data[2:,y,x] = 0.
        return data


    def get_sat_vals(self,file):
        """Read in SSB-format saturation reference file
        and deal with bad values
        """
        sat_model = SaturationModel(file)
        svals = sat_model.data

        #set pixels with saturation level of 0.
        #to 65535
        bad = svals == 0
        svals[bad] = 65535
        return svals


    def get_superbias(self,file):
        """Read in superbias and deal with bad pixels"""
        sb_model = SuperBiasModel(file)
        sb = sb_model.data
        bad = np.isnan(sb)
        sb[bad] = 0.
        return sb

        
    def make_file(self):
        """
        Create satruation values for grouped data

        Given an ungrouped saturation reference file
        (i.e. a standard JWST saturation reference file)
        and a description of the grouped data characteristics,
        calculate the grouped saturation values, which are 
        specific to the specified readout pattern

        Parameters
        ----------
        none

        Returns
        -------
        satramp
              3D numpy array of saturation values with
              format (groups x ydim x xdim)
        """
        
        #read in saturation file
        sat = self.get_sat_vals(self.satfile)
        
        #Read in linearity coefficients file
        lin = self.get_lin_coeffs(self.linfile)
        
        #superbias file
        superbias = self.get_superbias(self.superbiasfile)
        
        #if the output is a subarray, extract the proper
        #subarrays from the reference files
        ys,xs = sat.shape
        if self.xsize == 'full':
            self.xsize = xs
        if self.ysize == 'full':
            self.ysize = ys

        if ((self.xsize != xs) | (self.ysize != ys)):
            sat = sat[self.ystart:self.ystart+self.ysize,
                      self.xstart:self.xstart+self.xsize]
            lin = lin[:,self.ystart:self.ystart+self.ysize,
                      self.xstart:self.xstart+self.xsize]
            superbias = superbias[self.ystart:self.ystart+self.ysize,
                                  self.xstart:self.xstart+self.xsize]

        #t0 = d.datetime.now()
            
        #linearize the saturation values
        sat_lin = self.linearize(sat-superbias,lin)

        #read in coefficients to convert linear to non-linear signals
        if self.inverse_lin_file is not None:
            with fits.open(self.inverse_lin_file) as h:
                inv_lin_coeffs = h[1].data
        
        #loop over groups
        satramp = np.zeros((self.maxgroups,self.ysize,self.xsize))
        linsatramp = np.zeros((self.maxgroups,self.ysize,self.xsize))
        xfmeans = []
        for i in range(self.maxgroups):
            #exposure time to the final frame in the group
            exptime = self.frametime * ((self.nframe+self.nskip)
                                        * (i+1)) - 1
            satslope = sat_lin / exptime

            #now calculate signals for each frame within the
            #group, by reducing exposure time by one frametime
            #for each
            fsigs = np.zeros((self.nframe,ys,xs))
            fsigs[0,:,:] = sat_lin 
            for frame in range(1,self.nframe):
                fsigs[frame,:,:]  = satslope * (exptime
                                                -(self.frametime*frame))
                linsatramp[i,:,:] = np.mean(fsigs,axis=0) + superbias

            #non-linearize fsigs
            if self.inverse_lin_file is not None:
                fsigs_nl = self.unlinearize_polynomial(fsigs,
                                                       inv_lin_coeffs)
            else:
                fsigs_nl = self.unlinearize(fsigs,lin,sat-superbias)
                
            #add superbias back in
            fsigs_nl += superbias

            #Find the mean signal in frames. This is the group signal.
            satramp[i,:,:] = np.mean(fsigs_nl,axis=0)

        #t1 = d.datetime.now()
        #dt = t1-t0
        #print('Elspased time: {}'.format(dt))
            
        #save output file
        self.savefile(satramp,self.detector,self.readpatt)

        #return group saturation values
        return satramp
        

    def add_options(self,parser=None,usage=None):
        """Parse arguments"""
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,
                                             description=
                                             'Simulate JWST ramp')
        parser.add_argument("satfile",help='Satruation reference file from which to create the file for grouped data')
        parser.add_argument("linfile",help='Linearity reference file to use for linearizing/unlinearizing data')
        parser.add_argument("superbiasfile",help='Superbias reference file to use')
        parser.add_argument("readpatt",help='Readout pattern to use. Used only in output metadata')
        parser.add_argument("detector",help='Detector name. Used only in output metadata')
        parser.add_argument("nframe",help='Number of frames averaged to create a group',type=np.int)
        parser.add_argument("nskip",help='Number of frames skipped per group',type=np.int)
        parser.add_argument("frametime",help='Exposure time for a single frame.',type=np.float)
        #inputs for creating subarray output file from full frame inputs
        parser.add_argument("--xsize",help='x-size of array for output. "full" means keep full frame',default='full')
        parser.add_argument("--ysize",help='x-size of array for output. "full" means keep full frame',default='full')
        parser.add_argument("--xstart",help='x-coordinate of left side of output array in full-frame coords',default=0,type=np.int)
        parser.add_argument("--ystart",help='y-coordinate of left side of output array in full-frame coords',default=0,type=np.int)
        parser.add_argument("--maxiter",help='Maximum number of iterations to use in solving to add non-linearity to data',default=10)
        parser.add_argument("--accuracy",help='Accuracy limit at which to declare success when iterating to add non-linearity to data',default = 0.000001)
        parser.add_argument("--maxgroups",help='Number of groups to run the calculations on. (e.g. NIRCam allows obs of up to 10 groups in most cases.)',default=10,type=np.int)
        parser.add_argument("--inverse_lin_file",help="Name of file containing coeffs to translate between linear and non-linear data. If not provided, Newton's method will be used with the standard non-lin coefficients to convert to non-linear signals.",default=None)
        parser.add_argument("--outfile",help='Name of output file containing saturation data',default=None)
        return parser



if __name__ == '__main__':
    usagestring = 'group_saturation.py satfile.fits linfile.fits sb.fits SHALLOW4 NRCA1 4 1 10.73676'
    
    gs = GroupSat()
    parser = gs.add_options(usage=usagestring)
    args = parser.parse_args(namespace=gs)
    gs.make_file()
