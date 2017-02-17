#! /usr/bin/env python

'''
For a given detector, generate a map of RA and Dec across all x,y of the
detector, for a given roll angle.

These maps can help to define RA and Dec ranges to use for inputs 
to the ramp simulator

Inputs:
apername: Name of the aperture to simulate. The name must be present in the
          AperName column of the file that contains the distortion coefficients 
          (full_coeff_file input). Examples: NRCA1_FULL, NRCB4_SUB160

dist_model: Name of the distortion reference file that contains the models
            to translate between x,y and V2,V3. This is a CRDS-formatted
            ASDF file. 

ra: Right ascention at the reference location. For most apertures, the reference
    location is at the center of the aperture. Another way to think about it is
    the RA of the observation. RA must be given either in decimal degrees, or 
    as a string hh:mm:ss.ss

dec: Declination at the reference location. Dec must be given either as decimal
     degrees or as a string dd:mm:ss.ss

rotation: Rotation, in degrees, of the field of view. (i.e. roll angle of JWST)

Optional inputs:

refpixx - x coordinate of the reference location on the detector. 1023.5 for full 
         frame NIRCam arrays. If not provided, the value is read in from the
         full_coeff_file input file.

refpixy - y coordinate of the reference location on the detector. 1023.5 for full 
         frame NIRCam arrays. If not provided, the value is read in from the
         full_coeff_file input file.

refpixv2 - V2 value of the reference location on the detector. If not provided, 
         the value is read in from the full_coeff_file input file.

refpixv3 - V3 value of the reference location on the detector. If not provided, 
         the value is read in from the full_coeff_file input file.

full_coeff_file - Name of csv file that contains the full suite of coordinate
                  translation coefficients (sci->ideal->V2V3), along with 
                  reference location info for each aperture. Either this 
                  input must be provided, or the 4 refpix inputs above.

output - Name of file to save the contour plot of RA and Dec over the detector x,y

Output:

Image that shows a contour plot of RA and Dec values on top of the detector for the 
given RA, Dec, and rotation of the observation. The RA and Dec units in the output file
match those of the input RA and Dec. If RA is given as a string, then the contour labels
are strings. If the Dec is given in decimal degrees, then the contour labels are decimal
degrees.


Dependencies:
Uses several of Colin Cox's functions for translating coefficients. These are in his 
rotations.py and polynomial.py scripts.


Example call:

python RADec_vs_xy_map.py NRCA1_FULL NRCA1_FULL_distortion.asdf 64.5423 34.342 0.0 --full_coeff_file NIRCam_SIAF_2016-09-29.csv


'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import rotations  #Colin's functions
import polynomial #more of Colin's functions
import argparse, sys, os
from asdf import AsdfFile
from astropy.io import ascii

class RADecMap():

    def __init__(self):
        self.verbose = False

    def run(self):

        #check the entered RA and Dec. If entered as strings, then convert to decimal
        #degrees
        try:
            self.ra = float(self.ra)
            degraflag = True
        except:
            self.ra, junk = self.parseRADec(self.ra,"00:00:00.0")
            degraflag = False

        try:
            self.dec = float(self.dec)
            degdecflag = True
        except:
            junk,self.dec = self.parseRADec("00:00:00.0",self.dec)
            degdecflag = False
            
        #read in the transformation model
        coord_transform = self.get_coord_transform_model(self.dist_model)

        #if the reference location information is not provided by the user, then we
        #need to read it in from the full coefficient file
        if self.refpixx is None:
            x_sci2idl,y_sci2idl,v2_ref,v3_ref,x_ref,y_ref,xsize,ysize = self.read_coeff_file(self.full_coeff_file,self.apername)
            #reference pixel locations in the file are 1-based. Subtract 1.
            x_ref -= 1.
            y_ref -= 1.
        else:
            v2_ref = self.refpixv2
            v3_ref = self.refpixv3
            x_ref = self.refpixx
            y_ref = self.refpixy
            xsize = 2048
            ysize = 2048
            
        #Now create the attitude matrix needed for the coordinate transforms
        attitude_matrix = rotations.attitude(v2_ref,v3_ref,self.ra,self.dec,self.rotation)

        #ra,dec = self.XYToRADec(1502.45,1502.35,attitude_matrix,coord_transform,1023.5,1023.5,120.6714,-527.3877)
        #print(ra,dec)
        #sys.exit()

        
        #Generate a meshgrid of pixel values that cover the aperture
        x = np.arange(0,xsize,20)
        y = np.arange(0,ysize,20)
        #x = np.arange(0,2049,512)
        #y = np.arange(0,2049,512)
        xs,ys = np.meshgrid(x,y)


        #xs = np.zeros((3,2))
        #xs[:,0] = [0,1023.5,2047]
        #ys = xs
        yd,xd = xs.shape
        ra = np.zeros((yd,xd),dtype=np.float)
        dec = np.zeros((yd,xd),dtype=np.float)
        for i in xrange(xs.shape[1]):
            for j in xrange(xs.shape[0]):
        
                #calculate the RA and Dec
                ra[j,i],dec[j,i] = self.XYToRADec(xs[j,i],ys[j,i],attitude_matrix,coord_transform,x_ref,y_ref,v2_ref,v3_ref)
                #rr,dd = self.XYToRADec(1023.5,1024,attitude_matrix,coord_transform,x_ref,y_ref,v2_ref,v3_ref)



        #If the RA cycles past 360 and goes back to zero, we need to deal with that so the contours don't get screwy
        if np.max(ra) - np.min(ra) > 180:
            low = ra < 180
            ra[low] += 360.

                
        #plot results. use a contour plot to place the RA,Dec values on top of the x,y values
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        #f,a = plt.subplots()
        plt.figure()

        plt.grid(b=True, which='major', color='black', linestyle='--')
        plt.grid(b=True, which='minor', color='black', linestyle='--')
        cs = plt.contour(xs,ys,ra,6,colors='red')
        if degraflag == False:
            cs.levels = [self.makePos(val,0.0)[0][0:-2] for val in cs.levels]
        plt.clabel(cs, cs.levels, inline=True, fontsize=10)
        cs2 = plt.contour(xs,ys,dec,6,colors='blue')
        if degdecflag == False:
            cs2.levels = [self.makePos(0.0,val)[1][0:-3] for val in cs2.levels]
        plt.clabel(cs2, cs2.levels, inline=True, fontsize=10)
        plt.title(self.apername + ' Rotation = {}$^o$'.format(self.rotation))
        plt.text(xsize-248, ysize+20, 'RA', fontsize=14,color='red')
        plt.text(xsize-148, ysize+20, 'Dec', fontsize=14,color='blue')
        plt.xlabel('X Pixel')
        plt.ylabel('Y Pixel')
        if self.output is None:
            self.output = self.apername + '_rotation' + str(self.rotation) + 'deg_RADec_vs_XY_map.pdf'
        plt.savefig(self.output)

        
    def parseRADec(self,rastr,decstr):
        #convert the input RA and Dec strings to floats
        try:
            rastr=rastr.lower()
            rastr=rastr.replace("h",":")
            rastr=rastr.replace("m",":")
            rastr=rastr.replace("s","")
            decstr=decstr.lower()
            decstr=decstr.replace("d",":")
            decstr=decstr.replace("m",":")
            decstr=decstr.replace("s","")

            values=rastr.split(":")
            ra0=15.*(int(values[0])+int(values[1])/60.+float(values[2])/3600.)

            values=decstr.split(":")
            if "-" in values[0]:
                sign=-1
                values[0]=values[0].replace("-"," ")
            else:
                sign=+1
            dec0=sign*(int(values[0])+int(values[1])/60.+float(values[2])/3600.)
            return ra0,dec0
        except:
            print("Error parsing RA,Dec strings: {} {}".format(rastr,decstr))
            sys.exit()

    def getDistortionCoefficients(self,table,from_sys,to_sys,aperture):
        '''from the table of distortion coefficients, get the coeffs that correspond
        to the requested transformation and return as a list for x and another for y
        '''
        match = table['AperName'] == aperture
        if np.any(match) == False:
            print("Aperture name {} not found in input CSV file.".format(aperture))
            sys.exit()

        row = table[match]

        if ((from_sys == 'science') & (to_sys == 'ideal')):
            label = 'Sci2Idl'
        elif ((from_sys == 'ideal') & (to_sys == 'science')):
            label = 'Idl2Sci'
        else:
            print("WARNING: from_sys of {} and to_sys of {} not a valid transformation.".format(from_sys,to_sys))
            sys.exit()
        
        #get the coefficients, return as list
        X_cols = [c for c in row.colnames if label+'X' in c]
        Y_cols = [c for c in row.colnames if label+'Y' in c]
        x_coeffs = [row[c].data[0] for c in X_cols]
        y_coeffs = [row[c].data[0] for c in Y_cols]

        #get the V2,V3 values of the reference pixel
        v2ref = row['V2Ref'].data[0]
        v3ref = row['V3Ref'].data[0]

        #get the x,y values of the reference pixel
        xref = row['XSciRef'].data[0]
        yref = row['YSciRef'].data[0]

        #get the array size as well
        xsize = row['XSciSize'].data[0]
        ysize = row['YSciSize'].data[0]
        
        return x_coeffs,y_coeffs,v2ref,v3ref,xref,yref,xsize,ysize

        


    def read_coeff_file(self,filename,apername):
        #read in the distortion coefficients
        if os.path.isfile(filename):
            distortionTable = ascii.read(filename,header_start=1)
        else:
            print("WARNING: Input distortion coefficients file {} does not exist.".format(filename))
            sys.exit()

        #read in coefficients for the forward 'science' to 'ideal' coordinate transformation.
        #'science' is in units of distorted pixels, while 'ideal' is the undistorted
        #angular distance from the reference pixel

        x_sci2idl,y_sci2idl,v2_ref,v3_ref,x_ref,y_ref,xsize,ysize = self.getDistortionCoefficients(distortionTable,'science','ideal',apername)
            
        return x_sci2idl,y_sci2idl,v2_ref,v3_ref,x_ref,y_ref,xsize,ysize


    def get_coord_transform_model(self,filename):
        #read in asdf distortion transformation reference file
        #Read in the CRDS-format distortion reference file
        with AsdfFile.open(filename) as dist_file:
            coord_transform = dist_file.tree['model']

        return coord_transform

        

    def XYToRADec(self,pixelx,pixely,attitude_matrix,coord_transform,refpixx,refpixy,refpixv2,refpixv3):
        #Translate a given x,y location on the detector
        #to RA,Dec

        #Transform distorted pixels to V2,V3
        deltav2,deltav3 = coord_transform(pixelx-refpixx,pixely-refpixy)
        pixelv2 = deltav2 + refpixv2
        pixelv3 = deltav3 + refpixv3
        
        #Now translate V2,V3 to RA,Dec
        ra,dec = rotations.pointing(attitude_matrix,pixelv2,pixelv3)

        #Translate the RA/Dec floats to strings
        #ra_str,dec_str = self.makePos(ra,dec)

        return ra,dec#,ra_str,dec_str


    def makePos(self,alpha1,delta1):
        #given a numerical RA/Dec pair, convert to string
        #values hh:mm:ss
        if alpha1 < 0.: 
            alpha1=alpha1+360.
        if delta1 < 0.: 
            sign="-"
            d1=abs(delta1)
        else:
            sign="+"
            d1=delta1
        decd=int(d1)
        value=60.*(d1-float(decd))
        decm=int(value)
        decs=60.*(value-decm)
        a1=alpha1/15.0
        radeg=int(a1)
        value=60.*(a1-radeg)
        ramin=int(value)
        rasec=60.*(value-ramin)
        alpha2="%2.2d:%2.2d:%7.4f" % (radeg,ramin,rasec)
        delta2="%1s%2.2d:%2.2d:%7.4f" % (sign,decd,decm,decs)
        alpha2=alpha2.replace(" ","0")
        delta2=delta2.replace(" ","0")
        return alpha2,delta2




    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Generate RA,Dec map across detector')
        parser.add_argument("apername",help="Name of aperture to use. Follows AperName in Colin's spreadsheet")
        parser.add_argument("dist_model",help='Distortion model reference file (ASDF file)')
        parser.add_argument("ra",help="right ascention at the reference location of the aperture")
        parser.add_argument("dec",help="declination at the reference location of the aperture")
        parser.add_argument("rotation",help="telescope rotation of the field of view. (Degrees)",type=np.float)
        parser.add_argument("--refpixx",help="reference location x coordinate",default=None)
        parser.add_argument("--refpixy",help="reference location y coordinate",default=None)
        parser.add_argument("--refpixv2",help="reference location v2 coordinate",default=None)
        parser.add_argument("--refpixv3",help="reference location v3 coordinate",default=None)
        parser.add_argument("--full_coeff_file",help="csv file with full transformation coeffs.",default=None)
        parser.add_argument("--output",help="Name of file to output map to.")
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: RADec_vs_xy_map.py apername dist_model ra dec rotation'

    map = RADecMap()
    parser = map.add_options(usage = usagestring)
    args = parser.parse_args(namespace=map)

    if map.refpixx is None and map.full_coeff_file == None:
        print("WARNING: you need to define either the aperture's reference location in x,y,v2, and v3 manually, or")
        print("enter the name of the coefficient file that that information can be read from. Quitting.")
        sys.exit()
    
    map.run()
