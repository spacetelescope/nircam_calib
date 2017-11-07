#!/usr/bin/env python
'''
wrapper around recarray with convenience functions to ease handling of tables
A. Rest
'''
import sys,os,re,types
import numpy as np
from numpy import recarray
import numpy.lib.recfunctions as nprecf

from astropy.io import ascii
from astropy.time import Time
import astropy.io.fits as fits
import astropy

# def makepath(path,raiseError=True):
#     if path == '':
#         return(0)
#     if not os.path.isdir(path):
#         os.makedirs(path)
#         if not os.path.isdir(path):
#             if raiseError:
#                 raise RuntimeError, 'ERROR: Cannot create directory %s' % path
#             else:
#                 return(1)
#     return(0)
#
# def makepath4file(filename,raiseError=True):
#     path = os.path.dirname(filename)
#     if not os.path.isdir(path):
#         return(makepath(path,raiseError=raiseError))
#     else:
#         return(0)



#http://docs.astropy.org/en/stable/table/index.html#astropy-table
#table groups!! http://docs.astropy.org/en/stable/table/operations.html#table-operations
class astrotableclass:
    def __init__(self,**kwargs):
        self.t = astropy.table.Table(**kwargs)
        self.verbose = 0

    def load(self,filename,namesMapping={},formatMapping={},**kwargs):
        #self.t = ascii.read(filename,format='commented_header',delimiter='\s',fill_values=[('-',0),('--',0)])
        try:
            self.t = ascii.read(filename,delimiter='\s')
        except:
            print('ERROR: could not read %s!' % filename)
            return(1)

        if len(self.t.colnames)<1:
            print('ERROR: no data in %s?' % filename)
            return(1)

        if self.t.colnames[0]=='col1':
            if self.verbose: print('WARNING: it looks like there is no header info!')

        if len(namesMapping)>0:
            for name in self.t.colnames:
                if name in namesMapping:
                    self.t.rename_column(name,namesMapping[name])
        if len(formatMapping)>0:
            for name in formatMapping:
                if name in self.t.colnames:
                    self.t[name].format=formatMapping[name]
                else:
                    print('WARNING! col %s does not exist, so cannot format it!' % name)

    def write(self,filename,clobber=True,verbose=False,format='commented_header',makepathFlag=True,**kwargs):
        if verbose:
            print('Saving %s' % filename)

        # make the path to the file if necessary
        if makepathFlag:
            if makepath4file(filename,raiseError=False):
                print('ERROR: could not make directory for %s' % filename)
                return(1)

        # if clobber, then remove the old file first
        if clobber and os.path.isfile(filename):
            os.remove(filename)
            if os.path.isfile(filename):
                print('ERROR: could not save %s' % filename)
                return(2)

        ascii.write(self.t,filename,format=format,**kwargs)

    # def fitsheader2table(self,fitsfilecolname,rowindices=None,requiredfitskeys=None,optionalfitskey=None,raiseError=True,skipcolname=None,headercol=None):
    #     if rowindices==None:
    #         rowindices = xrange(len(self.t))
    #
    #     if len(rowindices)==0:
    #         print('no files!')
    #         return(0)
    #
    #     # initialize columns if necessary
    #     if requiredfitskeys!=None:
    #         for fitskey in requiredfitskeys:
    #             if not (fitskey in self.t.colnames):
    #                 self.t[fitskey]=None
    #     if optionalfitskey!=None:
    #         for fitskey in optionalfitskey:
    #             if not (fitskey in self.t.colnames):
    #                 self.t[fitskey]=None
    #
    #     if headercol!=None and (not (headercol in self.t.colnames)):
    #         self.t[headercol]=None
    #
    #     # loop through the images
    #     for rowindex in rowindices:
    #         header = fits.getheader(self.t[fitsfilecolname][rowindex])
    #         if headercol!=None:
    #             self.t[headercol]=header
    #
    #         if requiredfitskeys!=None:
    #             for fitskey in requiredfitskeys:
    #                 if fitskey in header:
    #                     self.t[fitskey][rowindex]=header[fitskey]
    #                 else:
    #                     self.t[fitskey][rowindex]=None
    #                     if raiseError:
    #                         raise RuntimeError,"fits key %s does not exist in file %s" % (fitskey,self.t[fitsfilecolname][rowindex])
    #                     else:
    #                         if skipcolname!=None:
    #                              self.t[skipcolname][rowindex]=True
    #         if optionalfitskey!=None:
    #             for fitskey in optionalfitskey:
    #                 if fitskey in header:
    #                     self.t[fitskey][rowindex]=header[fitskey]
    #                 else:
    #                     self.t[fitskey][rowindex]=None

    def dateobs2mjd(self,dateobscol,mjdcol,timeobscol=None):
        if not (mjdcol in self.t.colnames):
            self.t[mjdcol]=None

        if timeobscol!=None:
            dateobslist = list(self.t[dateobscol]+'T'+self.t[timeobscol])
        else:
            dateobslist = list(self.t[dateobscol])

        dateobjects = Time(dateobslist,  format='isot', scale='utc')
        mjds = dateobjects.mjd

        self.t[mjdcol]=mjds
