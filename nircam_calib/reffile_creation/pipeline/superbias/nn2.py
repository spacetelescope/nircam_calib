#!/usr/bin/env python

# File: nn2.py
# Author:  Michael Wood-Vasey <wmwood-vasey@cfa.harvard.edu>
# Date Crated:  2005 March 04
# Last Changed: 2009 Nov 9th
#
# CHANGELOG - GSN - modified to work with numpy instead of numarray
#
# Based on C code by John Tonry and Megan Novicki.
#
# Requires the 'numarray', 'numarray.linear_algebra', and 'pylab'/'matplotlib'
#  package(s) to be installed.
#
# Verified on
#   Python   version: 2.3-2.4
#   numarray version: 1.1-1.3
# matplotlib version: 0.80
#
# should work for higher versions of each
#
# Purpose:  This Python script reads in a lightcurve file of flux differences
#  between all (or an appropriate subset of all) N*(N-1)/2 difference fluxes
#  for a given object from a series of N observations.  The script then
#  creates the matrix of all known pairwise differences and solves for the
#  vector representing the lightcurve of the object.  The zero flux level
#  of the lightcurve can be with respect to a particular date or range
#  of dates.  By default the lightcurve is normalized to the last date
#  point.
#
#  The input lightcurve file must be in flux units as that is the natural space to do this calculation
#  The flux units can have an accompanying zeropoint in the last column
#  otherwise they are assumed to have been normalized to a zeropoint of 25.

import sys, re
from types import *

from math import *
from numpy import *
#from numarray import *
import numpy.linalg as la
#import numarray.linear_algebra as la
from optparse import OptionParser

#from numarray.ieeespecial import isfinite
def isfinite(val):
    import types
    """Returns True if value is finite, else returns False"""
    if type(val) is types.StringType:
        if val in ['nan', 'inf']:  return False

    # Use 'repr' because the following doesn't work due to a bug in Python 2.3
    if repr(val) in ['nan', 'inf']:  return False

    # Comparison with NaN is supposed to always return False
    # But under Python 2.3 it does not. This is unforgivable, but life.
    if not val < 0 and not val >= 0:
        return False
    else:
        return True


def readLightcurve(file):
    observation = []
    mjd = []
    flux = []
    dflux = []
    dfluxext = []
    dfluxint = []

    comment = re.compile('\s*#')

    header = []
    data = []

    try:
        for line in open(file).readlines():
            if comment.match(line):
                header.append(line)
                continue

            (m,f,df,dfe,dfi) = line.strip().split()

            if not isfinite(df):
                df = 1e6
                if isfinite(dfe):  df = dfe
                if isfinite(dfi):  df = dfi

            observation.append(o)
            mjd.append(float(m))
            flux.append(float(f))
            dflux.append(float(df))
            dfluxext.append(float(dfe))
            dfluxint.append(float(dfi))
    except:
        return (observation, mjd, flux, dflux, dfluxext, dfluxint)

    return (observation, mjd, flux, dflux, dfluxext, dfluxint)

def unique(seq):
    d = {}
    for x in seq:
       d[x] = 1
    return d.keys()

def normflux25(flux, zp, zpnorm=25):
    return flux * 10**(-0.4*(zp-zpnorm))

### 2005/03/07 - MWV:
# Ported from C 'antivec.c' written by John Tonry
#
# Fits an antisymmetric matrices as a vector term difference:
#        A[i][j] = V[i] - V[j]
# Inputs are a difference matrix and an error matrix.
# Outputs vector V and error vector for V.
# The first term of the output vector is set to zero.
#
def antivec(A, E, verbose=False):
    #import numarray.linear_algebra.lapack_lite as la

    floattype = 'Float64'

    # Number of terms
    n = len(A) # A is by definition square, so we just need one dimension
    # Vector we're going to return
    V = zeros(n,floattype)
    # Error vector we're going to return
    dVext = zeros(n,floattype)  # exterior
    dVint = zeros(n,floattype)  # interior

    # Degrees of freedom
    dof = (n*(n-1))/2 - (n-1)
    # Estimate an average weight
    wgtave = 0.0
    if verbose:  print "Calculating wgtave:"
    for j in range(n):
        for i in range(j):
#            print i, j, E[i,j]
            if E[i,j] > 0:
                wgtave += 1/(E[i,j]*E[i,j])
    # Normalize the weight
    wgtave /= (n*(n-1))/2;

    if verbose:  print "N: ", n
    if verbose:  print "wgtave: ", wgtave

    cov = zeros([n,n],floattype)     # Covariance matrix
    y   = zeros(   n ,floattype)     # work vector
    for j in range(n):
        for i in range(n):
#            wgt = (E[i,j] != 0 and i != j) ? 1/(E[i,j]*E[i,j]) : 0;
#            if i == j:  continue
            if E[i,j] > 0:
                wgt = 1.0/(E[i,j]*E[i,j])
            else:
                wgt = 0.0;

#            print i, j, E[i,j], wgt

            cov[j,j] += wgt*wgt
            if i != j:
               cov[i,j] = wgt*wgt
               y[j] += wgt

    if verbose:
       #print "External weight matrix -- Determinant: %f " % la.determinant(cov)
       print "External weight matrix -- Determinant: %f " % la.det(cov)
       print cov

    # Now we do a little inversion

    #GSN - 20091109 - modify to work with numpy
    #cov = la.inverse(cov)
    cov = la.inv(cov)
    #det = la.determinant(cov)
    det = la.det(cov)

    if verbose:
       print "External weight matrix after inversion -- Determinant: %f" % det
       print cov

    # Calculate something
    for j in range(n):
        for i in range(n):
            if verbose:  print "y[%2d]:  %f, cov[%2d,%2d]:  %f" % (i,y[i],i,j,cov[i,j])
            dVext[j] += y[i] * cov[i,j]

    for i in range(n):
        if verbose:
            print "dVext[%d]: %f" % (i, dVext[i])
        # 2006/02/25 - MWV:
        #   I don't entirely understand whether it's legitimate to
        #   take the sqrt of the absolute value, as Tonry does
        #   in his C code (as was pointed out to me by Drew Newman).
        dVext[i] = sqrt(dVext[i])
#        dVext[i] = sqrt(abs(dVext[i]))

    # Calculate the best-fit vector to the antisymmetric matrix
    #  cov[j][i] = 1/e[j+1][i+1]^2 off diagonal, and Sum(i!=j) -1/(e*e) on diag
    #  y[i] = Sum(j) a[j][i+1]/e[j][i+1]^2)
    # Add in a factor of (Sum(j) v(j))^2 to chi^2 to keep it non-singular (and
    #   force the Sum(j) v(j) = 0)
    cov[:,:] = wgtave  # This seems extraneous given the line below
    cov[:,:] = 1.0
    y[:] = 0.  # Zero out work vector
    for j in range(n):
        for i in range(n):
#           wgt = (E[i,j] != 0 and i != j) ? 1/(E[i,j]*E[i,j]) : 0
            if E[i,j] != 0 and i != j:
               wgt = 1/(E[i,j]*E[i,j])
            else:
               wgt = 0.0

            cov[i,j] += wgt
            cov[j,j] -= wgt    # Add and sub cancel when i==j
            y[j] += A[j,i] * wgt

    # Solve and multiply out to evaluate V
    if verbose:
        print "Hessian matrix: <w> = %f" % wgtave
        print cov

    # Another inversion of our now different covariance matrix
    #GSN - 20091109 - modify to work with numpy
    #cov = la.inverse(cov)
    #det = la.determinant(cov)
    cov = la.inv(cov)
    det = la.det(cov)

    for j in range(n):
        for i in range(n):
            V[j] += y[i] * cov[i,j]

    # Evaluate the consistency for V from the covariance matrix now in cov
    for j in range(n):
        dVint[j] = cov[j,j] = sqrt(-cov[j,j])
        for i in range(j):
            cov[i,j] = cov[j,i] = -cov[i,j] / (cov[i,i]*cov[j,j])

    # Collect chi^2 from the residual matrix
    chi = 0.0
    for j in range(n):
        for i in range(j+1,n):
            if E[i,j] > 0:
                chi += (A[i,j] - (V[j]-V[i])) * (A[i,j] - (V[j]-V[i])) / (E[i,j] * E[i,j])

 ### Tonry's comments:
 # Adjust the balance between external and internal error
 # We assume that some fraction of the error matrix comes from "external
 # error", i.e. errors in the vector values which are consistent between
 # all the difference measurements, hence invisible to the antisymmetric
 # matrix, and the rest comes from "internal error" which arose from the
 # process of getting the terms of the antisymmetric matrix.  The former
 # errors should be analyzed using the formalism for the "dvext" and the
 # latter using the "dvint" formalism.  Roughly speaking the dvext errors
 # are sqrt(2) smaller than the e terms and the dvint errors are sqrt(1/N)
 # smaller than the e terms.  We adjust the fraction so that chi/N = 1.
 # We expect normally that chi/N < 1, i.e. the difference vector matches
 # the antisymmetric matrix quite well, and the uncertainties in the e
 # vector are mainly arising from external uncertainties in v.  However,
 # there are at least two ways to partition things.
 #
 # Case 1: we take the e matrix at face value as telling us about the
 # external errors, use the full thing for the dvext, and scale it as much
 # as necessary for dvint:
 #
 #			dvext -> dvext
 #			dvint -> dvint * sqrt(chi/N)
 #
 # Case 2: we partition the e matrix into an internal piece and an external
 # piece, with the internal piece adjusted to make chi/N = 1 and the
 # external piece having a floor of zero.
 #
 #			dvext -> dvext * MAX(0, 1-sqrt(chi/N))
 #			dvint -> dvint * sqrt(chi/N)
 #
 # I'm dubious about Case 2, so I'm just going to report things as Case 1.
 # Users who may have error matrices which include a large piece of
 # internal matrix error (instead of external vector error) may want to
 # rescale the dvext and dvint before adding in quadrature to get a final
 # uncertainty.
 #

    for j in range(n):  dVint[j] *= sqrt(chi/dof)

    if verbose:
        print "Error matrix:"
        print cov

    return (V, dVext, dVint)


def NN2version():
    return "0.15"

# Container class for matrices
class nn2:
    def __init__(self, dates=[], A=[], E=[]):
        self.dates = dates
        self.A = A
        self.E = E

        self.V = []
        self.Verr = []
        self.Vext = []
        self.Vint = []

        # Type of floating-point values
        self.floattype = 'Float64'

    def solveMatrix(self, verbose=False, debug=False):
        #(V, S, WT) = la.singular_value_decomposition(self.A)
        (V, S, WT) = la.svd(self.A)
        (self.V, self.dVext, self.dVint) = antivec(self.A, self.E, verbose=verbose)
        if debug:
            print "V  vector: ", V
            print "S  matrix: ", S
            print "WT matrix: ", WT

        self.Verr = sqrt(self.dVint**2 + self.dVext**2)

    ### Use the last data point (in time) to set the zero flux
    ###   Should eventually make checks to decide about SN
    ###   or setting options
    def normToLastPoint(self):
        self.V -= self.V[-1]

    ### Use the data points after minMJD to set the zero flux
    ###   Should eventually make checks to decide about SN
    ###   or setting options.
    ###
    ### zerotype: "mean", "clipmean"
    def normToPoints(self, minMJD=None, maxMJD=None, zerotype="mean", clipsigma=3.):
        import numpy.core.ma as ma
        zerodates = None

        # Make sure things that should be floats if defined are floats
        if minMJD is not None:  minMJD = float(minMJD)
        if maxMJD is not None:  maxMJD = float(maxMJD)
        if clipsigma is not None:  clipsigma = float(clipsigma) # should always be defined

        if minMJD is None and maxMJD is None:
            self.normToLastPoint()
            return
        else:
            if minMJD is not None and maxMJD is not None:
                if minMJD < maxMJD:
                    zerodates = [i for (i, d) in enumerate(self.dates) if (d >= minMJD) and (d <= maxMJD)]
                else:
                    zerodates = [i for (i, d) in enumerate(self.dates) if (d >= minMJD) or (d <= maxMJD)]
            elif minMJD is not None:
                zerodates = [i for (i, d) in enumerate(self.dates) if d >= minMJD]
            elif maxMJD is not None:
                zerodates = [i for (i, d) in enumerate(self.dates) if d <= maxMJD]

        if zerodates is None or len(zerodates) <= 0:
            print "NN2: No points found in given MJD range to use as zero points for this filter."
            print "NN2: Leaving points unchanged."
            return

        # There's got to be a better way of doing the following masked array stuff
        zeroflux_mask    = [not isfinite(f ) for  f in self.V[zerodates]]
        zerofluxerr_mask = [not isfinite(df) for df in self.Verr[zerodates]]

        zeroflux    = ma.array(self.V[zerodates]   , mask = zeroflux_mask )
        zerofluxerr = ma.array(self.Verr[zerodates], mask = zerofluxerr_mask)

        # Now finally adjust the flux by the zeroflux calculated.
        #  If there's only one entry then we just take it
        if len(zeroflux) == 1:
            zerofluxnorm = zeroflux[0]

        elif zerotype == "mean":
#            zerofluxnorm = ma.average(zeroflux, weights=1./(zerofluxerr*zerofluxerr))
            zerofluxnorm = ma.average(zeroflux, weights=1./zerofluxerr*zerofluxerr)
        elif zerotype == "clipmean":
            zeroavg = ma.average(zeroflux, weights=1./(zerofluxerr*zerofluxerr))
            residual = zeroflux - zeroavg
            residual /= zerofluxerr # error-weighted

            # 2006/02/13 - MWV:
            #   Need the additional check to make sure the value isn't masked
            w = [i for (i,f) in enumerate(ma.fabs(residual)) \
                 if f is not ma.masked and f < clipsigma]

            if len(w) <= 0:
                print "No good zero flux points"
                print "No flux zero applied to points."

            try:
                zerofluxnorm  = ma.average(zeroflux[w], weights=1./(zerofluxerr[w]*zerofluxerr[w]))
                this_resid = residual[w] * residual[w]
                zerofluxchisq = this_resid.sum()
                print zerofluxchisq
                dof = len(residual) - 1
                zerofluxchisqnu = zerofluxchisq / dof
                print "Chi^2/DoF of zero points :  %f / %f  = %f " % ( zerofluxchisq, dof, zerofluxchisqnu )
            except 'DIE':
                zerofluxnorm = 0.
                zerofluxchisq = 0.
                dof = 0.

        else:
            print "zerotype == '%s' unknown or not yet supported." % (zerotype)
            print "No flux zero applied to points."
            zerofluxnorm = 0.

        # Subtract the fiducial flux
        self.V -= zerofluxnorm

class nn2plot(nn2):
    '''Plotting class for NN2 analysis.'''
    try:
       import pylab

       def num2date(self, mjd):
           import pylab
           return pylab.num2date(mjd)

    except:
       def num2date(self, mjd):
           return str(mjd)

    def mjd2date(self, mjd):
        '''Converts the input Modified Julian Day into YYYY/MM/DD format'''
        # There's probably some great general way to do this, but
        #  I just put in my own offset here because it seemed the simplest
        # Takes advantage of pylab's num2date
        # Offset from 0001-01-01 00:00:00 UTC and MJD
#        mjd_ce = +678575.0
        ### 2006/02/27 - MWV:
        ###   Updated according to the output of 'juldate' from IDL GSFC-Astro
        ###      IDL> juldate, [1,1,1], reduced_jd
        ###      IDL> mjd = reduced_jd - 0.5
        ###      IDL> print, mjd
        ###            -678577.00
        ### Note that the documentation for num2date is currently wrong.
        ### It says that num2date accepts an argument given the number of
        ### days since 0001/01/01 00:00:00, but it actually treats
        ### the input argument minutes 1 as the number of days since
        ###   0001/01/01 00:00:00.
        ### This simple but critical mistake in the pylab.num2date
        ###   documentation does not inspire confidence.
        ### In any event, that means that we subtract one from our MJD
        ###   base offset, which would otherwise be +678577.00
        ###
        mjd_ce = +678576.00

        #print [float(m)+mjd_ce for m in mjd]
        dateStrArr = [self.num2date(float(m)+mjd_ce).strftime("%Y/%m/%d") for m in mjd]
        return dateStrArr

    def __init__(self, *args):
        nn2.__init__(self, *args)

    def plotline(self, **args):
        import pylab
        mjd = [float(f) for f in self.dates]
        pylab.errorbar(mjd, self.V, self.Verr, **args)

    def visualize(self):
        import pylab
        '''Generates visualizationss of the A and E matrices.'''
        # Construct a matrix that is upper-triangle flux and lower-triangle error
        fluxAndError = self.A.copy()
        for i in range(len(self.A)):
            fluxAndError[i,0:i] = self.E[i,0:i]

#        pylab.clf()
        pylab.matshow(fluxAndError, cmap=pylab.cm.gray)
        l = len(self.A)
        pylab.plot([0,l],[0,l])
        pylab.xlim(0,l)
        pylab.ylim(0,1)
        pylab.gray()
        pylab.xlabel("Flux matrix: A(i,j)")
#        pylab.ylabel("Error matrix: E(i,j)")
        # I'm abusing the 'title' attribute here just to get the x2label where I want it
        pylab.title("Error matrix: E(i,j)")
        print self.dates
        datelabels = self.mjd2date(self.dates)
        datelabels.reverse()
        pylab.yticks(arange(len(self.A),0,-1),(datelabels))
#        pylab.xticks(arange(len(self.A))+0.5,(datelabels), rotation='vertical')
#        pylab.yticks(arange(len(self.A)),(self.gooddates))
#        pylab.xticks(arange(len(self.A)),(self.gooddates))
        pylab.show()


    # Remove dates with all 0 in A and/or E
    # Remove dates with one or fewer non-zero entries in a row/column
    # We also need something above and under the diagonal for each column
    #   for the matrix to be invertible
    def pruneMatrix(self, verbose=False):
        n = len(self.A)

        if verbose:  print "Pruning matrix"
        if verbose:  print "N: ", n

        # Construct our list of bad dates
        #    sum(self.A[i,:] == zeros(n))
        #  is the number of zero-element entries in the ith row of the A matrix
        #  If that quantity is greater than n-1, then all elements are 0
        #  I don't quite remember why I did it this way,
        #    but there must have been a reason
        baddates = []
        numbad = 0
        for i in range(n):
            baddates.append(False)
            if sum(self.A[i,:] == zeros(n)) >= n-1:
                baddates[i] = True
                numbad += 1

        if verbose:  print "Bad: ", baddates, "numbad: ", numbad

        numgooddates = n-numbad

        # Create new A, E matrices of only good values
        oldA = self.A.copy()
        oldE = self.E.copy()
        olddates = self.dates
        self.origA = oldA.copy()
        self.origE = oldE.copy()
        self.origdates = olddates
        newA = zeros([numgooddates,numgooddates], self.floattype)
        newE = zeros([numgooddates,numgooddates], self.floattype)
        self.gooddates = []
        newn = 0
        good = []
        for (oldi, bad, date) in zip(range(n), baddates, olddates):
            if bad: continue
            newn += 1
            if verbose:  print "newn: %d" % newn
            good.append(oldi)
            self.gooddates.append(date)
            for i in range(newn):  newA[i,newn-1] = oldA[good[i],newn-1]
            for j in range(newn):  newA[newn-1,j] = oldA[newn-1,good[j]]
            for i in range(newn):  newE[i,newn-1] = oldE[good[i],newn-1]
            for j in range(newn):  newE[newn-1,j] = oldE[newn-1,good[j]]



#         good = [0,1]
#         newn = 2
#         for i in range(newn):  newA[i,newn-1] = oldA[good[i],newn-1]
#         for j in range(newn):  newA[newn-1,j] = oldA[newn-1,good[j]]
#         for i in range(newn):  newE[i,newn-1] = oldE[good[i],newn-1]
#         for j in range(newn):  newE[newn-1,j] = oldE[newn-1,good[j]]

#         print "Checking determinants: %d, %d" % (newn, n)
#         for ind in range(newn,n):
#             if baddates[ind]:  continue

#             if self.verbose:  print "prunematrix: Checking row/column %d, newn %d" % (ind, newn)

#             good.append(ind)
#             newn += 1
#             print "Good: ", good
#             for i in range(newn):  print "%d  [%d,%d] %f " % (i, good[i], ind, oldA[good[i],ind])
#             for i in range(newn):  newA[i,newn-1] = oldA[good[i],ind]
#             for j in range(newn):  newA[newn-1,j] = oldA[ind,good[j]]
#             for i in range(newn):  newE[i,newn-1] = oldE[good[i],ind]
#             for j in range(newn):  newE[newn-1,j] = oldE[ind,good[j]]

#             if self.verbose:
#                 print "Latest Matrix A -- newn %d -- determinant %f: " % (newn, la.determinant(newA))
#                 print newA[0:newn,0:newn]

#             # Check that the determinant is non-zero
#             if la.determinant(newA[0:newn,0:newn]) == 0 and la.determinant(newE[0:newn,0:newn]):
#                 good.pop()
#                 newn -= 1
#                 continue

#             self.gooddates.append(self.alldates[ind])


        self.A = newA[0:newn,0:newn]
        self.E = newE[0:newn,0:newn]
        self.dates = self.gooddates




# This module reads in a difference lightcurve file,
# constructs the matrix A[i,j], where A[i,j] = DifferenceFluxFromSub[i,j],
# and then solves for V[i].
#
class nn2analyze:
    def __init__(self, args, usage=None):

        if usage is None:
            usage= """
    %prog lightcurvefile [options]

    Creates a lightcurve based on N*(N-1)/2 method of Tonry et al.

    Exit codes:
    101  --  Couldn't load lightcurve file
    111  --  Singular matrix encountered
    """

        self.usage = usage
        self.version = NN2version() # Version number

        self.defineoptions()

        self.date        = {}
        self.flux        = {}
        self.dflux       = {}
        self.zeropoint   = {}
        self.photcode    = {}
        self.refdate     = {}
        self.observation = {}
        self.gooddates   = {}  # The "good" date-pairs for the lightcurve

        self.matrix      = {}

        # Type of floating-point values
        self.floattype = 'Float64'

    def clear(self):
        self.date        = {}
        self.flux        = {}
        self.dflux       = {}
        self.zeropoint   = {}
        self.photcode    = {}
        self.refdate     = {}
        self.observation = {}
        self.gooddates   = {}

        self.matrix = {}

    def loadFile(self, lcfile=None):

        if lcfile==None:  lcfile=self.lcfile

        try:
            for line in open(lcfile).readlines():
               if re.search('#',line):  continue

               arr = line.strip().split()
               if len(arr) != 6:  continue

               # Dummy values
               obs   = arr[0]+arr[-1]
               templ = arr[-1]
               # Get variables and cast types from line
               d = float(arr[0])
               p = arr[1]
               f = float(arr[2])
               df= float(arr[3])
               zp= float(arr[4])
               rd= float(arr[5])

               # Skip if dflux or zeropoint is indicative of a bad point
               if df <= 0.0:  continue
               if zp <= 0.0:  continue
               # /*** HACK ALERT ***/
               # 2005/05/05 - Michael Wood-Vasey
               # Should break this out into some function that cleans up
               # but then I should generalize the stuff below as well
               # to make it easy to reject lines
               # 2005/05/05 - MWV.  Disabled because it's breaking the
               #   links between the different image pairs
 #              self.zp_min = 30.0
#               if zp <= zp_min: continue

               # Keep MJD dates as strings for now to use as hashes later if needed

               if not self.refdate.has_key(p):  self.refdate[p] = []
               if not self.date.has_key(p):  self.date[p] = []

               self.refdate[p].append(rd)
               self.date[p].append(d)
               self.flux[(rd,d,p)]      = f
               self.dflux[(rd,d,p)]     = df
               self.zeropoint[(rd,d,p)] = zp
               self.photcode[(rd,d,p)]  = p
               self.observation[(d,p)]    = obs
               self.observation[(rd,p)]   = templ

        except e:
            print e
            print "Unable to process file:  %s" % self.lcfile
            return False

        if self.options.verbose:  print self.date
        if self.options.verbose:  print "Flux: ", self.flux
        if self.options.verbose:  print "dFlux: ", self.dflux
        if self.options.verbose:  print "zeropoint: ", self.zeropoint

        # Create some derivative arrays
        self.alldates = {}
        self.uniqphotcodes = unique(self.photcode.values())
#        print "photcode: ", self.uniqphotcodes
        for p in self.uniqphotcodes:
            if self.options.verbose:  print "p: ", p
            self.alldates[p] = unique(self.date[p]+self.refdate[p])
            self.alldates[p].sort()
            if self.options.verbose:  print "sorted all dates: ", self.alldates[p]
            self.gooddates[p] = self.alldates[p]  # All dates good by default

        # For now just make observation the same as gooddates.  Eventually, it needs to be the filename or something
#        self.observation[p] = self.gooddates

        return True

    # Construct A_ij
    def createMatrix(self, p):
        # Build a dictionary
        index = {}
        for i in range(len(self.alldates[p])):
             index[self.alldates[p][i]] = i
        if self.options.verbose:  print len(self.alldates[p])

        dates   = [None for i in range(len(self.alldates[p]))] # Create blank date array
        rddates = [None for i in range(len(self.alldates[p]))] # Create blank date array
        A = zeros([len(self.alldates[p]),len(self.alldates[p])], self.floattype)
        E = zeros([len(self.alldates[p]),len(self.alldates[p])], self.floattype)

        #        for i in self.alldates:
        #            for j in self.alldates:
        #                if i!=j:
        #                    if self.options.verbose:  print i,j
        #                    self.A[i,j] = normflux25(self.flux[(i,j)] , self.zeropoint[(i,j)])
        #                    self.E[i,j] = normflux25(self.dflux[(i,j)], self.zeropoint[(i,j)])

        if self.options.verbose:  print "Filling A, E matrices:"
        for (rd,d) in zip(self.refdate[p],self.date[p]):
            i = index[rd]
            j = index[ d]
            zp = self.zeropoint[(rd,d,p)]
            if self.options.verbose:  print i,j,zp, d, rd
            if zp > 0:
                # We have to do both of these to make
                # sure we get the dates for both images in the subtraction
                ### 2005/06/21 - Michael Wood-Vasey :
                ####   Fixed bug here that was reversing the date
                ####   storage below.  This wasn't affecting
                ####   solution per se; it was causing the dates that
                ####   were printed with the output lightcurve to be wrong
                ####   in the case of incomplete matrix coverage.
                dates[i]   =  rd
                rddates[j] =   d

                A[i,j] =  normflux25(self.flux[(rd,d,p)] , zp)
                A[j,i] = -normflux25(self.flux[(rd,d,p)] , zp)
                E[i,j] =  normflux25(self.dflux[(rd,d,p)], zp)
                E[j,i] =  normflux25(self.dflux[(rd,d,p)], zp)


        for i in range(len(dates)):
            if dates[i] is None:
                if self.options.verbose:
                     print "setting dates[%d] = rddates[%d]" % (i, i), dates[i], rddates[i]
                dates[i] = rddates[i]

        if self.options.verbose:  print "createMatrix: dates: ", dates

        self.matrix[p] = nn2plot(dates, A, E)

        # Set error weights to very large for <= 0 values
        # Should be formally infinity
#        n = len(self.A)
#        for i in range(n):
#            for j in range(n):
#                if self.E[i,j] <= 0:
#                    self.E[i,j] = 1e30

        # Prune the matrix of bad dates
        if self.options.prune:  self.pruneMatrices()

#        if self.options.verbose:  print "origA: ", self.origA
#        if self.options.verbose:  print "origE: ", self.origE

#        if self.options.verbose:  print "goodates: ", self.gooddates
        if self.options.verbose:  print "A: ", self.matrix[p].A
        if self.options.verbose:  print "E: ", self.matrix[p].E


    ### Prune isolated dates from matrices
    def pruneMatrices(self):
        for k in self.matrix.keys():
                self.matrix[k].pruneMatrix(verbose=self.options.verbose)
        # We do this twice to see if there was a date left out
        #  that's now isolated or has only one.
        # Clearly this could go on forever, but we'll say two is good
        for k in self.matrix.keys():
                self.matrix[k].pruneMatrix(verbose=self.options.verbose)

    def LCheader(self):
        str = ''
        str = "# Flux normalized to zeropoint = 25\n"
        str += "# %-10s %-12s %12s %12s %12s %12s\n" % ("MJD", "photcode", "flux", "flux err (tot)", "flux err (ext)", "flux err (int)")
        return str


    def LCstring(self, photcode=None):
        if photcode is None:  photcode = self.uniqphotcodes
        if not isinstance(photcode, list):  photcode = [photcode]

        str = self.LCheader()

        for p in photcode:

            for (d,f,df,dfext,dfint) in zip(self.matrix[p].dates, self.matrix[p].V, self.matrix[p].Verr, self.matrix[p].dVext, self.matrix[p].dVint):
                if (d,p) in self.observation.keys():
                    o = self.observation[(d,p)]
                else:
                    o = d
                str += "%-12.4f %-30s %12.4f %12.4f %12.4f %12.4f\n" % (d,p,f,df,dfext,dfint)

        return str

    def printLightcurve(self, photcode=None):
#        print "Printing lightcurve [MJD, flux, flux err]"
#        for (d,f,df) in zip(self.gooddates, self.V, self.Verr):
#            print d,f,df
#        print "Printing lightcurve [MJD, flux, flux err (ext), flux err (int)]"
#        for (d,f,df,dfint) in zip(self.gooddates, self.V, self.dVext, self.dVint):
#            print d,f,df,dfint
        print "Printing lightcurve: "
        print self.LCstring(photcode)

    def outputFilename(self):
        self.outfile = self.lcfile.replace('.lc.dat','.nn2.dat')
        self.outfile = self.outfile.replace('.dat','.nn2.dat')
        self.outfile = self.outfile.replace('.lc','.nn2.dat')
        # If we get the same output name then add on an additional '.nn2.dat'
        #   just to make sure we don't overwrite it.
        if self.outfile == self.lcfile:  self.outfile += '.nn2.dat'

        return self.outfile

    def saveLightcurve(self):
        self.outputFilename()
        print "Saving lightcurve to '%s'" % self.outfile
        open(self.outfile,'w').write(self.LCstring())

    def plotLightcurve(self, photcode=None):
        import pylab
        if photcode is None:  photcode = self.uniqphotcodes
        if not isinstance(photcode, list):  photcode = [photcode]

        print photcode
        for p in photcode:
            print "Plotting ", p
            self.matrix[p].plotline()

        pylab.xlabel("MJD")
        pylab.ylabel("Flux (normalized to zp=25)")
        pylab.show()
        pylab.clf()

        ### Get options ###
    def defineoptions(self):
        pythonversion = (sys.version.split())[0]
        if pythonversion >= '2.4':
           defaultstr = " [default: %default]"
        else:
           defaultstr = ""

           self.usage += """
[If you want the default values of the options to show up in the list below
 then you should upgrade to Python 2.4 or later]"""

        parser = OptionParser(usage=self.usage, version=self.version)
        parser.add_option('-v', '--verbose', dest='verbose', default=False,
                          action="store_true",
                          help='Display extra verboseging output. '+defaultstr)
        parser.add_option('-s', '--show', dest='show', default=False,
                          action="store_true",
                          help='Plot the derived lightcurve. '+defaultstr)
        parser.add_option('--visualize',  default=False, action="store_true",
                          help='Plot the flux and error matrices. '+defaultstr)
        parser.add_option('--prune', dest='prune', default=False, action="store_true",
                          help='Prune out isolated dates or dates with only one point. '+defaultstr)
        parser.add_option('--noprune', dest='prune', default=True, action="store_false",
                          help='Do not prune out dates.'+defaultstr)

        parser.add_option('-n', '--norm', dest='norm', default=True,
                          action="store_true",
                          help="Use last point as lightcurve zero. "+defaultstr)
        parser.add_option('-r', '--raw', dest='norm',
                          action="store_false",
                          help="Don't use last point as lightcurve zero. "+defaultstr)
        parser.add_option('--zeroMJDmin', default=None, type='float',
                          help="Use average or median (see --zerotype) of flux values from dates >= given MJD as flux zero. "+defaultstr)
        parser.add_option('--zeroMJDmax', default=None, type='float',
                          help="Use average or median (see --zerotype) of flux values from dates <= given MJD as flux zero. "+defaultstr)
        parser.add_option('--zerotype', default='clipmean',
                          help="Specify 'mean', 'clipmean', or 'median' to calculate flux zero. "+defaultstr)
        parser.add_option('--clipsigma', default=3, type='float',
                          help="Specify sigma for zerotype == 'clipmean' option. "+defaultstr)
        parser.add_option('--LCversion', default=None,
                          help="Record the lightcurve version in the header information. "+defaultstr)

        self.parser = parser

    def getoptions(self, inargs):
        options, args = self.parser.parse_args(inargs)

        return (options, args)

    ### A hook for routines that wish to do processing after reading in
    ### args and options but before doing further processing
    def special_run(self):
        pass

    ### Do all of the steps to solve, output, and visualize the solution
    def run(self):
        # We get the options here so that the options can be added
        #  to before being processed, e.g. for the case of a derived class.
        (self.options, self.args) = self.getoptions(sys.argv[1:])

        if len(self.args) <= 0:
           print self.usage
           sys.exit()

        lcfile = self.args[0]
        self.lcfile = lcfile

        if not self.loadFile():
            print "Couldn't load lightcurve file: '%s'" % self.lcfile
            print "Giving up."
#            print "Usage:"
#            print self.usage
            sys.exit(101)

        ### Hook if an inherited class wants to run anything here
        self.special_run()

        for p in self.uniqphotcodes:
            self.createMatrix(p)

        for p in self.uniqphotcodes:
            if self.options.visualize:
                self.matrix[p].visualize()

            ### Actually solve the matrix
            try:
                self.matrix[p].solveMatrix(verbose=self.options.verbose)
            except 'DIE':
                print "Unable to solve matrix.  It's possible the matrix was singular."
                print "Try using the --visualize option "
                print "  to see how well the lightcurve difference matrix is filled in."
                sys.exit(111)

            if self.options.norm:
                print "Subtracting baseline flux for photcode: ", p
                self.matrix[p].normToPoints(minMJD=self.options.zeroMJDmin,
                                            maxMJD=self.options.zeroMJDmax,
                                            zerotype=self.options.zerotype,
                                            clipsigma=self.options.clipsigma)


#        self.printLightcurve()
        self.saveLightcurve()

        if self.options.show:  self.plotLightcurve()

#####################################
###  End of NN2 class definition  ###
#####################################




if __name__=='__main__':

    nn2LC = nn2analyze(sys.argv[1:])

    # run the main script
    nn2LC.run()

    # Clear
    nn2 = None
