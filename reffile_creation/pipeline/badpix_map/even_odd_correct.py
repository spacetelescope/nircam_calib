#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import os, argparse

qstart = [0,512,1024,1536,2048]

class Oddeven:
    def __init__(self):
        self.infile = None
        self.outfile = None
        self.data = None

    def calc_averages(self,data):
        allevenmean = np.zeros(4)
        alloddmean = np.zeros(4)
        template = np.zeros((2048,2048))
        for colnum in np.arange(2048):
            if colnum % 2 == 0:
                template[0:4,colnum] = 2
                template[2044:,colnum] = 2
            if colnum % 2 == 1:
                template[0:4,colnum] = 1
                template[2044:,colnum] = 1

        for amp in xrange(4):
            amptemplate = template[:,qstart[amp]:qstart[amp+1]]
            evencols = amptemplate == 2
            oddcols = amptemplate == 1
            ampdata = data[:,qstart[amp]:qstart[amp+1]]

            #sigma-clip once at 3-sigma, just to make the mean
            #a little more robust
            evendata = ampdata[evencols]
            evenmean = np.nanmean(evendata)
            evenstd = np.nanstd(evendata)
            #print(evendata)
            #print(qstart[amp],qstart[amp+1],evendata.shape,evenmean,evenstd)
            good = (evendata < (evenmean+3.*evenstd)) & (evendata > (evenmean-3.*evenstd))
            allevenmean[amp] = np.nanmean(evendata[good])
            
            odddata = ampdata[oddcols]
            oddmean = np.nanmean(odddata)
            oddstd = np.nanstd(odddata)
            good = (odddata < (oddmean+3.*oddstd)) & (odddata > (oddmean-3.*oddstd))
            alloddmean[amp] = np.nanmean(odddata[good])
        #print allevenmean
        #print alloddmean
        return allevenmean,alloddmean

    def submeans(self,data,allevenmean,alloddmean):
        template = np.zeros_like(data)
        for colnum in np.arange(2048):

            if colnum < qstart[1]:
                amp = 0
            elif colnum >= qstart[1] and colnum < qstart[2]:
                amp = 1
            elif colnum >= qstart[2] and colnum < qstart[3]:
                amp = 2
            elif colnum >= qstart[3]:
                amp = 3
            if colnum % 2 == 0:
                data[:,colnum] -= allevenmean[amp]
            if colnum % 2 == 1:
                data[:,colnum] -= alloddmean[amp]
        return data


    def run(self):
        #case where outfile was not set by user
        if self.outfile == 'unset':
            inbase = self.infile[0:-5]
            self.outfile = inbase + '_evenodd.fits'
        
        if os.path.isfile(self.outfile):
            print("WARNING! {} already exists. Quitting.".format(self.outfile))
            return 0

        h = fits.open(self.infile)
        data = h[1].data
        corrdata = np.zeros_like(data)

        #4-D data (i.e. SSB-format)
        if len(data.shape) == 4:
            for intnum in xrange(data.shape[0]):
                for groupnum in xrange(data.shape[1]):
                    evmn,odmn = oddeven.calc_averages(data[intnum,groupnum,:,:])
                    corrgroup = oddeven.submeans(data[intnum,groupnum,:,:],evmn,odmn)
                    corrdata[intnum,groupnum,:,:] = corrgroup
        #3-D data (i.e. UofA data, raw CV2 data)
        if len(data.shape) == 3:
            for groupnum in xrange(data.shape[1]):
                evmn,odmn = oddeven.calc_averages(data[groupnum,:,:])
                corrgroup = oddeven.submeans(data[groupnum,:,:],evmn,odmn)
                corrdata[groupnum,:,:] = corrgroup

        #2-D data
        if len(data.shape) == 2:
            evmn,odmn = oddeven.calc_averages(data)
            corrgroup = oddeven.submeans(data,evmn,odmn)
            corrdata = corrgroup
            

        h[1].data = corrdata
        h.writeto(self.outfile)


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument('file',help='file to be corrected')
        parser.add_argument('-o','--outfile',help='file to write result to',default='unset')
        return parser

if __name__ == '__main__':
    usagestr = 'python even_odd_correct.py filename'
        
    oddeven = Oddeven()
    parser = oddeven.add_options()
    args = parser.parse_args()
    
    oddeven.infile = args.file
    oddeven.outfile = args.outfile
    corrdata = oddeven.run()
    



