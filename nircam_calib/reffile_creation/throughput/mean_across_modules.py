#!/usr/bin/env python

'''
Input two throughput tables. One for mod A and one for mod B.
Calculate the average and print out into a new table.
'''

import numpy as np
import argparse,sys
from astropy.io import ascii
from astropy.table import Table
 
class ModMean:
    def __init__(self):
        self.verbose = False

    def read_tab(self,file):
        #read in a throughput table. assume ascii for the moment
        t = ascii.read(file)
        wave = t.columns[0].data
        thru = t.columns[1].data
        name1 = t.colnames[0]
        name2 = t.colnames[1]
        return wave,thru,name1,name2

    def mean_tab(self,tab1w,tab1t,tab2w,tab2t):
        #calculate the mean of the two tables
        
        #interpolate tab2 onto the wavelength entries of tab1
        tab2ti = np.interp(tab1w,tab2w,tab2t)

        #now calculate the mean
        tab1t = np.expand_dims(tab1t,axis=0)
        tab2ti = np.expand_dims(tab2ti,axis=0)
        concat = np.concatenate((tab1t,tab2ti))
        meantab = np.mean(concat,axis=0)
        return meantab


    def run(self):
        #read in table 1
        tab1w,tab1t,wname1,tname1 = self.read_tab(self.file1)
        
        #read in table 2
        tab2w,tab2t,wname2,tname2 = self.read_tab(self.file2)

        #put into ascending wavelength order
        srt = np.argsort(tab1w)
        tab1w = tab1w[srt]
        tab1t = tab1t[srt]
        srt = np.argsort(tab2w)
        tab2w = tab2w[srt]
        tab2t = tab2t[srt]

        #write out re-sorted individual tables
        outindtab1 = Table([tab1w,tab1t],names=(wname1,tname1))
        dot = self.file1.rfind('.')
        ascii.write(outindtab1,self.file1[0:dot]+'_sorted.txt')
        outindtab2 = Table([tab2w,tab2t],names=(wname2,tname2))
        dot = self.file1.rfind('.')
        ascii.write(outindtab2,self.file2[0:dot]+'_sorted.txt')

        #calculate the mean
        mtabt = self.mean_tab(tab1w,tab1t,tab2w,tab2t)

        #create table
        outtab = Table([tab1w,mtabt],names=('microns','throughput'))
        
        #write out mean table
        if self.outfile == None:
            self.outfile = 'mean_throughput_'+self.file1+'_'+self.file2+'.txt'
        ascii.write(outtab,self.outfile)


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description="make mean throughput curve across mod A and mod B")
        parser.add_argument("file1",help="Name of first table")
        parser.add_argument("file2",help="Name of second table")
        parser.add_argument("--outfile",help="Name of output file containing mean table.",default=None)
        return parser


if __name__ == '__main__':
    usagestring = "USAGE: filter_properties.py filters.list"


    mmean = ModMean()
    parser = mmean.add_options(usage=usagestring)
    args = parser.parse_args()

    mmean.file1 = args.file1
    mmean.file2 = args.file2
    mmean.outfile = args.outfile

    mmean.run()
