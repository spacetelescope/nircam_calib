#! /usr/bin/env python 

"""
Script name          : plot_pixel_integration.py
Author               : Bryan Hilbert
Created              : 23rd Feb 2016
Version              : 1.0

Description          : This program reads in a file that contains a list of 
                     : filenames. For each file in the list, it plots the signal 
                     : levels in the first integration for pixel (x,y), where the 
                     : user chooses x,y. The integrationss for all files in the list 
                     : are all plotted on a single plot.
#---------------------------------------------------------------------------


for a given (x,y) and list of files, plot the ramps for that pixel
from all files
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
import numpy as np


class Ramp_Plot:
    def __init__(self):
        #self.verbose == False
        #pass
        test = 0

    def run(self):
        
        #read in list of files to use
        files = []
        with open(self.infile) as f:
            for line in f:
                if len(line) > 2:
                    files.append(line.strip())


        #set up plot
        f,a = plt.subplots()

        #loop over files, read in, plot
        for file in files:
            with fits.open(file) as h:
                data = h[1].data

            if len(data.shape) == 4:
                data = data[0,:,:,:]
            else:
                print('File {} does not contain a ramp. Omitting.'.format(file))
                zd = -1
            
            zd = data.shape[0]

            if zd != -1:
                time = np.arange(zd) * 10.7
                a.plot(time,data[:,self.y,self.x],'bo')

        a.set_ylabel('Signal')
        a.set_xlabel('Time (sec)')
        a.set_title('Pixel ({},{})'.format(str(self.x),str(self.y)))
        if self.outfile == None:
            self.outfile = 'Ramp_plots_from_'+self.infile+'_pixel_'+str(self.x)+'_'+str(self.y)+'.pdf'

        f.savefig(self.outfile)
        plt.close(f)


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='Plot ramp for pixel (x,y) in all files in list.')
        parser.add_argument('infile',help='File listing ramps to use.')
        parser.add_argument('x',help='x-coordinate of pixel to plot')
        parser.add_argument('y',help='y-coordinate of pixel to plot')
        parser.add_argument('-o','--outfile',help='name of output file for plots.',default=None)
        return parser


if __name__ == '__main__':
    usagestring = 'python plot_ramps.py ramps.list x y'

    rplot = Ramp_Plot()
    parser = rplot.add_options(usage=usagestring)
    args = parser.parse_args()

    rplot.infile = args.infile
    rplot.x = args.x
    rplot.y = args.y
    rplot.outfile = args.outfile

    rplot.run()
