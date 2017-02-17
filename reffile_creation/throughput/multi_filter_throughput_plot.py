#!/usr/bin/env python

import argparse,sys,os
from astropy.io import ascii
from astropy.table import Table
import numpy as np
#from openpyxl import load_workbook
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

pwheel = {'F162M':'F150W2','F164N':'F150W2','F323N':'F322W2','F405N':'F444W','F466N':'F444W','F470N':'F444W'}

class Multiplot:
    def __init__(self):
        self.verbose = False

    def get_filterfilenames(self):
        #get the list of filter filenames and sort into lists by filter width
        w_list = []
        w2_list = []
        n_list = []
        m_list = []
        with open(self.filterlist) as f:
            for line in f:
                if len(line) > 3:
                    bp = line[4:6]
                    wav = line[1:4]
                    if bp[0] == 'N':
                        n_list.append(line.strip())
                    if bp[0] == 'M':
                        m_list.append(line.strip())
                    if bp == 'W_':
                        if (wav != '322'):
                            w_list.append(line.strip())
                        else:
                            w2_list.append(line.strip())
                    if bp == 'W2':
                        w2_list.append(line.strip())
        return w2_list,w_list,m_list,n_list

    def individual_plot(self,w,t,filter):
        #create plot of individal filter
        good = t > 0.01
    
        ftemp,axtemp = plt.subplots()
        axtemp.plot(w[good],t[good],'bo')
        axtemp.plot(w[good],t[good],color='blue')
        axtemp.set_xlabel('Wavelength (microns)')
        axtemp.set_ylabel('NIRCam Instrument-Only Throughput')
        indplotfile = filter+'_nrc_only_throughput_mod'+self.module+'.pdf'
        if self.ote != False:
            axtemp.set_ylabel('NIRCam + OTE System Throughput')
            indplotfile = filter+'_nrc_and_ote_throughput_mod'+self.module+'.pdf'
        if self.filteronly == True:
            axtemp.set_ylabel('Throughput: Filter Only')
            indplotfile = filter+'_filteronly_throughput_mod'+self.module+'.pdf'
            
        ftemp.savefig(indplotfile,clobber=True)
        plt.close(ftemp)


    def run(self):
        
        #get lists of filenames
        w2_list,w_list,m_list,n_list = self.get_filterfilenames()
        total_list = w2_list + w_list + m_list + n_list

        #create figure, assign size
        f,(ax1,ax2,ax3,ax4) = plt.subplots(4,sharex=True,sharey=True)
        f.set_size_inches(11.5,9)

        #set up colors
        spectral = cm = plt.get_cmap('spectral') 
        cNorm  = colors.Normalize(vmin=0.5, vmax=5.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=spectral)

        #set up minor tick marks
        XminorLocator = MultipleLocator(0.1)
        YminorLocator = MultipleLocator(0.1)

        #set up the delta for the labels
        delta = 0.
        if self.filteronly == True:
            delta = -0.35
        if self.nrc_optics == True:
            delta = 0.18
        if self.ote == True:
            delta = 0.20


        #loop through files and read in. Add to plot.
        for file in total_list:
            data = ascii.read(file)
            under = file.find('_')
            fname = file[0:under]
            print('Working on filter {}'.format(fname))
            waves_in = data.columns[0].data
            thru_in = data.columns[1].data

            #remove nans
            good = ~np.isnan(thru_in)
            waves = waves_in[good]
            thru = thru_in[good]

            #individual plot - if requested
            if self.ind_plots == True:
                    self.individual_plot(waves,thru,fname)

            #find the midpoint of the area with the transmission > 50%
            #put the filter label in that location
            #midpt_thresh = 0.5
            #if file in n_list:
            #    midpt_thresh = 0.4
            #if fname == 'F070W':
            #    midpt_thresh = 0.2
            #    print(thru)
                

            midpt_thresh = np.max(thru) - 0.15

            upthru = thru > midpt_thresh
            midpt = np.median(waves[upthru])
            colorVal = scalarMap.to_rgba(midpt)

            if file in w2_list:
                ax1.plot(waves,thru,color=colorVal)
                ax1.fill_between(waves,0,thru,color=colorVal,alpha=0.5)
                ax1.text(midpt-0.1,0.85-delta,fname,color='black')

            if file in w_list:
                ax2.plot(waves,thru,color=colorVal)
                ax2.fill_between(waves,0,thru,color=colorVal,alpha=0.5)
            
                #adjustment of label position for F090W, to prevent overlap
                if fname != 'F090W':
                    ax2.text(midpt-0.1,0.85-delta,fname,color='black')
                else:
                    ax2.text(midpt-0.1,0.7-delta,fname,color='black')
                
            if file in m_list:
                #ndelta = [-0.2,-0.14,-0.08,-0.04,-0.07,-0.1,-0.13,-0.1,-0.18,-0.08,-0.1,-0.05]
                ndelta = {'F140M':-0.2,'F162M':-0.14,'F182M':-0.08,'F210M':-0.04,'F250M':-0.07,'F300M':-0.1,'F335M':-0.13,'F360M':-0.1,'F410M':-0.18,'F430M':-0.08,'F460M':-0.1,'F480M':-0.05}

                ax3.plot(waves,thru,color=colorVal)
                ax3.fill_between(waves,0,thru,color=colorVal,alpha=0.5)
                if (fname != 'F480M') and (fname != 'F182M'):
                    ax3.text(midpt+ndelta[fname],0.85-delta,fname,color='black',fontsize=10)
                if fname == 'F480M':
                    ax3.text(midpt+ndelta[fname],0.73-delta,fname,color='black',fontsize=10)
                if fname == 'F182M':
                    ax3.text(midpt+ndelta[fname],0.9-delta,fname,color='black',fontsize=10)


            if file in n_list:
                ndelta = {'F164N':-0.15,'F187N':-0.1,'F212N':-0.05,'F323N':-0.1,'F405N':-0.1,'F466N':-0.15,'F470N':-0.05}
                ax4.plot(waves,thru,color=colorVal)
                ax4.fill_between(waves,0,thru,color=colorVal,alpha=0.5)
                if fname != 'F470N':
                    ax4.text(midpt+ndelta[fname],0.85-delta,fname,color='black',fontsize=10)
                else:
                    ax4.text(midpt+ndelta[fname],0.73-delta,fname,color='black',fontsize=10)

        #limit the number of y-ticks
        if self.ote == True:
            ax1.set_ylim(0,0.8)
            plotlim = [0,0.8]
        if (self.nrc_optics == True and self.ote == False):
            ax1.set_ylim(0,0.8)
            print("HERE!!!")
            plotlim = [0,0.8]
        if self.filteronly == True:
            ax1.set_ylim(0,1.4)
            plotlim = [0,1.4]
        #plotlim = ax1.get_ylim()
        #ax1.set_yticks(np.arange(0.,plotlim[1]+0.2,0.2))
        ax1.set_yticks(np.arange(0.,plotlim[1]+0.2,0.2))
        ax1.yaxis.set_minor_locator(YminorLocator)

        #create the band that shows the beamsplitter deadband
        ax1.fill_between([2.35,2.4],0,1-delta,facecolor='black',alpha=0.5)
        ax2.fill_between([2.35,2.4],0,1-delta,facecolor='black',alpha=0.5)
        ax3.fill_between([2.35,2.4],0,1-delta,facecolor='black',alpha=0.5)
        ax4.fill_between([2.35,2.4],0,1-delta,facecolor='black',alpha=0.5)
        
        #indicate filters mounted in the pupil wheel
        ax3.text(1.6,0.5-delta,'P',color='black')
        ax4.text(1.53,0.6-delta,'P',color='black')
        #ax4.text(2.18,0.6,'P',color='black')
        ax4.text(3.1,0.6-delta,'P',color='black')
        ax4.text(3.9,0.6-delta,'P',color='black')
        ax4.text(4.54,0.6-delta,'P',color='black')
        ax4.text(4.77,0.6-delta,'P',color='black')

        #x range and ticks
        ax1.set_xlim(0.5,5.25)
        ax1.set_xticks(np.arange(0.5,5.5,0.5))
        #ax1.set_xticks(np.arange(0.5,5.5,0.1))
        ax1.xaxis.set_minor_locator(XminorLocator)

        
        if self.filteronly == True:
            f.text(0.08, 0.5, 'Filter Only Throughput', ha='center', va='center', rotation='vertical',fontsize=18)
        if (self.nrc_optics == True) & (self.ote == False):
            f.text(0.08, 0.5, 'NIRCam Throughput', ha='center', va='center', rotation='vertical',fontsize=18)
        if self.ote == True:
            f.text(0.08, 0.5, 'NIRCam + OTE Throughput', ha='center', va='center', rotation='vertical',fontsize=18)
        
            
        f.text(0.5,0.045,'Wavelength (microns)',ha='center', va='center',fontsize=18)
        f.subplots_adjust(hspace=0.15)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

        if self.outfile == None:
            if self.filteronly == True:
                self.outfile = "filter_only_throughput_plot_for_"+self.filterlist+".pdf"
            if self.nrc_optics == True:
                self.outfile = "filter_and_nircam_optics_throughput_plot_for_"+self.filterlist+".pdf"
            if self.ote == True:
                self.outfile = "nircam_plus_ote_system_throughput_plot_for_"+self.filterlist+".pdf"
        
        f.savefig(self.outfile,clobber=True)
        print("Output plot saved to {}".format(self.outfile))


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description="Convert CSV files of filter transmissions from Marcia into ascii files")
        #parser.add_argument("opticsfile",help="Name of file containing optics throughputs")
        #parser.add_argument("dbsfile",help="Name of file containing dichroic beam splitter transmissions")
        parser.add_argument("--nrc_optics",help="NIRCam optics present in inputs?",action='store_true')
        parser.add_argument("--ote",help="OTE Present in inputs?",action='store_true')
        parser.add_argument("filterlist",help="File containing list of filter filenames, ascii files")
        parser.add_argument("-o","--outfile",help="Base name of output file.",default=None)
        parser.add_argument("-m","--module",help="Module of data. Choices are 'a' or 'b'.",default=None)
        parser.add_argument("-f","--filteronly",help="Are throughput plots for filter only?",action='store_true')
        parser.add_argument("--ind_plots",help="If set, create individual plots of filters in addition to mult-filter plot.",action='store_true')
        return parser


if __name__ == '__main__':
    usagestring = "USAGE: multi_filter_throughput_plot.py filterlist.list"

    multiplot = Multiplot()
    parser = multiplot.add_options(usage=usagestring)
    args = parser.parse_args()

    #set up variables
    multiplot.filterlist = args.filterlist
    multiplot.ote = args.ote
    multiplot.nrc_optics = args.nrc_optics
    if args.module != None:
        multiplot.module = args.module.lower()
        if (multiplot.module != 'a') & (multiplot.module != 'b') & (multiplot.module != 'mean'):
            print("Warning!! Incorrect module chosen. Pick 'a', 'b', or 'mean'.")
            sys.exit()
        if multiplot.module == 'mean':
            multiplot.module == 'meanModAB'
    else: 
        multiplot.module = None

    multiplot.outfile = args.outfile
    multiplot.ind_plots = args.ind_plots
    multiplot.filteronly = args.filteronly

    #extract info
    multiplot.run()
