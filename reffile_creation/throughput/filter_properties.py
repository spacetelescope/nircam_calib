#!/usr/bin/env python

'''
Calculate things like bandwidth, pivot wavelength, half-power wavelengths, etc
For a given set of filter bandpasses.
BNH - Dec 2015
'''

import numpy as np
import matplotlib.pyplot as plt
import argparse,sys
from astropy.io import ascii
from astropy.table import Table
from scipy.integrate import simps
from scipy.integrate import romb


#filter, location dictionary
location_dict = {'F070W':'Filter,1','F090W':'Filter,2','F115W':'Filter,3','F140M':'Filter,11','F150W':'Filter,4','F150W2':'Filter,12','F162M':'Pupil,6','F164N':'Pupil,11','F182M':'Filter,10','F187N':'Filter,8','F200W':'Filter,5','F210M':'Filter,9','F212N':'Filter,6','F250M':'Filter,10','F277W':'Filter,1','F300M':'Filter,4','F322W2':'Filter,12','F323N':'Pupil,6','F335M':'Filter,11','F356W':'Filter,2','F360M':'Filter,7','F405N':'Pupil,5','F410M':'Filter,6','F430M':'Filter,8','F444W':'Filter,3','F460M':'Filter,9','F466N':'Pupil,9','F470N':'Pupil,11','F480M':'Filter,5'}


class FiltProp:
    def __init__(self):
        self.verbose = False

    def get_filtfiles(self):
        filterfiles = []
        with open(self.listfile) as f:
            for line in f:
                if len(line) > 2:
                    filterfiles.append(line.strip())
        return filterfiles

    def get_ote(self):
        #read in the OTE throughput values
        ote = Table.read(self.otefile)
        ote_wave = ote['Wavelength'].data
        ote_thru = ote['Efficiency'].data
        return ote_wave,ote_thru

    def get_filter_data(self,file):
        #read in wavelength and throughput from a filter file
        filt = ascii.read(file,header_start=0,data_start=1)
        filt_wave = filt.columns[0].data
        filt_thru = filt.columns[1].data

        #some filters may be stored in wavelength-descending order. Re-order here.
        srt = np.argsort(filt_wave)
        filt_wave = filt_wave[srt]
        filt_thru = filt_thru[srt]
        return filt_wave,filt_thru

    def get_halfpower(self,wave,thru):
        #find the half-power wavelengths
        
        #scale the throughput to have a max of 1.0
        scalefactor = np.nanmax(thru)
        thru = thru / scalefactor

        #bandpass edges
        half_areas = np.where((thru > 0.3) & (thru < 0.7))[0]
        half_areas_shift = np.roll(half_areas,-1)
        diff = half_areas_shift - half_areas
        startofsegment = half_areas[0]
        jump = np.where(diff > 1)[0][0]
        endofsegment = half_areas[jump]
        segment = thru[startofsegment:endofsegment+1]
        segwave = wave[startofsegment:endofsegment+1]
        sw_half = np.interp(0.5,segment,segwave)

        r_half_areas = half_areas[::-1]
        r_half_areas_shift = np.roll(r_half_areas,-1)
        diff2 = r_half_areas - r_half_areas_shift
        endofsegment = r_half_areas[0]
        jump = np.where(diff2 > 1)[0][0]
        startofsegment = r_half_areas[jump]
        segment = thru[startofsegment:endofsegment+1]
        segwave = wave[startofsegment:endofsegment+1]
        #print("again, linear interpolation ok?")
        lw_half = np.interp(0.5,segment[::-1],segwave[::-1])
        #print("Found half-power wavelengths of {} and {}".format(sw_half,lw_half))
        return sw_half,lw_half


    def get_waveedges(self,wave,thru,thresh):
        #find the short and long wavelength edges of the passband, based on a throughput threshold
        #These values are NOT used to find the bandwidth.
        good = np.where(thru >= thresh)[0]
        sw = wave[good[0]]
        lw = wave[good[-1]]
        return sw,lw


    def get_bandwidth(self,wave,thru,thresh):
        #Calculate the bandwidth. Bandwidth is the integral of the normalized 
        #response. Defintion is in Rieke et al 2008, Appendix E eq 1
        good = np.where(thru >= thresh)[0]
        sw = wave[good[0]]
        lw = wave[good[-1]]
        bandwidth=simps(thru[good],x=wave[good],even='avg')/np.max(thru[good])
        return sw,lw,bandwidth

    def pivot_wavelength(self,wave,thru,thresh):
        #Calculate the pivot wavelength
        good = np.where(thru >= thresh)[0]        
        pivot = np.sqrt(simps(thru[good]*wave[good],x=wave[good],even='avg')/simps(thru[good]/wave[good],x=wave[good],even='avg'))
        return pivot

    def effective_response(self,wave,thru,pivot,bandwidth):
        #Calculate the effective response. This is defined as the mean response over the wavelength range of 
        #pivot wavelength +/- 0.5*bw
        good = (wave < (pivot+0.5*bandwidth)) & (wave > (pivot-0.5*bandwidth))
        eff = np.mean(thru[good])
        return eff

    def run(self):
        #read in the list of filter files
        filtfiles = self.get_filtfiles()
        
        #sort the files so that we have short to long wavelength files
        filtfiles.sort()

        #read in OTE file if requested
        if self.otefile == None:
            print("WARNING: no OTE throughput file specified. Continuing without including. (THIS COULD MEAN THAT THE INPUT")
            print("DATA ALREADY INCLUDE THE OTE EFFECTS, AND YOU JUST DON'T WANT TO MULTIPLY THE OTE IN A SECOND TIME.")
        else:
            ote_wave,ote_thru = self.get_ote()
        

        #create lists to save the eventual output values for all filters
        filter_name = []
        sw_halfpower = []
        lw_halfpower = []
        bandwidth = []
        pivot_wavelength = []
        bw_over_pivot = []
        eff_response = []
        location = []

        #loop over the filter files
        for file in filtfiles:

            #read in filter data
            filt_wave,filt_thru = self.get_filter_data(file)

            #interpolate the OTE data to match the filter data
            #and multiply that in, if requested
            if self.otefile != None:
                ote_interp = np.interp(filt_wave,ote_wave,ote_thru)
                filt_thru = filt_thru * ote_interp
            

            #find half-power response wavelengths
            sw_half, lw_half = self.get_halfpower(filt_wave,filt_thru)
            sw_halfpower.append('{:8.3f}'.format(sw_half))
            lw_halfpower.append('{:8.3f}'.format(lw_half))

            #get bandwidth, along with the short and long wavelength filter profile edges, as 
            #defined by a minimum throughput value
            sw_edge,lw_edge,bw = self.get_bandwidth(filt_wave,filt_thru,self.edge_thresh)
            bandwidth.append('{:8.3f}'.format(bw))

            #Calculate the pivot wavelength of the filter
            pivot = self.pivot_wavelength(filt_wave,filt_thru,self.edge_thresh)
            pivot_wavelength.append('{:8.3f}'.format(pivot))

            #Ratio of bandwidth/pivot wavelength
            bwp = bw / pivot * 100.
            bw_over_pivot.append('{:8.1f}'.format(bwp))

            #Effective response
            response = self.effective_response(filt_wave,filt_thru,pivot,bw)
            eff_response.append('{:8.3f}'.format(response))

            #Generate string listing the filter's location (wheel, slot)
            #At the moment we assume that the input file name begins with the filter name
            #followed by an underscore
            under = file.find('_')
            if under == -1:
                print("WARNING: Input file {}, does not have a name that begins with the filter name".format(file))
                print("followed by an underscore. Not sure where to grab the filter name from. Quitting.")
                sys.exit()
            fname = file[0:under]
            filter_name.append(fname)

            loc = location_dict[fname]
            location.append(loc)
            

        #combine output lists into a table
        outtab = Table()
        outtab['Filter'] = filter_name
        outtab['Pivot_Wavelength_(microns)'] = pivot_wavelength
        outtab['Bandwidth_(microns)'] = bandwidth
        outtab['BW/Pivot_%'] = bw_over_pivot
        outtab['Effective_Response'] = eff_response
        outtab['SW_Response_Half-Power_(microns)'] = sw_halfpower
        outtab['LW_Response_Half-Power_(microns)'] = lw_halfpower
        outtab['Wheel,_Slot'] = location

        #Save the results. Save as an ascii table with whitespace between entries, for easier readability,
        #and also save a csv version, so that it can be opened easily in excel
        whitesptab = self.outfile
        csvtab = self.outfile
        dot = self.outfile.rfind('.')
        suffix = self.outfile[dot+1:]
        if suffix == 'csv':
            whitesptab = self.outfile[0:dot+1] + 'txt'
        else:
            csvtab = self.outfile[0:dot+1] + 'csv'

        ascii.write(outtab,whitesptab)
        ascii.write(outtab,csvtab,format='csv')
        csvdot = csvtab.rfind('.')
        texname = csvtab[0:dot+1] + 'tex'
        ascii.write(outtab,texname,format='aastex')


        #html version for NIRCam webpage
        outtab.remove_column('Wheel,_Slot')
        htmlname = csvtab[0:dot+1] + 'html'
        ascii.write(outtab,htmlname,format='html')

    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description="Given filter transmission curves, calcualte bandwidth, pivot wavelength, effective response, half-power wavelengths, etc")
        parser.add_argument("listfile",help="Name of text file containing the list of filter bandpass files.")
        parser.add_argument("-p","--makeplots",help="If set, create and save plots of the filter bandpasses.",action='store_true')
        parser.add_argument("--outfile",help="Name of output file containing filter properties.",default=None)
        parser.add_argument("--otefile",help="File containing OTE transmission values.",default=None)
        parser.add_argument("--edge_thresh",help="Threshhold throughput value used to find the edges of the bandpass.",default=1.e-4)
        return parser


if __name__ == '__main__':
    usagestring = "USAGE: filter_properties.py filters.list"


    filtprop = FiltProp()
    parser = filtprop.add_options(usage=usagestring)
    args = parser.parse_args()

    filtprop.listfile = args.listfile
    filtprop.makeplots = args.makeplots
    filtprop.otefile = args.otefile
    filtprop.outfile  = args.outfile
    filtprop.edge_thresh = args.edge_thresh

    filtprop.run()
