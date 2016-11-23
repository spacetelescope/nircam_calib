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


class Throughput:
    def __init__(self):
        self.verbose = False

    def get_files(self):
        #get the list of filter filenames
        files = []
        with open(self.filterlist) as f:
            for line in f:
                if len(line) > 3:
                    files.append(line.strip())
        return files

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

    #Converts data to floats
    def tofloat(self,data):
        a = [x for x in data if x != '--']
        b = [float(i) for i in a]
        return b

    #replaces '---'' with nan
    def replace_nan(self,items):
        for index, item in enumerate(items):
            if (item == '---'):
                items[index] = float('nan')
        return items


    def read_optics(self):
        #read in the file containing the optics throughputs
        opttab = ascii.read(self.opticsfile,header_start=1,data_start=2,format='csv')

        #get info out of the table
        wave = opttab['Wavelength'].data.data
        nvr_thru = opttab['NVR_Transmission'].data.data
        nvr_wave = opttab['NVR_Wavelength'].data.data
        collimator = opttab['Collimator'].data.data
        sw_triplet = self.replace_nan(opttab['SW_Triplet'].data.data).astype('float')
        sw_mirrors = self.replace_nan(opttab['SW_Mirrors'].data.data).astype('float')
        lw_triplet = self.replace_nan(opttab['LW_triplet'].data.data).astype('float')
        lw_mirrors = self.replace_nan(opttab['LW_Mirrors'].data.data).astype('float')
        sw_particulates = self.replace_nan(opttab['SW_Particulates'].data.data).astype('float')
        lw_particulates = self.replace_nan(opttab['LW_Particulates'].data.data).astype('float')

        #manually extend the SW particulate column so that it is not the reason for F150W2 
        #to abruptly cut off at 2.44um
        sw_particulates = np.zeros(len(sw_triplet)) + sw_particulates[0]

        #do the same for LW particulates
        lw_particulates = np.zeros(len(lw_triplet)) + 0.984

        #remove extra entries in NVR columns
        good = np.where(nvr_wave != 0.)[0]
        nvr_thru = nvr_thru[good]
        nvr_wave = nvr_wave[good]

        #interpolate NVR to the same wavelength scale as the other columns
        nvr_interp = np.interp(wave,nvr_wave,nvr_thru)

        #combine the elements to produce a SW optics curve and a LW optics curve
        #The 0.98 factor is a 'reserve' that lockheed martin included in the throughput
        #budget, as the reported throughput values are supposed to represent those at the
        #END of the mission lifetime. --> from Marcia, via Stansberry, in 10 Feb 2016 email.
        #So, include the 0.98 factor here.
        sw_optics = nvr_interp * collimator * sw_triplet * sw_mirrors * sw_particulates * 0.98
        lw_optics = nvr_interp * collimator * lw_triplet * lw_mirrors * lw_particulates * 0.98

        fo,ao = plt.subplots()
        ao.plot(wave,collimator,color='red',label='collimator')
        ao.plot(wave,sw_triplet,color='blue',label='sw_trip')
        ao.plot(wave,sw_mirrors,color='black',label='sw_mirror')
        ao.plot(wave,lw_triplet,color='blue',linestyle='--',linewidth=2,label='lw_trip')
        ao.plot(wave,lw_mirrors,color='black',linestyle='--',linewidth=2,label='lw_mirror')
        ao.plot(wave,sw_particulates,color='green',label='sw_part')
        ao.plot(wave,lw_particulates,color='green',linestyle='--',linewidth=2,label='lw_part')
        ao.plot(nvr_wave,nvr_thru,color='orange',linewidth=2,label='NVR')
        ao.plot(wave,nvr_interp,color='magenta',linestyle='--',label='NVR_interp')
        ao.plot(wave,sw_optics,color='pink',label='Total SW')
        ao.plot(wave,lw_optics,color='pink',linestyle='--',linewidth=2,label='Total LW')
        ao.legend(loc='best')
        fo.savefig('optics_parts_throughputs.png')
        plt.close(fo)
        print("Throughputs of optics components saved to optics_parts_throughputs.png")

        return wave,sw_optics,lw_optics


    def get_DBS(self):
        #dbsdir = '/grp/jwst/wit/nircam/reference_files/SpectralResponse_2015-Final/DBS_QE_Optics/'
        dbsdir = ''
        #dbsfiles = ['DBS_SW_ModA.txt','DBS_SW_ModB.txt','DBS_LW_ModA.txt','DBS_LW_ModB.txt']
        #These files are the reflection and transmission columns from Marcia's spreadsheets. For SWA I have interpolated
        #from the low-res version over the anomalous 'knee' in the response curve from 2.3um onward.
        dbsfiles = ['DBS_SW_modA_highres_1_minus_trans_plus_absorb.txt','DBS_SW_ModB_highres.txt','DBS_LW_ModA_highres.txt','DBS_LW_ModB_highres.txt']
        swa = ascii.read(dbsdir+dbsfiles[0],header_start=0,data_start=1)
        swa_w = swa.columns[0].data
        swa_t = swa.columns[1].data

        swb = ascii.read(dbsdir+dbsfiles[1],header_start=0,data_start=1)
        swb_w = swb.columns[0].data
        swb_t = swb.columns[1].data

        lwa = ascii.read(dbsdir+dbsfiles[2],header_start=0,data_start=1)
        lwa_w = lwa.columns[0].data
        lwa_t = lwa.columns[1].data

        lwb = ascii.read(dbsdir+dbsfiles[3],header_start=0,data_start=1)
        lwb_w = lwb.columns[0].data
        lwb_t = lwb.columns[1].data

        return swa_w,swa_t,swb_w,swb_t,lwa_w,lwa_t,lwb_w,lwb_t


    def get_qe_curve(self,channel,module,wavelength):
        sw_coeffs = np.array([0.65830,-0.05668,0.25580,-0.08350])
        #from Stansberry:
        sw_exponential = 100.
        sw_wavecut = 2.38
        lw_coeffs_a = np.array([0.934871,0.051541,-0.281664,0.243867,-0.086009,0.014509,-0.001])
        lw_factor_a = 0.88
        lw_coeffs_b = np.array([2.9104951,-2.182822,0.7075635,-0.071767])

        if channel == 'SW':
            qe = sw_coeffs[0] + sw_coeffs[1]*wavelength + sw_coeffs[2]*wavelength**2 + sw_coeffs[3]*wavelength**3
            #insert the exponential decay for wavelengths redward of sw_wavecut
            red = wavelength > sw_wavecut
            qe[red] = qe[red] * np.exp((sw_wavecut-wavelength[red])*sw_exponential)
        if channel == 'LW':
            if module == 'a':
                qe = lw_factor_a * (lw_coeffs_a[0] + lw_coeffs_a[1]*wavelength + lw_coeffs_a[2]*wavelength**2 + lw_coeffs_a[3]*wavelength**3 + lw_coeffs_a[4]*wavelength**4 + lw_coeffs_a[5]*wavelength**5 + lw_coeffs_a[6]*wavelength**6)
            if module == 'b':
                qe = lw_coeffs_b[0] + lw_coeffs_b[1]*wavelength + lw_coeffs_b[2]*wavelength**2 + lw_coeffs_b[3]*wavelength**3

        return qe

    def get_ote(self):
        #read in the OTE throughput values
        ote = Table.read(self.otefile)
        ote_wave = ote['wavelength'].data
        ote_thru = ote['throughput'].data
        return ote_wave,ote_thru

    def individual_plot(self,wave,thru,filtname,ote=False):
        #Individual plot of filter throughput curve
        ylab = 'NIRCam Instrument-Only Throughput'
        fout = filtname+'_nrc_only_throughput_mod'+self.module+'.pdf'
        if self.otefile != None:
            ylab = 'NIRCam + OTE System Throughput'
            fout = filtname+'_nrc_and_ote_throughput_mod'+self.module+'.png'            
        ftemp,axtemp = plt.subplots()
        axtemp.plot(wave,thru,'bo')
        axtemp.set_xlabel('Wavelength (microns)')
        axtemp.set_ylabel(ylab)
        ftemp.savefig(fout,clobber=True)
        plt.close(ftemp)

    def load_blocking_filter(self,filtername,blockingfilter,filelist):
        #identify and load the throughput for the appropriate blocking filter
        blockingfile = [s for s in filelist if blockingfilter in s][0]

        data = ascii.read(blockingfile)
        under = blockingfile.find('_')
        bfname = blockingfile[0:under]
        print("Loading blocking filter {} for PW-based filter {}.".format(bfname,filtername))
        waves = data.columns[0].data
        through = data.columns[1].data

        #put wavelength in increasing order
        srt = np.argsort(waves)
        waves = waves[srt]
        through = through[srt]
        return waves,through


    def run(self):
        #read in optics file
        #multiply columns, interpolate nvr and multiply in to get total optics
        #where is 'scatter' from johns file? keep his fudge factor of 0.98??
        wave_optics, sw_optics, lw_optics = self.read_optics()

        #set up a function that will interpolate the optics transmission for each
        #filter later
        min_opt_wave = np.min(wave_optics)
        max_opt_wave = np.max(wave_optics)

        #get OTE if asked
        if self.otefile != None:
            ote_wave,ote_thru = self.get_ote()
            
        #get lists of the filter filenames
        #w2_list,w_list,m_list,n_list = self.get_filterfilenames()
        #allwide = np.concatenate((w2_list,w_list)) 

        #get list of filter filenames
        all_filts = self.get_files()

        #vertical offsets for labels, depending on the type of throughput
        #delta = 0.15
        #if self.otefile != None:
        #    delta = 0.20

        #read in 'scatter' throughput contribution, due to particulates.
        #All we have are filter-average values from Marica, via John Stansberry.
        #Read these values and multiply them in to the filter profiles before plotting
        scatter = ascii.read('nircam_optics_filter_average.dat',header_start=7,data_start=8)
        scatter_filt = scatter['Filter'].data
        scatter_val = scatter['Scatter'].data.astype('float')

        #get DBS curves
        dbs_swa_w,dbs_swa_t,dbs_swb_w,dbs_swb_t,dbs_lwa_w,dbs_lwa_t,dbs_lwb_w,dbs_lwb_t = self.get_DBS()


        #full_list = 
        #loop through files and read in. Add to plot.
        for file in all_filts:
            data = ascii.read(file)
            under = file.find('_')
            fname = file[0:under]
            print("Working on file for {}".format(fname))
            waves = data.columns[0].data
            through = data.columns[1].data

            #put wavelength in increasing order
            srt = np.argsort(waves)
            waves = waves[srt]
            through = through[srt]

            #restrict to the wavelength range covered by the optics
            goodwave = (waves >= min_opt_wave) & (waves <= max_opt_wave)
            through = through[goodwave]
            waves = waves[goodwave]

            #set up the plot that will show the throughput after each step
            fstep,astep = plt.subplots()
            astep.plot(waves,through,color='black',label='Filter Only')
            #also keep track of the throughput value at roughly the central wavelength after each step
            highthru = through > 0.5
            centralwave = np.median(waves[highthru])
            check = np.where(np.absolute(waves-centralwave) < 0.002)[0][0]
            #print(centralwave,check)
            #check = np.where(np.absolute((waves-3.5)) < 0.05)[0][0]
            checktable = Table(names=('Contributors','Element_Contribution','Throughput_at_'+str('%.4f' % centralwave)+'_microns'),dtype=['S55','f','f'])
            checktable.add_row(('Filter','%.4f' % through[check],'%.4f' % through[check]))


            #if the filter is in the pupil wheel, load the appropriate FW-based blocking
            #filter, and multiply its throughput in
            block_legend = ''
            if fname in pwheel:
                blockingname = pwheel[fname]
                blockwave,blockthru = self.load_blocking_filter(fname,blockingname,all_filts)

                #interpolate to match the wavelength scale of the primary filter
                blockingthru = np.interp(waves,blockwave,blockthru)

                #multiply the blocking filter's throughput in to the primary filter's throughput
                through = through * blockingthru
                
                #add to step-by-step plot
                block_legend = '*BlockingFilt'
                astep.plot(waves,through,color='silver',label='Filt'+block_legend)
                checktable.add_row(('Filter'+block_legend,'%.4f' % blockingthru[check],'%.4f' % through[check]))


            #multiply in scatter
            #fmatchname = file[0:6]
            smatch = np.where(scatter_filt == fname)[0][0]
            through = through * scatter_val[smatch]

            #add to step-by-step plot
            astep.plot(waves,through,color='red',label='Filt'+block_legend+'*Scatter')
            checktable.add_row(('Filter'+block_legend+'*Scatter','%.4f' % scatter_val[smatch],'%.4f' % through[check]))

            wname = float(fname[1:4])
            if wname < 249.:
                channel = 'SW'
                opt_match = np.interp(waves,wave_optics, sw_optics,left=np.nan,right=np.nan)
            else:
                channel = 'LW'
                opt_match = np.interp(waves,wave_optics, lw_optics,left=np.nan,right=np.nan)

            #fopt,aopt = plt.subplots()
            #aopt.plot(waves,opt_match,color='black')
            #fopt.savefig('optics_tmp.pdf')

            total_thru = through * opt_match

            #add to step-by-step plot
            astep.plot(waves,total_thru,color='orange',label='Filter'+block_legend+'*Scatter*Optics')
            checktable.add_row(('Filter'+block_legend+'*Scatter*Optics','%.4f' % opt_match[check],'%.4f' % total_thru[check]))

            #multiply the DBS in
            if self.module == 'a':
                if channel == 'SW':
                    dbs_val = np.interp(waves,dbs_swa_w,dbs_swa_t)
                if channel == 'LW':
                    dbs_val = np.interp(waves,dbs_lwa_w,dbs_lwa_t)
            else:
                if channel == 'SW':
                    dbs_val = np.interp(waves,dbs_swb_w,dbs_swb_t)
                if channel == 'LW':
                    dbs_val = np.interp(waves,dbs_lwb_w,dbs_lwb_t)

            total_thru = total_thru * dbs_val

            #add to step-by-step plot
            astep.plot(waves,total_thru,color='green',label='Filter'+block_legend+'*Scatter*Optics*DBS')
            checktable.add_row(('Filter'+block_legend+'*Scatter*Optics*DBS','%.4f' % dbs_val[check],'%.4f' % total_thru[check]))

            #multiply the QE in
            qe = self.get_qe_curve(channel,self.module,waves)
            fqe,aqe = plt.subplots()
            aqe.plot(waves,qe)
            aqe.set_title('QE, Module '+self.module.upper())
            aqe.set_xlabel('Wavelength (microns)')
            aqe.set_ylabel('QE')
            aqe.set_ylim(0.4,1.0)
            fqe.savefig('QE_mod'+self.module.upper()+'.pdf')
            plt.close(fqe)
            total_thru = total_thru * qe

            #check for negative throughput values and zero these out if found
            total_thru[total_thru < 0] = 0.

            #add to step-by-step plot
            astep.plot(waves,total_thru,color='blue',label='Filter'+block_legend+'*Scatter*Optics*DBS*QE')
            checktable.add_row(('Filter'+block_legend+'*Scatter*Optics*DBS*QE','%.4f' % qe[check],'%.4f' % total_thru[check]))

            #Save nircam-only table file  
            good = ~np.isnan(total_thru)
            outtab = Table([waves[good],total_thru[good]], names=['Wavelength_microns', 'Throughput'])
            tabfile = fname+'_nrc_only_throughput_mod'+self.module+'.txt'
            print("Tabulated output for NIRCam-only throughtputs saved in {}".format(tabfile))
            ascii.write(outtab, tabfile)

            #Individual plot
            self.individual_plot(waves,total_thru,fname,ote=False)

            #if requested, multiply the OTE in. Save a non-OTE version as well.
            if self.otefile != None:
                ote_thru_i = np.interp(waves,ote_wave,ote_thru)
                total_thru = total_thru * ote_thru_i

                #add to step-by-step plot, and set xrange of plot
                astep.plot(waves,total_thru,color='magenta',label='Filter'+block_legend+'*Scatter*Optics*DBS*QE*OTE')
                checktable.add_row(('Filter'+block_legend+'*Scatter*Optics*DBS*QE*OTE','%.4f' % ote_thru_i[check],'%.4f' % total_thru[check]))
                #save tabulated throughput to an ascii file
                good = ~np.isnan(total_thru)
                outtab = Table([waves[good], total_thru[good]], names=['Wavelength_microns', 'Throughput'])

                tabfile = fname+'_nircam_plus_ote_throughput_mod'+self.module+'.txt'
                print("Tabulated output for NIRCam+OTE throughtputs saved in {}".format(tabfile))
                ascii.write(outtab, tabfile)

                #individual plot of nircam+ote throughputs
                self.individual_plot(waves,total_thru,fname,ote=True)

            #save the step-by-step plot, as well as the step-by-step table
            astep.legend(loc='lower center',prop={'size':8})
            astep.set_ylabel('Throughput')
            astep.set_xlabel('Wavelength (microns)')
            astep.set_title(fname+' Module '+self.module.upper())
            astep.set_ylim(0,1)
            plotthru = np.where(total_thru > 0.005)[0]
            astep.set_xlim(waves[plotthru[0]],waves[plotthru[-1]])
            fstep.savefig('Step_by_step_plot_'+fname+'_mod'+self.module.upper()+'.pdf')
            plt.close(fstep)
            checktable['Contributors'].format = '<'
            checktable.write('Step_by_step_table_'+fname+'_mod'+self.module.upper()+'.txt',format='ascii.fixed_width',delimiter=None)


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description="Convert CSV files of filter transmissions from Marcia into ascii files")
        parser.add_argument("--opticsfile",help="Name of file containing optics throughputs")
        #parser.add_argument("dbsfile",help="Name of file containing dichroic beam splitter transmissions")
        parser.add_argument("--otefile",help="Name of file containing OTE throughputs",default=None)
        parser.add_argument("filterlist",help="File containing list of filter filenames, ascii files")
        parser.add_argument("-x","--opticsexcel",help="Set if the optics throughput file is an excel workbook.",action='store_true')
        parser.add_argument("-p","--makeplots",help="Create plots of filter transmissions",action='store_true')
        parser.add_argument("-o","--outfile",help="Base name of output file.",default=None)
        parser.add_argument("-m","--module",help="Module to use. Choices are 'a' or 'b'.",default='a')
        #parser.add_argument("-c","--csv",help="Set if the input file is a CSV file rather than an excel notebook",action='store_true')
        return parser


if __name__ == '__main__':
    usagestring = "USAGE: nircam_througputs.py opticsfile.csv dbsfile.csv filterfile_list.txt"

    throughput = Throughput()
    parser = throughput.add_options(usage=usagestring)
    args = parser.parse_args()

    #set up variables
    if args.otefile == None:
        args.otefile = '/itar/jwst/tel/share/Mirror_Reflectivity/jwst_telescope_ote_thruput.fits'

    throughput.opticsfile = args.opticsfile
    #throughput.dbsfile = args.dbsfile
    throughput.filterlist = args.filterlist
    throughput.otefile = args.otefile
    throughput.module = args.module.lower()
    if (throughput.module != 'a') & (throughput.module != 'b'):
        print("Warning!! Incorrect module chosen. Pick 'a' or 'b'.")
        sys.exit()

    throughput.outfile = args.outfile
    #if args.outfile == None:
    #    throughput.outfile = args.infile
    
    #extract info
    throughput.run()
