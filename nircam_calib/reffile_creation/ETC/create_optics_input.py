#! /usr/bin/env python

'''
Read in the optics + contaminants files for NIRCam and roll everything
together to produce a file for the ETC. Note that the ETC does not allow
for separate SW/LW files, so we'll have to choose a cutoff wavelength and 
chop input curves there. Also, at the moment the ETC cannot distinguish 
between module A and module B. For the filter throughputs, we have been
giving them the module mean. So to be consistent, it seems that we
should also provide the module mean for the optics and contaminants.

Produce some ancillary plots showing the mod A/B differences in each 
curve. The ETC is required to have an accuracy of 15%, with a goal of 
10%. The differences in the mod A/B LW QE curves is pretty large. If we
provide the module mean, the error introduced should be less than 15%...
I think.
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii

def replace_nan(items):
    for index, item in enumerate(items):
        if (item == '---'):
            items[index] = float('nan')
    return items


#DBS
lwafile = '../DBS_QE_Optics/DBS_LW_ModA_highres.txt'
lwbfile = '../DBS_QE_Optics/DBS_LW_ModB_highres.txt'
swagoodfile = '../DBS_QE_Optics/DBS_SW_modA_highres_1_minus_trans_plus_absorb.txt'
swbfile = '../DBS_QE_Optics/DBS_SW_ModB_highres.txt'

lwa_dbs = Table.read(lwafile,header_start=0,data_start=1,format='ascii')
lwb_dbs = Table.read(lwbfile,header_start=0,data_start=1,format='ascii')

swa_dbs = Table.read(swagoodfile,header_start=0,data_start=1,format='ascii')
swb_dbs = Table.read(swbfile,header_start=0,data_start=1,format='ascii')

#chop off the extra entries that are in the SWA file
lenb = len(swb_dbs['reflection'])
swa_dbs = swa_dbs[0:lenb]

lw_dbs_ratio = lwa_dbs['transmission'] / lwb_dbs['transmission']
good = np.where(swa_dbs['microns'] <= 2.7)[0]
sw_dbs_ratio = swa_dbs['reflection'][good] / swb_dbs['reflection'][good]

#create mean SW and LW DBS curves
#print(len(lwa_dbs['transmission']),len(lwb_dbs['transmission']))
#print(len(swa_dbs['reflection']),len(swb_dbs['reflection']))
lwdbs_mean = np.mean([lwa_dbs['transmission'],lwb_dbs['transmission']],axis=0)
swdbs_mean = np.mean([swa_dbs['reflection'],swb_dbs['reflection']],axis=0)

#plot of ratio across modules
f_dbs,a_dbs = plt.subplots()
a_dbs.plot(lwa_dbs['microns'],lw_dbs_ratio,color='black',linestyle='--',label='LW DBS')
a_dbs.plot(swa_dbs['microns'],sw_dbs_ratio,color='black',label='SW DBS')
a_dbs.set_xlabel('Wavlength (microns)')
a_dbs.set_ylabel('Mod A / Mod B Ratio')
a_dbs.legend(loc='best')
a_dbs.set_ylim(0.7,1.3)
f_dbs.savefig('DBS_modAB_ratio.pdf')
plt.close()

#plot of mean DBS throughput values
ff,aa = plt.subplots()
aa.plot(swa_dbs['microns'],swa_dbs['reflection'],color='red',label='SWA')
aa.plot(swb_dbs['microns'],swb_dbs['reflection'],color='blue',label='SWB')
aa.plot(swb_dbs['microns'][good],swdbs_mean,color='green',label='SWMean')
aa.plot(lwa_dbs['microns'],lwa_dbs['transmission'],linestyle='--',color='red',label='LWA')
aa.plot(lwb_dbs['microns'],lwb_dbs['transmission'],linestyle='--',color='blue',label='LWB')
aa.plot(lwb_dbs['microns'],lwdbs_mean,color='green',linestyle='--',label='LWMean')
aa.set_ylim(0,1.)
ff.savefig('DBS_curves.pdf')
plt.close()

#UPDATE: we are now providing separate LW and SW curves to the ETC
#so save the mean DBS curves in text files.
dbsswtab = Table()
dbsswtab['wavelength'] = swb_dbs['microns'][good]
dbsswtab['throughput'] = swdbs_mean
ascii.write(dbsswtab,'NIRCam_DBS_SW_modABmean.dat',overwrite=True)
dbslwtab = Table()
dbslwtab['wavelength'] = lwb_dbs['microns']
dbslwtab['throughput'] = lwdbs_mean
ascii.write(dbslwtab,'NIRCam_DBS_LW_modABmean.dat',overwrite=True)


#find a transition wavelength, where we'll need to cut all curves so that
#we can stitch together the SW and LW curves into a single curve

#first, cut out sections of the LW and SW curves around the wavelengths
#of interest
sw_good = np.where((swa_dbs['microns'] > 2.) & (swa_dbs['microns'] < 2.4))
lw_good = np.where((lwa_dbs['microns'] > 2.) & (lwa_dbs['microns'] < 2.4))

sw_cut = swdbs_mean[sw_good]
lw_cut = lwdbs_mean[lw_good]
sw_cut_wave = swa_dbs['microns'][sw_good]
lw_cut_wave = lwa_dbs['microns'][lw_good]

#now find the wavelength where the SW curve first dives lower than the LW
#curve. We'll call this our transition wavelength
transition_pt = np.where((sw_cut - lw_cut) < 0)[0][0]
transition_wave = sw_cut_wave[transition_pt]
print("In the mean DBS curves, the LW curve becomes higher than the SW curve at {} microns.".format(transition_wave))
print("Using this as the cutoff wavelength for all other optics/contaminants curves in order to combine SW and LW into a single curve, as required by the ETC.")
##Combine the SW and LW DBS cuves into a single curve, around the transition wavelength
#cutpt_sw = np.where(swa_dbs['microns'] == transition_wave)[0]
#cutpt_lw = np.where(lwa_dbs['microns'] == transition_wave)[0]
#final_dbs_wave = np.append(swa_dbs['microns'][0:cutpt_sw],lwa_dbs['microns'][cutpt_lw:])
#final_dbs_trans = np.append(swdbs_mean[0:cutpt_sw],lwdbs_mean[cutpt_lw:])

##f,a = plt.subplots()
##a.plot(final_dbs_wave,final_dbs_trans,color='black')
##f.savefig('DBS_ETC.pdf')
##plt.close()


#QE - coefficients are from Marcia's filter throughput spreadsheets
sw_qe_wavelength = np.arange(0.5,2.8,0.001)
lw_qe_wavelength = np.arange(2.25,5.5,0.001)

sw_qe_coeffs = np.array([0.65830,-0.05668,0.25580,-0.08350])
sw_qe_exponential = 100.
sw_qe_wavecut = 2.38
lw_qe_coeffs_a = np.array([0.934871,0.051541,-0.281664,0.243867,-0.086009,0.014509,-0.001])
lw_qe_factor_a = 0.88
lw_qe_coeffs_b = np.array([2.9104951,-2.182822,0.7075635,-0.071767])

sw_qe = sw_qe_coeffs[0] + sw_qe_coeffs[1]*sw_qe_wavelength + sw_qe_coeffs[2]*sw_qe_wavelength**2 + sw_qe_coeffs[3]*sw_qe_wavelength**3
red = sw_qe_wavelength > sw_qe_wavecut
sw_qe[red] = sw_qe[red] * np.exp((sw_qe_wavecut-sw_qe_wavelength[red])*sw_qe_exponential)


lw_qe_a = lw_qe_factor_a * (lw_qe_coeffs_a[0] + lw_qe_coeffs_a[1]*lw_qe_wavelength + lw_qe_coeffs_a[2]*lw_qe_wavelength**2 + lw_qe_coeffs_a[3]*lw_qe_wavelength**3 + lw_qe_coeffs_a[4]*lw_qe_wavelength**4 + lw_qe_coeffs_a[5]*lw_qe_wavelength**5 + lw_qe_coeffs_a[6]*lw_qe_wavelength**6)
lw_qe_b = lw_qe_coeffs_b[0] + lw_qe_coeffs_b[1]*lw_qe_wavelength + lw_qe_coeffs_b[2]*lw_qe_wavelength**2 + lw_qe_coeffs_b[3]*lw_qe_wavelength**3

#create a mean LW QE curve. NOTE THAT THE MOD A and B CURVES ARE QUITE DIFFERENT
lw_qe = np.mean([lw_qe_a,lw_qe_b],axis=0)


#-----------------------------------
#add F480M filter into the QE plot temporarily, to quantify the mod a/b difference
#filtfile_moda = '../Filters_ASCII-Tables/nircam_throughputs/modA/filters_only/F480M_FM.xlsx_filteronly_modA_sorted.txt'
#filtfile_modb = '../Filters_ASCII-Tables/nircam_throughputs/modB/filters_only/F480M_FM.xlsx_filteronly_modB_sorted.txt'
#filta = ascii.read(filtfile_moda)
#filtb = ascii.read(filtfile_modb)
#filta_interp = np.interp(lw_qe_wavelength,filta['microns'],filta['transmission'],right=0)
#filtb_interp = np.interp(lw_qe_wavelength,filtb['microns'],filtb['transmission'],right=0)

#filta_qe = filta_interp * lw_qe_a
#filtb_qe = filtb_interp * lw_qe_b
#filtab_qe_mean = np.mean([filta_qe,filtb_qe],axis=0)

#print("F480M, Mod A, Mod B, Mod AB Mean integrated throughput")
##cutoff_value = 0.1
##gooda = filta_qe >= cutoff_value
##goodb = filtb_qe >= cutoff_value
##goodab = filtab_qe_mean >=cutoff_value
#cutlambdalow = 4.662
#cutlambdahigh = 4.977
#gooda = ((lw_qe_wavelength >= cutlambdalow) & (lw_qe_wavelength <= cutlambdahigh))
#atot=np.sum(filta_qe[gooda])#,lw_qe_wavelength[gooda])
#btot=np.sum(filtb_qe[gooda])#,lw_qe_wavelength[gooda])
#meantot=np.sum(filtab_qe_mean[gooda])#,lw_qe_wavelength[gooda])
#print(atot,btot,meantot)
#print(atot/meantot,btot/meantot)

#ftmp,(atmp,btmp) = plt.subplots(2,sharex=True)
#atmp.plot(lw_qe_wavelength,filta_qe,color='blue',label='Mod A')
#atmp.plot(lw_qe_wavelength,filtb_qe,color='red',label='Mod B')
#atmp.plot(lw_qe_wavelength,filtab_qe_mean,color='green',label='Mod A/B mean')
#atmp.set_title('F480M * QE, Modules A and B')
#atmp.set_ylabel('Throughput')
#atmp.legend(loc='best',fontsize=10)
#atmp.set_xlim(4.6,5.1)

#btmp.plot(lw_qe_wavelength,(filta_qe/filtab_qe_mean),color='blue',label='Mod A / mean')
#btmp.plot(lw_qe_wavelength,(filtb_qe/filtab_qe_mean),color='red',label='Mod B / mean')
#btmp.plot([4.6,5.2],[1.1,1.1],color='lightgreen',linestyle='--',label='ETC Err Goal')
#btmp.plot([4.6,5.2],[0.9,0.9],color='lightgreen',linestyle='--')
#btmp.plot([4.6,5.2],[1.15,1.15],color='black',linestyle='--',label='ETC Err Req')
#btmp.plot([4.6,5.2],[0.85,0.85],color='black',linestyle='--')
#btmp.legend(loc='upper center',fontsize=10)
#btmp.set_ylabel('Throughput Ratio')
#btmp.set_xlabel('Wavelength (microns)')
#ftmp.savefig('F480M_and_qe_ModA_vs_ModB.pdf')
#-----------------------------------

fqe,aqe = plt.subplots()
aqe.plot(sw_qe_wavelength,sw_qe,color='green',linestyle='--',label='SW')
aqe.plot(lw_qe_wavelength,lw_qe_a,color='blue',label='ModA')
aqe.plot(lw_qe_wavelength,lw_qe_b,color='red',label='ModB')
aqe.plot(lw_qe_wavelength,lw_qe,color='green',label='Mean')
aqe.legend(loc='best')
aqe.set_ylabel('QE')
aqe.set_xlabel('Wavelegth (microns)')
fqe.savefig('QE_curves.pdf')
plt.close()

#save tabulated QE curves
swqetab = Table()
swqetab['wavelength'] = sw_qe_wavelength
swqetab['throughput'] = sw_qe
swqetab['wavelength'].format = '%.4g'
ascii.write(swqetab,'NIRCam_QE_SW.dat',overwrite=True)
lwqetab = Table()
lwqetab['wavelength'] = lw_qe_wavelength
lwqetab['throughput'] = lw_qe
lwqetab['wavelength'].format = '%.4g'
ascii.write(lwqetab,'NIRCam_QE_LW.dat',overwrite=True)

fqr,aqr = plt.subplots()
aqr.plot(lw_qe_wavelength,lw_qe_a/lw_qe,color='blue',label='ModA/Mean')
aqr.plot(lw_qe_wavelength,lw_qe_b/lw_qe,color='red',label='ModB/Mean')
aqr.set_xlabel('Wavelength (microns)')
aqr.set_ylabel('QE Ratios')
aqr.plot([2.25,5.5],[1.1,1.1],color='black',linestyle='--')
aqr.plot([2.25,5.5],[0.9,0.9],color='black',linestyle='--')
aqr.plot([2.25,5.5],[1.15,1.15],color='red',linestyle='--')
aqr.plot([2.25,5.5],[0.85,0.85],color='red',linestyle='--')
aqr.plot([4.662,4.662],[0,2],color='goldenrod',linestyle='--')
aqr.plot([4.977,4.977],[0,2],color='goldenrod',linestyle='--')
aqr.legend(loc='best')
fqr.savefig('QE_ratios_to_mean.pdf')
plt.close()

cutqe_sw = sw_qe_wavelength <= transition_wave
cutqe_lw = lw_qe_wavelength > transition_wave
final_qe_wave = np.append(sw_qe_wavelength[cutqe_sw],lw_qe_wavelength[cutqe_lw])
final_qe = np.append(sw_qe[cutqe_sw],lw_qe[cutqe_lw])


#Optics
ofile = '../DBS_QE_Optics/NIRCam_optics_transmission_29Oct2015.csv'

opttab = ascii.read(ofile,header_start=1,data_start=2,format='csv')

wave = opttab['Wavelength'].data.data
nvr_thru = opttab['NVR_Transmission'].data.data
nvr_wave = opttab['NVR_Wavelength'].data.data
collimator = opttab['Collimator'].data.data
sw_triplet = replace_nan(opttab['SW_Triplet'].data.data).astype('float')
sw_mirrors = replace_nan(opttab['SW_Mirrors'].data.data).astype('float')
lw_triplet = replace_nan(opttab['LW_triplet'].data.data).astype('float')
lw_mirrors = replace_nan(opttab['LW_Mirrors'].data.data).astype('float')
sw_particulates = replace_nan(opttab['SW_Particulates'].data.data).astype('float')
lw_particulates = replace_nan(opttab['LW_Particulates'].data.data).astype('float')

#remove extra entries in NVR columns
good = np.where(nvr_wave != 0.)[0]
nvr_thru = nvr_thru[good]
nvr_wave = nvr_wave[good]

#interpolate NVR to the same wavelength scale as the other columns
nvr_interp = np.interp(wave,nvr_wave,nvr_thru)

#combine the elements to produce a SW optics curve and a LW optics curve
#The 0.98 factor is a 'contingency factor' from Lockheed Martin that Marcia used
#for her calculations. She instructed us to keep it in for ETC testing. Now that
#ETC testing is over there's no need to keep it in, although it is less than the
#10% uncertainty in the ETC anyway, and also lower than the uncertainty in other
#contributing elements (such as QE), so keeping it in is fine (this was Marcia's
#reasoning. See emails from mid-Jan 2017). But since I'm remaking this file
#anyway, John agrees we can take the 0.98 out, and that'll be easier because we
#won't have to explain it.

#The 0.984 and 0.989 factors are for scatter, and are come from Marcia. This term is
#not present in the optics file above. The term is given as one value per filter, and
#is present in /grp/jwst/wit/nircam/hilbert/filter_throughputs/nircam_system_throughputs/final_version/no_nans/nircam_optics_filter_average.dat. 0.984 for all short-wave filters and
#0.989 for all long-wave filters.
sw_optics = collimator * sw_triplet * sw_mirrors * sw_particulates * nvr_interp * 0.984
lw_optics = collimator * lw_triplet * lw_mirrors * lw_particulates * nvr_interp * 0.989

#remove nan entries
keepsw = ~np.isnan(sw_optics)
sw_optics = sw_optics[keepsw]
sw_optics_wave = wave[keepsw]
keeplw = ~np.isnan(lw_optics)
lw_optics = lw_optics[keeplw]
lw_optics_wave = wave[keeplw]

#stitch together SW and LW curves into a single curve
cutopt_sw = sw_optics_wave <= transition_wave
cutopt_lw = lw_optics_wave > transition_wave
final_optics = np.append(sw_optics[cutopt_sw],lw_optics[cutopt_lw])
final_optics_wave = np.append(sw_optics_wave[cutopt_sw],lw_optics_wave[cutopt_lw])

#temp optics components plot
ff,aa = plt.subplots()
aa.plot(sw_optics_wave,sw_optics,color='blue',label='SW')
aa.plot(lw_optics_wave,lw_optics,color='red',label='LW')
aa.plot(final_optics_wave,final_optics,color='green',linestyle='--',label='Comb')
aa.legend(loc='best')
ff.savefig('combined_SWLW_opticsonly.pdf')
#t=Table()
#t['wave']=final_optics_wave
#t['final_optics']=final_optics
#ascii.write(t,'temp.dat')


#plot showing all stitched curves
f,a = plt.subplots()
#a.plot(final_dbs_wave,final_dbs_trans,color='black',label='DBS')
a.plot(final_qe_wave,final_qe,color='red',label='QE')
a.plot(final_optics_wave,final_optics,color='blue',label='Optics')
a.legend(loc='best')
a.set_ylim(0,1.05)
a.set_xlabel('Wavelength (microns)')
a.set_ylabel('Throughput')
f.savefig('ETC_component_curves.pdf')
plt.close()


#combine stitched curves

#interpolate to match wavelengths
qe_interp = np.interp(final_optics_wave,final_qe_wave,final_qe)
#optics_interp = np.interp(final_dbs_wave,final_optics_wave,final_optics)

#roll curves up together
#final_qeopt = qe_interp * final_optics

#plot of final curve
fff,aaa = plt.subplots()
aaa.plot(final_optics_wave,final_optics,color='black')
aaa.set_xlabel('Wavelength (microns)')
aaa.set_ylabel('Optics Throughput')
fff.savefig('final_opt_curve.pdf')
plt.close()

#write out file containing the final, rolled-up curve
fin = Table()
fin['wavelength'] = final_optics_wave
fin['throughput'] = final_optics
ascii.write(fin,'NIRCam_optics_contaminants_ETC_curve.dat',overwrite=True)

