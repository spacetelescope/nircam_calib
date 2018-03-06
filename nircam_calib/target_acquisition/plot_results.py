#! /usr/bin/env python

'''
Plot the final results of the centroiding with saturated pixels study
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.stats import sigmaclip

tabfile = 'Final_results_table.tab'
pixscale = 0.063
jitterarcsec = 0.007
jitterpix = jitterarcsec / pixscale

tab = Table.read(tabfile,format='ascii')

realx = tab['True_Centroid_x']
realy = tab['True_Centroid_y']
tablen = len(realx)

centerx = np.zeros_like(realx) + 15.
centery = np.zeros_like(realy) + 15.
y16 = realy == 16.5
print("Found {} entries where y is 16.5".format(np.sum(y16)))
centery[y16] = 16.
dist_from_pix_center = np.sqrt(np.abs(centerx-realx)**2 + np.abs(centery-realy)**2)

distances = np.unique(dist_from_pix_center)
print("Found {} unique distances.".format(len(distances)))
print(distances)
#colors = ['red','blue','green','orange','black','indigo','khaki','lightgray','purple','cyan','violet','lightsalmon']
psym = ['o','s','p','P','h','X','d','v','*']

f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='red',marker=psym[i],label=str(round(dist, 2)))
a.plot([0,11],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
a.set_xlabel('K_G2V Magnitude')
a.set_ylabel('Radial Centroid Error (arcsec)',color='red')
a.legend(loc='upper right',title='Dist. from pix center')
a.tick_params('y', colors='red')
aylim = a.get_ylim()
b = a.twinx()
b.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Pix'].data[good],color='red',marker=psym[i],label=str(dist))
b.tick_params('y', colors='blue')
b.set_ylabel('Radial Centroid Error (pixels)',color='blue')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['K_G2V_Mag'].data)-0.25,np.max(tab['K_G2V_Mag'].data)+0.25)
f.tight_layout()
f.savefig('radial_error_vs_kmag.pdf')


f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='red',marker=psym[i],label=str(round(dist, 2)))
a.plot([0,11],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
a.set_xlabel('K_G2V Magnitude')
a.set_ylabel('Radial Centroid Error (arcsec)',color='red')
a.legend(loc='upper right',title='Dist. from pix center')
a.tick_params('y', colors='red')
a.set_ylim(0,pixscale)
aylim = a.get_ylim()
b = a.twinx()
b.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Pix'].data[good],color='red',marker=psym[i],label=str(dist))
b.tick_params('y', colors='blue')
b.set_ylabel('Radial Centroid Error (pixels)',color='blue')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['K_G2V_Mag'].data)-0.25,np.max(tab['K_G2V_Mag'].data)+0.25)
f.tight_layout()
f.savefig('radial_error_vs_kmag_zoom.pdf')


f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='red',marker=psym[i],label=str(round(dist, 2)))
    #a[1].scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Pix'].data[good],color=colors[i],marker='o',label=str(round(dist, 2)))
a.plot([0,40],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
#a[1].plot([0,40],[jitterpix,jitterpix],color='black',linestyle='--',label='JWST jitter')
a.set_xlim(0,np.max(tab['Sat_in_all_groups'].data)+1)
a.set_xlabel('Num. Pixels Saturated in all 3 Groups')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.legend(loc='upper left',title='Dist. from pix center')
a.tick_params('y', colors='black')
aylim = a.get_ylim()
b = a.twinx()
b.scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Pix'].data[good],color='red',marker=psym[i],label=str(dist))
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['Sat_in_all_groups'].data)-0.25,np.max(tab['Sat_in_all_groups'].data)+0.25)
plt.show()
#f.tight_layout()
#f.savefig('radial_error_vs_fully_saturated_pix.pdf')


f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='red',marker=psym[i],label=str(round(dist, 2)))
    #a[1].scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Pix'].data[good],color=colors[i],marker='o',label=str(round(dist, 2)))
a.plot([0,40],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
#a[1].plot([0,40],[jitterpix,jitterpix],color='black',linestyle='--',label='JWST jitter')
a.set_xlim(0,np.max(tab['Sat_in_all_groups'].data)+1)
a.set_xlabel('Num. Pixels Saturated in all 3 Groups')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.set_ylim(0,pixscale)
a.legend(loc='upper left',title='Dist. from pix center')
a.tick_params('y', colors='red')
aylim = a.get_ylim()
b = a.twinx()
b.scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Pix'].data[good],color='red',marker=psym[i],label=str(dist))
b.tick_params('y', colors='blue')
b.set_ylabel('Radial Centroid Error (pixels)',color='blue')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['Sat_in_all_groups'].data)-0.25,10.25)
f.tight_layout()
f.savefig('radial_error_vs_fully_saturated_pix_zoom.pdf')



total_sat_23 = tab['Sat_in_all_groups'].data + tab['Sat_in_grps_2_3'].data
f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(total_sat_23[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='red',marker=psym[i],label=str(round(dist, 2)))
a.plot([0,100],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
a.set_xlim(0,np.max(total_sat_23)+1)
a.set_xlabel('Num. Pixels Saturated in Groups 2 and 3')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.legend(loc='upper left',title='Dist. from pix center')
a.tick_params('y', colors='red')
aylim = a.get_ylim()
b = a.twinx()
b.scatter(total_sat_23[good],tab['Centroid_Radial_Error_Pix'].data[good],color='red',marker=psym[i],label=str(dist))
b.tick_params('y', colors='blue')
b.set_ylabel('Radial Centroid Error (pixels)',color='blue')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(total_sat_23)-0.25,np.max(total_sat_23)+0.25)
f.tight_layout()
f.savefig('radial_error_vs_partially_saturated_pix.pdf')


f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(total_sat_23[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='red',marker=psym[i],label=str(round(dist, 2)))
a.plot([0,100],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
a.set_xlim(0,np.max(total_sat_23)+1)
a.set_xlabel('Num. Pixels Saturated in Groups 2 and 3')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.legend(loc='upper left',title='Dist. from pix center')
a.tick_params('y', colors='red')
a.set_ylim(0,pixscale)
aylim = a.get_ylim()
b = a.twinx()
b.scatter(total_sat_23[good],tab['Centroid_Radial_Error_Pix'].data[good],color='red',marker=psym[i],label=str(dist))
b.tick_params('y', colors='blue')
b.set_ylabel('Radial Centroid Error (pixels)',color='blue')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(total_sat_23)-0.25,19.25)
f.tight_layout()
f.savefig('radial_error_vs_partially_saturated_pix_zoom.pdf')








#Try plotting all points as back dots, then overplot the mean and stdev for
#each distance as a different color/symbol/etc
mns = {}
devs = {}
xs = {}
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    gooddata = tab['Centroid_Radial_Error_Arcsec'].data[good]
    goodx = tab['K_G2V_Mag'].data[good]
    xval = np.unique(goodx)
    mnlist = []
    devlist = []
    xlist = []
    for x in xval:
        goodmag = goodx == x
        ys = gooddata[goodmag]
        clipped,lower,upper = sigmaclip(ys,low=3,high=3)
        mn = np.mean(clipped)
        dev = np.std(clipped)
        mnlist.append(mn)
        devlist.append(dev)
        xlist.append(x)
    mns[i] = mnlist
    devs[i] = devlist
    xs[i] = xlist

f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    gooddata = tab['Centroid_Radial_Error_Arcsec'].data[good]
    a.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.02)#,label=str(round(dist, 2)))

colors = ['red','green','blue','orange','pink','skyblue']
for i,dist in enumerate(distances):
    key = i
    #a.errorbar(xs[key],mns[key],yerr=devs[key],fmt=psym[i],color=colors[i],label=str(round(dist, 2)))
    a.scatter(xs[key],mns[key],marker=psym[i],color=colors[i],label=str(round(dist, 2)))
a.plot([0,11],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
#a.set_xlabel('K_G2V Magnitude')
a.set_xlabel(r'$K_{Vega}$')
a.set_ylabel('Radial Centroid Error (arcsec)',color='black')
a.legend(loc='upper right',title='Dist. from pix center')
a.tick_params('y', colors='black')
a.set_xlim(2.8,7.6)
a.set_ylim(0,0.062)
aylim = a.get_ylim()
b = a.twinx()
#b.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Pix'].data[good],color='black',marker=psym[i],label=str(dist))
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam Longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['K_G2V_Mag'].data)-0.25,np.max(tab['K_G2V_Mag'].data)+0.25)
#plt.show()
f.tight_layout()
f.savefig('radial_error_vs_kmag_vs_distances.pdf')



#Try plotting all points as back dots, then overplot the overall mean and
#stdev as a red line with error bars
xval = np.unique(tab['K_G2V_Mag'].data)
mnlist = []
devlist = []
xlist = []
for x in xval:
    goodmag = tab['K_G2V_Mag'].data == x
    ys = tab['Centroid_Radial_Error_Arcsec'].data[goodmag]
    clipped,lower,upper = sigmaclip(ys,low=3,high=3)
    mn = np.mean(clipped)
    dev = np.std(clipped)
    mnlist.append(mn)
    devlist.append(dev)
    xlist.append(x)

f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    gooddata = tab['Centroid_Radial_Error_Arcsec'].data[good]
    a.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)

a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')#,label=str(round(dist, 2)))
#a.plot([0,11],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
a.set_xlabel(r'$K_{Vega}$')#'K_G2V Magnitude')
a.set_ylabel('Radial Centroid Error (arcsec)',color='black')
a.legend(loc='upper right')
a.tick_params('y', colors='black')
a.set_xlim(2.8,7.6)
a.set_ylim(0,0.062)
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam Longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['K_G2V_Mag'].data)-0.25,np.max(tab['K_G2V_Mag'].data)+0.25)
#plt.show()
f.tight_layout()
f.savefig('radial_error_vs_kmag_overallmean.pdf')


xval = np.unique(tab['K_G2V_Mag'].data)
mnlist = []
devlist = []
xlist = []
for x in xval:
    goodmag = tab['K_G2V_Mag'].data == x
    ys = tab['Centroid_Radial_Error_Arcsec'].data[goodmag]
    clipped,lower,upper = sigmaclip(ys,low=3,high=3)
    mn = np.mean(clipped)
    dev = np.std(clipped)
    mnlist.append(mn)
    devlist.append(dev)
    xlist.append(x)

f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    gooddata = tab['Centroid_Radial_Error_Arcsec'].data[good]
    a.scatter(tab['K_G2V_Mag'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)

a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')#,label=str(round(dist, 2)))
#a.plot([0,11],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
#a.set_xlabel('K_G2V Magnitude')
a.set_xlabel(r'$K_{Vega}$')
a.set_ylabel('Radial Centroid Error (arcsec)',color='black')
a.legend(loc='upper right')
a.tick_params('y', colors='black')
#a.set_xlim(2.8,7.6)
#a.set_ylim(0,0.062)
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam Longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
#a.set_xlim(np.min(tab['K_G2V_Mag'].data)-0.25,np.max(tab['K_G2V_Mag'].data)+0.25)
#plt.show()
f.tight_layout()
f.savefig('radial_error_vs_kmag_overallmean_expandedview.pdf')


total_sat_23 = tab['Sat_in_all_groups'].data + tab['Sat_in_grps_2_3'].data
xval = np.unique(total_sat_23)
mnlist = []
devlist = []
xlist = []
for x in xval:
    goodmag = total_sat_23 == x
    ys = tab['Centroid_Radial_Error_Arcsec'].data[goodmag]
    clipped,lower,upper = sigmaclip(ys,low=3,high=3)
    mn = np.mean(clipped)
    dev = np.std(clipped)
    mnlist.append(mn)
    devlist.append(dev)
    xlist.append(x)
f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(total_sat_23[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)
a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')
a.set_xlim(0,np.max(total_sat_23)+1)
a.set_xlabel('Number of Pixels Saturated in Group 2 of 3')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.tick_params('y', colors='black')
a.set_ylim(0,pixscale)
a.set_xticks(np.arange(0,17,2))
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam Longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(total_sat_23)-0.25,17.25)
f.tight_layout()
#plt.show()
f.savefig('radial_error_vs_partially_saturated_pix_zoom.pdf')


total_sat_23 = tab['Sat_in_all_groups'].data + tab['Sat_in_grps_2_3'].data
xval = np.unique(total_sat_23)
mnlist = []
devlist = []
xlist = []
for x in xval:
    goodmag = total_sat_23 == x
    ys = tab['Centroid_Radial_Error_Arcsec'].data[goodmag]
    clipped,lower,upper = sigmaclip(ys,low=3,high=3)
    mn = np.mean(clipped)
    dev = np.std(clipped)
    mnlist.append(mn)
    devlist.append(dev)
    xlist.append(x)
f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(total_sat_23[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)#label=str(round(dist, 2)))
a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')
a.set_xlabel('Number of Pixels Saturated in Group 2 of 3')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.tick_params('y', colors='black')
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam Longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(total_sat_23)-0.25,np.max(total_sat_23)+0.25)
f.tight_layout()
#plt.show()
f.savefig('radial_error_vs_partially_saturated_pix.pdf')



xval = np.unique(tab['Sat_in_all_groups'].data)
mnlist = []
devlist = []
xlist = []
for x in xval:
    goodmag = tab['Sat_in_all_groups'].data == x
    ys = tab['Centroid_Radial_Error_Arcsec'].data[goodmag]
    clipped,lower,upper = sigmaclip(ys,low=3,high=3)
    mn = np.mean(clipped)
    dev = np.std(clipped)
    mnlist.append(mn)
    devlist.append(dev)
    xlist.append(x)
f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)#label=str(round(dist, 2)))
a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')
a.set_xlabel('Number of Fully Saturated Pixels')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.set_ylim(0,pixscale)
a.tick_params('y', colors='black')
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['Sat_in_all_groups'].data)-0.25,10.25)
#plt.show()
f.tight_layout()
f.savefig('radial_error_vs_fully_saturated_pix_zoom.pdf')


xval = np.unique(tab['Sat_in_all_groups'].data)
mnlist = []
devlist = []
xlist = []
for x in xval:
    goodmag = tab['Sat_in_all_groups'].data == x
    ys = tab['Centroid_Radial_Error_Arcsec'].data[goodmag]
    clipped,lower,upper = sigmaclip(ys,low=3,high=3)
    mn = np.mean(clipped)
    dev = np.std(clipped)
    mnlist.append(mn)
    devlist.append(dev)
    xlist.append(x)
f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(tab['Sat_in_all_groups'].data[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)#label=str(round(dist, 2)))
a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')
a.set_xlabel('Number of Fully Satruated Pixels')
a.set_ylabel('Radial Centroid Error (arcsec)')
#a.set_ylim(0,pixscale)
a.tick_params('y', colors='black')
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(tab['Sat_in_all_groups'].data)-0.25,np.max(tab['Sat_in_all_groups'].data)+0.25)
#plt.show()
f.tight_layout()
f.savefig('radial_error_vs_fully_saturated_pix.pdf')




# plot of ETC's "partially" saturated pixels, meaning
# pixels that saturate anytime after group 1
#total_sat_123 = tab['Sat_in_all_groups'].data + tab['Sat_in_grps_2_3'].data + tab['Sat_in_grp_3_only'].data
total_sat_123 = tab['Sat_in_grps_2_3'].data + tab['Sat_in_grp_3_only'].data
xval = np.unique(total_sat_123)
mnlist = []
devlist = []
xlist = []
for x in xval:
    goodmag = total_sat_123 == x
    ys = tab['Centroid_Radial_Error_Arcsec'].data[goodmag]
    clipped,lower,upper = sigmaclip(ys,low=3,high=3)
    mn = np.mean(clipped)
    dev = np.std(clipped)
    mnlist.append(mn)
    devlist.append(dev)
    xlist.append(x)
f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(total_sat_123[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)#label=str(round(dist, 2)))
a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')
a.set_xlabel('Number of Partially Saturated Pixels')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.tick_params('y', colors='black')
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam Longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(total_sat_123)-0.25,np.max(total_sat_123)+0.25)
f.tight_layout()
#plt.show()
f.savefig('radial_error_vs_ETCs_partially_saturated_pix.png')

f,a = plt.subplots(figsize=(8,5))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    a.scatter(total_sat_123[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='black',marker='o',alpha=0.05)#label=str(round(dist, 2)))
a.errorbar(xlist,mnlist,yerr=devlist,fmt=psym[i],color='red')
a.set_xlabel('Number of Partially Saturated Pixels')
a.set_ylabel('Radial Centroid Error (arcsec)')
a.tick_params('y', colors='black')
a.set_ylim(0.,pixscale)
aylim = a.get_ylim()
b = a.twinx()
b.tick_params('y', colors='black')
b.set_ylabel('Radial Centroid Error (NIRCam Longwave pixels)',color='black')
b.set_ylim(aylim[0]/pixscale,aylim[1]/pixscale)
a.set_xlim(np.min(total_sat_123)-0.25,22.5)
f.tight_layout()
#plt.show()
f.savefig('radial_error_vs_ETCs_partially_saturated_pix_zoom.png')











sys.exit()

f,a = plt.subplots(2,sharex=True,figsize=(8,8))
for i,dist in enumerate(distances):
    good = dist_from_pix_center == dist
    #print("Plotting {} points for distance {}, {} {}".format(np.sum(good),dist,i,colors[i]))
    a[0].scatter(total_sat_23[good],tab['Centroid_Radial_Error_Arcsec'].data[good],color='red',marker='o',label=str(dist))
    a[1].scatter(total_sat_23[good],tab['Centroid_Radial_Error_Pix'].data[good],color=colors[i],marker='o',label=str(dist))
a[0].plot([0,100],[jitterarcsec,jitterarcsec],color='black',linestyle='--',label='JWST jitter')
a[1].plot([0,100],[jitterpix,jitterpix],color='black',linestyle='--',label='JWST jitter')
a[0].set_xlim(0,21)
a[0].set_ylim(0,0.1)
a[1].set_ylim(0,0.1/pixscale)
a[1].set_xlabel('Num. Pixels Saturated in Groups 2 and 3')
a[0].set_ylabel('Radial Centroid Error (arcsec)')
a[0].legend(loc='upper left',title='Dist. from pix center')
a[1].set_ylabel('Radial Centroid Error (pixels)')
f.savefig('radial_error_vs_partially_saturated_pix_zoom.pdf')
#plt.show()
