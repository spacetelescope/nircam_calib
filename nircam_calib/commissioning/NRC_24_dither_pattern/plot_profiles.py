#!/usr/bin/env python
import argparse
from astropy.io import fits
from astropy.table import Table,Column, MaskedColumn
#
from glob import glob
import inspect
#
from   matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
#
import os
from os.path import exists
from pathlib import Path
#
import re
import sys

#------------------------------------------------------------------------
#
def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno
#
#------------------------------------------------------------------------
#
#
def plot_profiles(sca, sim, filter, type, observation, car, debug):
    if(debug == 1) :
        print("sca is ", sca)
        print("sim is ", sim)
        print("filter is ", filter)
        print("type is ", type)
        print("car  is ", car)
        print("observation is ", observation)
    # Get files but first, what computer are we on ?
    host = os.environ.get('HOST')
    #root_dir = "car_"+car+"/data/"
    target_dir = "car_"+car+"/analysis/"
    plot_dir   = "car_"+car+"/plots/"
    sca_number = dict([
        ('a1','481'),
        ('a2','482'),
        ('a3','483'),
        ('a4','484'),
        ('a5','485'),
        ('b1','486'),
        ('b2','487'),
        ('b3','488'),
        ('b4','489'),
        ('b5','490')
    ])
    
    if(host == 'ema.as.arizona.edu'):
        plot_dir = "/home/cnaw/commissioning/car_24_apt_01073/plots/"
        if(sim == "mirage"):
            target_dir = '/home/cnaw/commissioning/car_24_apt_01073/mirage/analysis/'
            string   = target_dir+'*'+observation+'00*_*'+sca+'_'+type+'_'+filter+'*prof.fits'
        else:
            target_dir = '/home/cnaw/commissioning/car_24_apt_01073/guitarra/analysis/'
            if(type == "slp"):
                string   = target_dir+'*'+observation+'*0*_'+filter+'*'+sca_number[sca]+'*'+type+'*prof.fits'
            else:
                string   = target_dir+'*'+observation+'00*'+sca+'_'+filter+'*_'+type+'_*prof.fits'
        template_dir = '/home/cnaw/python/commissioning/templates/'
    if(host == 'orange.as.arizona.edu'):
        #    root_dir   = '/data1/car_24_apt_01073/mirage/reduced/'
        plot_dir = "/data1/car_24_apt_01073/plots/"
        if(sim == "mirage"):
            target_dir = '/data1//car_24_apt_01073/mirage/analysis/'
            string   = target_dir+'*'+observation+'00*_*'+sca+'_'+type+'_'+filter+'*prof.fits'
        else:
            target_dir = '/data1/car_24_apt_01073/guitarra/analysis/'
            if(type == "slp"):
                string   = target_dir+'apt*'+observation+'*_'+filter+'*'+sca_number[sca]+'*'+type+'*prof.fits'
            else:
                string   = target_dir+'*'+observation+'00*_*'+sca+'_'+type+'_'+filter+'*prof.fits'
            
        template_dir = './'

    status = exists(target_dir)

    colors =['black', 'red', 'green',  'blue','goldenrod', 'cyan', 'magenta',
             'teal','violet', 'grey', 'beige', 'olive','wheat', 'chartreuse',
             'orange',    'orchid', 'plum', 'purple', 'salmon','sienna',
             'silver', 'tan', 'tomato','turquoise', 'yellow','yellowgreen',
             'aqua', 'aquamarine', 'azure', 'chocolate', 'coral', 'crimson',
             'darkblue', 'darkgreen', 'fuchsia', 'gold', 'indigo', 'khaki',
             'lavender', 'lightblue', 'lightgreen', 'lime', 'maroon', 'navy',
             'orangered', 'pink', 'plum', 'purple','silver','tan']

    title = 'obs: '+observation+' '+filter+' sca: '+sca+' reduction:'+type
    title = sim+': '+title
    psf= dict([
        ('F070W', 'PSF_NIRCam_F070W_fov_386_os_2_prof.fits'),
        ('F115W', 'PSF_NIRCam_F115W_fov_508_os_2_prof.fits'),
        ('F150W', 'PSF_NIRCam_F150W_fov_671_os_2_prof.fits'),
        ('F277W', 'PSF_NIRCam_F277W_fov_595_os_2_prof.fits')
    ])

    template = template_dir+psf[filter]

    list =  sorted(glob(string))
    list.append(template)
    if(debug == 1) :
        print("at line : ", lineno(), " string is ", string)
        print("plot_profiles.py: ",len(list))
        print("plot_profiles.py:\n",list)
        
    if(len(list) == 1) :
        png_plot = None
        return png_plot

    
    png_plot = plot_dir+'nrc'+sca+'_obs_'+observation+'_'+filter+'_'+sim+'_'+type+'_profile.png'

    #   png_plot = '/home/cnaw/Pictures/car_24/nrc'+sca+'_'+filter+'_'+sim+'_'+type+'_profile.png'

    font = {'family' : 'DejaVu Sans',
            'weight' : 'bold',
            'size'   : 6}
    rc('font', **font)
    fig = plt.figure(figsize=(6,4.5))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(left=-0.01, right=1.1)
    ax.set_ylim(bottom=1.0e-5, top=1.5)
    plt.yscale("log")
    plt.xlabel("radius (arc sec)")
    plt.ylabel("Normalized surface brightess profile (flux/arcsec^2)")
    plt.title(title)

    for ii  in range(len(list)) :
        file = list[ii]
        hdulist = fits.open(file)
        table = hdulist[1].data
        hdulist.close()

        rt = table.field(0)
        # normalised by central flux
        #        encircled = table.field(4)
        # normalized by integral flux
        encircled = table.field(6)
        if(debug == 1):
            print("len(rt), len(encircled)",ii, len(rt), len(encircled))
        xx = 0.20
        if( ii == len(list)-1) :
            jj = 0
            plt.plot(rt, encircled,'+-', color='black')
            yy = 10.0**(-(ii*0.2))
            junk = file.split("/")
            text = junk[len(junk)-1]
            if(debug == 1):
                print("file ", xx, yy, text, colors[0])
            plt.text(xx,yy,text,color=colors[jj])
        else:
            jj = ii + 1
            if(jj >= len(colors)):
                jj= jj -50
            if(debug == 1):
                print("at line ", lineno(),"  file ", file, colors[jj])
            yy = 10.0**(-(ii*0.2))
            junk = file.split("/")
            text = junk[len(junk)-1]
            if(debug == 1):
                print("file ", ii, xx, yy, file, colors[jj])
            plt.text(xx,yy,text,color=colors[jj])
            plt.plot(rt, encircled,'.', color=colors[jj])

    
            
    plt.savefig(png_plot,bbox_inches='tight')
    #    plt.show()
#    print("plot saved as ", png_plot)
    return png_plot
#
#=======================================================================
#
if __name__ == "__main__":

#    filter      = "F150W"
#    observation = "00*"
#    sca         = "b3"
#    type        = "cal"
    # Use argument parser to pass keywords to program
    parser = argparse.ArgumentParser(description="Run multiple instances of average_psf.py")
    parser.add_argument("--sca", help="SCA name (a1...a5,b1...b5)", type=str, default="b3")
    parser.add_argument("--sim", help="guitarra or mirage", type=str, default="mirage")
    parser.add_argument("--filter", help="NIRCam filter", type=str, default="F150W")
    parser.add_argument("--type", help="file type (cal, i2d, slp)", type=str,default="cal")
    parser.add_argument("--obs", help="observation number", type=str, default="00*")
    parser.add_argument("--car", help="car number", type=str, default="24")
    parser.add_argument("--debug", help="to debug set to 1 ", type=int, default=0)

    #parser.add_argument("--test", help="Testing. Skips run os.system() command. Fake file names.", action="store_true")
    args = parser.parse_args()

    name = plot_profiles(args.sca, args.sim, args.filter, args.type, args.obs, args.car, args.debug)
    print(name)
