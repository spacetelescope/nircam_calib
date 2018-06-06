import copy,os,sys
import pytest
import numpy as np
from glob import glob
import json
from astropy.io import fits, ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from jwst.pipeline import Spec2Pipeline


def test_extract_2D(cases, direct_image_file):
    '''Test extract 2D.'''

    save_figs = True

    # load direct image
    direct_image = fits.getdata(direct_image_file,1)

    # get file names for catalog and dispersed image from association
    with open(cases) as json_data:
         d = json.load(json_data)
    members = d['products'][0]['members']
    for row in np.arange(0,len(members)):
        if members[row]['exptype'] == "sourcecat":
            catalog_file = members[row]['expname']
        elif members[row]['exptype'] == "science":
            dispersed_file = members[row]['expname']

    # get the catalog and dispersed image
    catalog_in = ascii.read(catalog_file)
    short_catalog = Table({'RA': catalog_in['icrs_centroid'].ra,
                              'Dec': catalog_in['icrs_centroid'].dec,
                              'PixelX': catalog_in['xcentroid'],
                              'PixelY': catalog_in['ycentroid'],
                              'Mag': catalog_in['abmag']},
                              names=['RA', 'Dec', 'PixelX', 'PixelY', 'Mag'])
    dispersed_image = fits.getdata(dispersed_file)

    # run the pipeline
    spec2 = Spec2Pipeline.call(cases,config_file="calwebb_spec2.cfg")
    extract_2D_output = glob(dispersed_file[:-10]+"*_cal.fits")[0]


    # save the plots for checking by eye (for now)
    if save_figs == True:

        # plot of direct image with sources in catalog circled
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=22)
        plt.xlabel('x pixels',fontsize=22)
        plt.imshow(direct_image, vmin=-0.02, vmax=0.08, cmap=plt.cm.gray, origin='lower')
        ax.scatter(short_catalog['PixelX'][0],short_catalog['PixelY'][0], s=5000,marker='o',facecolors='none', edgecolors='red',lw=3, label='Source 1')
        ax.scatter(short_catalog['PixelX'][1],short_catalog['PixelY'][1], s=1000,marker='o',facecolors='none', edgecolors='limegreen',lw=3, label='Source 2')
        ax.scatter(short_catalog['PixelX'][2],short_catalog['PixelY'][2], s=1000,marker='o',facecolors='none', edgecolors='lightblue',lw=3, label='Source 3')
        ax.text(short_catalog['PixelX'][0]-110,short_catalog['PixelY'][0]+130, "Source 1", fontsize=22,color='red')
        ax.text(short_catalog['PixelX'][1]-110,short_catalog['PixelY'][1]+90, "Source 2", fontsize=22,color='limegreen')
        ax.text(short_catalog['PixelX'][2]-110,short_catalog['PixelY'][2]+90, "Source 3", fontsize=22,color='lightblue')
        ax.set_title("Direct image and catalog sources",fontsize=22)
        plt.colorbar(orientation='horizontal',pad=0.05)
        plt.savefig(dispersed_file[:-5]+'_direct_img_sources.png')

        # source data and headers
        source1_order1 = fits.getdata(extract_2D_output,1)
        src1_hdr1 = fits.getheader(extract_2D_output,1)
        source1_order2 = fits.getdata(extract_2D_output,4)
        src1_hdr2 = fits.getheader(extract_2D_output,4)
        source2_order1 = fits.getdata(extract_2D_output,7)
        src2_hdr1 = fits.getheader(extract_2D_output,7)
        source2_order2 = fits.getdata(extract_2D_output,10)
        src2_hdr2 = fits.getheader(extract_2D_output,10)
        source3_order1 = fits.getdata(extract_2D_output,13)
        src3_hdr1 = fits.getheader(extract_2D_output,13)
        source3_order2 = fits.getdata(extract_2D_output,16)
        src3_hdr2 = fits.getheader(extract_2D_output,16)


        # plot showing dispersed image and extract 2D regions highlighted
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=22)
        plt.xlabel('x pixels',fontsize=22)
        plt.imshow(dispersed_image, vmin=0.2, vmax=1.5, cmap=plt.cm.gray, origin='lower')
        ax.scatter(src1_hdr1['SRCXPOS'],src1_hdr1['SRCYPOS'], s=1000,marker='o',facecolors='none', edgecolors='red',lw=3, label='Source 1')
        ax.scatter(src2_hdr1['SRCXPOS'],src2_hdr1['SRCYPOS'], s=1000,marker='o',facecolors='none', edgecolors='limegreen',lw=3, label='Source 2')
        ax.scatter(src3_hdr1['SRCXPOS'],src3_hdr1['SRCYPOS'], s=1000,marker='o',facecolors='none', edgecolors='lightblue',lw=3, label='Source 3')

        rect = patches.Rectangle((src1_hdr1['SLTSTRT1'],src1_hdr1['SLTSTRT2']),src1_hdr1['SLTSIZE1'],src1_hdr1['SLTSIZE2'],linewidth=1,edgecolor='r',lw=3,facecolor='none')
        ax.add_patch(rect)
        rect = patches.Rectangle((src1_hdr2['SLTSTRT1'],src1_hdr2['SLTSTRT2']),src1_hdr2['SLTSIZE1'],src1_hdr2['SLTSIZE2'],linewidth=1,edgecolor='r',lw=3,facecolor='none')
        ax.add_patch(rect)
        rect = patches.Rectangle((src2_hdr1['SLTSTRT1'],src2_hdr1['SLTSTRT2']),src2_hdr1['SLTSIZE1'],src2_hdr1['SLTSIZE2'],linewidth=1,edgecolor='limegreen',lw=3,facecolor='none')
        ax.add_patch(rect)
        rect = patches.Rectangle((src2_hdr2['SLTSTRT1'],src2_hdr2['SLTSTRT2']),src2_hdr2['SLTSIZE1'],src2_hdr2['SLTSIZE2'],linewidth=1,edgecolor='limegreen',lw=3,facecolor='none')
        ax.add_patch(rect)
        rect = patches.Rectangle((src3_hdr1['SLTSTRT1'],src3_hdr1['SLTSTRT2']),src3_hdr1['SLTSIZE1'],src3_hdr1['SLTSIZE2'],linewidth=1,edgecolor='lightblue',lw=3,facecolor='none')
        ax.add_patch(rect)
        rect = patches.Rectangle((src3_hdr2['SLTSTRT1'],src3_hdr2['SLTSTRT2']),src3_hdr2['SLTSIZE1'],src3_hdr2['SLTSIZE2'],linewidth=1,edgecolor='lightblue',lw=3,facecolor='none')
        ax.add_patch(rect)

        ax.text(src1_hdr1['SLTSTRT1']+10,src1_hdr1['SLTSTRT2']+10, "order 1", fontsize=16,color='red')
        ax.text(src1_hdr2['SLTSTRT1']+10,src1_hdr2['SLTSTRT2']+10, "order 2", fontsize=16,color='red')
        ax.text(src2_hdr1['SLTSTRT1']+10,src2_hdr1['SLTSTRT2']+10, "order 1", fontsize=16,color='limegreen')
        ax.text(src2_hdr2['SLTSTRT1']+10,src2_hdr2['SLTSTRT2']+10, "order 2", fontsize=16,color='limegreen')
        ax.text(src3_hdr1['SLTSTRT1']+10,src3_hdr1['SLTSTRT2']+10, "order 1", fontsize=16,color='lightblue')
        ax.text(src3_hdr2['SLTSTRT1']+10,src3_hdr2['SLTSTRT2']+10, "order 2", fontsize=16,color='lightblue')

        ax.set_title("Dispersed image showing source locations and the cutout bounding boxes",fontsize=22)
        leg = ax.legend(fontsize=15,facecolor='black',loc=2)
        for text in leg.get_texts():
            plt.setp(text, color = 'w')
        plt.colorbar(orientation='horizontal',pad=0.05)
        plt.savefig(dispersed_file[:-5]+'_dispersion_cutout_regions.png')

        # extraction for source 1 order 1
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=18)
        plt.xlabel('x pixels',fontsize=18)
        plt.imshow(source1_order1, vmin=0.2, vmax=1.7, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Source 1 (order 1) -- red left box above",fontsize=18)

        # extraction for source 1 order 2
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=18)
        plt.xlabel('x pixels',fontsize=18)
        plt.imshow(source1_order2, vmin=0.2, vmax=1.7, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Source 1 (order 2) -- red right box above",fontsize=18)
        plt.savefig(dispersed_file[:-5]+'_src1order2.png')

        # extraction for source 2 order 1
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=18)
        plt.xlabel('x pixels',fontsize=18)
        plt.imshow(source2_order1, vmin=0.2, vmax=1.7, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Source 2 (order 1) -- green left box above",fontsize=18)
        plt.savefig(dispersed_file[:-5]+'_src2order1.png')

        # extraction for source 2 order 2
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=18)
        plt.xlabel('x pixels',fontsize=18)
        plt.imshow(source2_order2, vmin=0.2, vmax=1.7, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Source 2 (order 2) -- green right box above",fontsize=18)
        plt.savefig(dispersed_file[:-5]+'_src2order2.png')

        # extraction for source 3 order 1
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=18)
        plt.xlabel('x pixels',fontsize=18)
        plt.imshow(source3_order1, vmin=0.2, vmax=1.7, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Source 3 (order 1) -- blue left box above",fontsize=18)
        plt.savefig(dispersed_file[:-5]+'_src3order1.png')

        # extraction for source 3 order 2
        fig, ax = plt.subplots(1,1,figsize=(20,20))
        plt.ylabel('y pixels',fontsize=18)
        plt.xlabel('x pixels',fontsize=18)
        plt.imshow(source3_order2, vmin=0.2, vmax=1.7, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Source 3 (order 2) -- blue right box above",fontsize=18)
        plt.savefig(dispersed_file[:-5]+'_src3order2.png')
