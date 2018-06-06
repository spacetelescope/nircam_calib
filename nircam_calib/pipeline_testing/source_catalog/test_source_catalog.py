import copy,os,sys
import pytest
import numpy as np
from glob import glob
import json
from astropy.io import fits, ascii
from astropy.coordinates import Angle
from astropy.table import Table, vstack, unique
from jwst.pipeline import calwebb_image3
import matplotlib.pyplot as plt

def get_input_table(num_source_files, sourcelists):
    '''Function to read in and access the simulator source input files.'''

    all_source_table = Table()

    # point source and galaxy source tables have different headers
    # change column headers to match for filtering later
    for i in np.arange(0, num_source_files):
        if "point" in sourcelists[i]:
            col_names = ["RA", "Dec", "RA_degrees", "Dec_degrees",
                         "PixelX", "PixelY", "Magnitude",
                         "counts_sec", "counts_frame"]
        elif "galaxy" in sourcelists[i]:
            col_names = ["PixelX", "PixelY", "RA", "Dec",
                         "RA_degrees", "Dec_degrees", "V2", "V3", "radius",
                         "ellipticity", "pos_angle", "sersic_index",
                         "Magnitude", "countrate_e_s", "counts_per_frame_e"]
        else:
            print('Error! Source list column names need to be defined.')
            sys.exit(0)

        # read in the tables
        input_source_table = ascii.read(sourcelists[i])
        orig_colnames = input_source_table.colnames

        for col, n in zip(orig_colnames, np.arange(0,len(orig_colnames))):
            input_source_table[col].name = col_names[n]

        # only grab values for source catalog analysis
        short_source_table = Table({'In_RA': input_source_table['RA_degrees'],
                                  'In_Dec': input_source_table['Dec_degrees'],
                                  'In_PixelX': input_source_table['PixelX'],
                                  'In_PixelY': input_source_table['PixelY'],
                                  'In_Magnitude': input_source_table['Magnitude']},
                                  names=['In_RA', 'In_Dec', 'In_PixelX', 'In_PixelY', 'In_Magnitude'])

        # combine source lists into one master list
        all_source_table = vstack([all_source_table, short_source_table])

    # set up columns to track which sources were detected by Photutils
    # mask off-detector sources, these should be accounted for in other images
    all_source_table['Detected'] = 'N'
    for row in np.arange(0,len(all_source_table)):
        pos = np.asfarray([all_source_table['In_PixelX'][row], all_source_table['In_PixelY'][row]])
        if np.any(pos < 0.) == True:
            all_source_table['Detected'][row] = 'O'
        if np.any(pos >= 2048.) == True:
            all_source_table['Detected'][row] = 'O'

    out_of_field = all_source_table['Detected'] == 'O'
    all_source_table = all_source_table[~out_of_field]
    all_source_table['Out_RA'] = np.nan
    all_source_table['Out_Dec'] = np.nan
    all_source_table['Out_mag'] = np.nan
    all_source_table['Src_sum'] = np.nan

    # filter by RA, Dec (for now)
    no_duplicates = unique(all_source_table,keys=['In_RA','In_Dec'])

    return no_duplicates


def  get_output_table(input_name_base):
    '''Function to read in and access the output source catalog.'''

    # use association file to read in output catalog from Image3 pipeline
    print('\nWorking on output table: '+input_name_base+'_cat.ecsv')
    cat_file = glob(input_name_base+'*cat.ecsv')[0]
    output_source_table = ascii.read(cat_file)
    try:
        short_output_table = Table({'RA': output_source_table['sky_centroid'].ra,
                                  'Dec': output_source_table['sky_centroid'].dec,
                                  'PixelX': output_source_table['xcentroid'],
                                  'PixelY': output_source_table['ycentroid'],
                                  'Magnitude': output_source_table['abmag'],
                                  'Source_sum': output_source_table['source_sum']},
                                  names=['RA', 'Dec', 'PixelX', 'PixelY', 'Magnitude', 'Source_sum'])
    except:
        short_output_table = Table({'RA': output_source_table['sky_centroid.ra'],
                                  'Dec': output_source_table['sky_centroid.dec'],
                                  'PixelX': output_source_table['xcentroid'],
                                  'PixelY': output_source_table['ycentroid'],
                                  'Magnitude': output_source_table['abmag'],
                                  'Source_sum': output_source_table['source_sum']},
                                  names=['RA', 'Dec', 'PixelX', 'PixelY', 'Magnitude', 'Source_sum'])

    return short_output_table


def generate_figure(image, x, y, vmin, vmax, title, filename):
    '''Function to generate images.'''


    # first look at input source positions overlayed on the Image2 calibrated exposure
    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    plt.ylabel('y pixels', fontsize=22)
    plt.xlabel('x pixels', fontsize=22)
    plt.imshow(image, vmin=vmin, vmax= vmax, cmap=plt.cm.gray, origin='lower')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.set_title(title, fontsize=22)
    for i in np.arange(0, len(x)):
        circle = plt.Circle((x[i], y[i]), 15., facecolor='None', edgecolor='red', linewidth=3)
        ax.add_artist(circle)
    plt.colorbar(orientation='horizontal', pad=0.05)
    plt.savefig(filename)


def generate_scatter(RA_in, Dec_in, RA_out, Dec_out, title, filename):
    '''Function to generate scatter plots of input and output sources.'''

    # look at input source and output catalog source RA and Dec
    fig, ax = plt.subplots(1,1,figsize=(20,20))
    plt.ylabel('Dec',fontsize=22)
    plt.xlabel('RA',fontsize=22)
    ax.set_facecolor('Black')

    ax.scatter(RA_in,Dec_in, s=170,marker='s',facecolors='none', edgecolors='lightblue',lw=5, label='Input sources')
    ax.scatter(RA_out,Dec_out, s=120,facecolors='none', edgecolors='r',lw=5, label='Source catalog')

    ax.set_title(title,fontsize=22)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    leg = ax.legend(loc='best',facecolor='None',edgecolor='black',fontsize=20)
    for text in leg.get_texts():
        plt.setp(text, color = 'w')
    plt.savefig(filename)
    plt.close()


def test_sources_detected(cases, RAtol, Dectol, kernel_fwhm, kernel_xsize, kernel_ysize, npixels, snr_threshold):
    '''Test how many sources were detected.'''

    # save source plots?
    save_figs = True

    # get association file and source file names
    # create complete list of input source files
    sourcelists = []
    with open(cases) as json_data:
         d = json.load(json_data)
    members = d['products'][0]['members']
    for row in np.arange(0,len(members)):
        idx = members[row]['expname'].find("W_")+1
        input_name_base = members[row]['expname'][:idx]
        s1 = glob(input_name_base+'*sources.list')
        s2 = glob(input_name_base+'*Sources.list')
        sourcelists = sourcelists + s1 + s2
    num_source_files = len(sourcelists)

    # get master table of input sources
    all_source_table = get_input_table(num_source_files,sourcelists)

    # run Image3 pipeline to get source catalog
    im3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image3.cfg')
    im3.tweakreg.skip = True
    im3.source_catalog.kernel_fwhm = kernel_fwhm
    im3.source_catalog.kernel_xsize = kernel_xsize
    im3.source_catalog.kernel_ysize = kernel_ysize
    im3.source_catalog.npixels = npixels
    im3.source_catalog.snr_threshold = snr_threshold
    im3.resample.blendheaders = False
    im3.run(cases)

    # read in output source catalog
    output_name_base = d['products'][0]['name']
    if '_cal' in output_name_base:
        cal_idx = output_name_base.find('_cal')
        output_name_base = output_name_base[:cal_idx]
    short_output_table = get_output_table(output_name_base)

    # load output data arrays
    final_output = output_name_base+'_i2d.fits'
    image3 = fits.getdata(final_output,1)

    # flag input sources with RA,Dec close to output source RA, Decs
    for i in np.arange(0,len(all_source_table['In_RA'])):
        for j in np.arange(0,len(short_output_table['RA'])):
            RA_diff = np.absolute(all_source_table['In_RA'][i]-short_output_table['RA'][j])
            Dec_diff = np.absolute(all_source_table['In_Dec'][i]-short_output_table['Dec'][j])
            if RA_diff < RAtol and Dec_diff < Dectol:
                all_source_table['Detected'][i] = 'Y'
                all_source_table['Out_RA'][i] = short_output_table['RA'][j]
                all_source_table['Out_Dec'][i] = short_output_table['Dec'][j]
                all_source_table['Out_mag'][i] = short_output_table['Magnitude'][j]
                all_source_table['Src_sum'][i] = short_output_table['Source_sum'][j]

    # count number of sources found by Photutils
    total_found = (all_source_table['Detected']=='Y').sum()
    total_missed = (all_source_table['Detected']=='N').sum()
    total_percent_found = (total_found/len(all_source_table))*100
    mask = (all_source_table['Detected']=='Y')

    print('\n')
    print('Total found:',total_found)
    print('Total number of sources:',len(all_source_table))
    print('Total percent found:',total_percent_found)
    print('Number of false detections: ',len(short_output_table)-len(all_source_table))
    print('\n')

    out_ext = '_RA'+str(RAtol)+'Dec'+str(Dectol)+'_kern'+str(kernel_fwhm)+str(kernel_xsize)+str(kernel_ysize)+'_npix'+str(npixels)+'_snr'+str(snr_threshold)
    ascii.write(all_source_table, output_name_base+'_input_comparison'+out_ext+'.dat',format='fixed_width_two_line')

    if save_figs == True:

        generate_scatter(all_source_table[mask]['In_RA'],all_source_table[mask]['In_Dec'],all_source_table[mask]['Out_RA'],all_source_table[mask]['Out_Dec'],
                         'Input sources and corresponding catalog sources (SNR = '+str(snr_threshold)+')',
                         output_name_base+'_detected_sources'+out_ext+'.png')
        generate_scatter(all_source_table['In_RA'],all_source_table['In_Dec'],short_output_table['RA'],short_output_table['Dec'],
                         'Input source RA,Dec overlayed with source catalog RA,Dec (SNR = '+str(snr_threshold)+')',
                         output_name_base+'_input_output_RADec_overplot'+out_ext+'.png')
        generate_figure(image3, short_output_table['PixelX'], short_output_table['PixelY'],np.nanmean(image3)-0.08,np.nanmean(image3)+0.001,
                        'Output image with source catalog overlayed (SNR = '+str(snr_threshold)+')',
                        output_name_base+'_output_img_source_cat'+out_ext+'.png')

    assert (total_percent_found > 90.0) == True




# #pytest test_source_catalog.py --html=report.html --self-contained-html
