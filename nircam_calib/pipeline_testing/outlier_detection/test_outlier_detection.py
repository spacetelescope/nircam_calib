import copy,os,sys
import pytest
import numpy as np
from glob import glob
import json
from astropy.io import fits, ascii
from astropy.coordinates import Angle
from astropy.table import Table, vstack, unique
from astropy.stats import sigma_clip
from jwst.pipeline import calwebb_image3
import matplotlib.pyplot as plt



def test_median(median):
    '''Test median image.'''

    save_figs = True

    # create fake data for subpixel dither 1 to test median
    with fits.open("V54321001002P0000000001101_A1_F150W_cal.fits") as h:
        h['SCI'].data[:,:] = 3.0
        h['DQ'].data[:,:] = 0
        h.writeto("V54321001002P0000000001101_A1_F150W_cal_mediantest.fits",overwrite=True)


    # create fake data for subpixel dither 2 to test median
    with fits.open("V54321001002P0000000001102_A1_F150W_cal.fits") as h:
        h['SCI'].data[:,:] = 5.0
        h['DQ'].data[:,:] = 0
        h.writeto("V54321001002P0000000001102_A1_F150W_cal_mediantest.fits",overwrite=True)


    # run Image3 pipeline to get outlier detection outputs
    im3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image3.cfg')
    im3.tweakreg.skip = True
    im3.resample.blendheaders = False
    im3.outlier_detection.save_intermediate_results = True
    im3.run(median)

    # get file names and bases to load output data
    input_files = []
    with open(median) as json_data:
         d = json.load(json_data)
    members = d['products'][0]['members']
    base = members[0]['expname'][1:8]
    output_files = glob("jw"+base+"*01101_00001_outlier_i2d.fits") + \
                   glob("jw"+base+"*01102_00001_outlier_i2d.fits")

    input_file_base = members[0]['expname']

    # put results into array
    with fits.open(output_files[0]) as h:
        shape = h[1].data.shape
    all_dithers = np.zeros((len(output_files),shape[0],shape[1]),dtype='float32')
    for i in np.arange(0,len(all_dithers)):
        print(np.shape(fits.getdata(output_files[i],1)))
        all_dithers[i,:,:] = fits.getdata(output_files[i],1)

    # get output median image
    first = input_file_base
    cal = first.find("_F")
    median_file = first[:cal]+"_median.fits"
    median_image = fits.getdata(median_file,1)

    # manually calculate the median
    calc_median_array = np.nanmedian(all_dithers,axis=0)
    diff_array = np.abs(calc_median_array - median_image)

    # sigma-clipping
    clip = sigma_clip(diff_array)
    clip.data[clip.mask] = np.nan
    diff_mean = np.nanmean(clip.data)

    if save_figs == True:

        # save figure to show median
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        plt.ylabel('y pixels',fontsize=15)
        plt.xlabel('x pixels',fontsize=15)
        plt.imshow(median_image, vmin=3.0, vmax=5.0, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Median image",fontsize=12)
        plt.colorbar(orientation='horizontal',pad=0.09)
        plt.savefig(median[:5]+"_median_image.png")

        # save figure to show difference between calculated and output medians
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        plt.ylabel('y pixels',fontsize=15)
        plt.xlabel('x pixels',fontsize=15)
        plt.imshow(diff_array, vmin=-0.5, vmax=0.5, cmap=plt.cm.gray, origin='lower')
        ax.set_title("Difference between median image and manually calculated median",fontsize=12)
        plt.colorbar(orientation='horizontal',pad=0.09)
        plt.savefig(median[:5]+"_median_diff_image.png")

    assert np.isclose(diff_mean,0.0) == True


def test_outlier(outlier):
    '''Test outlier detection.'''

    save_figs = True

    # create fake outliers for subpixel dither 1 to test outlier detection
    pix1loc = [300,300]
    pix2loc = [700,700]
    pix3loc = [1500,1500]
    pix4loc = [2000,2000]

    with fits.open("V54321001002P0000000001101_A1_F150W_cal.fits") as h:
        h['SCI'].data[pix1loc[0],pix1loc[1]] = 5.0
        h['SCI'].data[pix2loc[0],pix2loc[1]] = 10.0
        h['SCI'].data[pix3loc[0],pix3loc[1]] = 25.0
        h['SCI'].data[pix4loc[0],pix4loc[1]] = 50.0
        h['DQ'].data[:,:] = 0
        h.writeto("V54321001002P0000000001101_A1_F150W_cal_outliertest.fits",overwrite=True)


    # create fake outliers for subpixel dither 2 to test outlier detection
    pix5loc = [1300,1300]
    pix6loc = [1700,1700]
    pix7loc = [600,600]
    pix8loc = [50,50]

    with fits.open("V54321001002P0000000001102_A1_F150W_cal.fits") as h:
        h['SCI'].data[pix5loc[0],pix5loc[1]] = 5.0
        h['SCI'].data[pix6loc[0],pix6loc[1]] = 10.0
        h['SCI'].data[pix7loc[0],pix7loc[1]] = 25.0
        h['SCI'].data[pix8loc[0],pix8loc[1]] = 50.0
        h['DQ'].data[:,:] = 0
        h.writeto("V54321001002P0000000001102_A1_F150W_cal_outliertest.fits",overwrite=True)

    # run Image3 pipeline
    im3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image3.cfg')
    im3.tweakreg.skip = True
    im3.resample.blendheaders = False
    im3.outlier_detection.save_intermediate_results = True
    im3.run(outlier)

    # get filenames and outlier detection outputs
    output_files = []
    input_files = []
    with open(outlier) as json_data:
         d = json.load(json_data)
    members = d['products'][0]['members']
    for item in np.arange(0,len(members)):
        input_files.append(members[item]['expname'])
        expname = members[item]['expname'][:-5]+"_a3001_crf.fits"
        output_files.append(expname)
    output_files.sort()

    all_out_dqs = []
    dq_before_dith1 = []
    dq_after_dith1 = []

    # get pixel DQ values before (should be 0.0) for dither 1
    with fits.open(input_files[0]) as h:
        pixels = [pix1loc,pix2loc,pix3loc,pix4loc]
        for pix in pixels:
            dq_before_dith1.append([pix,h['DQ'].data[pix[0],pix[1]]])

            if save_figs == True:

                # save figure of input dq vals
                fig, ax = plt.subplots(1,1,figsize=(10,10))
                plt.ylabel('y pixels',fontsize=15)
                plt.xlabel('x pixels',fontsize=15)
                plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
                ax.set_title("Dither 1 input DQ",fontsize=15)
                plt.colorbar(orientation='horizontal',pad=0.09)
                plt.savefig(outlier[:5]+"_dither1_inputDQ.png")

    # get pixel DQ values after (should be 4.0) for dither 1
    with fits.open(output_files[0]) as h:
        pixels = [pix1loc,pix2loc,pix3loc,pix4loc]
        for pix in pixels:
            dq_after_dith1.append([pix,h['DQ'].data[pix[0],pix[1]]])
            all_out_dqs.append((h['DQ'].data[pix[0],pix[1]] == 4.0))

            if save_figs == True:

                # save figure of output dq vals
                fig, ax = plt.subplots(1,1,figsize=(10,10))
                plt.ylabel('y pixels',fontsize=15)
                plt.xlabel('x pixels',fontsize=15)
                plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
                ax.set_title("Dither 1 output DQ",fontsize=15)
                plt.colorbar(orientation='horizontal',pad=0.09)
                plt.savefig(outlier[:5]+"_dither1_outputDQ.png")

    dq_before_dith2 = []
    dq_after_dith2 = []

    # get pixel DQ values before (should be 4.0) for dither 2
    with fits.open(input_files[1]) as h:
        pixels = [pix5loc,pix6loc,pix7loc,pix8loc]
        for pix in pixels:
            dq_before_dith2.append([pix,h['DQ'].data[pix[0],pix[1]]])

            if save_figs == True:

                # save figure of output dq vals
                fig, ax = plt.subplots(1,1,figsize=(10,10))
                plt.ylabel('y pixels',fontsize=15)
                plt.xlabel('x pixels',fontsize=15)
                plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
                ax.set_title("Dither 2 input DQ",fontsize=15)
                plt.colorbar(orientation='horizontal',pad=0.09)
                plt.savefig(outlier[:5]+"_dither2_inputDQ.png")

    # get pixel DQ values after (should be 4.0) for dither 2
    with fits.open(output_files[1]) as h:
        pixels = [pix5loc,pix6loc,pix7loc,pix8loc]
        for pix in pixels:
            dq_after_dith2.append([pix,h['DQ'].data[pix[0],pix[1]]])
            all_out_dqs.append((h['DQ'].data[pix[0],pix[1]] == 4.0))

            if save_figs == True:

                # save figure of output dq vals
                fig, ax = plt.subplots(1,1,figsize=(10,10))
                plt.ylabel('y pixels',fontsize=15)
                plt.xlabel('x pixels',fontsize=15)
                plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
                ax.set_title("Dither 2 output DQ",fontsize=15)
                plt.colorbar(orientation='horizontal',pad=0.09)
                plt.savefig(outlier[:5]+"_dither2_outputDQ.png")

    print('Output DQ values: ',all_out_dqs)

    assert np.alltrue(all_out_dqs) == True


def test_image(cases):
    '''Test outlier detection on simulated images.'''

    save_figs = True

    # run Image3 pipeline on regular image file to make sure it works
    im3 = calwebb_image3.Image3Pipeline(config_file='calwebb_image3.cfg')
    im3.tweakreg.skip = True
    im3.resample.blendheaders = False
    im3.outlier_detection.save_intermediate_results = True
    im3.run(cases)

    # get filenames and outputs
    output_files = []
    input_files = []
    with open(cases) as json_data:
         d = json.load(json_data)
    members = d['products'][0]['members']
    for item in np.arange(0,len(members)):
        input_files.append(members[item]['expname'])
        cal = members[item]['expname'].find("_cal")
        expname = members[item]['expname'][:cal]+"_a3001_crf.fits"
        output_files.append(expname)
    output_files.sort()

    # count number of flagged pixels before for dither 1
    with fits.open(input_files[0]) as h:
        dith1_before_counts = np.count_nonzero(h['DQ'].data == 4.0)

        if save_figs == True:

            # save figure of output dq vals
            fig, ax = plt.subplots(1,1,figsize=(10,10))
            plt.ylabel('y pixels',fontsize=15)
            plt.xlabel('x pixels',fontsize=15)
            plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
            ax.set_title("Dither 1 input DQ",fontsize=15)
            plt.colorbar(orientation='horizontal',pad=0.09)
            plt.savefig(cases[:5]+"_sim_dither1_inputDQ.png")

    # count number of flagged pixels after for dither 1
    with fits.open(output_files[0]) as h:
        dith1_after_counts = np.count_nonzero(h['DQ'].data == 4.0)

        if save_figs == True:

            # save figure of output dq vals
            fig, ax = plt.subplots(1,1,figsize=(10,10))
            plt.ylabel('y pixels',fontsize=15)
            plt.xlabel('x pixels',fontsize=15)
            plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
            ax.set_title("Dither 1 output DQ",fontsize=15)
            plt.colorbar(orientation='horizontal',pad=0.09)
            plt.savefig(cases[:5]+"_sim_dither1_outputDQ.png")

    # count number of flagged pixels before for dither 2
    with fits.open(input_files[1]) as h:
        dith2_before_counts = np.count_nonzero(h['DQ'].data == 4.0)

        if save_figs == True:

            # save figure of output dq vals
            fig, ax = plt.subplots(1,1,figsize=(10,10))
            plt.ylabel('y pixels',fontsize=15)
            plt.xlabel('x pixels',fontsize=15)
            plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
            ax.set_title("Dither 2 input DQ",fontsize=15)
            plt.colorbar(orientation='horizontal',pad=0.09)
            plt.savefig(cases[:5]+"_sim_dither2_inputDQ.png")

    # count number of flagged pixels after for dither 2
    with fits.open(output_files[1]) as h:
        dith2_after_counts = np.count_nonzero(h['DQ'].data == 4.0)

        if save_figs == True:

            # save figure of output dq vals
            fig, ax = plt.subplots(1,1,figsize=(10,10))
            plt.ylabel('y pixels',fontsize=15)
            plt.xlabel('x pixels',fontsize=15)
            plt.imshow((h['DQ'].data == 4.0), vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
            ax.set_title("Dither 2 output DQ",fontsize=15)
            plt.colorbar(orientation='horizontal',pad=0.09)
            plt.savefig(cases[:5]+"_sim_dither2_outputDQ.png")

    before_counts = dith1_before_counts + dith2_before_counts
    after_counts = dith1_after_counts + dith2_after_counts

    assert (before_counts < after_counts) == True

# pytest test_outlier_detection.py --html=report.html --self-contained-html
