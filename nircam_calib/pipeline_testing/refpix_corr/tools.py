#! /usr/bin/env python

'''
This is a list of functions used in the main pipeline test scripts.
'''

def extract(array,xstart,xend,ystart,yend):
    '''Extract the described subarray from a full frame array.'''

    # (Designed to be used for extracting subarrays from reference files for
    # manual comparisons to pipeline outputs
    dimensions = array.shape
    if len(dimensions) == 2:
        return array[ystart:yend,xstart:xend]
    elif len(dimensions) == 3:
        return array[:,ystart:yend,xstart:xend]
    elif len(dimensions) == 4:
        return array[:,:,ystart:yend,xstart:xend]


def get_coords(file):
    '''Get coords of subarray to be extracted.'''

    # Subtract 1 because coords in the header are indexed to 1
    # but python indexes to zero
    with fits.open(file) as h:
        xstart = h[0].header['SUBSTRT1'] - 1
        ystart = h[0].header['SUBSTRT2'] - 1
        xend = xstart + h[0].header['SUBSIZE1']
        yend = ystart + h[0].header['SUBSIZE2']
    return xstart,xend,ystart,yend


def get_coords_rampmodel(model):
    '''Get coords of subarray from a rampmodel.'''

    # Subtract 1 because coords in the header are indexed to 1
    # but python indexes to zero
    xstart = model.meta.subarray.xstart - 1
    ystart = model.meta.subarray.ystart - 1
    xend = xstart + model.meta.subarray.xsize
    yend = ystart + model.meta.subarray.ysize
    return xstart,xend,ystart,yend


def display_img(image,title, vmin, vmax):
    '''Display image of data (need 2D array).'''

    import matplotlib.pyplot as plt

    plt.figure(figsize=(20,20))
    plt.ylabel('y pixels',fontsize=22)
    plt.xlabel('x pixels',fontsize=22)
    plt.title('\n\n'+title+'\n',fontsize=22)
    plt.imshow(image, vmin = vmin, vmax=vmax, cmap=plt.cm.gray, origin='lower')
    plt.colorbar(orientation='horizontal',pad=0.05)
    return plt


def save_df_table(table,title):
    '''Save out table of results to compare groups and amps.'''

    import pandas as pd
    from astropy.table import Table
    from astropy.io import ascii

    columns=['Group','Amp','Data_Mean','Refpix_Mean']
    df = pd.DataFrame(table, columns=columns)
    table = Table.from_pandas(df)
    ascii.write(table,title,format='fixed_width_two_line',overwrite=True)
