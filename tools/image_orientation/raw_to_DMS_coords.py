#! /usr/bin/env python

'''
Translate detector coordinates from those in the raw orientation to 
DMS orientation. Assumption is that raw coordinates are 1-indexed, while
DMS coordinates are 0-indexed. Also assume that coordinates are full-frame
coordinates.
'''

import argparse,sys
from astropy.io import ascii

hflip = ['A1','A3','B2','B4','A5']
vflip = ['A2','A4','B1','B3','B5']
maxindex = {'A1':2048,'A2':2048,'A3':2048,'A4':2048,'A5':2048,'B1':2048,'B2':2048,'B3':2048,'B4':2048,'B5':2048}

class rawToDMS:
    def __init__(self):
        self.verbose = False

    def run(self):

        #read in ascii file of raw subarray coordinates
        tab = ascii.read(self.infile)

        #add columns for xend and yend
        tab['xend'] = 0
        tab['yend'] = 0

        #One by one, translate the subarray coords
        for i in xrange(len(tab)):
            
            #translate any subarray entry, skip FULL 
            if (tab['Name'][i] != 'FULL'):
                #also skip any entries with 0 or lower in coordinate entries
                if ((tab['xstart'][i] > 0) & (tab['ystart'][i] > 0)):

                    #first get the detector
                    det = tab['Name'][i][3:5]

                    #calcualte the ending coordinates of the subarray
                    tab['xend'][i] = tab['xstart'][i] + tab['xlen'][i] - 1
                    tab['yend'][i] = tab['ystart'][i] + tab['ylen'][i] - 1

                    #translate
                    newx1,newy1 = self.translate(det,tab['xstart'].data[i],tab['ystart'].data[i])
                    newx2,newy2 = self.translate(det,tab['xend'].data[i],tab['yend'].data[i])
                
                    newx = [newx1,newx2]
                    newy = [newy1,newy2]

                    #insert the translated values back into the table
                    tab['xstart'][i] = min(newx)
                    tab['ystart'][i] = min(newy)
                    tab['xend'][i] = max(newx)
                    tab['yend'][i] = max(newy)

            else:
                #reduce the coords in the FULL entry by 1, to get into 0-indexed system
                tab['xstart'][i] = tab['xstart'][i] - 1
                tab['xend'][i] = tab['xstart'][i] + tab['xlen'][i] - 1
                tab['ystart'][i] = tab['ystart'][i] - 1
                tab['yend'][i] = tab['ystart'][i] + tab['ylen'][i] - 1
                

        #remove the x and y length columns
        tab.remove_columns(['xlen', 'ylen'])

        #save the updated table
        if self.outfile == None:
            self.outfile = self.infile + '_translated_to_DMS.tab'
        ascii.write(tab,self.outfile)




    def translate(self,detector,x,y):

        print('within raw_to_DMS_coords, detector is {}'.format(detector))
        
        #find whether we want a vertical or horizontal flip
        flip_dir = self.findFlipDir(detector)
        
        #what is the maximum index value (2048 for SW, 1024 for LW)
        mxval = maxindex[detector]
        
        #subtract 1 from the coordinate in the non-flipping
        #direction in order to go from 1-indexed to 0-indexed
        #Similarly, subtract x or y from mxval, rather than
        #(mxval+1) in order to go to 0-indexed coords.
        if flip_dir == 'horizontal':
            newx = mxval - x
            newy = y - 1

        if flip_dir == 'vertical':
            newx = x - 1
            newy = mxval - y

        return newx,newy


    def findFlipDir(self,detector):
        if detector in hflip:
            return 'horizontal'
        if detector in vflip:
            return 'vertical'
        print('Error: unrecognized detector')
        sys.exit()


    def add_options(self,parser=None,usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,description='Translate raw to DMS coordinates for NIRCam.')
        parser.add_argument("infile",help="Name of table file with coordinates to translate")
        parser.add_argument("--outfile",help="Name of file to save table of translated values.",default=None)
        return parser
        

if __name__ == '__main__':

    usagestring = 'raw_to_DMS_coords.py inputfile.tab --outfile translated.tab'

    translate = rawToDMS()
    parser = translate.add_options(usage=usagestring)
    args = parser.parse_args(namespace=translate)
    translate.run()

    
