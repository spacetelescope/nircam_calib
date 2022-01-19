#! /usr/bin/env python

"""For a given WFSS image and list of coordinates, calculate the orientation of the trace and report
how far out of alignment with rows/columns it is.

Call the script as such:

python trace_orientation.py coordinates_horiz_obs4visit1.txt --image_file jw01076004001_01101_00013_nrcb5_rate.fits --direction horizontal --n_collapse 20 --cross_disp_width 20
python trace_orientation.py coordinates_vert_obs4visit1.txt --image_file jw01076004001_01101_00009_nrcb5_rate.fits --direction vertical --n_collapse 20 --cross_disp_width 20

The coordinate*txt files should be ascii files that look like this:

# Trace coords
x1    y1    x2    y2
943   780  943   1748
1589  425  1591  1428
1512  475  1512  1004

x1 and y1 in a given row should contain the x and y pixel values associated with the left (or bottom) of a source's trace.
x2 and y2 in a given row should contain the x and y pixel values for the corresponding right (or top) or the trace. Since the
code will stack a certain number of rows/columns, make sure these points are not too close to the edge of the detector.

The direction argument indicates which way the traces fall. 'horizontal' means GRISMR. 'vertical' means GRISMC.

n_collapse is the number of rows/columns to collapse down into a mean 1D profile of the trace. These collapsed
rows/columns will be centered around the input coordinates.

cross_disp_width is the number of rows/columns in the cross dispersed direction to use in order to contain the 1D profile.

For GRISMR data, the code will use 'n_collapse' columns surrounding 'x1' and calculate the median column (keeping only 'cross_disp_width'
pixels, surrounding y1). It will then fit a 1D Gaussian to this median profile, and save the mean value as the y-coordinate
of the trace. It will then repeat this process using (x2, y2). Using these results for the two input coordinates, the angle of the
trace is then calculated as arctan(dy/dx). Coordinates, x/y differences, and angles are saved in a table named:
trace_angles_<name of coord file>_<name of image file>.dat. It also prints the median value of all angles to the screen.

The process is the same for GRISMC data, but with rows/columns switched.
"""
import argparse
from astropy.io import fits, ascii
from astropy.modeling import models, fitting
from astropy.table import Table
import numpy as np
import os


class Source:

    def __init__(self, pt1, arr, horizontal=True, num_collapse=10, cross_disp_width=8):
        """
        Paramters
        ---------
        pt1 : tuple
            2-tuple containing the coordinates of the source in the image

        arr : numpy.ndarray
            The 2d numpy array containing the image

        horizontal : bool
            If True, the source is assumed to be dispersed along rows. If False,
            along columns

        num_collapse : int
            Number of columns (if the trance is horizontal) or rows (if the source is vertical)
            to collapse in order to generate a profile

        cross_disp_width : int
            Number of rows or columns in the cross-dispersion direction to use
        """
        self.pt1 = pt1
        self.array = arr
        self.horizontal = horizontal
        self.num_collapse = num_collapse
        self.cross_disp_width = cross_disp_width

        self.generate_profile()
        self.fit_gaussian()

    def generate_profile(self):
        """Collapse rows/columns to produce a 1D profile of the trace. pt1 will be
        centered within the rows to be collapsed
        """
        #For trace aligned along rows
        if self.horizontal:
            minxpt = int(self.pt1[0] - self.num_collapse/2)
            maxxpt = int(minxpt + self.num_collapse)
            minypt = int(self.pt1[1] - self.cross_disp_width/2)
            maxypt = int(minypt + self.cross_disp_width)
            self.mean_profile = np.median(self.array[minypt:maxypt+1, minxpt:maxxpt+1], axis=1)
            self.coords = np.arange(minypt, maxypt+1)
        else:
            # For a trace aligned along columns
            minypt = int(self.pt1[1] - self.num_collapse/2)
            maxypt = int(minypt + self.num_collapse)
            minxpt = int(self.pt1[0] - self.cross_disp_width/2)
            maxxpt = int(minxpt + self.cross_disp_width)
            self.mean_profile = np.median(self.array[minypt:maxypt+1, minxpt:maxxpt+1], axis=0)
            self.coords = np.arange(minxpt, maxxpt+1)

    def fit_gaussian(self):
        """Fit Gaussian to the mean profile"""
        g_init = models.Gaussian1D(amplitude=np.max(self.mean_profile), mean=np.min(self.coords) + self.cross_disp_width/2, stddev=1.)
        fit_g = fitting.LevMarLSQFitter()
        self.g = fit_g(g_init, self.coords, self.mean_profile)




def define_options(parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument('coords', default=None,
                            help='Filename continaing coordinates to be measured')
        parser.add_argument('--image_file', default=None,
                            help='Fits file containing data')
        parser.add_argument('--direction', default=None,
                            help='Dispersion direction: horizontal or vertical')
        parser.add_argument('--n_collapse', default=None, type=int,
                            help='Number of rows/cols to collapse when making mean profile')
        parser.add_argument('--cross_disp_width', default=None, type=int,
                            help='Number of rows/cols in cross dispersion direction to use')
        return parser


def read_coord_file(filename):
    """Read in ASCII file with coordinates
    """
    tab = ascii.read(filename)
    return tab


if __name__ == '__main__':
    parser = define_options()
    args = parser.parse_args()

    if args.direction.lower() == 'horizontal':
        horiz = True
    elif args.direction.lower() == 'vertical':
        horiz = False
    else:
        raise ValueError(f'Unrecognized direction {args.direction}. Must be "horizontal" or "vertical".')

    # Read in coordinates file
    tab = read_coord_file(args.coords)

    # Read in fits file
    image = fits.getdata(args.image_file)

    src1_list = []
    src2_list = []
    coord2_src1 = []
    coord2_src2 = []
    diff_list = []
    dx_list = []
    angle_list = []
    for i, row in enumerate(tab):
        src1 = Source((row['x1'], row['y1']), image, horizontal=horiz, num_collapse=args.n_collapse, cross_disp_width=args.cross_disp_width)
        coord1_src1 = src1.g.mean.value
        src1_list.append(coord1_src1)
        src2 = Source((row['x2'], row['y2']), image, horizontal=horiz, num_collapse=args.n_collapse, cross_disp_width=args.cross_disp_width)
        coord1_src2 = src2.g.mean.value
        src2_list.append(coord1_src2)

        if horiz:
            coord2_src1.append(row['x1'])
            coord2_src2.append(row['x2'])
        else:
            coord2_src1.append(row['y1'])
            coord2_src2.append(row['y2'])

        diff = coord1_src2 - coord1_src1
        diff_list.append(diff)
        if horiz:
            dx = row['x2'] - row['x1']
        else:
            dx = row['y2'] - row['y1']
        dx_list.append(dx)
        angle = np.arctan2(diff, dx) * 180./np.pi
        angle_list.append(angle)

    median_angle = np.median(angle_list)
    stdev_angle = np.stdev(angle_list)
    print(f'Median and stdev of trance angles relative to rows/cols: {median_angle} +/- {stdev_angle} degrees.')

    if horiz:
        results = Table([coord2_src1, src1_list, coord2_src2, src2_list, diff_list, dx_list, angle_list],
                        names = ('src1_x', 'src1_y', 'src2_x', 'src2_y', 'diff', 'trace_length', 'angle'))
    else:
        results = Table([src1_list, coord2_src1, src2_list, coord2_src2,  diff_list, dx_list, angle_list],
                        names = ('src1_x', 'src1_y', 'src2_x', 'src2_y', 'diff', 'trace_length', 'angle'))

    base_coords = os.path.basename(args.coords)
    base_img = os.path.basename(args.image_file)
    outfile = f'trace_angles_{base_coords}_{base_img}.dat'
    ascii.write(results, outfile, overwrite=True)
    print("Saved results to: {outfile}")


