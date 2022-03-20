#! /usr/bin/env python

"""This module is designed to search for snowballs in JWST data

Add more details later

High level method description:

Input: list of jump.fits filenames
1. For a given file output from the jump step:
   A. Extract the DQ array, strip away all flags other than JUMP
   B. Transform the DQ array into a 2D map for each integration, where the pixel value is the
      group number within which a JUMP flag is present. In the case of multiple CR hits in a single
      pixel, just keep the first group number
   C. Run photutils source detection on the DQ map. Only return sources with more than (?20?) adjaent
      pixels, (and impose a roundness constraint as well?).
   D. For each found source:
      i. Return a cutout of the CDS image from the appropriate integration/group in the data extension.
      ii. Run photometry in order to measure the total signal added by the snowball
      iii. Find and keep track of the abosolute time of the group
2. Calculate the rate of snowball appearances (total #snowballs / total exposure time)
3. Stitch together the CRs into a gallery for easy viewing
4. Create histogram of total signal added
5. Create a summary table that lists the number of snowballs per file
"""
import os
from photutils.segmentation import detect_sources, detect_threshold, SourceCatalog
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from jwst.datamodels import dqflags
import matplotlib.pyplot as plt

class SnowballDetector():

    def __init__(self, filenames, output_dir='./', histogram_filebase='snowball_histogram',
                 gallery_filename='snowball_gallery.png', summary_file='snowball_summary.ecsv',
                 snowball_table='snowball_table.ecsv', min_snowball_area=100, max_ellipticity=0.9,
                 min_total_signal=100000., save_segment_map=False, make_gallery_view=True, log_histograms=True):

        self.filenames = filenames
        self.output_dir = output_dir
        self.histogram_filebase = histogram_filebase
        self.summary_file = summary_file
        self.snowball_table = snowball_table
        self.snowball_info = None
        self.gallery_filename = gallery_filename
        self.min_snowball_area = min_snowball_area
        self.max_ellipticity = max_ellipticity
        self.min_total_signal = min_total_signal
        self.save_segment_map = save_segment_map
        self.make_gallery_view = make_gallery_view
        self.log_histograms = log_histograms
        self.run()

    def run(self):
        self.cutouts = []
        self.gallery_info = Table(names=('filename', 'integration', 'group', 'x', 'y'), dtype=('S80', 'i4', 'i4', 'f4', 'f4'))
        self.total_exposure_time = 0.
        self.number_of_snowballs = 0
        self._snowballs_per_exp = []
        self._median_intensity = []
        self._median_area = []
        self._exposure_time = []
        for i, filename in enumerate(self.filenames):
            print(f"Working on {filename}. File {i+1} out of {len(self.filenames)}")
            sb = ExpSnowballs(filename, min_snowball_area=self.min_snowball_area, max_ellipticity=self.max_ellipticity,
                              min_total_signal=self.min_total_signal, create_segment_map=self.save_segment_map)

            if self.save_segment_map:
                seg_filename = os.path.join(self.output_dir, f'segmap_{os.path.basename(self.filename)}.fits')
                h0 = fits.PrimaryHDU(sb.segmap)
                hl = fits.HDUList([h0])
                hl.writeto(seg_filename, overwrite=True)

            self._exposure_time.append(sb.exposure_time)
            if sb.source_catalog is not None:
                if self.snowball_info is None:
                    self.snowball_info = sb.source_catalog
                else:
                    self.snowball_info = vstack([self.snowball_info, sb.source_catalog])

                # Add the file's exposure time to the total
                self.total_exposure_time += sb.exposure_time

                # Add the number of snowballs to the total
                self.number_of_snowballs += sb.num_snowballs

                # File-specific results, to be placed in the final summary table
                self._snowballs_per_exp.append(sb.num_snowballs)
                med_intensity = np.median(sb.source_catalog["segment_flux"].data - sb.source_catalog["local_background"].data)
                med_area = np.median(sb.source_catalog["area"].data)
                self._median_intensity.append(med_intensity)
                self._median_area.append(med_area)

                # Add the cutouts to the master list
                # Impose a maximum number of cutouts in order to
                # prevent the later gallery image from becoming too
                # large
                if len(self.cutouts) < 400:
                    self.cutouts.extend(sb.cutouts)

                    # Keep track of cutout information
                    self.gallery_info = vstack([self.gallery_info, sb.cutout_info])
            else:
                self._median_intensity.append(0.)
                self._median_area.append(0.)
                self._snowballs_per_exp.append(0)

        if self.snowball_info is not None:
            self.snowball_info.write(self.snowball_table, overwrite=True)

            # Calculate the overall snowball rate
            self.snowball_rate = self.number_of_snowballs / self.total_exposure_time

            # Make histogram of snowball intensity
            self.intensity_histogram()

            # Make a histogram of snowball areas
            self.area_histogram()

            # Make a histogram of snowball ellipticities
            self.ellipticity_histogram()

            # Make a gallery that shows all cutouts
            if self.make_gallery_view:
                self.make_gallery_image()
        else:
            print("No snowballs found.")
            self.snowball_rate = 0.

        # Summary table
        self.create_summary_table()


    def area_histogram(self):
        n_area, bins_area, patches_area = plt.hist(self.snowball_info["area"].value, 50, facecolor='blue', alpha=0.5)
        plt.xlabel('Snowball Area (pixels^2)')
        plt.ylabel('Number of snowballs')
        if self.log_histograms:
            plt.gca().set_yscale("log")
        plt.gca().ticklabel_format(axis='x', style='plain')
        plt.savefig(os.path.join(self.output_dir, f'{self.histogram_filebase}_area.png'))
        plt.close()


    def ellipticity_histogram(self):
        n_area, bins_area, patches_area = plt.hist(self.snowball_info["ellipticity"].value, 50, facecolor='blue', alpha=0.5)
        plt.xlabel('Snowball Ellipticity')
        plt.ylabel('Number of snowballs')
        plt.gca().ticklabel_format(axis='x', style='plain')
        plt.savefig(os.path.join(self.output_dir, f'{self.histogram_filebase}_ellipticity.png'))
        plt.close()


    def create_summary_table(self):
        self.table = Table()
        self.table["filename"] = self.filenames
        self.table["exposure_time"] = self._exposure_time
        self.table["number_of_snowballs"] = self._snowballs_per_exp
        self.table["median_signal"] = self._median_intensity
        self.table["median_area"] = self._median_area
        comment = f"Overall snowball rate: {self.snowball_rate} snowballs per detector per second"
        if 'comments' in self.table.meta:
            self.table.meta['comments'].append(comment)
        else:
            self.table.meta['comments'] = comment
        self.table.write(os.path.join(self.output_dir, self.summary_file), overwrite=True)


    def intensity_histogram(self):
        intensities = self.snowball_info["segment_flux"] - self.snowball_info["local_background"]
        n, bins, patches = plt.hist(intensities, 50, facecolor='blue', alpha=0.5)
        plt.xlabel('Signal (DN)')
        plt.ylabel('Number of snowballs')
        if self.log_histograms:
            plt.gca().set_yscale("log")
        plt.gca().ticklabel_format(axis='x', style='plain')
        plt.savefig(os.path.join(self.output_dir, f'{self.histogram_filebase}_signal.png'))
        plt.close()


    def make_gallery_image(self):
        """Create a gallery of cutout images. If cutouts have inconsistent
        dimensions, padding will be added to all of those that are smaller
        than the largest, in order to form a grid.
        """
        num = len(self.cutouts)
        ncols = int(np.sqrt(num))
        nrows = int(num / ncols)
        extras = num % nrows
        if extras > 0:
            nrows += 1

        ysizes = [cutout.shape[0] for cutout in self.cutouts]
        xsizes = [cutout.shape[1] for cutout in self.cutouts]
        maxy = np.max(ysizes)
        maxx = np.max(xsizes)

        self.gallery = np.zeros((nrows*maxy, ncols*maxx))
        i = 0
        j = 0
        ilist = []
        jlist = []
        minxlist = []
        minylist = []
        for cutout in self.cutouts:
            ydim, xdim = cutout.shape
            xpad = int((maxx - xdim) / 2)
            ypad = int((maxy - ydim) / 2)
            miny = (j * maxy) + ypad
            minx = (i * maxx) + xpad
            self.gallery[miny:miny+ydim, minx:minx+xdim] = cutout
            ilist.append(i)
            jlist.append(j)
            minxlist.append(minx)
            minylist.append(miny)
            if i < (ncols-1):
                i += 1
            else:
                i = 0
                j += 1

        # Update the cutout_info table with the location of the cutout within
        # the gallery
        self.gallery_info['gallery_xindex'] = ilist
        self.gallery_info['gallery_yindex'] = jlist
        self.gallery_info['gallery_xloc'] = minxlist
        self.gallery_info['gallery_yloc'] = minylist

        self.save_gallery_image()

    def save_gallery_image(self):
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, ax = plt.subplots(figsize=(12, 12))  # Adapt this to the number of rows/cols
        ax.imshow(self.gallery, origin='lower', cmap='Greys_r', norm=norm)
        ax.set_title('Snowball Gallery')
        fig.savefig(self.gallery_filename)

        h0 = fits.PrimaryHDU(self.gallery)
        hlist = fits.HDUList([h0])
        hlist.writeto(f'{os.path.splitext(self.gallery_filename)[0]}.fits', overwrite=True)

        # Also save the cutout information table in order to tell where in the gallery
        # each cutout lands
        full_filename = os.path.join(self.output_dir, f'{os.path.splitext(self.gallery_filename)[0]}.ecsv')
        self.gallery_info.write(full_filename, overwrite=True)


class ExpSnowballs():
    def __init__(self, filename, min_snowball_area=20, padding=10, max_ellipticity=0.9, min_total_signal=100000., create_segment_map=False):
        self.filename = filename
        self.source_catalog = None
        self.num_snowballs = 0
        self.cutouts = []
        self.cutout_info = Table(names=('filename', 'integration', 'group', 'x', 'y'), dtype=('S80', 'i4', 'i4', 'f4', 'f4'))
        self.max_ellipticity = max_ellipticity
        self.min_total_signal = min_total_signal
        self.create_segment_map = create_segment_map

        # This is the minimum number of flagged pixels which will
        # result in a source being considered a snowball. Below this
        # limit, it will simply be considered a CR
        self.min_snowball_area = min_snowball_area

        # This is the number of rows/columns to add around the photutils
        # bounding box when making the cutout image
        self.padding = padding

        self.run()

    def run(self):
        self.get_data()
        self.find_sources()
        self.exposure_time = self.header['EFFEXPTM']
        self.make_jump_map()
        self.find_jump_sources()

    def add_info_cols(self, catalog, integration, group):
        """Add the filename, integration number, and group number to the input catalog
        """
        catalog_len = len(catalog)
        fname_col = [self.filename] * catalog_len
        integ_col = [integration] * catalog_len
        grp_col = [group] * catalog_len

        # Calculate the absolute time associated with the group
        expstart = self.header['EXPSTART']
        tgroup = self.header['TGROUP']
        tframe = self.header['TFRAME']
        tint = tgroup * integration
        # Calculated group time is exposure start + integration number * int_time +
        # 1 reset frame per int after the first + group time * num_groups
        grouptime = expstart + (tint*integration + tframe*(integration-1) + tgroup*group) / 86400.
        grptime_col = [grouptime] * catalog_len

        catalog['filename'] = fname_col
        catalog['integration'] = integ_col
        catalog['group'] = grp_col
        catalog['absolute_time'] = grptime_col
        return catalog

    def create_cutouts(self, catalog):
        """Create a cutout of the snowball. Use the bounding box defined by photutils
        and pad with 10 rows/cols on each side
        """
        integration = catalog["integration"][0]
        group = catalog["group"][0]
        for row in catalog:
            xmin = row["bbox_xmin"] - self.padding
            xmax = row["bbox_xmax"] + self.padding + 1
            ymin = row["bbox_ymin"] - self.padding
            ymax = row["bbox_ymax"] + self.padding + 1

            if xmin < 0:
                xmin = 0
            if ymin < 0:
                ymin = 0

            if xmin == xmax:
                print('0 x length for cutout', self.filename, xmin, xmax, ymin, ymax)
            if ymin == ymax:
                print('0 y length for cutout', self.filename, xmin, xmax, ymin, ymax)


            if group != 0:
                self.cutouts.append(self.cds_data[integration, group-1, ymin:ymax, xmin:xmax])
                if np.max(self.cds_data[integration, group-1, ymin:ymax, xmin:xmax]) == 0:
                    print('Zero signal in cutout.', self.filename, xmin, xmax, ymin, ymax)
            else:
                # If the snowball is in group 0 of the data, then the best we can do for a cutout
                # is to look at group 0 of the original data, rather than a CDS image
                self.cutouts.append(self.data[integration, group, ymin:ymax, xmin:xmax])
                if np.max(self.cds_data[integration, group, ymin:ymax, xmin:xmax]) == 0:
                    print('Zero signal in cutout.', self.filename, xmin, xmax, ymin, ymax)

            self.cutout_info.add_row([self.filename, integration, group, (xmax+xmin)/2, (ymax+ymin)/2])

    def create_snowball_catalog(self, segmentation_map, integ_num, grp_num):
        # Create a catalog of all snowballs. In this case, we want to supply
        # the CDS data, rather than the dq values
        if grp_num == 0:
            print("Found a snowball in group 0. I didn't think the jump step could do that.")

        #cat = SourceCatalog(self.data[integ_num, grp_num, :, :], segmentation_map)
        cat = SourceCatalog(self.cds_data[integ_num, grp_num-1, :, :], segmentation_map, localbkg_width=50)
        cols_for_table = ['label', 'xcentroid', 'ycentroid', 'bbox_xmin', 'bbox_xmax', 'bbox_ymin', 'bbox_ymax', 'area',
                          'orientation', 'ellipticity', 'elongation', 'eccentricity', 'min_value',
                          'max_value', 'local_background', 'segment_flux', 'segment_fluxerr']
        cat_table = cat.to_table(columns=cols_for_table)

        # Snowballs are round. Filter out sources that are too elliptical
        good = ((cat_table['ellipticity'] <= self.max_ellipticity) & (np.abs(cat_table['segment_flux']) > self.min_total_signal))
        if np.any(good):
            cat_table = cat_table[good]
        else:
            cat_table = None

        return cat_table

    def find_jump_sources(self):
        nint, ngrp, ny, nx = self.jump_map.shape

        if self.create_segment_map:
            self.segmap = np.zeros((nint, ngrp, ny, nx))

        for integ in range(nint):
            for grp in range(ngrp):
                print(f'Working on int {integ+1}, group {grp}...')
                frame = self.jump_map[integ, grp, :, :]
                threshold = detect_threshold(self.jump_map[integ, grp, :, :], nsigma=2.)

                segm = detect_sources(self.jump_map[integ, grp, :, :], threshold, npixels=self.min_snowball_area,
                                      mask=self.data_segmap.data.astype(bool)) #, kernel=kernel)
                if segm is not None:

                    if self.create_segment_map:
                        self.segmap[integ, grp, :, :] = segm

                    cat = self.create_snowball_catalog(segm, integ, grp)

                    if cat is not None:
                        cat = self.add_info_cols(cat, integ, grp)
                        if self.source_catalog is None:
                            self.source_catalog = cat
                        else:
                            self.source_catalog = vstack([self.source_catalog, cat])
                        print(f'Found {len(cat)} snowballs!')

                        # Create cutouts of the found snowballs
                        self.create_cutouts(cat)

        if self.source_catalog is not None:
            self.num_snowballs = len(self.source_catalog)
        else:
            self.num_snowballs = 0

    def find_sources(self):
        """Testing shows that sometimes, highly saturated sources are being flagged
        as snowball-like CRs. So let's find the sources in the data, and exclude any
        pixels that contain real sources from the source detection on the jump maps
        """
        snapshot = self.data[0, 0, :, :] #- self.data[0, 0, :, :]
        threshold = detect_threshold(snapshot, nsigma=2.)
        self.data_segmap = detect_sources(snapshot, threshold, npixels=9)

    def get_data(self):
        with fits.open(self.filename) as hdulist:
            self.data = hdulist['SCI'].data
            self.dq = hdulist['GROUPDQ'].data
            self.header = hdulist[0].header

        self.cds_data = self.data[:, 1:, :, :] - self.data[:, 0:-1, :, :]

        # In order to fill in holes in snowballs (since the cores are often flagged as
        # saturated, and therefore not as jumps), add a jump flag to any pixel that
        # is flagged as saturated. Since we do source detection on the 0th group and
        # don't look for snowballs in pixels with sources, we should be filtering out
        # bright sources that will saturate. Do not adjust flags in the 0th group, since
        # we don't look for snowballs there.

        #We need to add the CR flag to a saturated pixel only in the first group where it saturates.
        #Once its saturated, it wont register any more jumps...
        box_hw = 25
        jump_flag_threshold = 30

        # Make a map of saturation flags
        sat_flag_map = flag_map(self.dq, 'SATURATED')   #  #1

        # Make a map of the group number in which each pixel is first saturated
        first_sat_group = collapse_dq_flag_map(sat_flag_map) #    #2
        jump_flag_map = flag_map(self.dq, 'JUMP_DET')

        # For each saturated pixel above, look at dq frame. if there are >threshold number of
        # jump flags nearby, then add a jump flag to the saturation flag
        nints, ngrps, ny, nx = self.dq.shape
        for intnum in range(nints):
            first_sat_int_dq = first_sat_group[intnum, :, :]
            for grpnum in range(ngrps):
                first_sat_y, first_sat_x = np.where(first_sat_int_dq == grpnum)

                # Loop over each saturated pixel. If there are a lot of jump flags nearby, then add
                # the jump flag to this pixel in this group and integration
                for satx, saty in zip(first_sat_x, first_sat_y):

                    # Sanity check
                    if first_sat_int_dq[saty, satx] != grpnum:
                        print(self.dq[intnum, grpnum, saty, satx])
                        raise ValueError('Something is wrong. Maybe x/y index switch')

                    box = jump_flag_map[intnum, grpnum, saty-box_hw:saty+box_hw, satx-box_hw:satx+box_hw]
                    njumps = np.sum(box)

                    if (jump_flag_map[intnum, grpnum, saty, satx] & dqflags.pixel["JUMP_DET"] == 0):
                        if njumps > jump_flag_threshold:
                            self.dq[intnum, grpnum, saty, satx] += dqflags.pixel["JUMP_DET"]

    def make_jump_map(self):
        # Keep only jump flags
        self.jump_map = (self.dq & dqflags.pixel['JUMP_DET'] > 0)


def collapse_dq_flag_map(dq_map):
    nints, ngroups, ny, nx = dq_map.shape

    # Create an array containing all group indexes
    all_groups = np.zeros((1, ngroups, 1, 1), dtype=np.int)
    all_groups[0, :, 0, 0] = np.arange(ngroups)
    intermediate1 = np.repeat(all_groups, nints, axis=0)
    intermediate2 = np.repeat(intermediate1, ny, axis=2)
    all_indexes = np.repeat(intermediate2, nx, axis=3)

    # Array to contain only group numbers of CR hits
    hit_indexes = np.zeros_like(all_indexes) + ngroups #np.nan

    # Find the CR flag locations
    hits = np.where(dq_map != 0)

    # All elements are NaN except the groups with CR hits
    hit_indexes[hits] = all_indexes[hits]

    # Find the minimum group number of the CR-hit groups for each pixel
    index_map = np.nanmin(hit_indexes, axis=1)
    return index_map


def flag_map(dq_array, mnemonic):
    """Return a map of cosmic ray flags given a data quality array
    Parameters
    ----------
    dq_array : numpy.ndarray
        GROUP_DQ extension of an integration
    mnemonic : str
        Name of a type of bad pixel. This must be recognized by the
        calibration pipeline (i.e. it must come from jwst.datamodels.dqflags)
    Returns
    -------
    cr_map : numpy.ndarray
        GROUP_DQ extension with all flags other than jumps removed
    """
    cr_map = (dq_array & dqflags.pixel[mnemonic.upper()] > 0)
    return cr_map
