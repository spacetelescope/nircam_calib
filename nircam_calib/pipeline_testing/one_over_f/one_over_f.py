#! /usr/bin/env python

"""
Functions dealing with 1/f noise in JWST data.
"""

import numpy as np

import nn2 # need to put this in nircam_calib, then update this line

class Oof:
    def __init__(self):
        self.outdir = './'

        # Number of rows and columns of reference pixels
        self.refpix_rows = 4
        self.refpix_cols = 4
        
        # Fast read direction. 1/f noise appears as banding
        # along the fast read direction. Value can be
        # 'horizontal' or 'vertical'. NIRCam is horizontal.
        self.fast_read_dir = 'horizontal'
        
        # Coordinates of the 4 amps. Need to work on
        # each amp separately, and different instruments
        # have different orientations of the 4 amps.
        # Dictionary where for each intrument there are
        # 4 entries: lower left, upper left, upper right,
        # and lower right for each of the 4 amps. (x,y)
        # pairs for each point.
        self.amps = {}
        self.amps["nircam"] = [[4, 4], [512, 4], [1024, 4], [1536, 4]]
                          

        
    def remove_with_nn2(self, data, boxwidth, appwidth):
        """
        Remove 1/f noise from input data using NN2 method.
        This should remove ~100% of the 1/f noise in a dark
        current integration. The removed 1/f noise is saved
        separately from the corrected input file
        """
        print('Beginning 1/f correction, using NN2 method.')

        # Read in data. Assume JWST file format
        with fits.open(filename) as h:
            exposure = h['SCI'].data
            instrument = h[0].header['INSTRUME'].lower()            
            
        # Amp boundaries - make sure to subtract the reference
        # pixel rows and columns
        inst_amps = self.amps[instrument]
        inst_amps = [[c[0]-self.refpix_cols, c[1]-self.refpix_rows] for c in inst_amps]

        # If the fast read direction is verical, transpose
        # the data
        if self.fast_read_dir == 'vertical':
            exposure = np.transpose(exposure, [0, 1, 3, 2])
            inst_amps = [[c[1], c[0]] for c in inst_amps]

        # Add an entry to the amp_boundary list specifying the
        # far edge of amp 4
        expshape = exposure.shape
        inst_amps.append([expshape[-1], 0])
            
        # Define output filenames:
        indir, infile = os.path.split(filename)
        # File to contain the removed 1/f noise
        voutname = os.path.join(self.outdir, "NN2_V_for_" + infile)
        # File to contain the 1/f-corrected data
        coutname = os.path.join(self.outdir,
                                "one_over_f_corrected_data_from_" + infile)

        # Send only the science pixels into the NN2 engine
        yd, xd = expshape[-2:]
        xmax = np.int(xd - self.refpix_cols)
        ymax = np.int(yd - self.refpix_cols)
        corrdata, corrmean, corrdataerr, corrmeanunc = self.run_nn2(exposure[:, :, \
                                                    self.refpix_rows:xmax, self.refpix_cols:ymax],
                                                                    inst_amps, boxwidth, appwidth)

        # If the fast readout direction is vertical, transpose the NN2 signals back
        if self.fast_read_dir == 'vertical':
            corrdata = np.transpose(corrdata, [0, 1, 3, 2])
            corrmean = np.transpose(corrmean, [0, 2, 1])
            corrdataerr = np.transpose(corrdataerr, [0, 1, 3, 2])
            corrmeanunc = np.transpose(corrmeanunc, [0, 2, 1])

        # Save 1/f signals - put back into full frame
        corr = np.zeros(expshape)
        corr[:, :, self.refpix_rows: expshape[-2] - self.refpix_rows,
             self.refpix_cols: expshape[-1] - self.refpix_cols] = corrdata
        corr_err = np.zeros(expshape)
        corr_err[:, :, self.refpix_rows: expshape[-2] - self.refpix_rows,
             self.refpix_cols: expshape[-1] - self.refpix_cols] = corrdataerr
        corr_mean = np.zeros((expshape[0], expshape[2], expshape[3]))
        corr_mean[:, self.refpix_rows: expshape[-2] - self.refpix_rows,
             self.refpix_cols: expshape[-1] - self.refpix_cols] = corrmean
        corr_meanerr = np.zeros((expshape[0], expshape[2], expshape[3]))
        corr_meanerr[:, self.refpix_rows: expshape[-2] - self.refpix_rows,
             self.refpix_cols: expshape[-1] - self.refpix_cols] = corrmeanunc
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(corr)
        h2 = fits.ImageHDU(corr_err)
        h3 = fits.ImageHDU(corr_mean)
        h4 = fits.ImageHDU(corr_meanerr)
        hlist = fits.HDUList([h0, h1, h2, h3, h4])
        hlist.writeto(voutfile, overwrite=True)
        
        # Subtract calculated 1/f signal from integration
        exposure[:,:,4:yd-4,4:xd-4] += (corrdata - corrmean)

        # Save corrected exposure
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(exposure)
        hlist = fits.HDUList([h0, h1])
        hlist.writeto(coutname, overwrite=True)

    def run_nn2(self, data, amp_bounds, box, app):
        """
        Run NN2
        """
        shape = data.shape
        
        # Output variables
        nn2corrdata = zeros_like(data)
        nn2corrdataerr = zeros_like(data)
        nn2corrmean = zeros((shape[0], shape[-2], shape[-1]))
        nn2corrmeanunc = zeros_like(nn2corrmean)
        
        # Keep each integration separate
        for integration in range(shape[0]):
        
            # We don't want to mix pixels from different amplifiers,
            # so work on each quad separately
            for amp in xrange(4):
                qdata = data[integration, :,
                             amp_bounds[amp][1]: amp_bounds[amp+1][1],
                             amp_bounds[amp][0]: amp_bounds[amp+1][0]]

                # Calculate mean values here. One mean value per box of boxwidth
                # in each row

                # Work one row at a time to save memory
                zd, yd, xd = qdata.shape
                for rownum in range(yd):
                    if rownum % 100 == 0:
                        print("Amp: {}, Row: {}".format(amp+1, rownum))
                    rowdata = qdata[:, rownum, :]

                    # Calculate the matrix of differences
                    diffrowdata = np.zeros((zd * (zd-1) / 2, xd))
                    outidx = 0
                    for idx in range(zd):
                        for idx2 in range(idx+1, zd):
                            diffrowdata[outidx, :] = rowdata[idx, :] - rowdata[idx2, :]
                            outidx = outidx + 1

                # Calculate the means in overlapping boxes of size boxwidth
                # Output means,uncs are 2d arrays, groups x number of boxes in the row
                means, uncs, targx = self.calculate_means(diffrowdata, boxwidth)

                # Now we need to work on each meanbox of the row. For each box, construct
                # the NxN matrix of all possible difference values. This matrix is then fed
                # into the NN2 machinery.
                for box in range(means.shape[1]):
                    allmeans = zeros((zd, zd))
                    alluncs = zeros_like(allmeans)
                    idx = 0
                    for y in range(0, zd):
                        for x in range(y+1, zd):
                            allmeans[y, x] = means[idx, box]
                            allmeans[x, y] = 0. - means[idx, box]
                            alluncs[y, x] = uncs[idx, box]
                            alluncs[x, y] = uncs[idx, box]
                            idx = idx + 1

                    # Pass the means and uncertainties to the NN2 engine
                    nn2inst = nn2.nn2()
                    nn2inst.A = allmeans
                    nn2inst.E = alluncs
                    nn2inst.solveMatrix()

                    # Save the NN2 output vector and uncertainty vector for each
                    # pixel. We also need the mean vector, for later subtraction
                    #for xpt in range(qstart[amp]+targx[box]-(appwidth/2),qstart[amp]+targx[box]+(appwidth/2)):
                    xstart = amp_bound[amp][0] + targx[b] - appwidth/2
                    xend = amp_bound[amp][0] + targx[box] + appwidth/2 + 1
                    for xpt in range(xstart, xend):
                        if xpt >= (amp_bound[amp][0]) and (xpt < amp_bound[amp+1][0]):
                            nn2corrdata[integration, :, rownum, xpt] = nn2inst.V
                            nn2corrdataerr[integration, :, rownum, xpt] = nn2inst.Verr
                            nn2corrmean[integration, rownum, xpt] = np.mean(nn2inst.V)
                            nn2corrmeanunc[integration, rownum, xpt] = np.std(nn2inst.V) / np.sqrt(len(nn2inst.V))
        return nn2corrdata, nn2corrmean, nn2corrdataerr, nn2corrmeanunc

    def calculate_means(self, data, boxwidth):
        """
        Given a 2D array of many differences for a single row,
        calulate and return the means and uncertainties for every
        boxwidth pixel. 
        NOTE: if boxwidth=100, then we will calculate the mean in a
        box that is +/-50 pixels around the central pixel. Later, the
        solution will be applied within a box +/-25 pixels around the
        central pixel. This implies that we need overlap between 
        boxes. So for boxwidth of 100, we perform the calculation
        centered around every 50th pixel.

        Parameters:
        -----------
        data : obj
            numpy ndarray
        boxwidth : int
            With of box in pixels to use when calculating means
        
        Returns:
        --------
        returns : tuple
            means - 2d numpy array of mean values
            uncs - 2d numpy array of standard deviations
            targets - array of integers giving box centers
        """
        targets = np.arange(0, data.shape[1], boxwidth/2)
        means = np.zeros((data.shape[0], len(targets)))
        uncs = np.zeros_like(means)
        for group in range(data.shape[0]):
            for i, center in enumerate(targets):
                if center < boxwidth/2:
                    pix = data[group, 0:boxwidth/2+1]
                elif data.shape[1]-boxwidth/2 < (boxwidth/2):
                    pix = data[group, data.shape[1] - boxwidth/2:]
                else:
                    pix = data[group, center - (boxwidth/2.): center + (boxwidth/2.) + 1]
                    
                clipped, tlow, thigh = sigmaclip(pix, low=3., high=3.)
                means[group, i] = np.nanmean(clipped)
                uncs[group, i] = np.nanstd(clipped)
        return means, uncs, targets



def remove_with_pipeline(filename):
    """
    Remove 1/f noise from input data using the JWST
    Calibration pipeline's refpix step
    """
    from jwst.refpix import RefpixStep
    m = RefpixStep.call(filename, set_options, output_data=outname)
    

def measure_one_over_f():
    Is this an NN2 run, or some other way to measure?

def pipeline_efficacy(filename):
    from jwst.refpix import RefpixStep

    # Measure 1/f in input data
    orig = remove_with_nn2(filename)

    # Run pipeline refpix step
    pipeout = remove_with_pipeline(filename)
    
    # Measure 1/f in output data
    fixed = remove_with_nn2(outname)

    # Make meaningful plots of results
    something
