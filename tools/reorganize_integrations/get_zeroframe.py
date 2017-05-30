from astropy.io import fits
from jwst.datamodels import RampModel
import argparse
import numpy as np

class GetZeroframe:
    """This class extracts the simulated zeroframe extension with proper headers.

    This function extracts the zeroframe extension from a simulated exposure
    and gives it the same headers as the simulation.

    Attributes:
        get_zeroframe:           Extracts the simulated zeroframe extension.

    Usage (in JWST pipeline environment):
        $: python fix_headers.py simulated_file.fits

    """


    def __init__(self):
        self.verbose = False


    def add_options(self, parser=None, usage=None):
        '''Adds in command line arguments.'''

        if parser is None:
            parser = argparse.ArgumentParser(usage=usage)
            parser.add_argument("simfile",
                                help="File that needs correct headers.")
        return parser



    def get_zeroframe(self,simfile,model):
        '''Extracts the zeroframe.'''

        # grab zeroframe data and put it into model instance
        empty = np.zeros((1,1,2048,2048))
        with fits.open(simfile) as h:
            zero = h[2].data

        empty[0,:,:,:] = zero
        model.data = empty

        # for now, just give zeros for other extensions
        model.pixeldq = np.zeros((2048,2048))
        model.groupdq = np.zeros((1,1,2048,2048))
        model.err = np.zeros((1,1,2048,2048))

        # save the unformatted zeroframe file
        outname = 'zeroframe_delme_'+simfile
        model.save(outname,clobber=True)

        return outname



    def run(self):
        '''Main function.'''

        # initialize ramp model instance and get zeroframe file
        model = RampModel(self.simfile)
        outname = self.get_zeroframe(self.simfile,model)

        bad_hdulist = fits.open(outname)
        print('\nUnformatted file: '+str(outname))
        dat = bad_hdulist[1].data
        shape = dat.shape

        # get headers from original simulated exposure
        good_hdulist = fits.open(self.simfile)
        print('\nFile with correct formatting:  '+str(self.simfile))
        headers = good_hdulist[0].header

        for k,v in zip(headers.keys(), headers.values()):
            #print(k,v)
            bad_hdulist[0].header[k]=v

        # fix NGROUPS for zeroframe size
        bad_hdulist[0].header['NGROUP'] = shape[1]
        bad_hdulist[0].header['NGROUPS'] = shape[1]

        print('\nTesting header values...')
        print('TFRAME: '+str(bad_hdulist[0].header['TFRAME']))
        print('DETECTOR: '+str(bad_hdulist[0].header['DETECTOR']))

        # save out the formatted file
        final_outname = outname[:-5]+"_properHeaders_delme.fits"
        print('\nSaving final file: '+str(final_outname))
        bad_hdulist.writeto(final_outname,clobber=True)



if __name__ == '__main__':
    usagestring = 'USAGE: get_zeroframe.py good_file.fits'

    fix = GetZeroframe()
    parser = fix.add_options(usage=usagestring)
    args = parser.parse_args(namespace=fix)

    fix.run()
