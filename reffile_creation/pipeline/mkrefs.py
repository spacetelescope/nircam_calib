#!/usr/bin/env python
'''                                                                                                                                                                                            
create reference files for JWST instruments
A. Rest
'''

import sys, os,re,types,glob
import scipy,argparse
import numpy as np
import astropy.io.fits as fits

# get the root dir of the code. This is needed only if the scripts are not installed as a module!
if 'JWST_MKREFS_SRCDIR' in os.environ:
    rootdir = os.environ['JWST_MKREFS_SRCDIR']
    sys.path.extend([rootdir,
                     '%s/gain' % rootdir,
                     '%s/badpix_map' % rootdir])
#else:
#    rootdir = os.path.dirname(os.path.realpath(__file__))

from tools import astrotableclass,yamlcfgclass

class mkrefsclass(astrotableclass):
    def __init__(self):
        astrotableclass.__init__(self)

        #config file 
        self.cfg = None

        self.imtable = astrotableclass()
        self.darks = None
        self.flats = None
        
        #
        self.DDtable = astrotableclass()
        self.FFtable = astrotableclass()
        self.DDFFtable = astrotableclass()
    
    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument('--verbose', '-v', action='count')
        parser.add_argument('-d','--debug',help="debug: do lou's phot", action='count')

        parser.add_argument('-o','--outbasename', default=None,
                            help='full output basename, e.g. for photometry table. Supercedes all other output options, e.g. --outrootdir and --outsubdir (default=%default)')
        parser.add_argument('--outrootdir', default=None,
                            help='output root directory. If not specified, then the directory of input fits file is used. (default=%default)')
        parser.add_argument('--outsubdir', default=None,
                            help='subdir added to the output root directory (default=%default)')
        parser.add_argument('--addsuffix', default=None,
                            help='suffix added to the output basename (default=%default)')
        parser.add_argument('--detector2subdir',help="add instrument_detector as a subdir to output basename",action="store_true",default=False)
        parser.add_argument('--propID2subdir',help="add propID as a subdir to output basename",action="store_true",default=False)
        parser.add_argument('--visit2subdir',help="add visit as a subdir to output basename",action="store_true",default=False)

        # options for config file
        if 'JWST_MKREFS_CONFIGFILE' in os.environ and os.environ['JWST_MKREFS_CONFIGFILE']!='':
            cfgfile = os.environ['JWST_MKREFS_CONFIGFILE']
        else:
            cfgfile = None
        parser.add_argument('-c','--cfgfile', default=cfgfile,
                            help='main config file. (default=%default)')
        parser.add_argument('-e','--extracfgfile', action='append', default=None, 
                            help='additional config file. These cfg files do not need to have all parameters. They overwrite the parameters in the main cfg file. (default=%default)')
        parser.add_argument('-p', '--params', action='append', default=None, nargs=2,
                            help='"param val": change parameters in config file (section independent) (default=%default)')
        parser.add_argument('--pp', action='append', default=None, nargs=3,
                            help='"section param val". change parameters in section of config file (default=%default)')
        
        return parser

    def loadcfgfiles(self,maincfgfile,extracfgfiles=None,params=None,params4sections=None,requireParamExists=True):
        if self.cfg == None:
            self.cfg = yamlcfgclass()
        if self.cfg.loadcfgfiles(maincfgfile,extracfgfiles=extracfgfiles,
                                 params=params,params4sections=params4sections,
                                 requireParamExists=requireParamExists,verbose=self.verbose):
            raise RuntimeError,"Something went wrong when loading config files!"
        return(0)
    
    def parse_reftypes_images(self,reftypes_and_imagelist):
        reftypelist = []
        imagelist = []
        for s in reftypes_and_imagelist:
            if s in self.cfg.params['reftypes']:
                reftypelist.append(s)
            else:
                if not os.path.isfile(s):
                    raise RuntimeError,"ERROR: file %s does not exist, thus not a viable input file" % s
                imagelist.append(s)

        return(reftypelist,imagelist)

    def getimtypes(self):
        if not ('imtype' in self.imtable.t.colnames):
            self.imtable.t['imtype']=None

        darkpattern = re.compile(self.cfg.params['dark_pattern'])
        flatpattern = re.compile(self.cfg.params['flat_pattern'])
            
        for i in xrange(len(self.imtable.t)):
            shortfilename = os.path.basename(self.imtable.t['fitsfile'][i])
            if darkpattern.search(shortfilename):
                self.imtable.t['imtype'][i]='dark'
            elif flatpattern.search(shortfilename):
                self.imtable.t['imtype'][i]='flat'
            else:
                print 'WARNING: image type of image %s is unknown!'
        
    def getimageinfo(self,imagelist):
        
        #self.imtable['fitsfile'].format('%s')
        self.imtable.t['fitsfile']=imagelist
        self.imtable.t['imtype']=None
        self.imtable.t['skip']=False
        self.imtable.fitsheader2table('fitsfile',
                                      requiredfitskeys=self.cfg.params['requiredfitskeys'],
                                      optionalfitskey=self.cfg.params['optionalfitskeys'],
                                      raiseError=False,skipcolname='skip')
        self.imtable.dateobs2mjd('DATE-OBS','MJD-OBS',timeobscol='TIME-OBS')

        self.getimtypes()

        self.darks = self.imtable.t[np.where(self.imtable.t['imtype']=='dark')]
        self.flats = self.imtable.t[np.where(self.imtable.t['imtype']=='flat')]

        return(0)
        
    def organize_inputfiles(self,reftypes_and_imagelist):

        # parse teh command line arguments for reftypes and images
        (reftypelist,imagelist) = self.parse_reftypes_images(reftypes_and_imagelist)
        self.reftypelist = reftypelist

        # make the image table and populate it with info. Also get teh darks and flats table
        self.getimageinfo(imagelist)

        self.detectors = set(self.imtable.t['DETECTOR'])
        
        if self.verbose:
            print '#################\n### %d images found!' % len(self.imtable.t)
            print '### %d darks, %d flats' % (len(self.darks),len(self.flats))
            print '### %d detectors:' % (len(self.detectors)),", ".join(self.detectors)
            if self.verbose>1:
                print self.imtable.t
        

if __name__=='__main__':

    mkrefs=mkrefsclass()
    parser = mkrefs.add_options()
    parser.add_argument("reftypes_and_imagelist",nargs='+',help="list of ref types to be done and image (or file patterns) lists")
    args = parser.parse_args()

    # set verbose level
    mkrefs.verbose = args.verbose
    mkrefs.debug = args.debug
    
    # Load config files
    mkrefs.loadcfgfiles(args.cfgfile,
                        extracfgfiles=args.extracfgfile,
                        params=args.params,
                        params4sections=args.pp)

    
    mkrefs.organize_inputfiles(args.reftypes_and_imagelist)

    
