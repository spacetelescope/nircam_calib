#!/usr/bin/env python
'''                                                                                                                                                                                            
create reference files for JWST instruments
A. Rest
'''

import sys, os,re,types,glob
import scipy,argparse
import numpy as np
import astropy.io.fits as fits


from jwst_lib.models import RampModel

# get the root dir of the code. This is needed only if the scripts are not installed as a module!
if 'JWST_MKREFS_SRCDIR' in os.environ:
    rootdir = os.environ['JWST_MKREFS_SRCDIR']
    sys.path.extend([rootdir,
                     '%s/gain' % rootdir,
                     '%s/badpix_map' % rootdir])
#else:
#    rootdir = os.path.dirname(os.path.realpath(__file__))

from tools import scitableclass,yamlcfgclass

def defaultrefoptions(parser=None, usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')
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


class mkrefclass:
    def __init__(self):
        self.ref = RampModel()

if __name__=='__main__':

    mkrefs=mkrefsclass()
    parser = mkrefs.add_options()
    parser.add_argument("filepattern",help="file pattern list of input files")
    args = parser.parse_args()

    # set verbose level
    mkrefs.verbose = args.verbose
    mkrefs.debug = args.debug
