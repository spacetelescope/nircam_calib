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
from mkref import mkrefclass

class cmdsclass(astrotableclass):
    def __init__(self):
        astrotableclass.__init__(self)
        
        
        
class mkrefsclass(astrotableclass):
    def __init__(self):
        astrotableclass.__init__(self)

        #config file 
        self.cfg = None

        self.imtable = astrotableclass()
        self.images4ssb = astrotableclass()
        self.darks = None
        self.flats = None
        
        #
        self.DDtable = astrotableclass()
        self.FFtable = astrotableclass()
        self.DDFFtable = astrotableclass()

        self.cmdtable = cmdsclass()
        
    def define_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')
        parser.add_argument("reftypes_and_imagelist",nargs='+',help="list of ref types to be done and image (or file patterns) lists")
        parser.add_argument('-b','--batchmode',help="run the commands in batch mode (default=%(default)s)",action="store_true",default=False)
        parser.add_argument('--onlyshow',help="only show what would be done, but don't do it... (default=%(default)s)",action="store_true",default=False)

        mkref = mkrefclass()
        parser = mkref.refoptions4mkrefs(parser=parser)
        return(parser)
        
    def loadcfgfiles(self,maincfgfile,extracfgfiles=None,params=None,params4all=None,params4sections=None,requireParamExists=True):
        if self.cfg == None:
            self.cfg = yamlcfgclass()
        if self.cfg.loadcfgfiles(maincfgfile,extracfgfiles=extracfgfiles,
                                 params=params,params4all=params4all,params4sections=params4sections,
                                 requireParamExists=requireParamExists,verbose=self.verbose):
            raise RuntimeError,"Something went wrong when loading config files!"
        return(0)
    
    def collapse_imagelist(self,imagelist):
        '''
        remove the .fits ending, and then look for the unique basenames
        '''

        print imagelist
        if imagelist==None or len(imagelist)==0:
            print 'Nothing to do, no images!!!'
            return(imagelist)

        subpattern = re.compile('\.fits$')
        skipdict = {}
    
        # get rid of the 'fits' at the end of the filenames
        for i in xrange(len(imagelist)):
            imagelist[i]=subpattern.sub('.',imagelist[i])
            skipdict[imagelist[i]]=False


        for i in xrange(len(imagelist)):
            
            # if skipdict[imagelist[i]], then this image is already determined to be not one of the basenames, so skip it!
            if skipdict[imagelist[i]]:
                continue
            
            for s in xrange(len(imagelist)):
                if i==s: continue
                
                # if skipdict[imagelist[i]], then this image is already determined to be not one of the basenames, so skip it!
                if skipdict[imagelist[s]]:
                    #print 'SKIPPED FROM BEFORE',imagelist[s]
                    continue

                commonrootfilename = os.path.commonprefix([imagelist[i],imagelist[s]])
                if commonrootfilename == imagelist[i]:
                    if commonrootfilename == imagelist[s]:
                        raise RuntimeError,"%s=%s=%s, that should not happen!!!" % (commonrootfilename,imagelist[i],imagelist[s])
                    # if the commonrootfilename is equal to imagelist[i], then imagelist[s] must be a derivative!
                    skipdict[imagelist[s]]=True
                    #print 'iiiiiiiiiii skip',imagelist[s]
                if commonrootfilename == imagelist[s]:
                    if commonrootfilename == imagelist[i]:
                        raise RuntimeError,"%s=%s=%s, that should not happen!!!" % (commonrootfilename,imagelist[i],imagelist[s])
                    # if the commonrootfilename is equal to imagelist[s], then imagelist[i] must be a derivative!
                    skipdict[imagelist[i]]=True
                    #print 'sssssssssss skip',imagelist[i]
                
        # only keep the ones that are not skipped, and add the 'fits' back to it!!
        collapsed_imagelist=[]            
        for i in xrange(len(imagelist)):
            #print 'BBB',imagelist[i]+'fits',skipdict[imagelist[i]]
            if not skipdict[imagelist[i]]:
                collapsed_imagelist.append(imagelist[i]+'fits')
            
        #print collapsed_imagelist

        return(collapsed_imagelist)
        
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

        imagelist = self.collapse_imagelist(imagelist)
              
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
        self.imtable.t['fitsID']=range(len(imagelist))
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

    def get_optional_arguments(self,args,sysargv):
        fitspattern = re.compile('\.fits$')
        
        opt_arg_list = []
        for i in xrange(1,len(sysargv)):
            if sysargv[i] in args.reftypes_and_imagelist:
                if sysargv[i] in self.reftypelist:
                    print 'this is a reftype',sysargv[i]
                else:
                    # test if it is an image
                    if not fitspattern.search(sysargv[i]):
                        print 'not a fits file! filepattern?',sysargv[i]
            else:
                opt_arg_list.append(sysargv[i])
        print 'optional arguments:',opt_arg_list
        return(opt_arg_list)
        
    def getDlist(self,detector):
        Dlist = list(self.darks[np.where(np.logical_and(self.darks['DETECTOR']==detector,np.logical_not(self.darks['skip'])))]['fitsfile'])
        return(Dlist)
    
        
    def get_inputimage_sets(self,reftype,detector):
        imtypes = self.cfg.params[reftype]['imtypes']
        imagesets = []
        if imtypes == 'D':
            imagesets = self.getDlist(detector)
        else:
            raise RuntimeError,"ERROR: imtypes=%s not yet implemented!" % imtypes
        return(imagesets)
    
    def cmds4mkref(self,optionalargs):

        if self.verbose:
            print '##################################\n### Constructing commands'
            
        for reftype in self.reftypelist:

            if self.verbose:
                print '### Constructing %s commands' % reftype
            counter=0
            for detector in self.detectors:
                inputimagesets = self.get_inputimage_sets(reftype,detector)
                for inputimageset in inputimagesets:
                    cmdargs = '%s' % reftype
                    
                    if type(inputimageset) is types.ListType:
                        cmdargs+=' %s' % ' '.join(inputimageset)
                    else:
                        cmdargs+=' %s' % inputimageset
                        
                    cmdargs += ' '
                    cmdargs += ' '.join(optionalargs)

                    mkref = mkrefclass()
                    mkref.mkref(cmdargs.split(),onlyinit=True)
                    
                    if len(self.cmdtable.t)==0:
                        self.cmdtable.t['reftype']=np.array([reftype])
                        self.cmdtable.t['detector']=detector
                        self.cmdtable.t['outbasename']=None
                        self.cmdtable.t['cmdargs']=cmdargs
                    else:
                        self.cmdtable.t.add_row({'reftype':reftype,'detector':detector,'outbasename':None,'cmdargs':cmdargs})
                    counter+=1
                print '%d %s commands for detector %s' % (counter,reftype,detector)
        print '### in total, %d commands' % (len(self.cmdtable.t))
        if self.verbose>1:
            print 'Commands constructed:'
            print self.cmdtable.t
        

    def submitbatch(self):
        print "### submitbatch: NOT YET IMPLEMENTED!!!"
        
    def mkrefloop(self):
        for i in xrange(len(self.cmdtable.t)):
            
            print '### running mkref.py %s' % self.cmdtable.t['cmdargs'][i]
            mkref = mkrefclass()
            
            mkref.mkref(self.cmdtable.t['cmdargs'][i].split())
            
    def combinerefs(self):
        print "### combinerefs: NOT YET IMPLEMENTED!!!"
        
            
    def overview(self):
        print "### overview: NOT YET IMPLEMENTED!!!"
        
if __name__=='__main__':

    mkrefs=mkrefsclass()
    parser = mkrefs.define_options()
    args = parser.parse_args()

    # set verbose level
    mkrefs.verbose = args.verbose
    mkrefs.debug = args.debug
    
    # Load config files
    mkrefs.loadcfgfiles(args.cfgfile,
                        extracfgfiles=args.extracfgfile,
                        params=args.params,
                        params4all=args.pall,
                        params4sections=args.pp)

    print mkrefs.cfg.params
    mkrefs.organize_inputfiles(args.reftypes_and_imagelist)

    optionalargs = mkrefs.get_optional_arguments(args,sys.argv)

    mkrefs.cmds4mkref(optionalargs)

    
    
    if args.batchmode:
        mkrefs.submitbatch()
    else:
        mkrefs.mkrefloop()

    mkrefs.combinerefs()

    mkrefs.overview()
