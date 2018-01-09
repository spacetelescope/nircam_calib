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
        #self.darks = None
        #self.flats = None
        
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
    
    def trim_imagelist(self,imagelist,basenamepattern=None):
        '''
        only keep fits file that match basenamepattern
        '''

        #print imagelist
        if imagelist==None or len(imagelist)==0:
            print 'Nothing to do, no images!!!'
            return(imagelist)

        if basenamepattern!=None:
            newimagelist = []
            basenamepattern_compiled = re.compile(basenamepattern)
            # make sure the input files are all ok!
            for i in xrange(len(imagelist)):
            
                m = basenamepattern_compiled.search(imagelist[i])
                if m == None:
                    print 'SKIPPING',imagelist[i]
                else:
                    #if len(m.groups())==0:
                    #    raise RuntimeError,"%s is matching %s, but does not return basename" % (basenamepattern_compiled.pattern,imagelist[i])
                    #basename = m.groups()[0]
                    #print 'VVV',imagelist[i]
                    newimagelist.append(imagelist[i])
            if len(imagelist)!=len(newimagelist):
                if self.verbose>1:
                    print 'skipping %d out of %d images, %d left' % (len(imagelist)-len(newimagelist),len(imagelist),len(newimagelist))
            imagelist=newimagelist

        return(imagelist)
    
    def parse_reftypes_images(self,reftypes_and_imagelist,basenamepattern=None):
        reftypelist = []
        imagelist = []
        for s in reftypes_and_imagelist:
            if s in self.cfg.params['reftypes']:
                reftypelist.append(s)
            else:
                if not os.path.isfile(s):
                    raise RuntimeError,"ERROR: file %s does not exist, thus not a viable input file" % s
                imagelist.append(os.path.abspath(s))

        imagelist = self.trim_imagelist(imagelist,basenamepattern)
              
        return(reftypelist,imagelist)

    def getimtypes(self):
        if not ('imtype' in self.imtable.t.colnames):
            self.imtable.t['imtype']=None

        darkpattern = re.compile(self.cfg.params['inputfiles']['dark_pattern'])
        flatpattern = re.compile(self.cfg.params['inputfiles']['flat_pattern'])
            
        for i in xrange(len(self.imtable.t)):
            shortfilename = os.path.basename(self.imtable.t['fitsfile'][i])
            if darkpattern.search(shortfilename):
                self.imtable.t['imtype'][i]='dark'
            elif flatpattern.search(shortfilename):
                self.imtable.t['imtype'][i]='flat'
            else:
                print 'ERROR: image type of image %s is unknown!'
        

    def getimageinfo(self,imagelist,dateobsfitskey=None,timeobsfitskey=None,mjdobsfitskey=None):
        
        #self.imtable['fitsfile'].format('%s')
        self.imtable.t['fitsfile']=imagelist
        self.imtable.t['fitsID']=range(len(imagelist))
        self.imtable.t['imtype']=None
        self.imtable.t['skip']=False

        requiredfitskeys = self.cfg.params['inputfiles']['requiredfitskeys']
        if requiredfitskeys==None: requiredfitskeys=[]
        if type(requiredfitskeys) == types.StringType: requiredfitskeys=[requiredfitskeys]
        if dateobsfitskey!=None: requiredfitskeys.append(dateobsfitskey)
        if timeobsfitskey!=None: requiredfitskeys.append(timeobsfitskey)
        if mjdobsfitskey!=None:
            requiredfitskeys.append(mjdobsfitskey)
            mjdcol = mjdobsfitskey
        else:
            mjdcol = 'MJD'

        
        self.imtable.fitsheader2table('fitsfile',
                                      requiredfitskeys=requiredfitskeys,
                                      optionalfitskey=self.cfg.params['inputfiles']['optionalfitskeys'],
                                      raiseError=False,skipcolname='skip')
        self.imtable.dateobs2mjd(dateobsfitskey,mjdcol,mjdobscol=mjdobsfitskey,timeobscol=timeobsfitskey)
            
        self.getimtypes()
        # sort by MJD
        self.imtable.t.sort('MJD')

        #self.darks = self.imtable.t[np.where(self.imtable.t['imtype']=='dark')]
        #self.flats = self.imtable.t[np.where(self.imtable.t['imtype']=='flat')]

        return(0)
        
    def organize_inputfiles(self,reftypes_and_imagelist):

        # parse teh command line arguments for reftypes and images
        (reftypelist,imagelist) = self.parse_reftypes_images(reftypes_and_imagelist,basenamepattern=self.cfg.params['inputfiles']['basenamepattern'])
        self.reftypelist = reftypelist

        # make the image table and populate it with info. Also get teh darks and flats table
        self.getimageinfo(imagelist,
                          dateobsfitskey=self.cfg.params['inputfiles']['dateobs_fitskey'],
                          timeobsfitskey=self.cfg.params['inputfiles']['timeobs_fitskey'],
                          mjdobsfitskey=self.cfg.params['inputfiles']['mjdobs_fitskey'])

        self.detectors = set(self.imtable.t['DETECTOR'])
        
        if self.verbose:
            print '#################\n### %d images found!' % len(self.imtable.t)
            print '### %d darks, %d flats' % (len(np.where(self.imtable.t['imtype']=='dark')[0]),len(np.where(self.imtable.t['imtype']=='flat')[0]))
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
        '''
        returns list of Dark indeces, where the indices refer to the self.darks table
        '''

        self.imtable.t['skip'][7]=True
        # indices for dark frames
        dindex, = np.where(self.imtable.t['imtype']=='dark')
        # indices for detector and not skipped
        dindex4detector = dindex[np.where(np.logical_and(self.imtable.t['DETECTOR'][dindex]==detector,np.logical_not(self.imtable.t['skip'][dindex])))]
        if self.verbose>2:
            print 'Possible %d Darks for detector %s' % (len(dindex4detector),detector)
            print self.imtable.t[dindex4detector]

        #D = astrotableclass(names=('D1index','D1fitsID'),dtype=('i4', 'i4'))
        D = astrotableclass()
        D.t['D1index']=dindex4detector
        D.t['D1fitsID']=self.imtable.t['fitsID'][dindex4detector]
        D.t['D1index','D1fitsID'].format='%d'
        #print D.t
        #sys.exit(0)
        return(D)
            
    
    def getDDlist(self,detector,max_Delta_MJD=None):
        '''
        returns list of Dark-Dark pair indeces,  where the indices refer to the self.imtable table
        '''
        if self.verbose>1: print '# Getting DD list'
        #self.imtable.t['skip'][7]=True
        #print  self.imtable.t[6:11]
        
        # indices for dark frames
        dindex, = np.where(self.imtable.t['imtype']=='dark')
        # indices for detector and not skipped
        dindex4detector = dindex[np.where(np.logical_and(self.imtable.t['DETECTOR'][dindex]==detector,np.logical_not(self.imtable.t['skip'][dindex])))]
        if self.verbose>2:
            print 'Possible %d Darks for detector %s' % (len(dindex4detector),detector)
            print self.imtable.t[dindex4detector]

        DD = astrotableclass(names=('D1index','D2index','D1fitsID','D2fitsID'),dtype=('i4', 'i4', 'i4', 'i4'))
        i=0
        while i<len(dindex4detector)-1:
            if max_Delta_MJD!=None:
                if self.verbose>2: print 'Checking if fitsID=%d and %d can be DD' % (self.imtable.t['fitsID'][dindex4detector[i]],self.imtable.t['fitsID'][dindex4detector[i+1]])
                dMJD =self.imtable.t['MJD'][dindex4detector[i+1]]-self.imtable.t['MJD'][dindex4detector[i]]
                print 'dMJD:',dMJD
                if dMJD>max_Delta_MJD:
                    if self.verbose>2:
                        print 'Skipping fitsID=%d (MJD=%f) since fitsID=%d is not within timelimit (Delta MJD = %f>%f)!' % (self.imtable.t['fitsID'][dindex4detector[i]],self.imtable.t['MJD'][i],self.imtable.t['fitsID'][dindex4detector[i+1]],dMJD,max_Delta_MJD)
                    i+=1
                    continue
            if self.verbose>1: print 'Adding DD pair with fitsID=%d and %d' % (self.imtable.t['fitsID'][dindex4detector[i]],self.imtable.t['fitsID'][dindex4detector[i+1]])
            DD.t.add_row({'D1index':dindex4detector[i],'D2index':dindex4detector[i+1],'D1fitsID':self.imtable.t['fitsID'][dindex4detector[i]],'D2fitsID':self.imtable.t['fitsID'][dindex4detector[i+1]]})
            i+=2

        if self.verbose>2:
            print DD.t
        return(DD)
    
    def getFFlist(self,detector,max_Delta_MJD=None):
        '''
        returns list of Flat-Flat pair indeces, where the indices refer to the self.imtable table
        '''
        if self.verbose>1: print '# Getting FF list'
        #self.imtable.t['skip'][7]=True
        #print  self.imtable.t[6:11]
        
        # indices for dark frames
        findex, = np.where(self.imtable.t['imtype']=='flat')
        # indices for detector and not skipped
        findex4detector = findex[np.where(np.logical_and(self.imtable.t['DETECTOR'][findex]==detector,np.logical_not(self.imtable.t['skip'][findex])))]
        if self.verbose>2:
            print 'Possible %d Flats for detector %s' % (len(findex4detector),detector)
            print self.imtable.t[findex4detector]
        
        FF = astrotableclass(names=('F1index','F2index','F1fitsID','F2fitsID'),dtype=('i4', 'i4', 'i4', 'i4'))
        i=0
        while i<len(findex4detector)-1:
            if max_Delta_MJD!=None:
                if self.verbose>2: print 'Checking if fitsID=%d and %d can be FF' % (self.imtable.t['fitsID'][findex4detector[i]],self.imtable.t['fitsID'][findex4detector[i+1]])
                dMJD =self.imtable.t['MJD'][findex4detector[i+1]]-self.imtable.t['MJD'][findex4detector[i]]
                print 'dMJD:',dMJD
                if dMJD>max_Delta_MJD:
                    if self.verbose>2:
                        print 'Skipping fitsID=%d (MJD=%f) since fitsID=%d is not within timelimit (Delta MJD = %f>%f)!' % (self.imtable.t['fitsID'][findex4detector[i]],self.imtable.t['MJD'][i],self.imtable.t['fitsID'][findex4detector[i+1]],dMJD,max_Delta_MJD)
                    i+=1
                    continue
            if self.verbose>1: print 'Adding FF pair with fitsID=%d and %d' % (self.imtable.t['fitsID'][findex4detector[i]],self.imtable.t['fitsID'][findex4detector[i+1]])
            FF.t.add_row({'F1index':findex4detector[i],'F2index':findex4detector[i+1],'F1fitsID':self.imtable.t['fitsID'][findex4detector[i]],'F2fitsID':self.imtable.t['fitsID'][findex4detector[i+1]]})
            i+=2
        
        if self.verbose>2:
            print FF.t
        return(FF)
    
    def getDDFFlist(self,detector,DD_max_Delta_MJD=None,FF_max_Delta_MJD=None,DDFF_max_Delta_MJD=None):
        '''
        returns list of Flat-Flat pair indeces, where the indices refer to the self.imtable table
        '''
        if self.verbose>0: print '\n### Getting DDFF list'
        DD = self.getDDlist(detector,max_Delta_MJD=DD_max_Delta_MJD)
        FF = self.getFFlist(detector,max_Delta_MJD=FF_max_Delta_MJD)
        if self.verbose>0: print '### DD and FF lists created, no matching them!!!'
        
        DDFF = astrotableclass(names=('F1index','F2index','F1fitsID','F2fitsID','D1index','D2index','D1fitsID','D2fitsID'),dtype=('i4', 'i4', 'i4', 'i4','i4', 'i4', 'i4', 'i4'))
        ddcount = np.zeros(len(DD.t))

        for f in xrange(len(FF.t)):
            if self.verbose>2: print '# Finding DD pair for FF pair with fitsID=%d and %d' % (FF.t['F1fitsID'][f],FF.t['F2fitsID'][f])
            if DDFF_max_Delta_MJD!=None:
                FF_MJDmin = np.amin(np.array((self.imtable.t['MJD'][FF.t['F1index'][f]],self.imtable.t['MJD'][FF.t['F2index'][f]])))
                FF_MJDmax = np.amax(np.array((self.imtable.t['MJD'][FF.t['F1index'][f]],self.imtable.t['MJD'][FF.t['F2index'][f]])))
                if self.verbose>3: print 'FF MJDs:',FF_MJDmin,FF_MJDmax

            d_best = None
            ddcountmin=None
            for d in  xrange(len(DD.t)):
                if ddcountmin==None or ddcount[d]<ddcountmin:
                    if DDFF_max_Delta_MJD!=None:
                        if self.verbose>2: print '# Testing DD pair with fitsID=%d and %d' % (DD.t['D1fitsID'][d],DD.t['D2fitsID'][d])
                        DD_MJDmin = np.amin(np.array((self.imtable.t['MJD'][DD.t['D1index'][d]],self.imtable.t['MJD'][DD.t['D2index'][d]])))
                        DD_MJDmax = np.amax(np.array((self.imtable.t['MJD'][DD.t['D1index'][d]],self.imtable.t['MJD'][DD.t['D2index'][d]])))
                        print 'DD MJD range',DD_MJDmin,DD_MJDmax
                        dMJD = np.fabs([DD_MJDmax-FF_MJDmin,DD_MJDmin-FF_MJDmax,DD_MJDmax-FF_MJDmax,DD_MJDmin-FF_MJDmin])
                        max_dMJD = np.amax(dMJD)
                        if max_dMJD > DDFF_max_Delta_MJD:
                            print 'DD pair with fitsID=%d and %d cannot be used, dMJD=%f>%f' % (DD.t['D1fitsID'][d],DD.t['D2fitsID'][d],max_dMJD,DDFF_max_Delta_MJD)
                            continue
                    ddcountmin=ddcount[d]
                    d_best = d

            if d_best == None:
                print 'SKIPPING following FF pair since there is no matching DD pair!'
                print FF.t[f]
            else:
                if self.verbose>2: print 'SUCCESS! DD pair with fitsID %d and %d found for FF pair with fitsID %d and %d' % (DD.t['D1fitsID'][d_best],DD.t['D2fitsID'][d_best],FF.t['F1fitsID'][f],FF.t['F2fitsID'][f])
                ddcount[d_best]+=1
                DDFF.t.add_row({'D1index':DD.t['D1index'][d_best],
                                'D2index':DD.t['D2index'][d_best],
                                'D1fitsID':DD.t['D1fitsID'][d_best],
                                'D2fitsID':DD.t['D2fitsID'][d_best],
                                'F1index':FF.t['F1index'][f],
                                'F2index':FF.t['F2index'][f],
                                'F1fitsID':FF.t['F1fitsID'][f],
                                'F2fitsID':FF.t['F2fitsID'][f]})
                    
                  
        print DDFF.t
        return(DDFF)
        
    def get_inputimage_sets(self,reftype,detector,DD_max_Delta_MJD=None,FF_max_Delta_MJD=None,DDFF_max_Delta_MJD=None):
        imtypes = self.cfg.params[reftype]['imtypes']
        imagesets = []
        if imtypes == 'D':
            imagesets = self.getDlist(detector)
        elif imtypes == 'DD':
            imagesets = self.getDDlist(detector,max_Delta_MJD=DD_max_Delta_MJD)
        elif imtypes == 'FF':
            imagesets = self.getFFlist(detector,max_Delta_MJD=FF_max_Delta_MJD)
        elif imtypes == 'DDFF' or imtypes == 'FFDD':
            imagesets = self.getDDFFlist(detector,DD_max_Delta_MJD=DD_max_Delta_MJD,FF_max_Delta_MJD=FF_max_Delta_MJD,DDFF_max_Delta_MJD=DDFF_max_Delta_MJD)
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
                inputimagesets = self.get_inputimage_sets(reftype,detector,
                                                          DD_max_Delta_MJD=self.cfg.params['DD']['max_Delta_MJD'],
                                                          FF_max_Delta_MJD=self.cfg.params['FF']['max_Delta_MJD'],
                                                          DDFF_max_Delta_MJD=self.cfg.params['DDFF']['max_Delta_MJD'])
                continue
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
        

        sys.exit(0)
            
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
