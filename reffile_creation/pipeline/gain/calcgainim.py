#!/usr/bin/env python


import sys, os, re, math, scipy, glob
import astropy.io.fits as pyfits
import numpy as np
import optparse
from refpixcorr import refpixcorrclass,frameresetclass

# put the tools directory into the path
# add pythonmodules to PATH
if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
else:
    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
    sys.exit(0)

from sigmacut import calcaverageclass
from texttable import txttableclass
from tools import rmfile

class gainimclass:
    def __init__(self):
        self.filename = None
        self.Nx=None
        self.Ny=None

        self.bpmmask = None

        # array list with readnoise info 
        self.rdnoise={}
        # set the extensions that are required to None
        for e in ['im','err','flag']:self.rdnoise[e]=None

        self.readnoisemethod=None

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('-v', '--verbose', action='count', dest='verbose',default=0)
        parser.add_option('--agressivemasking', default=False, action="store_true",
                          help='make aggressive bpm masks, up to 20-30% of masking')
        parser.add_option('-t','--test', default=False, action="store_true",
                          help='test')
        parser.add_option('--outsubdir'  , default=None , type="string",
                          help='Add this subdir to path (default=%default)')
        parser.add_option('-b','--boxsize'  , default=64 , type="int",
                          help='boxsize for stats (default=%default)')
        parser.add_option('--darks', default=(None,None), nargs=2, type="string",
                          help='specify the two darks to be used for the readnoise determination')
        parser.add_option('-a','--autodark4readnoise', default=False, action="store_true",
                          help='search for two images in workspace subdir with *DARK* fits and determine the readnoise from these images.')
        parser.add_option('--gmax4dark'  , default=None , type="int",
                          help='maximum group for dark frame (default=%default)')
        parser.add_option('--xmin'  , default=0 , type="int",
                          help='xmin for stats (default=%default)')
        parser.add_option('--xmax'  , default=2048 , type="int",
                          help='xmax for stats (default=%default)')
        parser.add_option('--ymin'  , default=0 , type="int",
                          help='ymin for stats (default=%default)')
        parser.add_option('--ymax'  , default=2048 , type="int",
                          help='ymax for stats (default=%default)')
        parser.add_option('--Ngroupskeyword'  , default='NGROUP' , type="string",
                          help='keyword for # of groups (default=%default)')
#        parser.add_option('-r','--readnoise'  , default=7.084/math.sqrt(2) , type="float",
#                          help='readnoise (default=%default)')
        parser.add_option('--gmin'  , default=0 , type="int",
                          help='minimum group (default=%default)')
        parser.add_option('--gmax'  , default=None , type="int",
                          help='maximum group (default=%default)')
        parser.add_option('--Smeanmax4fit'  , default=20000.0 , type="float",
                          help='all groups with a mean signal>Smeanmax4fit are not used for stdev and gain linear fit (default=%default)')
        parser.add_option('--skipcorrect4framereset' , default=False, action="store_true",
                          help='do not correct for frame reset (default=%default)')
        parser.add_option('--refpixboxsize'  , default=12, type="int",
                          help='boxsize for boxcar smoothing of biasdrift vector (default=%default)')
        parser.add_option('--skiprefpixcorr' , default=False, action="store_true",
                          help='do not correct 1/f with side refpix (default=%default)')
        return(parser)

    def calcgain(self,xinfo,yinfo,tsub,imtype,Smean_max=None):

        gainim = scipy.zeros((self.Ny,self.Nx),dtype=float)
        gainim_err = scipy.zeros((self.Ny,self.Nx),dtype=float)

        keys_all  =  tsub.CUT_inrange('S_mean',None,Smean_max)
        
        x=xinfo[0]
        while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            keys_xcut  = tsub.CUT_inrange('xmin',x,x,keys=keys_all)
            y=yinfo[0]
            while y<yinfo[1] and y+yinfo[2]<=self.Ny:
                keys  = tsub.CUT_inrange('ymin',y,y,keys=keys_xcut)
                keys  = tsub.CUT_none_vals('gain',keys=keys)
                
                # fits a slope through stdev2 versus g for dsub2
                (slope_b,slope_b_err,offset_a,offset_a_err) = tsub.straightline(keys,'geff','gain',dyval=None)
                if self.options.verbose:
                    print '%4d,%4d: %s fit to gain versus geff: slope=%.4f(%.4f) offset=%.4f(%.4f)' % (x,y,imtype,slope_b,slope_b_err,offset_a,offset_a_err) 
                    
                # the true gain at t=0 is at g=-1, since at g=0, already delta t time has passed
                gain = -1.0*slope_b+offset_a
                gain_err = math.sqrt((-1)*(-1)*slope_b_err*slope_b_err + offset_a_err*offset_a_err)
                gainim[y:y+yinfo[2],x:x+xinfo[2]]=gain
                gainim_err[y:y+yinfo[2],x:x+xinfo[2]]=gain_err

                y+=yinfo[2]
            x+=xinfo[2]

        gainfilename = self.outfilebasename+'.gain.%s.fits' % imtype
        phdu = pyfits.PrimaryHDU(data=None)
        hdu_gain = pyfits.ImageHDU(gainim,name='gain')
        hdu_gain_err = pyfits.ImageHDU(gainim_err,name='gain_err')
        hdulist = pyfits.HDUList([phdu,hdu_gain,hdu_gain_err])
        print 'Saving ',gainfilename
        rmfile(gainfilename)
        hdulist.writeto(gainfilename)

        del gainim, gainim_err, hdulist, hdu_gain, hdu_gain_err


    def correct_stdev2_with_rdnoiseterms(self,xinfo,yinfo,tsubik,tsubg0,tdsub1,tdsub2,tdsubij,Smean_max=None,Pmax=20.0,Nmin=3):

        print '### Correcting for readnoise'

        if self.readnoisemethod=='dsub12':
            if self.options.verbose: print 'Determining readnoise from dsub1 and dsub2 for correction, readnoisemethod="dsub12"!!'
        else:
            if self.options.verbose: print 'readnoisemethod="%s"' % self.readnoisemethod
            

        # arrays for the readnoise determined from dsub1 and dsub2. This is only to save the readnoise image
        rdnoiseim = scipy.zeros((self.Ny,self.Nx),dtype=float)
        rdnoiseim_err = scipy.zeros((self.Ny,self.Nx),dtype=float)
        rdnoiseim_Pskip = scipy.zeros((self.Ny,self.Nx),dtype=float)
        rdnoiseim_Nused = scipy.zeros((self.Ny,self.Nx),dtype=np.int32)
        rdnoiseim_flag = scipy.zeros((self.Ny,self.Nx),dtype=np.int16)

        sigmacut = calcaverageclass()
        keyssubik_all   =  tsubik.CUT_inrange('S_mean',None,Smean_max)
        keyssubg0_all =  tsubg0.CUT_inrange('S_mean',None,Smean_max)
        keyssub1_all  =  tdsub1.CUT_inrange('S_mean',None,Smean_max)
        keyssub2_all  =  tdsub2.CUT_inrange('S_mean',None,Smean_max)
        keyssubij_all  =  tdsubij.CUT_inrange('S_mean',None,Smean_max)

        x=xinfo[0]
        while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            keyssub1_xcut  = tdsub1.CUT_inrange('xmin',x,x,keys=keyssub1_all)
            keyssub2_xcut  = tdsub2.CUT_inrange('xmin',x,x,keys=keyssub2_all)
            keyssubij_xcut  = tdsubij.CUT_inrange('xmin',x,x,keys=keyssubij_all)
            keyssubg0_xcut = tsubg0.CUT_inrange('xmin',x,x,keys=keyssubg0_all)
            keyssubik_xcut   = tsubik.CUT_inrange('xmin',x,x,keys=keyssubik_all)
            y=yinfo[0]
            while y<yinfo[1] and y+yinfo[2]<=self.Ny:
                keyssub1  = tdsub1.CUT_inrange('ymin',y,y,keys=keyssub1_xcut)
                keyssub2  = tdsub2.CUT_inrange('ymin',y,y,keys=keyssub2_xcut)
                keyssubij  = tdsubij.CUT_inrange('ymin',y,y,keys=keyssubij_xcut)
                keyssubg0 = tsubg0.CUT_inrange('ymin',y,y,keys=keyssubg0_xcut)
                keyssubik   = tsubik.CUT_inrange('ymin',y,y,keys=keyssubik_xcut)
                
                tdsub2.printtxttable(keys=keyssub2,cols=['geff','stdev2'])

                # fits a slope through stdev2 versus g for dsub2
                (slope_b,slope_b_err,offset_a,offset_a_err) = tdsub2.straightline(keyssub2,'geff','stdev2',dyval=None)
                if self.options.verbose:
                    print '%4d,%4d: dsub1 fit to stdev2 versus geff: slope=%.4f(%.4f) offset=%.4f(%.4f)' % (x,y,slope_b,slope_b_err,offset_a,offset_a_err) 

                # subtract stdev2 from dsub2 slope fit from stdev2 of dsub1. This difference should be 2*readnoise^2
                stdev2array = scipy.zeros((len(keyssub1)),dtype=float)
                counter=0
                for key1 in keyssub1:
                    stdev2diff = tdsub1.getentry(key1,'stdev2')-(tdsub1.getentry(key1,'geff')*slope_b+offset_a)
                    stdev2array[counter]=stdev2diff
                    counter+=1
                #calculate the readnoise^2 with dsub12 method
                sigmacut.calcaverage_sigmacutloop(stdev2array,Nsigma=3.0,verbose=0)
                if sigmacut.mean!=None and sigmacut.mean>0.0:
                    rdnoise2_dsub12 = sigmacut.mean*0.5
                    rdnoise_dsub12 = math.sqrt(rdnoise2_dsub12) 
                    rdnoise_dsub12_err = math.sqrt(0.5/(4.0*sigmacut.mean)*sigmacut.mean_err*sigmacut.mean_err) 
                    rdnoiseim[y:y+yinfo[2],x:x+xinfo[2]]=rdnoise_dsub12
                    rdnoiseim_err[y:y+yinfo[2],x:x+xinfo[2]]=rdnoise_dsub12_err
                    Pskip=100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)
                    rdnoiseim_Nused[y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.Nused
                    rdnoiseim_Pskip[y:y+yinfo[2],x:x+xinfo[2]]=Pskip
                    if  Pskip<=Pmax and sigmacut.Nused>=Nmin:
                        rdnoiseim_flag[y:y+yinfo[2],x:x+xinfo[2]]=0
                    else:
                        rdnoiseim_flag[y:y+yinfo[2],x:x+xinfo[2]]=1
                else:
                    rdnoise_dsub12 = np.nan

                # Which readnoise to take for correction?
                if self.rdnoise['im'] is None:
                    rdnoise = rdnoise_dsub12
                else:
                    rdnoise = self.rdnoise['im'][int(y+0.5*yinfo[2]),int(x+0.5*xinfo[2])]
                rdnoise2=rdnoise*rdnoise

                if self.options.verbose:
                    print 'rdnoise=%.3f' % (rdnoise)

                # get the correction for the default diffim stdev2
                keyssubik_g0 = tsubik.CUT_inrange('g',0,0,keys=keyssubik)
                if len(keyssubik_g0)!=1:
                    print 'ERROR keyssubik_g0!',keyssubik_g0
                    tsubik.printtxttable(keys=keyssubik)
                    print 'g0'
                    tsubik.printtxttable(keys=keyssubik_g0)
                    sys.exit(0)
                correction4subik = tsubik.getentry(keyssubik_g0[0],'stdev2')

                # correct all the tables with the readnoise^2
                for key in keyssubik:
                    tsubik.setentry(key,'stdev2_rdnoisecorr',tsubik.getentry(key,'stdev2')-correction4subik)
                for key in keyssubg0:
                    tsubg0.setentry(key,'stdev2_rdnoisecorr',tsubg0.getentry(key,'stdev2')-4.0*rdnoise2)
                for key in keyssub1:
                    tdsub1.setentry(key,'stdev2_rdnoisecorr',tdsub1.getentry(key,'stdev2')-6.0*rdnoise2)
                for key in keyssub2:
                    tdsub2.setentry(key,'stdev2_rdnoisecorr',tdsub2.getentry(key,'stdev2')-4.0*rdnoise2)
                for key in keyssubij:
                    tdsubij.setentry(key,'stdev2_rdnoisecorr',tdsubij.getentry(key,'stdev2')-4.0*rdnoise2)

                y+=yinfo[2]
            x+=xinfo[2]
            
        rdnoisefilename = self.outfilebasename+'.rdnoise_dsub.fits'
        readnoisehdulist = pyfits.HDUList([pyfits.PrimaryHDU(data=None),
                                           pyfits.ImageHDU(rdnoiseim,name='rdnoise'),
                                           pyfits.ImageHDU(rdnoiseim_err,name='rdnoiseerror'),
                                           pyfits.ImageHDU(rdnoiseim_flag,name='flag'),
                                           pyfits.ImageHDU(rdnoiseim_Nused,name='Nused'),
                                           pyfits.ImageHDU(rdnoiseim_Pskip,name='Pskip')
                                           ])
        print 'Saving ',rdnoisefilename
        rmfile(rdnoisefilename)
        readnoisehdulist.writeto(rdnoisefilename)

    def calcgain_grid(self,diffim,xinfo,yinfo,t,g,mask=None,S=None,dS=None,outfile=None,fixmean=None,imtype=None):
        keys=[]
        sigmacut = calcaverageclass()
        x=xinfo[0]

        print '### calculating gain for imtype',imtype
        while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            y=yinfo[0]
            while y<yinfo[1] and y+yinfo[2]<=self.Ny:
                mask2use = None
                if not (mask is None):
                    mask2use = mask[y:y+yinfo[2],x:x+xinfo[2]]

                #print 'x:%4d %4d %4d y:%4d %4d %4d ' % (x,x+xinfo[2],xinfo[2],y,y+yinfo[2],yinfo[2])
                sigmacut.calcaverage_sigmacutloop(diffim[y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,saveused=True,verbose=0,fixmean=fixmean)
                key = sigmacut.results2texttable(t)
                keys.append(key)

                # For S and dS calculation: use exactly same pixels used for the stdev calculations!
                if not (mask is None):
                    mask4S = np.logical_or(mask2use,sigmacut.clipped)
                else:
                    mask4S = sigmacut.clipped

                if imtype!=None:
                    stdev2 = t.getentry(key,'stdev')*t.getentry(key,'stdev')
                    t.setentry(key,'stdev2',stdev2)

                if imtype in ['dsub2','dsubij']:
                    t.add2row(key,{'g':g,'geff':g-0.5,'xmin':x,'xdelta':xinfo[2],'ymin':y,'ydelta':yinfo[2]})
                else:
                    t.add2row(key,{'g':g,'geff':g,'xmin':x,'xdelta':xinfo[2],'ymin':y,'ydelta':yinfo[2]})
                if not (S is None):
                    #sigmacut.calcaverage_sigmacutloop(S[y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,verbose=0)
                    sigmacut.calcaverage_sigmacutloop(S[y:y+yinfo[2],x:x+xinfo[2]],mask=mask4S,Nsigma=0.0,verbose=0)
                    sigmacut.results2texttable(t,key=key,meancol='S_mean',meanerrcol='S_mean_err',stdevcol='S_stdev',Nusedcol='S_Nused',Nskippedcol='S_Nskipped',Pskippedcol='S_Pskipped',convergedcol='S_converged',iterationcol='S_i')
                if not (dS is None):
                    #sigmacut.calcaverage_sigmacutloop(dS[y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,verbose=0)
                    sigmacut.calcaverage_sigmacutloop(dS[y:y+yinfo[2],x:x+xinfo[2]],mask=mask4S,Nsigma=0.0,verbose=0)
                    sigmacut.results2texttable(t,key=key,meancol='dS_mean',meanerrcol='dS_mean_err',stdevcol='dS_stdev')

                y+=yinfo[2]
            x+=xinfo[2]

        if outfile != None:
            t.save2file(outfile,keys=keys,verbose=True)
        return(keys)

    def calc_flattenfactor(self,S1,S2,mask=None):
        mask2use = scipy.where(S2<=0.0,True,False)
        if not (mask is None): mask2use |= mask
        ratio = S1/S2
        sigmacut = calcaverageclass()
        sigmacut.calcaverage_sigmacutloop(ratio,mask=mask2use,Nsigma=3.0,verbose=0)
        print sigmacut.__str__()
        return(sigmacut.mean)

    def calc_cum_stdev2(self,t,xinfo,yinfo,gmin,gmax):
        t.configcols(['cum_stdev2_rdnoisecorr','ave_Smean','cum_avestdev2'],'f','%.3f',visible=1)
        x=xinfo[0]
        #print 'CCC',xinfo,yinfo
        while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            y=yinfo[0]
            while y<yinfo[1] and y+yinfo[2]<=self.Ny:
                keys = t.CUT_inrange('xmin',x,x)
                keys = t.CUT_inrange('ymin',y,y,keys=keys)
                keys = t.sortkeysbycols(keys,'g',asstring=0)
                sum=0.0
                for key in keys:
                    if t.getentry(key,'stdev2_rdnoisecorr')==None:
                        break
                    sum+=t.getentry(key,'stdev2_rdnoisecorr')
                    t.setentry(key,'cum_stdev2_rdnoisecorr',sum)
                y+=yinfo[2]
            x+=xinfo[2]

        for g in xrange(gmin,gmax):
            keys = t.CUT_inrange('g',g,g)
            keys = t.CUT_none_vals('cum_stdev2_rdnoisecorr',keys=keys)
            if len(keys)>0:
                cum_avestdev2 = t.calcaverage(keys,'cum_stdev2_rdnoisecorr')
                ave_Smean = t.calcaverage(keys,'S_mean')
                for key in keys:
                    t.setentry(key,'ave_Smean',ave_Smean)
                    t.setentry(key,'cum_avestdev2',cum_avestdev2)

        #t.printtxttable(keys=keys,cols=['g','xmin','xdelta','ymin','ydelta','stdev2_rdnoisecorr','cum_stdev2_rdnoisecorr'])

    #def readnoise_sq_4imtype(self,imtype):
    #    if imtype!=None:
    #        rdnoise2=self.options.readnoise*self.options.readnoise
    #        if imtype=='default':rdnoise2_tot = 2.0*rdnoise2
    #        elif imtype=='subg0':rdnoise2_tot = 4.0*rdnoise2
    #        #elif imtype=='subgm1':rdnoise2_tot = 2.0*rdnoise2
    #        elif imtype=='dsub1':rdnoise2_tot = 6.0*rdnoise2
    #        elif imtype=='dsub2':rdnoise2_tot = 4.0*rdnoise2
    #        else: raise RuntimeError,"Wrong imtype %s" % imtype
    #        return(rdnoise2_tot)
    #    else:
    #        return(None)

    #def calc_varying_readnoise(self,t,imtype,g0=2):
    #    t.configcols(['rdnoise2_var'],'f','%.3f',visible=1)
    #    keys = t.CUT_inrange('g',g0,g0)
    #    dS0  = t.calcaverage(keys,'dS_mean')
    #    print 'dS0:',dS0
    #    for key in t.allrowkeys:
    #        if t.getentry(key,'dS_mean')>=10.0:
    #            rdnoise2_tot = self.readnoise_sq_4imtype(imtype)
    #            rdnoise2_var = rdnoise2_tot * (t.getentry(key,'dS_mean')/dS0*t.getentry(key,'dS_mean')/dS0)
    #            stdev2 = t.getentry(key,'stdev2')
    #            t.setentry(key,'rdnoise2_var',rdnoise2_var)
    #            t.setentry(key,'stdev2_rdnoisecorrX',stdev2-rdnoise2_var)

    def calcstdev2_model(self,t,gmin,gmax,g0=2):
        t.configcols(['stdev2_rdnc_model','ave_Smean','ave_stdev2_rdnc_model'],'f','%.3f',visible=1)
        keys  = t.CUT_inrange('g',g0,g0)
        stdev2_g0= t.calcaverage(keys,'stdev2_rdnoisecorr')
        dS0   = t.calcaverage(keys,'dS_mean')
        print 'stdev2_g0:',stdev2_g0
        print 'dS0:',dS0
        for key in t.allrowkeys:
            if t.getentry(key,'dS_mean')>=10.0:
                dSratio = t.getentry(key,'dS_mean')/dS0
                t.setentry(key,'stdev2_rdnc_model',stdev2_g0*dSratio*dSratio)

        for g in xrange(gmin,gmax):
            keys = t.CUT_inrange('g',g,g)
            keys = t.CUT_none_vals('stdev2_rdnc_model',keys=keys)
            if len(keys)>0:
                ave_stdev2_rdnc_model = t.calcaverage(keys,'stdev2_rdnc_model')
                ave_Smean = t.calcaverage(keys,'S_mean')
                for key in keys:
                    t.setentry(key,'ave_Smean',ave_Smean)
                    t.setentry(key,'ave_stdev2_rdnc_model',ave_stdev2_rdnc_model)

    def calcgainmodel(self,t,xinfo,yinfo,g0=2):
        t.configcols(['gain_theo','gain0'],'f','%.3f',visible=1)
        keys = t.CUT_inrange('g',g0,g0)
        gain0= t.calcaverage(keys,'gain')
        dS0  = t.calcaverage(keys,'dS_mean')
        print 'gain0:',gain0
        print 'dS0:',dS0
        for key in t.allrowkeys:
            if t.getentry(key,'dS_mean')>=10.0:
                t.setentry(key,'gain_theo',gain0*dS0/t.getentry(key,'dS_mean'))
            if t.getentry(key,'gain')!=None:
                t.setentry(key,'gain0',t.getentry(key,'gain')*t.getentry(key,'dS_mean')/dS0)


    def calceffgain_dS(self,t):
        t.configcols(['gain'],'f','%.3f',visible=1)
        for key in t.allrowkeys:
            #if t.getentry(key,'S_mean')<=35000.0:
            if t.getentry(key,'stdev2_rdnoisecorr')>=3.0:
                gain = 2.0*t.getentry(key,'dS_mean')/t.getentry(key,'stdev2_rdnoisecorr')
                t.setentry(key,'gain',gain)

    def subtract_stdev2_g0(self,t):
        keys0 = t.CUT_inrange('g',0,0)
        stdev2_g0 = t.calcaverage(keys0,'stdev2_rdnoisecorr')
        keys = t.CUT_none_vals('stdev2_rdnoisecorr')
        for key in keys:
            t.setentry(key,'stdev2_rdnoisecorr',t.getentry(key,'stdev2_rdnoisecorr')-stdev2_g0)

    def getfitsfile(self,filename):
        if self.options.verbose:
            print '# Loading ',filename
        self.filename=filename
        data, hdr = pyfits.getdata(filename, 0, header=True)

        if not ('NGROUP' in hdr):
            hdr = pyfits.getheader(filename,0)
        if not ('NGROUP' in hdr):
            print 'ERROR, could not find NGROUP in header!!'
            sys.exit(0)

        if len(data.shape)<3:
            print 'ERROR: data cube needs to have at least 3 dimensions for %s, only %d!!!' %(filename,len(data.shape))
            sys.exit(0)

        Ngroups = hdr['NGROUP']
        if self.options.verbose>1:
             print 'Ngroups ',Ngroups

        Nint = hdr['NINT']
        if self.options.verbose>1:
             print 'Nint ',Nint

        if len(data.shape)==3:
            if self.options.verbose>1:
                print 'CV data format'
                # make sure data is in float!
                data=data*1.0
        elif len(data.shape)==4:
            if self.options.verbose>1:
                print 'SSB data format'
            scinew=scipy.zeros((Nint*Ngroups,data.shape[2],data.shape[3]), dtype=float)
            for i in xrange(Nint):
                scinew[i:i+Ngroups,:,:]=data[i,:,:,:]
            data = scinew
        else:
            sys.exit(0)

        # Make sure the NAXIS1,2 are consistent for all images
        if self.Nx==None:
            self.Nx = data.shape[2]
        else:
            if data.shape[2] !=self.Nx:
                print 'ERROR! Nx %d!=%d' % (data.shape[2],self.Nx)
                sys.exit(0)

        if self.Ny==None:
            self.Ny = data.shape[1] 
        else:
            if data.shape[1]!=self.Ny:
                print 'ERROR! Nx %d!=%d' % (data.shape[1],self.Ny)
                sys.exit(0)

        return(data,hdr)

    def darks4readnoise(self,darks4readnoiselist,xinfo,yinfo,mask=None,gmax=None,Pmax=20.0,Nmin=40,saveramps=False,
                        skiprefpixcorr=False,refpixboxsize=10,skipcorrect4framereset=False):

        print '\n###### Calculating the readnoise from the double difference of two dark frames'
        if self.readnoisemethod!='doublediff':
            raise RuntimeError,'ERROR: attempting to calculate the readnoise from the dark frames, but readnoise methos = %s!' % self.readnoisemethod

        if len(darks4readnoiselist)<2:
            print 'Not enough darks, at least 2 needed!',darks4readnoiselist
            sys.exit(0)

        data1,hdr1 = self.getfitsfile(darks4readnoiselist[0])
        data2,hdr2 = self.getfitsfile(darks4readnoiselist[1])
        
        if gmax == None:
            gmax = data1.shape[0]
        else:
            if data1.shape[0]<gmax:
                gmax = data1.shape[0]   

        if self.options.verbose>1:
            print 'gmax for dark frames:',gmax

        #print 'VVVVVVVVV11'
        #pyfits.writeto('test.diffim.fits',data1[1:,:,:]-data1[0,:,:],clobber=True)
        # make the frame reset correcton and side ref pix correction if wanted
        self.make_refpixcorr(data=data1,gmax=gmax,skiprefpixcorr=skiprefpixcorr,refpixboxsize=refpixboxsize,skipcorrect4framereset=skipcorrect4framereset)
        #print 'VVVVVVVVV'
        #pyfits.writeto('test.finaldiffim.fits',data1[1:,:,:]-data1[0,:,:],clobber=True)
        #sys.exit(0)
        self.make_refpixcorr(data=data2,gmax=gmax,skiprefpixcorr=skiprefpixcorr,refpixboxsize=refpixboxsize,skipcorrect4framereset=skipcorrect4framereset)


        # double difference. stepsize is two groups so that the differences are completely independent!
        diffimdoublediff = (data1[1:gmax:2,:,:]-data1[0:gmax-1:2,:,:]) - (data2[1:gmax:2,:,:]-data2[0:gmax-1:2,:,:])
        #pyfits.writeto('diffimdoublediff.fits',data=diffimdoublediff)

        # define image array ramps for readnoise determined from the double difference of dark groups
        rdnoiseramp={}
        rdnoiseramp['im']= scipy.zeros(diffimdoublediff.shape,dtype=float)
        rdnoiseramp['err']= scipy.zeros(diffimdoublediff.shape,dtype=float)
        rdnoiseramp['flag']= scipy.ones(diffimdoublediff.shape,dtype=np.int16)
        rdnoiseramp['meanval']=scipy.zeros(diffimdoublediff.shape,dtype=float) 
        rdnoiseramp['Pskip']= scipy.zeros(diffimdoublediff.shape,dtype=float)
        rdnoiseramp['Nused']=scipy.zeros(diffimdoublediff.shape,dtype=np.int32) 

        # define image arrays for average readnoise determined from the different groups of the double difference of dark groups
        self.rdnoise['im']= scipy.zeros((self.Ny,self.Nx),dtype=float)
        self.rdnoise['err']= scipy.zeros((self.Ny,self.Nx),dtype=float)
        self.rdnoise['flag']= scipy.ones((self.Ny,self.Nx),dtype=np.int16)
        self.rdnoise['stdev']=scipy.zeros((self.Ny,self.Nx),dtype=float) 
        self.rdnoise['Pskip']= scipy.zeros((self.Ny,self.Nx),dtype=float)
        self.rdnoise['Nused']=scipy.zeros((self.Ny,self.Nx),dtype=np.int32) 

        # for consistency check: readnoise determined from the single difference of dark groups
        diffimsinglediff = (data1[1:gmax:2,:,:]-data1[0:gmax-1:2,:,:])
        rdnoiseramp['singlediff']= scipy.zeros(diffimsinglediff.shape,dtype=float)
        self.rdnoise['singlediff']= scipy.zeros((self.Ny,self.Nx),dtype=float)

        sigmacut = calcaverageclass()
        
        x=xinfo[0]
        # loop through all boxes
        while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            y=yinfo[0]
            print 'x=%d' % x
            while y<yinfo[1] and y+yinfo[2]<=self.Ny:
                mask2use = None
                if not (mask is None):
                    mask2use = mask[y:y+yinfo[2],x:x+xinfo[2]]

                # loop through all boxes and groups
                for g in xrange(0,diffimdoublediff.shape[0]):
                    # get statistics for double diff
                    sigmacut.calcaverage_sigmacutloop(diffimdoublediff[g,y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,verbose=3)
                    #sigmacut.calcaverage_sigmacutloop(diffimdoublediff[g,y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,verbose=3,percentclip_firstiteration=20)

                    #print sigmacut.__str__()
                    if sigmacut.mean!=None:
                        rdnoiseramp['im'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev*0.5
                        rdnoiseramp['Nused'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.Nused
                        rdnoiseramp['err'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev_err*0.5
                        Pskip=100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)
                        rdnoiseramp['Pskip'][g,y:y+yinfo[2],x:x+xinfo[2]]=Pskip
                        rdnoiseramp['meanval'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean
                        if Pskip<=Pmax and sigmacut.Nused>=Nmin:
                            rdnoiseramp['flag'][g,y:y+yinfo[2],x:x+xinfo[2]]=0
                        else:
                            rdnoiseramp['flag'][g,y:y+yinfo[2],x:x+xinfo[2]]=1

                    # get statistics for single diff
                    sigmacut.calcaverage_sigmacutloop(diffimsinglediff[g,y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,verbose=0)
                    rdnoiseramp['singlediff'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev*0.70710678118655

                (ymid,xmid) = (int(y+0.5*yinfo[2]),int(x+0.5*xinfo[2]))
                sigmacut.calcaverage_sigmacutloop(rdnoiseramp['im'][:,ymid,xmid],mask=rdnoiseramp['flag'][:,ymid,xmid],Nsigma=3.0,verbose=0)
                if sigmacut.mean!=None:
                    self.rdnoise['im'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean
                    self.rdnoise['Nused'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.Nused
                    self.rdnoise['err'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean_err
                    Pskip=100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)
                    self.rdnoise['Pskip'][y:y+yinfo[2],x:x+xinfo[2]]=Pskip
                    self.rdnoise['stdev'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev
                    if sigmacut.Nused>=Nmin:
                        self.rdnoise['flag'][:y+yinfo[2],x:x+xinfo[2]]=0
                    else:
                        self.rdnoise['flag'][y:y+yinfo[2],x:x+xinfo[2]]=1
                sigmacut.calcaverage_sigmacutloop(rdnoiseramp['singlediff'][:,ymid,xmid],mask=rdnoiseramp['flag'][:,ymid,xmid],Nsigma=3.0,verbose=0)
                self.rdnoise['singlediff'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean

                y+=yinfo[2]
            x+=xinfo[2]
            
        # save the results in an hdulist
        rdnoisefilename = self.outfilebasename+'.doublediffrdnoise.fits'
        readnoisehdulist = pyfits.HDUList([pyfits.PrimaryHDU(data=None),
                                           pyfits.ImageHDU(self.rdnoise['im'],name='rdnoise'),
                                           pyfits.ImageHDU(self.rdnoise['err'],name='rdnoiseerror'),
                                           pyfits.ImageHDU(self.rdnoise['stdev'],name='stdev'),
                                           pyfits.ImageHDU(self.rdnoise['flag'],name='flag'),
                                           pyfits.ImageHDU(self.rdnoise['Nused'],name='Nused'),
                                           pyfits.ImageHDU(self.rdnoise['Pskip'],name='Pskip'),
                                           pyfits.ImageHDU(self.rdnoise['singlediff'],name='rdnoiseimsinglediff')
                                           ])
        print 'Saving rdnoise',rdnoisefilename
        rmfile(rdnoisefilename)
        readnoisehdulist.writeto(rdnoisefilename)

        if saveramps:
            # save the results in an hdulist
            rdnoiserampfilename = self.outfilebasename+'.doublediffrdnoiseramp.fits'
            readnoiseramphdulist = pyfits.HDUList([pyfits.PrimaryHDU(data=None),
                                                   pyfits.ImageHDU(rdnoiseramp['im'],name='rdnoise'),
                                                   pyfits.ImageHDU(rdnoiseramp['err'],name='rdnoiseerror'),
                                                   pyfits.ImageHDU(rdnoiseramp['flag'],name='flag'),
                                                   pyfits.ImageHDU(rdnoiseramp['Nused'],name='Nused'),
                                                   pyfits.ImageHDU(rdnoiseramp['Pskip'],name='Pskip'),
                                                   pyfits.ImageHDU(rdnoiseramp['meanval'],name='meanval'),
                                                   pyfits.ImageHDU(rdnoiseramp['singlediff'],name='rdnoiseimsinglediff')
                                                   ])
            print 'Saving rdnoise ramp',rdnoiserampfilename
            rmfile(rdnoiserampfilename)
            readnoiseramphdulist.writeto(rdnoiserampfilename)
            del readnoiseramphdulist

        del data1,data2,hdr1,hdr2,rdnoiseramp,diffimdoublediff,diffimsinglediff,readnoisehdulist
      

    def mkbpmmask(self,flatdata,hdr,xinfo,yinfo,fluxcutfraction=0.7,minmedian=35000.0,maxmedian=40000.0):
        print '\n##### Making bpm mask',xinfo,yinfo

        bpm = scipy.zeros((self.Ny,self.Nx),dtype=np.int16)
        flatdata_corr = flatdata[1:,:,:]-flatdata[0,:,:]

        g0 = None
        for g in xrange(flatdata_corr.shape[0]-1,0,-1):
            if self.options.verbose>1: print 'g=%d' % g
            median = np.median(flatdata_corr[g,:,:])
            if (median>minmedian) and (median<maxmedian):
                g0=g
            if (median<=minmedian):
                if g0==None:g0=g
                break
        if g0==None:
            raise RuntimeError,'Could not determine a good g for which the median flux is between %.0f and %.0f' % (minmedian,maxmedian) 

        flatdata_corr0=flatdata_corr[g0,:,:]

        flatdata4bpmfilename= self.outfilebasename+'.flatdata4bpm.fits'
        print 'Saving flatdata4bpm:',flatdata4bpmfilename
        rmfile(flatdata4bpmfilename)
        pyfits.writeto(flatdata4bpmfilename,data=flatdata_corr0)

        # crude cut: get rid of overscan and the most blatant bad pixels
        median0 = np.median(flatdata_corr0)
        bpm[scipy.where(flatdata_corr0<fluxcutfraction*median0)]=1
        bpm[scipy.where(np.isnan(flatdata_corr0))]=1
        bpm[:4,:]=1
        bpm[-4:,:]=1
        bpm[:,:4]=1
        bpm[:,-4:]=1

        #bpmfilename= self.outfilebasename+'.bpmmaskcrude.fits'
        #print 'Saving bpm:',bpmfilename
        #rmfile(bpmfilename)
        #pyfits.writeto(bpmfilename,data=bpm)

        sigmacut = calcaverageclass()
        x=xinfo[0]
        # loop through all boxes
        print '# looping through boxes to identify pixels with similar fluxes'
        while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            y=yinfo[0]
            while y<yinfo[1] and y+yinfo[2]<=self.Ny:

                flatbox = flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]

                #medianflux = np.median(flatdata_corr[g0,y:y+yinfo[2],x:x+xinfo[2]])
                #medianflux = np.median(flatbox[scipy.where(bpm[y:y+yinfo[2],x:x+xinfo[2]]==0)])
                #if self.options.verbose>1:
                #    print 'x=%d, y=%d, median flux=%f' % (x,y,medianflux)
                #print 'TEST',np.median(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]])

                flatvals = sorted(flatbox[scipy.where(bpm[y:y+yinfo[2],x:x+xinfo[2]]==0)])
                N = len(flatvals)
                medianflux = flatvals[int(N*0.5)]
                fluxplus = flatvals[int(N*(0.5+0.5*0.6827))]
                fluxminus = flatvals[int(N*(0.5-0.5*0.6827))]
                if self.options.verbose>1:
                    print 'median flux: %.1f, sigma+=%.1f, sigma-=%.1f' % (medianflux,np.abs(medianflux-fluxplus),np.abs(medianflux-fluxminus))
                sigma=min([np.abs(medianflux-fluxplus),np.abs(medianflux-fluxminus)])
                Nsigma=3.0
                if self.options.agressivemasking:
                    bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]<medianflux-Nsigma*sigma,2,0)
                    bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]>medianflux+Nsigma*sigma,4,0)
                else:
                    bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]<medianflux-Nsigma*np.abs(medianflux-fluxminus),2,0)
                    bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]>medianflux+Nsigma*np.abs(medianflux-fluxplus),4,0)

                if self.options.verbose>1:
                    print '%d' % (100.0*np.count_nonzero(bpm[y:y+yinfo[2],x:x+xinfo[2]])/(xinfo[2]*yinfo[2])),'% of pixels masked in this box'

                # old way
                # get statistics for double diff
                #sigmacut.calcaverage_sigmacutloop(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]],mask=bpm[y:y+yinfo[2],x:x+xinfo[2]],fixmean=medianflux,Nsigma=3.0,saveused=True,verbose=0)
                #print sigmacut.__str__()
                #if sigmacut.mean!=None:
                #    bpm[y:y+yinfo[2],x:x+xinfo[2]] |= sigmacut.clipped*2

                y+=yinfo[2]
            x+=xinfo[2]

        bpmfilename= self.outfilebasename+'.bpmmask.fits'
        print 'Saving bpm:',bpmfilename
        rmfile(bpmfilename)
        pyfits.writeto(bpmfilename,data=bpm)

        return(bpm)

    def make_refpixcorr(self,data,gmax=None,skiprefpixcorr=False,refpixboxsize=10,skipcorrect4framereset=False):
        test=False
        if gmax == None:
            gmax = data.shape[0]

        framereset = frameresetclass()
        refpixcorr = refpixcorrclass()

        if not skipcorrect4framereset:
            framereset.correct4framereset(data,method='first')
           
        if not skiprefpixcorr:
            refpixcorr.outfilebasename='test'
            refpixcorr.verbose=0
            refpixcorr.mkframeresetcorrection(data,boxsize=refpixboxsize,testflag=False)

        if not skipcorrect4framereset:
            framereset.correct4framereset(data,method='next')

            
    def make_refpixcorr_delme(self,data,gmax=None,skiprefpixcorr=False,refpixboxsize=10,skipcorrect4framereset=False):
        test=False
        if gmax == None:
            gmax = data.shape[0]

        if not skiprefpixcorr:
            # refpix correction
            refpixcorr = refpixcorrclass()
            refpixcorr.outfilebasename='test'
            refpixcorr.verbose=3
            refpixcorr.data = data
            refpixcorr.mkdiffim(savediffimflag=True)
            refpixcorr.calcbiasdriftvec_allgroups(boxsize=refpixboxsize)
            refpixcorr.savebiasdriftvecs_as_txt(boxsize=refpixboxsize)
            refpixcorr.savebiasdriftvecs(boxsize=refpixboxsize)
            refpixcorr.correctimage4leftrightrefpix(boxsize=refpixboxsize)
            refpixcorr.mkdiffimcorr(savediffimflag=True)

        if not skipcorrect4framereset:
            sigmacutrefpixbottom = calcaverageclass()
            sigmacutrefpixtop = calcaverageclass()
            print 'Frame reset correction!'
            for g in xrange(self.options.gmin+1,gmax):
                if self.options.verbose:
                    print 'g',g

                for amp in xrange(1,5):
                    print 'amp',amp
                    xmin4amp = (amp-1)*512
                    xmax4amp = amp*512
                    
                    # image 1
                    refpixbottom = data[g,:4,xmin4amp:xmax4amp]-data[g-1,:4,xmin4amp:xmax4amp]
                    sigmacutrefpixbottom.calcaverage_sigmacutloop(refpixbottom,Nsigma=3,verbose=3)
                    if sigmacutrefpixbottom.mean==None:
                        raise RuntimeError,"Could not correct for frame reset!"

                    refpixtop = data[g,-4:,xmin4amp:xmax4amp]-data[g-1,-4:,xmin4amp:xmax4amp]
                    if test: refpixtop = data[g,:4,xmin4amp:xmax4amp]-data[g-1,:4,xmin4amp:xmax4amp]
                    sigmacutrefpixtop.calcaverage_sigmacutloop(refpixtop,Nsigma=3,verbose=3)
                    if sigmacutrefpixtop.mean==None:
                        raise RuntimeError,"Could not correct for frame reset!"

                    frameresetcorrection = 0.5*(sigmacutrefpixbottom.mean+sigmacutrefpixtop.mean)
                    data[g:,:,xmin4amp:xmax4amp] -= frameresetcorrection
        

    def gainim(self,infile1,infile2,darks4readnoiselist=None,pattern='.fits',format='g%03d',imtype=None,
               skiprefpixcorr=False,
               refpixboxsize=10,
               skipcorrect4framereset=False):

        if darks4readnoiselist!=None:
            self.readnoisemethod='doublediff'
        else:
            self.readnoisemethod='dsub12'
            
        self.outfilebasename = re.sub('\.fits$','',infile1)+'.%s' % self.readnoisemethod
        if self.options.outsubdir!=None:
            (head,tail)=os.path.split(os.path.abspath(self.outfilebasename))
            newoutdir = '%s/%s' % (head,self.options.outsubdir)
            if not os.path.isdir(newoutdir):
                os.makedirs(newoutdir)
            if not os.path.isdir(newoutdir):
                raise RuntimeError, 'ERROR: Cannot create directory %s' % newoutdir

            self.outfilebasename = '%s/%s' % (newoutdir,tail)

        print 'outfilebasename:',self.outfilebasename

        def initcols(ts):
            for t in ts:
                t.configcols(['g'],'d','%d',visible=1)
                t.configcols(['geff'],'f','%.1f',visible=1)
                t.configcols(['xmin','xdelta','ymin','ydelta'],'d','%d',visible=1)
                t.configcols(['stdev2','stdev2_rdnoisecorr'],'f','%.3f',visible=1)

        #if self.options.test:
        #    xmin = 1700
        #    ymin = 1700
        #    xmax = int(xmin+2*128)
        #    ymax = int(ymin+2*128)
        #    boxsize = 128
        #    if self.options.gmax==None:
        #        self.options.gmax = 5
        #    self.options.gmax4dark=10
        if self.options.test:
            boxsize = 32
            xmin = 1024
            ymin = 0
            xmax = int(xmin+2*boxsize)
            ymax = int(ymin+2*boxsize)
            #xmin = 0
            #xmax = 256
            #ymin = 0
            #ymax = 256
            if self.options.gmax==None:
                self.options.gmax = 5
            self.options.gmax4dark=10
        else:
            xmin = self.options.xmin
            xmax = self.options.xmax
            ymin = self.options.ymin
            ymax = self.options.ymax
            boxsize = self.options.boxsize

        print 'x:%d-%d  y:%d-%d' % (xmin,xmax,ymin,ymax)


        # load the flats
        self.data1,self.hdr1 = self.getfitsfile(infile1)
        self.data2,self.hdr2 = self.getfitsfile(infile2)

        # make the mask files
        self.bpmmask = self.mkbpmmask(self.data1,self.hdr1,(xmin,xmax,boxsize),(ymin,ymax,boxsize))

        # get the readnoise from the dark frames if wanted...
        if darks4readnoiselist!=None:
            self.darks4readnoise(darks4readnoiselist,(xmin,xmax,boxsize),(ymin,ymax,boxsize),mask=self.bpmmask,gmax=self.options.gmax4dark,
                                 skiprefpixcorr=skiprefpixcorr,refpixboxsize=refpixboxsize )

        Ngroups = self.hdr2[self.options.Ngroupskeyword]
        if self.options.verbose:
             print 'Ngroups ',Ngroups
        tsubik = txttableclass()
        #tsubik_subgm1 = txttableclass()
        tsubg0 = txttableclass()
        tdsub1 = txttableclass()
        tdsub2 = txttableclass()
        tdsubij = txttableclass()
        initcols([tsubik,tsubg0,tdsub1,tdsub2,tdsubij])

        if self.options.gmax == None:
            gmax = Ngroups
        else:
            gmax = min(self.options.gmax,Ngroups)

        self.make_refpixcorr(data=self.data1,gmax=gmax,skiprefpixcorr=skiprefpixcorr,refpixboxsize=refpixboxsize,skipcorrect4framereset=skipcorrect4framereset)
        self.make_refpixcorr(data=self.data2,gmax=gmax,skiprefpixcorr=skiprefpixcorr,refpixboxsize=refpixboxsize,skipcorrect4framereset=skipcorrect4framereset)
        #print 'VVVVVVVVV'
        #pyfits.writeto('test.finaldiffim.fits',self.data1[1:,:,:]-self.data1[0,:,:],clobber=True)
        #sys.exit(0)

        #if not skipcorrect4framereset:
        #    sigmacutrefpixbottom = calcaverageclass()
        #    sigmacutrefpixtop = calcaverageclass()
        #    print 'Frame reset correction!'
        #    for g in xrange(self.options.gmin,gmax):
        #        if self.options.verbose:
        #            print 'g',g

        #        for amp in xrange(1,5):
        #            xmin4amp = (amp-1)*512
        #            xmax4amp = amp*512
        #            
        #            # image 1
        #            refpixbottom = self.data1[g,:4,xmin4amp:xmax4amp]-self.data1[g-1,:4,xmin4amp:xmax4amp]
        #            sigmacutrefpixbottom.calcaverage_sigmacutloop(refpixbottom,Nsigma=3,verbose=3)
        #            if sigmacutrefpixbottom.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"

        #            #refpixtop = self.data1[g,:4,xmin4amp:xmax4amp]-self.data1[g-1,:4,xmin4amp:xmax4amp]
        #            refpixtop = self.data1[g,-4:,xmin4amp:xmax4amp]-self.data1[g-1,-4:,xmin4amp:xmax4amp]
        #            sigmacutrefpixtop.calcaverage_sigmacutloop(refpixtop,Nsigma=3,verbose=3)
        #            if sigmacutrefpixtop.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"

        #            frameresetcorrection = 0.5*(sigmacutrefpixbottom.mean+sigmacutrefpixtop.mean)
        #            self.data1[g:,:,xmin4amp:xmax4amp] -= frameresetcorrection
                    
        #            # image 2
        #            refpixbottom = self.data2[g,:4,xmin4amp:xmax4amp]-self.data2[g-1,:4,xmin4amp:xmax4amp]
        #            sigmacutrefpixbottom.calcaverage_sigmacutloop(refpixbottom,Nsigma=3,verbose=3)
        #            if sigmacutrefpixbottom.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"
                    
        #            #refpixtop = self.data2[g,:4,xmin4amp:xmax4amp]-self.data2[g-1,:4,xmin4amp:xmax4amp]
        #            refpixtop = self.data2[g,-4:,xmin4amp:xmax4amp]-self.data2[g-1,-4:,xmin4amp:xmax4amp]
        #            sigmacutrefpixtop.calcaverage_sigmacutloop(refpixtop,Nsigma=3,verbose=3)
        #            if sigmacutrefpixtop.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"
                    
        #            frameresetcorrection = 0.5*(sigmacutrefpixbottom.mean+sigmacutrefpixtop.mean)
        #            self.data2[g:,:,xmin4amp:xmax4amp] -= frameresetcorrection

                
        print '\n##### Looping through boxes, calculating stdevs....'
        for g in xrange(self.options.gmin,gmax):
            if self.options.verbose:
                print 'g',g

            print 'calculating diffims'
            #diffim =  self.data2[g,:,:]-self.data1[g,:,:]
            diffim =  self.data1[g,:,:]-self.data2[g,:,:]

            if g>0:
                diffim_subg0 =  (self.data1[g,:,:]-self.data1[0,:,:])-(self.data2[g,:,:]-self.data2[0,:,:])
                # calculate the average count between g and g-1
                diffim_dsubij =  (self.data1[g,:,:]-self.data1[g-1,:,:])-(self.data2[g,:,:]-self.data2[g-1,:,:])
                dS_dsubij = 0.5*(self.data1[g,:,:]-self.data1[g-1,:,:]+self.data2[g,:,:]-self.data2[g-1,:,:])
                S_dsubij = 0.25*(self.data1[g,:,:]+self.data1[g-1,:,:]+self.data2[g,:,:]+self.data2[g-1,:,:]-2.0*self.data1[0,:,:]-2.0*self.data2[0,:,:])

            if g>0 and g<Ngroups-1:
                diffim_dsub1 =  (self.data1[g+1,:,:]-self.data1[g,:,:])-(self.data1[g,:,:]-self.data1[g-1,:,:])
                dS_dsub1 = 0.5*(self.data1[g+1,:,:]-self.data1[g-1,:,:])

            if g>1 and g<Ngroups-1:
                diffim_dsub2 =  (self.data1[g+1,:,:]-self.data1[g,:,:])-(self.data1[g-1,:,:]-self.data1[g-2,:,:])
                dS_dsub2 = 1.0/3.0*(self.data1[g+1,:,:]-self.data1[g-2,:,:])
                S_dsub2 = 0.5*(self.data1[g,:,:]-self.data1[0,:,:]+self.data1[g-1,:,:]-self.data1[0,:,:])

            #im = 0.5*(self.data1[g,:,:]-self.data1[0,:,:]+self.data2[g,:,:]-self.data2[0,:,:])
            #im_subgm1 = 0.5*(self.data1[g,:,:]-self.data1[0,:,:]+self.data1[g-1,:,:]-self.data1[0,:,:])
            S = self.data1[g,:,:]-self.data1[0,:,:]
            #dS = self.data1[g,:,:]-self.data1[g-1,:,:]

            self.calcgain_grid(diffim,(xmin,xmax,boxsize),(ymin,ymax,boxsize),tsubik,g,mask=self.bpmmask,S=S,imtype="default")

            if g>0:
                self.calcgain_grid(diffim_subg0,(xmin,xmax,boxsize),(ymin,ymax,boxsize),tsubg0,g,mask=self.bpmmask,S=S,imtype="subg0")
                self.calcgain_grid(diffim_dsubij,(xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsubij,g,mask=self.bpmmask,S=S_dsubij,dS=dS_dsubij,imtype="dsubij")
            if g>0 and g<Ngroups-1:
                self.calcgain_grid(diffim_dsub1,(xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsub1,g,mask=self.bpmmask,S=S,dS=dS_dsub1,imtype="dsub1")
            if g>1 and g<Ngroups-1:
                self.calcgain_grid(diffim_dsub2,(xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsub2,g,mask=self.bpmmask,S=S_dsub2,dS=dS_dsub2,imtype="dsub2")

            if 0==1:
                outname = self.outfilebasename+'.%02d.diff.fits' % g
                print 'Saving ',outname
                rmfile(outname)
                pyfits.writeto(outname,diffim,self.hdr1)

                outname = self.outfilebasename+'.%02d.diff.subg0.fits' % g
                print 'Saving ',outname
                rmfile(outname)
                pyfits.writeto(outname,diffim_subg0,self.hdr1)

                #outname = self.outfilebasename+'.%02d.diff.subgm1.fits' % g
                #print 'Saving ',outname
                #rmfile(outname)
                #pyfits.writeto(outname,diffim_subgm1,self.hdr1)

            print 'cleaning up'
            del(diffim,S)
            if g>0:del(diffim_subg0,diffim_dsubij,dS_dsubij,S_dsubij)
            if g>0 and g<Ngroups-1:del(diffim_dsub1,dS_dsub1)
            if g>1 and g<Ngroups-1:del(diffim_dsub2,dS_dsub2,S_dsub2)

        self.correct_stdev2_with_rdnoiseterms((xmin,xmax,boxsize),(ymin,ymax,boxsize),tsubik,tsubg0,tdsub1,tdsub2,tdsubij,Smean_max=self.options.Smeanmax4fit)
            
        # calculate the effective gain for each group
        self.calceffgain_dS(tdsub1)
        self.calceffgain_dS(tdsub2)
        self.calceffgain_dS(tdsubij)

        self.calcgain((xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsub1,'dsub1',Smean_max=self.options.Smeanmax4fit)
        self.calcgain((xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsub2,'dsub2',Smean_max=self.options.Smeanmax4fit)
        self.calcgain((xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsubij,'dsubij',Smean_max=self.options.Smeanmax4fit)

        self.calcstdev2_model(tdsub1,self.options.gmin,gmax,g0=1)
        self.calcstdev2_model(tdsub2,self.options.gmin,gmax,g0=2)
        self.calcstdev2_model(tdsubij,self.options.gmin,gmax,g0=1)

        self.calcgainmodel(tdsub1,(xmin,xmax,boxsize),(ymin,ymax,boxsize))
        self.calcgainmodel(tdsub2,(xmin,xmax,boxsize),(ymin,ymax,boxsize),g0=2)
        self.calcgainmodel(tdsubij,(xmin,xmax,boxsize),(ymin,ymax,boxsize))

        self.calc_cum_stdev2(tdsub1,(xmin,xmax,boxsize),(ymin,ymax,boxsize),self.options.gmin,gmax)
        self.calc_cum_stdev2(tdsub2,(xmin,xmax,boxsize),(ymin,ymax,boxsize),self.options.gmin,gmax)
        self.calc_cum_stdev2(tdsubij,(xmin,xmax,boxsize),(ymin,ymax,boxsize),self.options.gmin,gmax)

        tsubik.save2file('%s.all.subij' % (self.outfilebasename),verbose=True)
        tsubg0.save2file('%s.all.subg0' % (self.outfilebasename),verbose=True)
        tdsub1.save2file('%s.all.dsub1' % (self.outfilebasename),verbose=True)
        tdsub2.save2file('%s.all.dsub2' % (self.outfilebasename),verbose=True)
        tdsubij.save2file('%s.all.dsubij' % (self.outfilebasename),verbose=True)

if __name__=='__main__':

    usagestring='USAGE: calcgainim.py flatfilename1 flatfilename2'

    gainim=gainimclass()
    parser = gainim.add_options(usage=usagestring)
    gainim.options,  args = parser.parse_args()

    if len(args)!=2:
        sys.argv.append('--help')
        options,  args = parser.parse_args()
        sys.exit(0)

    (infile1,infile2) = args
    if re.search('DARK',infile1) or re.search('DARK',infile2):
        print 'DARK frame in infiles, exiting!!'
        print "SUCCESS calcgainim" # Say success here, makes the pipeline not barf...
        sys.exit(0)

    if gainim.options.darks!=(None,None):
        darkfilelist = [os.path.abspath(gainim.options.darks[0]),os.path.abspath(gainim.options.darks[1])]
        print 'Using the following dark frames:',darkfilelist
    elif gainim.options.autodark4readnoise:
        m = re.search('(\d+)_SE',os.path.basename(infile1))
        if m==None:
            raise RuntimeError,'file %s' % filename
        SCA = int(m.groups()[0])
        path = os.path.dirname(os.path.abspath(infile1))
        darkfilelist = glob.glob(path+'/*DARK*_%d_SE*.fits' % SCA)
        if len(darkfilelist)==0:
            raise RuntimeError,'Could not find dark frames *DARK*.fits in %s' % (path)
        print 'Found the following dark frames:',darkfilelist
    else:
        darkfilelist=None

    gainim.gainim(infile1,infile2,darks4readnoiselist=darkfilelist,
                  skiprefpixcorr=gainim.options.skiprefpixcorr,
                  refpixboxsize=gainim.options.refpixboxsize,
                  skipcorrect4framereset=gainim.options.skipcorrect4framereset)

    print "SUCCESS calcgainim"
