#!/usr/bin/env python


import sys, os,re,math,copy
import astropy.io.fits as pyfits
import optparse
import numpy as np
import scipy
from scipy.signal import convolve, boxcar
import scipy.ndimage as ndi
from tools import rmfile,makepath4file
from sigmacut import calcaverageclass
from texttable import txttableclass

class frameresetclass:
    def __init__(self, instrument='nircam', imagesize='fullimage'):

        self.debug = False
        self.verbose = 0
        self.sigmacut = calcaverageclass()

        if instrument=='nircam':
            if imagesize=='fullimage':
                self.nx = 2048
                self.ny = 2048
                self.nxamp = 512
                self.refpixdxleft=4
                self.refpixdxright=4
                self.refpixdytop=4
                self.refpixdybottom=4
                self.Namps = 4
                self.refleftindex = 0
                self.refrightindex = self.Namps+1
                self.refbottomindex = 0
                self.reftopindex=2
            else:
                raise RuntimeError,'imagesize!=fullimage not yet implemented, check __init__ of frameresetclass'
        else:
            raise RuntimeError,'instrument!=nircam not yet implemented, check __init__ of frameresetclass'

        # get the min/max for the different regions (amplifiers, overscan etc)
        self.getxyminmax4regions()

    def getxyminmax4regions(self):
        ### overscan left: index=self.refleftindex=0
        ### amp 1-4 is index 1-4
        ### overscan right: index=self.refrightindex=self.Namps+1=5
        ### overscan bottom: index=0
        ### amps: index=1
        ### overscan bottom: index=2

        self.xmins={}
        self.xmaxs={}

        # left overscan region
        if self.refpixdxleft>0:
            self.xmins[self.refleftindex]=0
            self.xmaxs[self.refleftindex]=self.refpixdxleft

        # right overscan region
        if self.refpixdxright>0:
            self.xmins[self.refrightindex]=self.nx-self.refpixdxright
            self.xmaxs[self.refrightindex]=self.nx

        for amp in xrange(1,self.Namps+1):

            xmin = (amp-1)*self.nxamp
            xmax = (amp)*self.nxamp

            # exclude overscan
            if xmin<self.refpixdxleft: xmin=self.refpixdxleft
            if xmax>self.nx-self.refpixdxright: xmax = self.nx-self.refpixdxright

            self.xmins[amp]=xmin
            self.xmaxs[amp]=xmax

        ### y limits

        self.ymins={}
        self.ymaxs={}

        # bottom overscan region
        if self.refpixdybottom>0:
            self.ymins[self.refbottomindex]=0
            self.ymaxs[self.refbottomindex]=self.refpixdybottom

        # regular region
        self.ymins[1]=0+self.refpixdybottom
        self.ymaxs[1]=self.ny-self.refpixdybottom

        # top overscan region
        if self.refpixdytop>0:
            self.ymins[self.reftopindex]=self.ny-self.refpixdytop
            self.ymaxs[self.reftopindex]=self.ny

    def save2Dvec_as_txt(self,vec,vecfilename,vecnoise=None,ymin=None,ymax=None,ycolname='y',xcolnameprefix='amp'):
        if ymin==None: ymin=0
        if ymax==None: ymax = vec.shape[0]
        rmfile(vecfilename)
        print 'Saving vector into',vecfilename
        f=open(vecfilename,'w')
        s = '#   %s' % ycolname
        for x in xrange(vec.shape[1]):
            s+=' %6s%d' % (xcolnameprefix,x)
            if vecnoise!=None:
                s+=' %3s%derr' % (xcolnameprefix,x)
        f.write(s+'\n')

        for y in xrange(ymin,ymax):
            s = ' %4d' % y
            for x in xrange(vec.shape[1]):
                s+= ' %7.3f' % vec[y,x]
                if vecnoise!=None:
                    s+=' %7.3f' %  vecnoise[y,x]
            f.write(s+'\n')
        f.close()

    def save2Dvec_as_txt(self,vec,vecfilename,vecnoise=None,ymin=None,ymax=None,ycolname='y',xcolnameprefix='amp'):
        if ymin==None: ymin=0
        if ymax==None: ymax = vec.shape[0]
        rmfile(vecfilename)
        print 'Saving vector into',vecfilename
        f=open(vecfilename,'w')
        s = '#   %s' % ycolname
        for x in xrange(vec.shape[1]):
            s+=' %6s%d' % (xcolnameprefix,x)
            if vecnoise!=None:
                s+=' %3s%derr' % (xcolnameprefix,x)
        f.write(s+'\n')

        for y in xrange(ymin,ymax):
            s = ' %4d' % y
            for x in xrange(vec.shape[1]):
                s+= ' %7.3f' % vec[y,x]
                if vecnoise!=None:
                    s+=' %7.3f' %  vecnoise[y,x]
            f.write(s+'\n')
        f.close()

    def calcframereset(self,im,xmin,xmax,method='first'):
        ### calculate the offset between the different groups of a ramp
        if len(im.shape)!=3:
            print 'Shape: ',im.shape
            raise RuntimeError,'3-dim image required!'

        itop    = self.reftopindex
        ibottom = self.refbottomindex
        iaverage=1

        # vec[g,i]:
        # bottom: i=0=self.reftopindex
        # top: i=2=self.refbottomindex
        # average of top,bottom: i=1
        vec = scipy.zeros((im.shape[0],3),dtype=float)
        vecerr = scipy.zeros((im.shape[0],3),dtype=float)

        framereset = scipy.zeros((im.shape[0],3),dtype=float)
        dframereset = scipy.zeros((im.shape[0],3),dtype=float)

        # calculate the average differences for top and bottom
        for i in [ibottom,itop]:
            if self.verbose>1:
                if i==ibottom: print '### bottom ref pixels'
                else:  print '### top ref pixels'
            if method == 'next':
                gmax = im.shape[0]
                diffim = im[1:gmax,self.ymins[i]:self.ymaxs[i],xmin:xmax]-im[0:gmax-1,self.ymins[i]:self.ymaxs[i],xmin:xmax]
            elif method == 'first':
                gmax = im.shape[0]
                #for g in xrange(1,im.shape[0]):
                diffim = im[1:gmax,self.ymins[i]:self.ymaxs[i],xmin:xmax]-im[0,self.ymins[i]:self.ymaxs[i],xmin:xmax]
            else:
                raise RuntimeError, 'method %s not yet implemented!' % method

            for g in xrange(diffim.shape[0]):
                self.sigmacut.calcaverage_sigmacutloop(diffim[g,:,:],median_firstiteration=True,verbose=(self.verbose>2),Nsigma=3.0)
                vec[g+1,i]=self.sigmacut.mean
                vecerr[g+1,i]=self.sigmacut.mean_err

        # calculate the average bias offset
        for i in [ibottom,itop]:
            if method == 'next':
                vecsum = 0
                vecerrsum2 = 0
                for g in xrange(1,vec.shape[0]):
                    framereset[g,i] = (vec[g,i] + vecsum)
                    dframereset[g,i] = math.sqrt(vecerr[g,i]*vecerr[g,i]+vecerrsum2)
                    vecsum += vec[g,i]
                    vecerrsum2 += vecerr[g,i]*vecerr[g,i]
            elif method == 'first':
                for g in xrange(1,vec.shape[0]):
                    framereset[g,i] = vec[g,i]
                    dframereset[g,i] = vecerr[g,i]
            else:
                raise RuntimeError, 'method %s not yet implemented!' % method

        # calculate average offset
        framereset[1:vec.shape[0],1] = 0.5*(framereset[1:vec.shape[0],ibottom]+framereset[1:vec.shape[0],itop])
        # make that the offsets are on average zero
        mean = scipy.mean(framereset[:,1])
        framereset[:,:] -= mean

        #for g in xrange(1,vec.shape[0]):
        #    dframereset[g,1] = 0.5 * math.sqrt(dframereset[g,ibottom]*dframereset[g,ibottom]+dframereset[g,itop]*dframereset[g,itop])
        # faster and as good:
        meanerror_bottom = scipy.mean(dframereset[:,ibottom])
        meanerror_top = scipy.mean(dframereset[:,itop])
        dframereset[1:vec.shape[0],1] = 0.5 * math.sqrt(meanerror_bottom**2+meanerror_top**2)

        return(framereset,dframereset)

    def calcframereset_allamps(self,im,method='first'):
        ### framereset(g,a), a=1-4 is amp 1-4: average of bottom and top reference pixels
        ### bottom_deltaoffset,reftop_deltaoffset:difference between the ref bottom/top to average

        framereset = scipy.zeros((im.shape[0],self.Namps+1), dtype=float)
        dframereset = scipy.zeros((im.shape[0],self.Namps+1), dtype=float)

        # difference between the ref bottom to average
        refbottom_deltaoffset = scipy.zeros((im.shape[0],self.Namps+1), dtype=float)
        # difference between the ref top to average
        reftop_deltaoffset = scipy.zeros((im.shape[0],self.Namps+1), dtype=float)

        for amp in xrange(1,self.Namps+1):
            if self.verbose:
                print '########### Calculating frame reset for amp %d, method %s' % (amp,method)
            (framereset_amp,dframereset_amp) = self.calcframereset(im,self.xmins[amp],self.xmaxs[amp],method=method)
            #framereset_amp[:,0]: ref bottom
            #framereset_amp[:,1]: ref average(top,bottom)
            #framereset_amp[:,2]: ref top
            framereset[:,amp]=framereset_amp[:,1]
            dframereset[:,amp]=dframereset_amp[:,1]
            refbottom_deltaoffset[:,amp] = framereset_amp[:,self.refbottomindex]-framereset_amp[:,1]
            reftop_deltaoffset[:,amp] = framereset_amp[:,self.reftopindex]-framereset_amp[:,1]
            if self.debug: self.save2Dvec_as_txt(framereset_amp,'framereset.amp%d.debug.txt' % amp,vecnoise=dframereset_amp,ycolname='g',xcolnameprefix='')
        if self.debug:
            self.save2Dvec_as_txt(framereset,'framereset.debug.txt',vecnoise=dframereset,ycolname='g',xcolnameprefix='')
            self.save2Dvec_as_txt(reftop_deltaoffset,'framereset.top.debug.txt',vecnoise=dframereset,ycolname='g',xcolnameprefix='')
            self.save2Dvec_as_txt(refbottom_deltaoffset,'framereset.bottom.debug.txt',vecnoise=dframereset,ycolname='g',xcolnameprefix='')

        return(framereset,dframereset,refbottom_deltaoffset,reftop_deltaoffset)

    def subtractframeresetoffset(self,im,framereset):
        for amp in xrange(1,self.Namps+1):
            xmin = self.xmins[amp]
            xmax = self.xmaxs[amp]
            if amp==1: # include left overscan
                xmin = self.xmins[amp-1]
            if amp==self.Namps:
                xmax = self.xmaxs[amp+1]
            ymin = 0
            ymax = self.ny
            if self.verbose: print "subtracting frame reset for amp %d: x=%d-%d, y=%d-%d" % (amp,xmin,xmax,ymin,ymax)
            for g in xrange(im.shape[0]):
                im[g,ymin:ymax,xmin:xmax] -= framereset[g,amp]
                #im[0:im.shape[0],ymin:ymax,xmin:xmax] -= framereset[0:im.shape[0],amp]

    def correct4topbottomdifference(self,im,refbottom_deltaoffset,reftop_deltaoffset):
        for amp in xrange(1,self.Namps+1):
            xmin = self.xmins[amp]
            xmax = self.xmaxs[amp]
            if amp==1: # include left overscan
                xmin = self.xmins[amp-1]
            if amp==self.Namps:
                xmax = self.xmaxs[amp+1]
            ymin = 0
            ymax = self.ny
            if self.verbose: print "Correcting for differences between top and bottom ref pixels for amp %d: x=%d-%d, y=%d-%d" % (amp,xmin,xmax,ymin,ymax)
            yrange = xrange(ymin,ymax)
            for g in xrange(im.shape[0]):
                print 'nnn',reftop_deltaoffset[g,amp],refbottom_deltaoffset[g,amp]
                if self.verbose>1 and (g % 10) == 0: print "group %d" % g
                line_slope = (reftop_deltaoffset[g,amp]-refbottom_deltaoffset[g,amp])/(ymax-ymin)
                line_offset = refbottom_deltaoffset[g,amp] - line_slope*ymin
                for y in yrange:
                    im[g,y:y+1,xmin:xmax] -= line_slope*y+line_offset
                    #print 'VVV',y,line_slope*y+line_offset
                #im[0:im.shape[0],ymin:ymax,xmin:xmax] -= framereset[0:im.shape[0],amp]
                #if g==1:sys.exit(0)

    def correct4framereset(self,im,method='first',correct4topbottomdifference=False):
        # get the framereset vector
        (framereset,dframereset,refbottom_deltaoffset,reftop_deltaoffset) = self.calcframereset_allamps(im,method=method)

        # subtract the framereset vector
        self.subtractframeresetoffset(im,framereset)

        # make a linear correction for the differences between top and bottom
        if correct4topbottomdifference:
            self.correct4topbottomdifference(im,refbottom_deltaoffset,reftop_deltaoffset)

        return(im)

class refpixcorrclass(frameresetclass):
    def __init__(self):
        frameresetclass.__init__(self)
        self.filename = None
        self.sigmacut = calcaverageclass()

        self.data=None
        self.hdr=None
        self.data_is_refpixcorrected = False

        self.noisedata=None
        self.noisehdr=None

        self.maskdata=None
        self.maskhdr=None

        self.verbose=False

        self.gmin=None
        self.gmax=None
        self.ymin=None
        self.ymax=None

        self.biasdriftvec={}
        self.biasdriftvec_smoothed={}
        self.biasdriftvec_refpixcorrleft={}
        self.biasdriftvec_refpixcorrright={}

        self.cols4vectxt = ['sci','err','stdev','X2norm']

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_option('--debug', default=False, action="store_true",
                          help='debug mode: more output and debug files')
        parser.add_option('-f','--force', default=False, action="store_true",
                          help='if output file already exists, still force to do it...')
        parser.add_option('--gmin'  , default=None, type="int",
                          help='minimum group (default=%default)')
        parser.add_option('--gmax'  , default=None, type="int",
                          help='maximum group (default=%default)')
        parser.add_option('--ymin'  , default=None, type="int",
                          help='minimum y (default=%default)')
        parser.add_option('--ymax'  , default=None, type="int",
                          help='maximum y (default=%default)')
        parser.add_option('--bpm'  , default=None, type="string",
                          help='bad pixel mask filename (default=%default)')
        parser.add_option('--noise'  , default=None, type="string",
                          help='noise filename (default=%default)')
        parser.add_option('-b','--boxsize'  , default=None, type="int",
                          help='boxsize for boxcar smoothing of biasdrift vector (default=%default)')
        parser.add_option('-s','--save', default=False, action="store_true",
                          help='Save the corrected image')
        parser.add_option('-o','--outfilebasename'  , default='auto' , type="string",
                          help='file basename of output file. If \'auto\', then basename is input filename with fits removed (default=%default)')
        parser.add_option('--outsubdir'  , default=None , type="string",
                          help='if specified gets added to output directory (default=%default)')
        parser.add_option('--outsuffix'  , default=None , type="string",
                          help='if specified gets added to output filename (default=%default)')
        parser.add_option('--skipcorrect4framereset', default=False, action="store_true",
                          help='Skip the correction for the bias offset (mainly for debu to see the effect of the bias offset!)')
        parser.add_option('--savediffim', default=False, action="store_true",
                          help='Save the diffim (for debugging) as *.diffim.fits')
        parser.add_option('--savebiasdriftvec', default=False, action="store_true",
                          help='Save the biasdrift vec as fits files')
        parser.add_option('--savebiasdriftvec_as_txt', default=False, action="store_true",
                          help='Save the biasdrift vec as *.g?.biasdrift.txt ASCII file')
        parser.add_option('-p','--showplot', default=False, action="store_true",
                          help='show the linear fits')
        parser.add_option('-t','--test', default=False, action="store_true",
                          help='test: instead of left ref pix, take the 4 pixels just next to it')
        return(parser)

    def mkoutfilebasename(self,filename, outfilebasename='auto',outsuffix=None,outsubdir=None,testflag=False):

        if  outfilebasename.lower() == 'auto':
            outfilebasename = re.sub('\.fits$','',filename)
            if outfilebasename ==filename:
                raise RuntimeError,'BUG!!! %s=%s' % (outfilebasename,filename)
        if outsuffix!=None:
            outfilebasename += '.'+outsuffix
        if outsubdir!=None:
            (d,f)=os.path.split(outfilebasename)
            outfilebasename = os.path.join(d,outsubdir,f)
        if testflag:
            outfilebasename += '.test'
        self.outfilebasename = outfilebasename

        #makepath4file(self.outfilebasename)
        return(0)

    def biasdriftvecfilename(self,filename=None):
        if filename == None:
            filename = self.outfilebasename + '.biasdrift.fits'
        return(filename)

    def loadimage(self, filename, bpmfilename=None, noisefilename=None):
        self.data_is_refpixcorrected = False

        # load image
        if self.verbose: print 'Loading ',filename
        self.data, self.hdr = pyfits.getdata(filename, 0, header=True)
        if len(self.data.shape)==4:
            print 'SSB-conform image'
            newdata = scipy.zeros(self.data.shape[1:], dtype=float)
            newdata[:,:,:] = self.data[0,:,:,:]
            self.data = newdata

        self.filename=filename

        #load mask
        if bpmfilename!=None:
            if self.verbose: print 'Loading bpm',bpmfilename
            self.maskdata, self.maskhdr = pyfits.getdata(bpmfilename, 0, header=True)

        # load noise
        if noisefilename!=None:
            if self.verbose: print 'Loading noise',noisefilename
            self.noisedata, self.noisehdr = pyfits.getdata(noisefilename, 0, header=True)

    def mkdiffimcorr(self, savediffimflag=False):
        # calculate diffim
        print 'corrected image shape:',self.datacorr.shape

        gmax = self.datacorr.shape[0]
        self.diffimcorr = self.datacorr[1:gmax,:,:]-self.datacorr[0,:,:]

        # save diffim if wanted
        if savediffimflag:
            diffimoutfilename = self.outfilebasename + '.diffimcorr.fits'
            rmfile(diffimoutfilename)
            print 'Saving corrected diffim',diffimoutfilename
            pyfits.writeto(diffimoutfilename,self.diffimcorr,self.hdr)

    def mkdiffim(self, savediffimflag=False, outsuffix='auto'):
        # calculate diffim
        print 'image shape:',self.data.shape

        gmax = self.data.shape[0]
        self.diffim = self.data[1:gmax,:,:]-self.data[0,:,:]

        # save diffim if wanted
        if savediffimflag:
            if outsuffix=='auto':
                diffimoutfilename = self.outfilebasename + '.diffim.fits'
            else:
                diffimoutfilename= self.outfilebasename + outsuffix
            rmfile(diffimoutfilename)
            print 'Saving diffim',diffimoutfilename
            pyfits.writeto(diffimoutfilename,self.diffim,self.hdr)

    def initbiasdriftvec(self,diffimshape):
        vec = {}
        for exttype in ['sci','err','stdev','X2norm']:
            vec[exttype] = scipy.zeros((diffimshape[0],diffimshape[1],len(self.xmins)), dtype=float)
        for exttype in ['Nused','Nclipped']:
            vec[exttype] = scipy.zeros((diffimshape[0],diffimshape[1],len(self.xmins)), dtype=np.uint16)

        return(vec)

    def calcbiasdriftvec(self,im,xmin,xmax,ymin,ymax,noise=None,mask=None,median_firstiteration=True,Nsigma=3.0):
        ### calculate the average for each row
        if len(im.shape)!=2:
            print 'Shape: ',im.shape
            raise RuntimeError,'2-dim image required!'

        vec = {}
        for exttype in ['sci','err','stdev','X2norm']:
            vec[exttype] = scipy.zeros((ymax-ymin),dtype=float)
        for exttype in ['Nused','Nclipped']:
            vec[exttype] = scipy.zeros((ymax-ymin), dtype=np.uint16)
 
        for y in xrange(ymin,ymax):
            noiserow = maskrow = None
            if noise!=None:
                noiserow = noise[y:y+1,xmin:xmax]
            if mask!=None:
                maskrow = mask[y:y+1,xmin:xmax]
                if noise!=None:
                    # skip pixels without noise
                    maskrow = scipy.where(noiserow == 0,4,maskrow)
            #print noiserow
            self.sigmacut.calcaverage_sigmacutloop(im[y:y+1,xmin:xmax],mask=maskrow,noise=noiserow,median_firstiteration=median_firstiteration,verbose=0,Nsigma=Nsigma)
            vec['sci'][y-ymin]=self.sigmacut.mean
            vec['err'][y-ymin]=self.sigmacut.mean_err
            vec['stdev'][y-ymin]=self.sigmacut.stdev
            vec['X2norm'][y-ymin]=self.sigmacut.X2norm
            vec['Nused'][y-ymin]=self.sigmacut.Nused
            vec['Nclipped'][y-ymin]=self.sigmacut.Nskipped
        return(vec)

    def getbiasdriftvecfilename(self,filename=None,txtflag=False,extrasuffix=None,g=None):
        if filename == None:
            filename = self.outfilebasename
            if g!=None:
                filename += '.g%d' % g
            if extrasuffix!=None:
                if not re.match('^\.',extrasuffix):
                    filename += '.'
                filename += extrasuffix
            filename += '.biasdrift'
            if txtflag:
                filename += '.txt'
            else:
                filename += '.fits'
        return(filename)

    def savebiasdriftvec_as_txt(self,vec,vecfilename,gmin=None,gmax=None,ymin=None,ymax=None,zcolname='g',ycolname='y',xcolnameprefix='amp',cols4vectxt=None):
        def colname(x,coltype,Nx):
            if x in [0,Nx-1]:
                if x == 0:
                    prefix = 'ref0'
                else:
                    prefix = 'ref1'
            else:
                prefix='amp%d' % (x-1)
            colname = '%4s%s' % (prefix,coltype)
            return(colname)

        if cols4vectxt == None:
            cols4vectxt=self.cols4vectxt

        if len(vec['sci'].shape)!=3:
            raise RuntimeError,"wrong shape of %d, shape=3 expected" % (vec['sci'].shape)
        if ymin==None: ymin=0
        if ymax==None: ymax = vec['sci'].shape[1]
        if gmin==None: gmin=0
        if gmax==None: gmax = vec['sci'].shape[0]

        rmfile(vecfilename)
        print 'Saving vector into',vecfilename
        f=open(vecfilename,'w')

        Nx = vec['sci'].shape[2]

        # make header
        s = '#  %s    %s' % (zcolname,ycolname)
        for x in xrange(Nx):
            for coltype in cols4vectxt:
                s+=' %10s' % colname(x,coltype,Nx)
        f.write(s+'\n')

        for g in xrange(gmin,gmax):
            for y in xrange(ymin,ymax):
                s = ' %3d %4d' % (g,y)
                for x in xrange(Nx):
                    for coltype in cols4vectxt:
                        if np.isnan(vec[coltype][g,y,x]):
                            s+= ' %10s' % '-'
                        else:
                            s+= ' %10.4f' % vec[coltype][g,y,x]
                f.write(s+'\n')
        f.close()

    def savebiasdriftvecs_as_txt(self,boxsize=None,save_g_separately=False):
        
        self.savebiasdriftvec_as_txt(self.biasdriftvec,self.getbiasdriftvecfilename(txtflag=True))
        if boxsize!=None:
            self.savebiasdriftvec_as_txt(self.biasdriftvec_smoothed,self.getbiasdriftvecfilename(txtflag=True,extrasuffix='smoothed'),cols4vectxt=['sci'])
            self.savebiasdriftvec_as_txt(self.biasdriftvec_refpixcorrleft,self.getbiasdriftvecfilename(txtflag=True,extrasuffix='refpixcorrleft'),cols4vectxt=['sci'])
            self.savebiasdriftvec_as_txt(self.biasdriftvec_refpixcorrright,self.getbiasdriftvecfilename(txtflag=True,extrasuffix='refpixcorrright'),cols4vectxt=['sci'])

        if save_g_separately:
            for g in xrange(gmin,gmax):
                self.savebiasdriftvec_as_txt(self.biasdriftvec,self.getbiasdriftvecfilename(txtflag=True,g=g),gmin=g,gmax=g+1)
                if boxsize!=None:
                    self.savebiasdriftvec_as_txt(self.biasdriftvec_smoothed,self.getbiasdriftvecfilename(txtflag=True,extrasuffix='smoothed',g=g),gmin=g,gmax=g+1,cols4vectxt=['sci'])
                    self.savebiasdriftvec_as_txt(self.biasdriftvec_refpixcorrleft,self.getbiasdriftvecfilename(txtflag=True,extrasuffix='refpixcorrleft',g=g),gmin=g,gmax=g+1,cols4vectxt=['sci'])
                    self.savebiasdriftvec_as_txt(self.biasdriftvec_refpixcorrright,self.getbiasdriftvecfilename(txtflag=True,extrasuffix='refpixcorrright',g=g),gmin=g,gmax=g+1,cols4vectxt=['sci'])
                    
        return(0)

    def savebiasdriftvec(self,vec,biasdriftvecfilename):
        # make the hdu with all extensions
        phdu = pyfits.PrimaryHDU(data=None,header=self.hdr)
        hdulist = [phdu]
        for exttype in ['sci','err','stdev','X2norm','Nused','Nclipped']:
            hdulist.append(pyfits.ImageHDU(vec[exttype], name=exttype))
        output = pyfits.HDUList(hdulist)

        # save it....
        print 'Saving biasdriftvec',biasdriftvecfilename
        output.writeto(biasdriftvecfilename,clobber=True)

    def savebiasdriftvecs(self,boxsize=None):
        self.savebiasdriftvec(self.biasdriftvec,self.biasdriftvecfilename())
        if boxsize!=None:
            self.savebiasdriftvec(self.biasdriftvec_smoothed,self.getbiasdriftvecfilename(extrasuffix='smoothed'))
            self.savebiasdriftvec(self.biasdriftvec_refpixcorrleft,self.getbiasdriftvecfilename(extrasuffix='refpixcorrleft'))
            self.savebiasdriftvec(self.biasdriftvec_refpixcorrright,self.getbiasdriftvecfilename(extrasuffix='refpixcorrright'))

    def smoothboxcar(self,vec,boxsize,mode='nearest'):
        b = boxcar(boxsize)
        vec_smooth = ndi.convolve(vec['sci'],b,mode=mode)/boxsize
#        print b,vec['sci'],vec_smooth
        return vec_smooth

    def calcbiasdriftvec_allgroups(self,gmin=None,gmax=None,ymin=None,ymax=None,boxsize=None,savebiasdriftvec_as_txt=False):
        print 'SHAPE:',self.diffim.shape

        if gmin == None: gmin=0
        if gmax == None: gmax=self.diffim.shape[0]

        #if ymin == None: ymin=0+self.refpixdybottom
        #if ymax == None: ymax=  self.diffim.shape[1]-self.refpixdytop
        #include the top/bottom refpix in this!!

        if ymin == None: ymin=0
        if ymax == None: ymax=self.diffim.shape[1]

        self.gmin=gmin
        self.gmax=gmax
        self.ymin=ymin
        self.ymax=ymax

        # initialize vec
        self.biasdriftvec = self.initbiasdriftvec(self.diffim.shape)
        
        if boxsize!=None:
            self.biasdriftvec_smoothed=self.initbiasdriftvec(self.diffim.shape)
            self.biasdriftvec_refpixcorrleft=self.initbiasdriftvec(self.diffim.shape)
            self.biasdriftvec_refpixcorrright=self.initbiasdriftvec(self.diffim.shape)

        if self.verbose:
            '### calculating refpix vector ###'

        # loop through groups
        for g in xrange(gmin,gmax):
            if self.verbose:
                print 'g: ',g
            # loop through regions: i=0: left overscan, 1-4: amps 1-4, i=5: right overscan
            for i in [self.refleftindex,self.refrightindex]:
                if self.verbose>1:
                    print 'x=%d-%d' % (self.xmins[i],self.xmaxs[i])

                median_firstiteration = True
                Nsigma=3.0
                # if overscan region, don't do Nsigma cut.
                if i in [self.refleftindex,self.refrightindex]:
                    median_firstiteration = False
                    Nsigma=0.0

                (vec) = self.calcbiasdriftvec(self.diffim[g,:,:],self.xmins[i],self.xmaxs[i],ymin,ymax,
                                              noise=self.noisedata,mask=self.maskdata,
                                              median_firstiteration=median_firstiteration, Nsigma=Nsigma)

                # save the result in the different extensions
                for exttype in ['sci','err','stdev','X2norm','Nused','Nclipped']:
                    self.biasdriftvec[exttype][g,ymin:ymax,i]=vec[exttype]

                # smooth the bioasdriftvec if wanted
                if boxsize!=None:
                    self.biasdriftvec_smoothed['sci'][g,ymin:ymax,i]= self.smoothboxcar(vec,boxsize)
                    #vec_smoothed = self.smoothboxcar(vec,boxsize)
                    #for exttype in ['sci',rr','stdev','X2norm','Nused','Nclipped']:
                    #    self.biasdriftvec_smoothed[exttype][g,ymin:ymax,i]=vec_smoothed[exttype]

            if boxsize!=None:
                for i in [self.refleftindex,self.refrightindex]:
                    self.biasdriftvec_refpixcorrleft['sci'][g,ymin:ymax,i]  = self.biasdriftvec['sci'][g,ymin:ymax,i] - self.biasdriftvec_smoothed['sci'][g,ymin:ymax,self.refleftindex]
                    self.biasdriftvec_refpixcorrright['sci'][g,ymin:ymax,i] = self.biasdriftvec['sci'][g,ymin:ymax,i] - self.biasdriftvec_smoothed['sci'][g,ymin:ymax,self.refrightindex]

    def correctimage4leftrightrefpix(self,boxsize=None,testflag=False):
        self.data_is_refpixcorrected = False
        if self.verbose:
            print '### Appplying ref pix correction ###'

        #self.datacorr = copy.deepcopy(self.data)
        corrim = np.zeros((self.data.shape[0]-1,self.data.shape[1],self.data.shape[2]),dtype=np.double)

        if boxsize!=None:
            vec = self.biasdriftvec_smoothed
        else:
            vec = self.biasdriftvec


        vecleft = vec['sci'][:,self.ymin:self.ymax,self.refleftindex:self.refleftindex+1]
        vecright = vec['sci'][:,self.ymin:self.ymax,self.refrightindex:self.refrightindex+1]
        # middle amps: calculate the average between left and right
        vecmiddle = 0.5*(vecleft+vecright)

        # create the correction image
        corrim[:,:,0:512]= vecleft.repeat(512).reshape((vecleft.shape[0],vecleft.shape[1],512))
        corrim[:,:,512:1536]= vecmiddle.repeat(1024).reshape((vecleft.shape[0],vecleft.shape[1],1024))
        corrim[:,:,1536:2048]= vecright.repeat(512).reshape((vecleft.shape[0],vecleft.shape[1],512))

        if testflag:
            pyfits.writeto('test.diffimbeforerefpix.subnext.fits',self.data[1:,:,:]-self.data[:-1,:,:],clobber=True)
            pyfits.writeto('test.diffimbeforerefpix.subg0.fits',self.data[1:,:,:]-self.data[0,:,:],clobber=True)
            pyfits.writeto('test.imbeforerefpix.fits',self.data,clobber=True)

        # subtract it!
        self.data[1:,:,:]-=corrim

        if testflag:
            pyfits.writeto('test.corrim.fits',corrim,clobber=True)
            pyfits.writeto('test.diffimafterrefpix.subnext.fits',self.data[1:,:,:]-self.data[:-1,:,:],clobber=True)
            pyfits.writeto('test.diffimafterrefpix.subg0.fits',self.data[1:,:,:]-self.data[0,:,:],clobber=True)
            pyfits.writeto('test.imafterrefpix.fits',self.data,clobber=True)

        self.data_is_refpixcorrected = True
            
        return self.data      

    def savecorrectedimage(self,filename=None):
        if not self.data_is_refpixcorrected:
            raise RuntimeError,'Cannot save corrected image, doesn\'t exist yet!'
        if filename==None:
            filename = self.outfilebasename + '.refpixcorr.fits'
        rmfile(filename)
        print 'Saving corrected image',filename
        pyfits.writeto(filename,self.data,self.hdr,clobber=True)
        
    def mkframeresetcorrection(self,data,boxsize=10,savediffimflag=False,testflag=False):
        self.data_is_refpixcorrected = False
        self.data = data
        self.mkdiffim(savediffimflag=savediffimflag)
        self.calcbiasdriftvec_allgroups(boxsize=boxsize)
        self.correctimage4leftrightrefpix(boxsize=boxsize,testflag=testflag)
        return(self.data)

if __name__=='__main__':

    usagestring='USAGE: refpixcorr.py infile'

    refpixcorr=refpixcorrclass()
    parser = refpixcorr.add_options(usage=usagestring)
    options,  args = parser.parse_args()

    refpixcorr.verbose=options.verbose
    refpixcorr.debug=options.debug


    if options.test:
        refpixcorr.xmins[0]=4
        refpixcorr.xmaxs[0]=8


    for infile in args:

        refpixcorr.mkoutfilebasename(infile, outfilebasename=options.outfilebasename,outsuffix=options.outsuffix,outsubdir=options.outsubdir,testflag=options.test)

        refpixcorr_outputfilename = refpixcorr.outfilebasename + '.refpixcorr.fits'
        if not options.force and os.path.isfile(refpixcorr_outputfilename):
            print '%s already exist, skipping recreating it!! (use -f to force recreation...)' % refpixcorr_outputfilename
            continue

        biasdriftvecfilename = refpixcorr.biasdriftvecfilename()
        
        refpixcorr.loadimage(infile,bpmfilename=options.bpm,noisefilename=options.noise)

        refpixcorr.mkdiffim(savediffimflag=options.savediffim, outsuffix='.diffim0.fits')

        if not options.skipcorrect4framereset:
            refpixcorr.correct4framereset(refpixcorr.data)
        else:
            print 'WARNING: The correction for the bias offset is skipped!!'

        refpixcorr.mkdiffim(savediffimflag=options.savediffim)

        refpixcorr.calcbiasdriftvec_allgroups(gmin=options.gmin,gmax=options.gmax,
                                              ymin=options.ymin,ymax=options.ymax,
                                              boxsize=options.boxsize
                                              )

        if options.savebiasdriftvec_as_txt:
            refpixcorr.savebiasdriftvecs_as_txt(boxsize=options.boxsize)
        if options.savebiasdriftvec:
            refpixcorr.savebiasdriftvecs(boxsize=options.boxsize)

        refpixcorr.correctimage4leftrightrefpix(boxsize=options.boxsize)

        # redo frame reset correction
        if not options.skipcorrect4framereset:
            refpixcorr.correct4framereset(refpixcorr.datacorr)
        else:
            print 'WARNING: The correction for the bias offset is skipped!!'

        if options.savediffim:
            refpixcorr.mkdiffimcorr(savediffimflag=options.savediffim)

        if options.save:
            refpixcorr.savecorrectedimage()
            
