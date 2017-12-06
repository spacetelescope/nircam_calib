#!/usr/bin/env python

import re,sys,string,os,types

#from ConfigParser import SafeConfigParser


def executecommand(cmd,successword,errorlog=None,cmdlog=None,verbose=1):
    if verbose: print('executing: ',cmd)
    (cmd_in,cmd_out)=os.popen4(cmd)
    output = cmd_out.readlines()
    if successword=='':
        successflag = 1
    else:
        m = re.compile(successword)
        successflag = 0
        for line in output:
            if m.search(line):
                successflag = 1                       
    errorflag = not successflag
    if errorflag:
        print('error executing:',cmd)
        if errorlog != None:
            append2file(errorlog,['\n error executing: '+cmd+'\n'])
            append2file(errorlog,output)
        if cmdlog != None:
            append2file(cmdlog,['\n error executing:',cmd])
    else:
        if cmdlog != None:
            append2file(cmdlog,[cmd])        
    return errorflag, output

def makepath(path,raiseError=1):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=1):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)

def rmfile(filename,raiseError=1):
    " if file exists, remove it "
    if os.path.isfile(filename): 
        os.remove(filename)
        if os.path.isfile(filename): 
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    return(0)

# some little defs to compare lists and tuples
def unique(seq):
    if seq == None or seq == []: return []
    d = {}
    for x in seq:
        d[x] = 1
    return d.keys()

def AnotB(A,B):
    "returns elements that are in A, but not in B"
    if A == None: return []
    if not (type(A) in [types.ListType,types.TupleType]): A = [A,]
    if B == None: return A
    if not (type(B) in [types.ListType,types.TupleType]): B = [B,]
    c = {}
    d = {}
    for x in A:
        c[x] = 1
        d[x] = 1
    for x in B:
        if c.has_key(x):
            if d.has_key(x):            
                del(d[x])
    del c
    return d.keys()

# Extend generic ConfigParser so that environment variables are substituted
'''
class ConfigParser_env(SafeConfigParser):
    def __init__(self):
        SafeConfigParser.__init__(self)
        self.envvarpattern=re.compile('\$(\w+)')
        
    def getstring(self,section,paramname):
        s = self.subenvvarplaceholder(self.get(section,paramname))
        return(s)

    def subenvvarplaceholder(self,s):
        envvarnames=self.envvarpattern.findall(s)
        if envvarnames:
            for name in envvarnames:
                envval=os.environ[name]
                subpattern='\$%s' % (name)
                s = re.sub(subpattern,envval,s)
        return(s)

    def setval_nosection(self,option,val,allflag=False,throwerror=True):
        sections = self.sections()        
        foundflag=False
        errorflag=False
        for section in sections:
            if self.has_option(section,option):
                errorflag = self.set(section,option,val)
                if errorflag!=None and throwerror:
                    raise RuntimeError("ERROR %s %s %s" % (section,option,val))
                foundflag=True
                if not allflag:
                    break
        if (not foundflag) and throwerror:
            raise RuntimeError("ERROR: could not find section=%s, parameter=%s!" % (section,option))
            
        return(not foundflag)

    def setvals_nosection(self,option_val_list,allflag=False,throwerror=True):
        if option_val_list == None: return(0)
        errorflag = False
        for (option,val) in option_val_list:
            errorflag |= self.setval_nosection(option,val,allflag=allflag,throwerror=throwerror)
        return(errorflag)

    def setvals(self,section_option_val_list,throwerror=True):
        if section_option_val_list == None: return(0)
        errorflagall = False
        for (section,option,val) in section_option_val_list:
            if not self.has_option(section,option):
                errorflagall = True
                if throwerror:
                    raise RuntimeError("ERROR: section=%s, parameter=%s does not exist!" % (section,option))
                continue
            errorflag = self.set(section,option,val)
            if errorflag != None:
                errorflagall = True
                if throwerror:
                    raise RuntimeError("ERROR: could not set section=%s, parameter=%s to value %s!" % (section,option,val))
        return(errorflagall)
'''
