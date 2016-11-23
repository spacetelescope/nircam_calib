#! /usr/bin/env python

from multiprocessing import Pool
import superbias_create_ngroups_parallelized_Build5version as sb
import argparse

def setup(listfile):

    listfiles = []
    with open(listfile) as f:
        for line in f:
            if len(line) > 3:
                listfiles.append(line.strip())

    tuplist = []
    for file in listfiles:
        tup = (file,1024,1024,False,False,'./',None,False,False,False,0,20,0,False,False)
        tuplist.append(tup)
    return tuplist

def add_options(parser=None,usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage,description='calibrate cv3 data')
    parser.add_argument('listfile',help='list file that contains a list of files to process.')
    return parser

if __name__ == '__main__':
    usagestring = 'USAGE: run_pipeline.py'
    
    parser = add_options(usage=usagestring)
    args = parser.parse_args()

    input_tuples = setup(args.listfile)
    n_cores = 2
    pool = Pool(n_cores)
    pool.map(sb.create_superbias,input_tuples)
