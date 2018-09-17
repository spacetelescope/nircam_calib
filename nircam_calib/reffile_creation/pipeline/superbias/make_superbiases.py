#! /usr/bin/env python

'''
starting with a list of bias files created from individual
ramps, create mean superbias files to be used in testing
'''

import average_biases
from astropy.io import ascii

infiles = []
with open('individual_biases.list') as f:
    for line in f:
        if len(line) > 2:
            infiles.append(line.strip())


#in order to use average_bias.py, we need to create listfiles
#for each combination of bias files that we want to combine
for i in range(2,len(infiles)):
    ofiles = infiles[0:i]
    lfilename = 'ind_biases_1_to_'+str(i)+'.list'
    fout = open(lfilename,'wb')
    for indfile in ofiles:
        fout.write(indfile+'\n')
    fout.close()

    sbname = 'NRCBLONG_superbias_from_ramps_1_through_'+str(i)+'.fits'
    print()
    print()
    print(ofiles)
    print(('saving superbias file as {}'.format(sbname)))
    print()
    print()
    average_biases.run(lfilename,sbname)
    
