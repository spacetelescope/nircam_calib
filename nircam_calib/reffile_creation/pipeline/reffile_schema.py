#! /usr/bin/env python

'''
Return a list of required metadata for a given reference file 
datamodel.

This class looks in the yaml files that define the schema for 
each type of reference file, and returns a dictionary of these
values. The dictionary contains the schema entry:fits keyword.

For example:
{'meta.exposure.readpatt': 'READPATT'}

DEPENDENCIES:
The JWST Calibration pipeline must be installed.

INPUTS: 
reftype: (string) Type of reference file to be investigated.
         e.g 'readnoise', 'gain', 'dark'

OUTPUTS:
metadict: dictionary of applicable, required metadata for
          the given reference file data model

HISTORY:
1 Nov 2017: created. B. Hilbert
'''

import os
import sys
import imp
import argparse
import yaml


class Schema():
    def __init__(self):
        self.reftype = ''

        
    def find_schema(self):
        # Return a list of schema for a given datamodel
        path = imp.find_module('jwst')
        yamlname = self.reftype.lower() + '.schema.yaml'
        file = os.path.join(path[1],'datamodels/schemas/',yamlname)

        if os.path.isfile(file) == False:
            print("Unable to locate schema file:")
            print("{}".format(file))
            sys.exit(0)
        else:
            # Read in yaml file
            with open(file,'r') as infile:
                refdict = yaml.load(infile)

            # Grab names of yaml files that together
            # comprise the full set of metadata
            schema_files = []
            for schtype in refdict['allOf']:
                for key in schtype:
                    if key == '$ref':
                        schema_files.append(schtype[key])

                # Save a dictionary with the metadata path names
                # and corresponding fits header keyword names
                self.metadict = {}
                for sfile in schema_files:
                    fullpath = os.path.join(path[1],'datamodels/schemas/',sfile)
                    with open(fullpath,'r') as infile:
                        newdict = yaml.load(infile)
        
                    subdict = newdict['properties']['meta']['properties']
                    for key in subdict:
                        if 'properties' not in subdict[key].keys():
                            self.metadict['meta.'+key] = subdict[key]['fits_keyword']
                        else:
                            l2dict = subdict[key]['properties']
                            for l2key in l2dict:
                                self.metadict['meta.'+key+'.'+l2key] = l2dict[l2key]['fits_keyword']

            return self.metadict


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Find required metadata')
        parser.add_argument("reftype",help='Type of datamodel to investigate')
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: reffile_schema.py readnoise'

    metadata = Schema()
    parser = metadata.add_options(usage = usagestring)
    args = parser.parse_args(namespace = metadata)
    metadata.find_schema()
    
