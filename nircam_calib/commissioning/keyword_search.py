#! /usr/bin/env python
"""
Search the JWST Keyword Dictionary for header keyword information. The Keyword
Dictionary is here: https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/
Right now, there is no API so the code downloads a JSON file version of the
keyword database and searches the JSON file.

Authors:
--------
    - Alicia Canipe

Dependencies:
-------------
    - astroquery

Use:
----
1. To pull up all information about a keyword from all schemas:

    python keyword_search.py keyword

2. To pull up details about a particular keyword attribute (e.g., source):

    python keyword_search.py keyword --attrib='source'

3. To pull up details about a keyword for a specific instrument:

    python keyword_search.py keyword --attrib='source' --inst='nircam'

4. To pull up keyword information from a different keyword dictionary version:

    python keyword_search.py keyword --version='JWSTDP-2018_2b-190328'

5. To turn off the verbosity and/or save out a text file with the results:

    python keyword_search.py keyword --noverb --save

"""

# Required packages
import sys
import argparse
import pprint
from time import time
from astroquery.mast import Mast


def mast_request(version):
    """Access MAST and get JSON version of the keyword dictionary.

    Parameters
    ----------
    version : string
        Keyword dictionary version

    Returns
    -------
    kwd_dict : json object
        Keyword dictionary as a JSON nested dictionary object
    """
    service = "Mast.Jwst.Keywords"

    params = {"filename": "all.json",
              "version": str(version),
              "cacheBreaker": time()}

    resp = Mast.service_request_async(service, params)
    kwd_dict = resp[0].json()

    return kwd_dict


def getpath(kwd_dict, keyword, prepath=()):
    """Get dictionary location of keywords.

    Parameters
    ----------
    kwd_dict : json object
        Keyword dictionary as a JSON nested dictionary object

    keyword : string
        Keyword of interest

    Returns
    -------
    path : string
        Dictionary location of keyword
    """

    for k, v in kwd_dict.items():
        path = prepath + (k,)

        # if v is a dictionary
        if hasattr(v, 'items'):

            # recursive call
            p = getpath(v, keyword, path)
            if p is not None:
                return p

        # if keyword is found
        elif v == keyword:
            return path


def recursive_keyword_search(kwd_dict, keyword):
    """Do recursive search of nested keyword dictionary to find keyword info.

    Parameters
    ----------
    kwd_dict : json object
        Keyword dictionary as a JSON nested dictionary object

    keyword : string
        Keyword of interest

    Returns
    -------
    results : list
        List of all the schemas that contain the keyword
    """

    # check that input is a dictionary
    if not isinstance(kwd_dict, dict):
        return False

    results = []
    for k, v in kwd_dict.items():

        # recursively search for the schemas containing the keyword
        if hasattr(v, 'items') and getpath(v, keyword) is not None:
            tmp = []
            for i in (getpath(v, keyword), k):
                if type(i) == (tuple) or type(i) == (list) or type(i) == (set):
                    for j in i:
                        tmp.append(j)
                else:
                    tmp.append(i)
            results.append(tmp)
        else:

            # skip the schema if it doesn't contain the keyword
            pass

    return results


def get_keyword_info(kwd_dict, kd_search, keyword, attrib, inst):
    """Grab information on the header keyword (e.g., source).

    Parameters
    ----------
    kwd_dict : json object
        Keyword dictionary as a JSON nested dictionary object

    kd_search : list
        List of keyword dictionary paths to access information about a keyword.

    keyword : string
        Keyword of interest

    attrib : string
        KD attribute to pull from dictionary (e.g., source)

    inst : string
        Instrument of interest

    Returns
    -------
    details : list
        Info from all search results.
    """

    if args.save is True:
        log = open("kd_search_"+keyword+".txt", "w")
    else:
        log = None

    for result in kd_search:
        if inst not in ['all', 'ALL']:
            # save out keyword information according to arguments
            if inst in result[-1]:
                schema = result[-1]
                kd_loc = result[0]
                info = '\nKeyword: '+keyword+'\n'+'Schema: '+schema+'\n'
                output_options(info, args.noverb, args.save, log, outst='\n')

                if attrib not in ['all', 'ALL']:
                    details = kwd_dict[schema][kd_loc][attrib]
                    detst = '{} : {}'.format(attrib, details)
                    output_options(detst, args.noverb, args.save, log)

                else:
                    details = kwd_dict[schema][kd_loc]
                    detst = 'Details: '+pprint.pformat(details)
                    output_options(detst, args.noverb, args.save, log)

            if inst not in result[-1] and 'core' in result[-1]:
                schema = result[-1]
                kd_loc = result[0]
                info = '\nKeyword: '+keyword+'\n'+'Schema: '+schema+'\n'
                output_options(info, args.noverb, args.save, log, outst='\n')

                if attrib not in ['all', 'ALL']:
                    details = kwd_dict[schema][kd_loc][attrib]
                    detst = '{} : {}'.format(attrib, details)
                    output_options(detst, args.noverb, args.save, log)

                else:
                    details = kwd_dict[schema][kd_loc]
                    detst = 'Details: '+pprint.pformat(details)
                    output_options(detst, args.noverb, args.save, log)

            else:
                pass

        else:
            schema = result[-1]
            kd_loc = result[0]
            info = '\nKeyword: '+keyword+'\n'+'Schema: '+schema+'\n'
            output_options(info, args.noverb, args.save, log, outst='\n')

            if attrib not in ['all', 'ALL']:
                details = kwd_dict[schema][kd_loc][attrib]
                detst = '{} : {}'.format(attrib, details)
                output_options(detst, args.noverb, args.save, log)
            else:
                details = kwd_dict[schema][kd_loc]
                detst = 'Details: '+pprint.pformat(details)
                output_options(detst, args.noverb, args.save, log)

    # for prettier formatting
    output_options('\n\n', args.noverb, args.save, log)


def output_options(string, verbose, save, log, outst=''):
    """Control terminal output and log file with a convenience function.

    Parameters
    ----------
    string : string
        Text to print to terminal or save in output file

    verbose : bool
        Turn verbosity on or off

    save : bool
        Save out text file with search results

    log : file
        Log file to write out search results

    Returns
    -------

    """

    if verbose is True:
        sys.stdout.write(outst)
        sys.stdout.write(string)

    if save is True:
        log.write(outst)
        log.write(string)


def main(args):
    """Run the main function.
    """

    # Access the JWST Keyword Dictionary
    kwd_dict = mast_request(args.kd_version)

    # Do a search for the keyword of interest
    kwd = str(args.kwd).upper()
    kd_search = recursive_keyword_search(kwd_dict, kwd)

    # Throw an error and exit if keyword not found
    if not kd_search:
        print('\nERROR! Keyword not found in any schemas. Exiting.\n')
        sys.exit(0)

    elif kd_search:
        # Grab the information from the dictionary
        attrib = str(args.attrib).lower()
        inst = str(args.inst).lower()
        get_keyword_info(kwd_dict, kd_search, kwd, attrib, inst)


if __name__ == '__main__':
    # Command line argument handler.
    parser = argparse.ArgumentParser(
        description='Programmatically search the JWST Keyword Dictionary',
        epilog='example: python keyword_search.py EXP_TYPE --attrib="source" --save')

    parser.add_argument('kwd',
                        help='keyword to search for in the dictionary.')
    parser.add_argument('--attrib', default="ALL",
                        help='attribute of interest (e.g., source)')
    parser.add_argument('--inst', default="ALL",
                        help='instrument to search')
    parser.add_argument('--kd_version', default="JWSTDP-2019.1.0-95~3f0d4fc0",
                        help='keyword dictionary version')
    parser.add_argument('--save', action='store_true',
                        help='save output to a text file')
    parser.add_argument('--noverb', action='store_false',
                        help='print output to terminal')
    args = parser.parse_args()

main(args)
