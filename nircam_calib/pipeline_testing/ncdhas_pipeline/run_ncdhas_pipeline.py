#! /usr/bin/env python

"""
Wrapper to run NCDHAS pipeline. Must be in server environment (e.g., witserv1)
under bash shell with access to numpy packages. Input file must be in
FITSwriter format (i.e., must have data in the primary extension). You can use
ncdhas_format.py to do rough and quick file reformat.

Environment variables to set:
    export NCDHAS_PATH=/grp/software/Linux/RH6/x86_64/ncdhas
    export PATH=$PATH:$NCDHAS_PATH

To access help file for NCDHAS, type: ncdhas -h

Authors:
--------
    - Alicia Canipe
    - Brian Brooks

Dependencies:
-------------
    - numpy
    - NCDHAS
    - ncdhas_dict (dictionary to map cal steps to dhas commands)


Inputs:
-------------
    - file to be calibrated (NCDHAS expects FITSwriter format)
    - calibration steps to be performed (dq sat sup ref lin jump ramp)

Outputs:
--------
    - calibrated file (*.red.fits, renamed to reflect steps performed)
    - slope image (*.slp.fits)
    - diagnostic file (*.dia.fits)
    (all outputs are moved to individual directories)

Issues:
-------
    - NCDHAS doesnt have the ability to output intermediate steps, but only the
      combination of all performed steps to avoid taking up disk space. This
      script outputs intermediate steps by "turning off" the steps that are not
      requested. I wrote the script to run calibration steps up through those
      that are requested, for example: running ncdhas for refpix completes
      steps dq_init, saturation, superbias, and then refpix.

    - There is also no way to change the output file name for NCDHAS, so this
      is hard coded into the script. I did it by creating a directory for each
      step and moving the pipeline outputs into their respective directories.

Use:
----
    python RunNCDHASpipeline.py input_file.fits dq sat sup ref
"""


# Required packages
import os, sys
import argparse
import numpy as np
from ncdhas_dict import *


def terminal_cmd(calibstep, filebase, outdir):
    """Builds NCDHAS pipeline terminal command.

    Parameters
    ----------
    calibstep : string
        Pipeline steps to run: dq sat sup ref lin jump flat

    outdir : string
        Directory that pipeline outputs will be stored in for each step

    Returns
    -------
    name: string
        Name of calibration step to perform

    cmd : string
        Terminal command to execute to run NCDHAS pipeline

    output : string
        Default output from NCDHAS

    rename : string
        New filename for NCDHAS outputs to match SSB outputs

    """

    # get step name and terminal command, along with output extensions to add
    # to processsed file
    name = ext[calibstep]
    indivcmd = step[name]
    outname = outext[calibstep]

    # put individual calibration commands into final pipeline command
    # move outputs into the appropriate directory
    mv_cmd = 'mv '+str(filebase)+'.*.fits '+outdir+'/'
    cmd = 'ncdhas '+args.infile+' '+str(indivcmd)+' && '+mv_cmd

    # rename output files to avoid confusion
    # (all ncdhas outputs are in the format: 'fileBase.slp.fits')
    output = outdir+'/'+str(filebase)+'.red.fits'
    rename = outdir+'/'+str(filebase)+'_NCDHAS_'+str(outname)+'.fits'

    return name, cmd, output, rename


def main(args):
    """Main function."""

    # check that environment variables are set
    if os.environ.get('NCDHAS_PATH') is None:
        print('\n###################')
        print('Error! Need to set the following environment variables:\n')
        print('export NCDHAS_PATH=/grp/software/Linux/RH6/x86_64/ncdhas')
        print('export PATH=$PATH:$NCDHAS_PATH\n')
        print('###################\n')

        sys.exit(0)

    # collect calibration steps to be performed into an array for iteration
    allsteps = np.asarray(args.steps)
    print('Steps are: ', allsteps)

    filebase = args.infile.split('.fits')[0]

    # for each calibration step, build the appropriate command for the pipeline
    # make a directory for each calibration step to hold outputs
    for calib in allsteps:
        outdir = 'temp_'+calib
        os.mkdir(outdir)

        name, cmd, output, rename = terminal_cmd(calib, filebase, outdir)

        print('\n\nRunning pipeline with command: ', cmd)
        print('\nCalibration steps completed up through: ', name)
        os.system(cmd)

    print('\nMoving ', output, 'to ', rename)
    if os.path.isfile(output):
        os.rename(output, rename)


if __name__ == '__main__':

    # Command line argument handler.
    parser = argparse.ArgumentParser(
        description='Wrapper to run NCDHAS pipeline on an exposure.',
        epilog='example: python run_ncdhas_pipeilne.py infile.fits dq sat')

    parser.add_argument("infile",
                        help="File to run through pipeline.", type=str)
    parser.add_argument("steps", nargs='+',
                        help='Calibration steps to perform: dq, sat, sup, ref,\
                        lin, jump, flat.', type=str)

    args = parser.parse_args()
    main(args)
