# Notes on distortion reference file generation code

The code contained in this directory:
nircam_distortion_reffiles_from_pysiaf.py
make_all_imaging_distortion_reffiles_from_pysiaf.py

can be used to create distortion reference files for NIRCam, and in fact is the
code used to create the reference files delievered to CRDS in October of 2019.

The metadata in those delivered files references the spacetelescope/jwreftools
repository, as historically that repo contained the code for distortion reference files.
We still consider the copy of the code in jwreftools to be the official copy of the
code. We have placed a copy of the code in this repo as a convenience, since most
NIRCam team members come to this repo for NIRCam reference file code. If you wish
to modify the code used to generate the distortion reference files, you should
either modify the version in jwreftools, or modify this copy and then add a note
to the jwreftools repo that you have moved the official copy of the code to
nircam_calib (and be sure to update the referenece file metadata to point to
nircam_calib).
