Readme file for:

	NRC24_dither_offsets.py

	Author: Anton M. Koekemoer, STScI

	v1.1	2022-01-31	Added input args and env variables, tidied some up.

	v1.2	2022-01-31	Fixed some bugs so that it should at least ru for analysis_type='absolute'



Overview:

 - This script is designed to run on a set of calibrated, detector-frame (unrectified) Pipeline Stage 2 output files (*_cal.fits).

 - Running the script:
 
	- Shell:

		bash (mostly for consistency with conda install of JWST pipeline and MIRAGE)

 	- Usage:

		python NRC24_dither_offsets.py

The current version has no input args but this will change in the next version,
which will have these inputs (currently these are just hardwirsed):

	- Input args:

		-x	--xmlfile',        Input xml file from APT	(default='nrc24-1073_same_PA_141deg31.xml')
		-p	--pointing',       Input pointing file from APT	(default='nrc24-1073_same_PA_141deg31.pointing') 
		-a	--analysis_type',  Type of analysis		(default="absolute", other allowed choice is "relative")
		-r	--refcat',         Reference astrometric cat	(default='lmc_catalog_flag1.cat')


	- Additional environment variables (optional, default to './' if not provided)

		analysis_dir		(eg 'nrc24_analysis/')

		pipeline_outputs_stage2 (eg 'pipeline_outputs_stage2/')

	- Note that it will currently just look for all *_cal.fits files in path "pipeline_outputs_stage2"

	- Also, the default "Analysis_Type" is :absolute", for which it needs this catalog:

		lmc_catalog_flag1.cat

		(which is available from our NRC24 page  https://outerspace.stsci.edu/display/JN/CAP%3A+NIRCam-24



Script execution:

	- based on the input files, construct dicts of propid, observation number, filter, exposure number, and SCA fits files

	- detect all sources in all SCAs, and convert them all to RA,Dec

	- for the 'absolute' option, transform all RA,Dec coords in reference cat to Ideal X,Y  for NIRCAM_FULL, for a (0,0) dither offset
	  (currently the (0,0) dither file needs to be hardcoded, since not all the dither patterns include a (0,0) dither)

	- for each combination of propid / observation / filter:

		- loop through all the exposures, and all the SCA files for each exposure

		- match the RA,Dec of all detected sources to those in the reference catalog

		- transform all sources to IDeal X,Y for NIRCAM_FULL (via v2,v3 first)

		- calculate dx_IDL, dy_IDL offsets w.r.t. the reference catatalog Ideal X,Y

		- compare the resulting measured dx_IDL, dy_IDL offsets wirh the commanded dithers

		- produce diagnostic plots and print out tables comparing commanded vs measured offsets.




