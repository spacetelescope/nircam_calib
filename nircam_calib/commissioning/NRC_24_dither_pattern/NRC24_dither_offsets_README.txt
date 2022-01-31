Readme file for:

	NRC24_dither_offsets.py

Anton M. Koekemoer, STScI


Overview:

 - This script is designed to run on a set of calibrated, detector-frame (unrectified) Pipeline Stage 2 output files (*_cal.fits).

 - Running the script:
 
	- Environment variables needed:

		MIRAGE_DATA
		PYSYN_CDBS
		WEBBPSF_PATH
		CRDS_CONTEXT
		CRDS_PATH
		CRDS_SERVER_URL

	- Shell:

		bash (mostly for consistency with conda install of JWST pipeline and MIRAGE)

 	- Usage:

		python NRC24_dither_offsets.py

The current version has no input args but this will change in the next version,
which will have these inputs (currently these are just hardwirsed):

	- Input args:

		xml_file	(eg 'nrc24-1073_same_PA_141deg31.xml')

		pointing_file	(eg 'nrc24-1073_same_PA_141deg31.pointing')

		analysis_type	(either 'absolute', which needs external ref catalog, or 'relative', w.r.t. first exposure)

	- Additional environment variables:

		catfile_lmc		(eg 'lmc_catalog_flag1.cat')

		analysis_dir		(eg 'nrc24_analysis/')

		pipeline_outputs_stage2 (eg 'pipeline_outputs_stage2/')

	And whatever other inputs / env variables are useful.


Script execution:

	- based on the input files, construct dicts of propid, observation number, filter, exposure number, and SCA fits files

	- detect all sources in all SCAs, and convert them all to RA,Dec

	- for the 'absolute' option, transform all RA,Dec coords in reference cat to Ideal X,Y  for NIRCAM_FULL, for a (0,0) dither offset
	  (currently the (0,0) dither file needs to be hardcoded, since not all the dither patterns include a (0,0) dither)

	- for each combination of propid / observation / filter:

		- match the RA,Dec of all detected sources to those in the reference catalog

		- transform all sources to IDeal X,Y for NIRCAM_FULL (via v2,v3 first)

		- calculate dx_IDL, dy_IDL offsets w.r.t. the reference catatalog Ideal X,Y

		- compare the resulting measured dx_IDL, dy_IDL offsets wirh the commanded dithers

		- produce diagnostic plots and print out tables comparing commanded vs measured offsets.




