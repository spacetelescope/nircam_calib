python test_saturation.py ../exposures/nrca1_47Tuc_subpix_dither1_newpos_uncal.fits nrca1_modified_saturation_reffile.fits  --run_dq --maskfile ../reffiles/NRCA1_17004_BPM_ISIMCV3_2016-01-21_ssbspmask_DMSorient.fits --outfile nrca1_47Tuc_subpix_dither1_newpos_dq_init_saturation.fits

python test_saturation.py ../dq_init/nrca1_47Tuc_subpix_dither1_newpos_dq_init.fits nrca1_modified_saturation_reffile.fits  --outfile nrca1_47Tuc_subpix_dither1_newpos_dq_init_saturation.fits

