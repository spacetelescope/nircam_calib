NRC-24 Dither Pattern Verification  (NGAS-81, APT-1073)

This github area is for our simulation scripts and other materials for this CAR.

The simulations use this catalog from Matteo Correnti, based on Jay Anderson's astrometric catalog of the LMC (selecting objects with flag=1):
  https://stsci.box.com/s/dtdo5k9qd13w8b34wyt62rw2u7otqgnb      (lmc_catalog_flag1.cat)
  
psf_template_mp.py  This is a calling programme for stellar_photometry.py for the calculation of the average PSF for each image. This uses
multi-processing to carry out the task.

stellar_photometry.py  is the code that calculates the average PSF.
