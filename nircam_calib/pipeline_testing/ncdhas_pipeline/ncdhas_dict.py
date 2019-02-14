step = {}
step['dq_init'] = '-dr -ipc -cbs -cs -cl'
step['saturation'] = '-dr -ipc -cbs -cl'
step['superbias'] = '-dr -ipc -cl'
step['refpix'] = '-ipc -cl'
step['linearity'] = '-ipc'
step['jump'] = '-ipc +cr'

# flat field step has a problem -- have to manually point to calibration file directory for now
step['flat_field'] = '-ipc +cr cf +FFf /grp/software/Linux/RH6/x86_64/ncdhas/cal/Flat/NRCA1_16989_PFlat_F090W_CLEAR_2015-09-11.fits'

ext = {}
ext['dq'] = 'dq_init'
ext['sat'] = 'saturation'
ext['sup'] = 'superbias'
ext['ref'] = 'refpix'
ext['lin'] = 'linearity'
ext['jump'] = 'jump'
ext['flat'] = 'flat_field'

outext = {}
outext['dq'] = 'dq_init'
outext['sat'] = 'dq_init_saturation'
outext['sup'] = 'dq_init_saturation_superbias'
outext['ref'] = 'dq_init_saturation_superbias_refpix'
outext['lin'] = 'dq_init_saturation_superbias_refpix_linearity'
outext['jump'] = 'dq_init_saturation_superbias_refpix_linearity_jump'
outext['flat'] = 'dq_init_saturation_superbias_refpix_linearity_jump_flat_field'

readsteps = {}
readsteps['step'] = step
readsteps['ext'] = ext
readsteps['outext'] = outext
