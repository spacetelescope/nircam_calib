#! /usr/bin/env python

'''find filters where the effective response is more than 10% different
between module A and module B'''

import numpy as np

filts = np.array(['F070W','F090W','F115W','F140M','F150W2','F150W','F162M','F164N','F182M','F187N','F200W','F210M','F212N','F250M','F277W','F300M','F322W2','F323N','F335M','F356W','F360M','F405N','F410M','F430M','F444W','F460M','F466N','F470N','F480M'])
moda = np.array([.200,.290,.327,.399,.416,.422,.415,.354,.455,.365,.475,.463,.400,.401,.409,.380,.450,.284,.424,.485,.459,.355,.458,.453,.453,.379,.276,.255,.319])
modb = np.array([.201,.294,.323,.398,.414,.422,.416,.366,.454,.393,.470,.462,.388,.365,.373,.347,.437,.281,.412,.475,.458,.394,.486,.486,.497,.412,.313,.291,.387])

diff = (moda - modb)/moda
print(diff)
bad = np.where((diff > 0.05) | (diff < -0.05))[0]
print('Inconsistent effective response: ',filts[bad])
print(diff[bad])

bwa = np.array([.132,.189,.224,.141,1.169,.319,.168,.02,.238,.024,.457,.203,.027,.18,.68,.311,1.392,.038,.364,.775,.377,.044,.449,.223,.993,.228,.054,.052,.276])

bwb = np.array([.133,.199,.224,.142,1.178,.317,.168,.02,.237,.023,.456,.203,.028,.179,.67,.319,1.317,.038,.34,.782,.363,.046,.425,.231,1.060,.227,.052,.05,.329])

diffbw = (bwa-bwb)/bwa
badbw = np.where((diffbw > 0.05) | (diffbw < -0.05))[0]
print('Inconsistent bandwidth: ',filts[badbw]) 
print(diffbw[badbw])
