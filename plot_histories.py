#!/usr/bin/env python
import numpy as np
import sys, os, time
from dragons import meraxes

name      = 'tuned_reion_T125'
redshift  = float(sys.argv[1])
cut_property = 'Mvir'
prop_low     = int(sys.argv[2])
prop_high    = int(sys.argv[3])

feature   = sys.argv[4]
statistic = sys.argv[5]
Nboot     = int(sys.argv[6])
fmeraxes  = '/home/mmarshal/data_dragons/%s/output/meraxes.hdf5'%(name)
snapshot  = meraxes.io.check_for_redshift(fmeraxes,redshift)[0]
results   = np.zeros([snapshot+1,3])
save_dir  = '/home/mmarshal/data_dragons/histories/%s/history_by%s/%03d/%s%02d-%02d/'%(name, cut_property,snapshot, cut_property,prop_low, prop_high)


f=np.fromfile(save_dir+'StellarMass_%s_Nboot%s.bin'%(statistic,Nboot))
f=f.reshape(snapshot+1,3)


