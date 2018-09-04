#!/usr/bin/env python
import numpy as np
import sys, os, time
from dragons import meraxes
sys.path.append('Yuxiang/')
from _function import _bootstrap

t0        = time.time()
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

if statistic == 'mean':
    statistic_input = np.mean
    package         = False
elif statistic == 'median':
    statistic_input = np.median
    package         = False
elif statistic == 'std':
    statistic_input = np.std
    package         = False
else:
    print("unrecognized statistic (%s), choose from mean, median and std!"%statistic)
    exit(-1)

Types     = np.fromfile(save_dir+'Type.bin').reshape(-1,(snapshot+1))
Ghosts    = np.fromfile(save_dir+'GhostFlag.bin').reshape(-1,(snapshot+1))
Features  = np.fromfile(save_dir+feature+'.bin').reshape(-1,(snapshot+1))

if len(Features) > 0:
    Features[np.where(Ghosts!=0)] = np.nan #no ghosts
    Features[np.where(Types!=0)]  = np.nan #only central
    if feature in ["Mvir", "StellarMass"]:
        Features  = Features * 1e10
    for snap in range(snapshot+1):
        t0 = time.time()
        results[snap] = _bootstrap(Features[:,snap], num_samples = Nboot, statistic=statistic_input)
        print("[Snap%03d]%d s"%(snap, time.time() - t0))
results.tofile(save_dir+'%s_%s_Nboot%d.bin'%(feature,statistic,Nboot))
