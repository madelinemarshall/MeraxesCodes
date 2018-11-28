#!/usr/bin/env python
import numpy as np
import sys, os, time
from dragons import meraxes
import matplotlib.pyplot as plt

name      = 'tuned_reion'
redshift  = float(sys.argv[1])
cut_property = 'BlackHoleMass'
prop_low     = int(sys.argv[2])
prop_high    = int(sys.argv[3])

statistic = sys.argv[4]
Nboot     = int(sys.argv[5])
fmeraxes  = '/home/mmarshal/data_dragons/%s/output/meraxes.hdf5'%(name)
snapshot  = meraxes.io.check_for_redshift(fmeraxes,redshift)[0]
results   = np.zeros([snapshot+1,3])
save_dir  = '/home/mmarshal/data_dragons/histories/%s/history_by%s/%03d/%s%02d-%02d/'%(name, cut_property,snapshot, cut_property,prop_low, prop_high)
snapshots=range(37,159)
EXPANSION_FACTOR_PATH = "/fred/oz013/simulations/Tiamat/a_list.txt"
a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
z_list = np.array(1.0/a_list - 1.0)
z_list = z_list[snapshots]

fig,axes=plt.subplots(1,3)

#________
feature = 'BlackHoleMass'
bh=np.fromfile(save_dir+'%s_%s_Nboot%s.bin'%(feature,statistic,Nboot)).reshape(len(snapshots),3)
axes[0].fill_between(z_list, bh[:,0]-bh[:,1], bh[:,0]+bh[:,2], alpha=0.5,label='__nolegend__')
axes[0].plot(z_list,bh[:,0])
#axes[0].plot(z_list,1e3*np.exp(14/(1+z_list)**1.5/0.45/0.2))
axes[0].set_xlabel('Redshift')
axes[0].set_ylabel('%s'%(feature))
axes[0].set_yscale('log')
axes[0].invert_xaxis()
axes[0].set_ylim([1e4,1e9])
#_________
feature = 'StellarMass'
sm=np.fromfile(save_dir+'%s_%s_Nboot%s.bin'%(feature,statistic,Nboot)).reshape(len(snapshots),3)

axes[1].fill_between(z_list, sm[:,0]-sm[:,1],sm[:,0]+sm[:,2], alpha=0.5,label='__nolegend__')
axes[1].plot(z_list,sm[:,0])
axes[1].set_xlabel('Redshift')
axes[1].set_ylabel('%s'%(feature))
axes[1].set_yscale('log')
axes[1].invert_xaxis()
#_________
feature = 'BulgeStellarMass'
bsm=np.fromfile(save_dir+'%s_%s_Nboot%s.bin'%(feature,statistic,Nboot)).reshape(len(snapshots),3)

axes[2].fill_between(z_list, bsm[:,0]/sm[:,0], bsm[:,0]/sm[:,0], alpha=0.5,label='__nolegend__')
axes[2].plot(z_list,bsm[:,0]/sm[:,0])
axes[2].set_xlabel('Redshift')
axes[2].set_ylabel('%s'%(feature))
axes[2].set_yscale('log')
axes[2].invert_xaxis()

plt.tight_layout()
plt.show()
