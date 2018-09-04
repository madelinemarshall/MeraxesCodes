import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,6)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,meraxes_loc,snapshot,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[(gals['BlackHoleMass']*1e10>1e6)]
  return gals


if __name__=="__main__":
  #Setup
  cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  }
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,192:1,250:0}
  snapshots=[52, 63,78,100,116,134,158,192,250]
  fig, axes = plt.subplots(1, 3)
  
  filenames=['tuned_reion_T125','tuned_T125_MR','tuned_T125_NR','tuned_reion']
  for filename in filenames:
    ii=-1
    if filename=='tuned_reion':
      snapshots=[52, 63,78,100,116,134,158]
      redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
    BHgrowth=np.zeros(len(snapshots))
    BHgrowth_MD=np.zeros(len(snapshots))
    BHgrowth_ID=np.zeros(len(snapshots))
    for snapshot in snapshots:
      ii+=1
      gals_bulges=load_data(filename,meraxes_loc,snapshot,cosmo)
      gals_bulges_past=load_data(filename,meraxes_loc,snapshot-10,cosmo)
      BHgrowth[ii]=(np.sum(gals_bulges['BlackHoleMass'])-np.sum(gals_bulges_past['BlackHoleMass']))/len(gals_bulges)
      BHgrowth_MD[ii]=(np.sum(gals_bulges['BlackHoleMass_MD'])-np.sum(gals_bulges_past['BlackHoleMass_MD']))/len(gals_bulges)
      BHgrowth_ID[ii]=(np.sum(gals_bulges['BlackHoleMass_ID'])-np.sum(gals_bulges_past['BlackHoleMass_ID']))/len(gals_bulges)
    axes[0].plot(np.array(list(redshift.values())),np.log10(BHgrowth))
    axes[1].plot(np.array(list(redshift.values())),np.log10(BHgrowth_ID/BHgrowth))
    axes[2].plot(np.array(list(redshift.values())),np.log10(BHgrowth_MD/BHgrowth))
  axes[0].set_xlabel('Average Total BH Growth')
  axes[1].set_xlabel('Average ID BH Growth')
  axes[2].set_xlabel('Average MD BH Growth')
  axes[2].legend(['HR','MR','NR','Tiamat'])
  plt.show()

