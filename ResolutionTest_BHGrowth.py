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
      snapshot=snapshot,props=['BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag','BlackHoleMass_Coalescence','BlackHoleMass_Radio'],\
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
  #redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,192:1,250:0}
  #snapshots=[52, 63,78,100,116,134,158,192,250]
  snapshots=range(30,158,5)
  redshift={}
  #fig, axes = plt.subplots(1, 3)
  fig, axes = plt.subplots(1, 1)
  
  filenames=['tuned_reion','tuned_reion_T125']#,'tuned_T125_MR','tuned_T125_NR','tuned_reion']
  for filename in filenames:
    ii=-1
    #if filename=='tuned_reion':
    #  snapshots=range(30,158,10)#[52, 63,78,100,116,134,158]
    #  redshift={}#{52:8,63:7,78:6,100:5,116:4,134:3,158:2}
    BHgrowth=np.zeros((len(snapshots),3))
    BHgrowth_MD=np.zeros((len(snapshots),3))
    BHgrowth_ID=np.zeros((len(snapshots),3))
    BHgrowth_Coalescence=np.zeros((len(snapshots),3))
    BHgrowth_Radio=np.zeros((len(snapshots),3))
    for snapshot in snapshots:
      ii+=1
      redshift[snapshot] = meraxes.io.grab_redshift(data_folder+filename+meraxes_loc,snapshot)
      gals_bulges=load_data(filename,meraxes_loc,snapshot,cosmo)
      gals_bulges_past=load_data(filename,meraxes_loc,snapshot-2,cosmo)
      jj=0
      for BHlimit in [6]:#,7,8]:
        gals_bulges=gals_bulges[gals_bulges['BlackHoleMass']*1e10>10**BHlimit]
        gals_bulges_past=gals_bulges_past[gals_bulges_past['BlackHoleMass']*1e10>10**BHlimit]
        BHgrowth[ii,jj]=(np.sum(gals_bulges['BlackHoleMass'])-np.sum(gals_bulges_past['BlackHoleMass']))/len(gals_bulges)
        BHgrowth_MD[ii,jj]=(np.sum(gals_bulges['BlackHoleMass_MD'])-np.sum(gals_bulges_past['BlackHoleMass_MD']))/len(gals_bulges)
        BHgrowth_Coalescence[ii,jj]=(np.sum(gals_bulges['BlackHoleMass_Coalescence'])-np.sum(gals_bulges_past['BlackHoleMass_Coalescence']))/len(gals_bulges)
        BHgrowth_ID[ii,jj]=(np.sum(gals_bulges['BlackHoleMass_ID'])-np.sum(gals_bulges_past['BlackHoleMass_ID']))/len(gals_bulges)
        BHgrowth_Radio[ii,jj]=(np.sum(gals_bulges['BlackHoleMass_Radio'])-np.sum(gals_bulges_past['BlackHoleMass_Radio']))/len(gals_bulges)
        jj+=1
    #axes[0].plot(np.array(list(redshift.values())),BHgrowth_ID[:,0]/BHgrowth_MD[:,0])
    #axes[1].plot(np.array(list(redshift.values())),BHgrowth_ID[:,0]/BHgrowth[:,0])
    #axes[2].plot(np.array(list(redshift.values())),BHgrowth_MD[:,0]/BHgrowth[:,0])
    
    #axes[0].plot(np.array(list(redshift.values())),BHgrowth_ID[:,1]/BHgrowth_MD[:,1])
    #axes[1].plot(np.array(list(redshift.values())),BHgrowth_ID[:,1]/BHgrowth[:,1])
    #axes[2].plot(np.array(list(redshift.values())),BHgrowth_MD[:,1]/BHgrowth[:,1])
    
    #axes[0].plot(np.array(list(redshift.values())),BHgrowth_ID[:,2]/BHgrowth_MD[:,2])
    #axes[1].plot(np.array(list(redshift.values())),BHgrowth_ID[:,2]/BHgrowth[:,2])
    #axes[2].plot(np.array(list(redshift.values())),BHgrowth_MD[:,2]/BHgrowth[:,2])
    axes.plot(np.array(list(redshift.values())),BHgrowth_ID[:,0]/BHgrowth[:,0],label='Instability-Driven Quasar-Mode')
    axes.plot(np.array(list(redshift.values())),BHgrowth_MD[:,0]/BHgrowth[:,0],label='Merger-Driven Quasar-Mode')
    axes.plot(np.array(list(redshift.values())),BHgrowth_Coalescence[:,0]/BHgrowth[:,0],label='BH-BH Coalescence')
    axes.plot(np.array(list(redshift.values())),BHgrowth_Radio[:,0]/BHgrowth[:,0],label='Radio-Mode')
    
  #axes[0].set_ylabel('Average BH Growth: ID/MD')
  #axes[1].set_ylabel('Average BH Growth: ID/Total')
  #axes[2].set_ylabel('Average BH Growth: MD/Total')
  axes.set_ylabel('Average BH Growth Relative To Total')
  axes.set_xlabel('Redshift')
  axes.set_yscale('log')
  plt.legend()
  #axes[2].legend(['HR','MR','NR','Tiamat'])
  #axes[2].legend(['Tiamat-125','Tiamat'])
  #axes[2].legend(['Tiamat'])
  plt.savefig('/home/mmarshal/results/plots/BHGrowthMode_redshift.pdf',format='pdf')
  plt.show()

