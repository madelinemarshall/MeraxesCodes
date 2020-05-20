import numpy as np
from dragons import meraxes
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
from _load_data import load_data

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.6,4)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4

data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'

if __name__=="__main__":
  snapshots=range(63,158,8)
  snapshots_125=range(63,194,8) #to z=1
  props=['BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag','BlackHoleMass_Coalescence','BlackHoleMass_Radio']
  fig, axes = plt.subplots(1, 1)
  
  filenames=['paper2','paper2_T125']#'tuned_reion_T125','tuned_T125_MR','tuned_T125_NR','tuned_reion']
  lines=['-','--',':']
  linewidth=[2.5,1.5]
  jj=0
  lgd=''
 
  for filename in filenames:
    if jj>0:
      snapshots=snapshots_125
    redshift={}
    ii=-1
    BHgrowth=np.zeros((len(snapshots),1))
    BHgrowth_MD=np.zeros((len(snapshots),1))
    BHgrowth_ID=np.zeros((len(snapshots),1))
    BHgrowth_Coalescence=np.zeros((len(snapshots),1))
    BHgrowth_Radio=np.zeros((len(snapshots),1))
    for snapshot in snapshots:
      ii+=1
      redshift[snapshot] = meraxes.io.grab_redshift(data_folder+filename+meraxes_loc,snapshot)
      gals_bulges=load_data(filename,snapshot,props)
      gals_bulges_past=load_data(filename,snapshot-1,props)
      for BHlimit in [6]:
        gals_bulges=gals_bulges[gals_bulges['BlackHoleMass']*1e10>10**BHlimit]
        gals_bulges_past=gals_bulges_past[gals_bulges_past['BlackHoleMass']*1e10>10**BHlimit]
        BHgrowth[ii]=(np.sum(gals_bulges['BlackHoleMass'])-np.sum(gals_bulges_past['BlackHoleMass']))/len(gals_bulges)
        BHgrowth_MD[ii]=(np.sum(gals_bulges['BlackHoleMass_MD'])-np.sum(gals_bulges_past['BlackHoleMass_MD']))/len(gals_bulges)
        BHgrowth_Coalescence[ii]=(np.sum(gals_bulges['BlackHoleMass_Coalescence'])-np.sum(gals_bulges_past['BlackHoleMass_Coalescence']))/len(gals_bulges)
        BHgrowth_ID[ii]=(np.sum(gals_bulges['BlackHoleMass_ID'])-np.sum(gals_bulges_past['BlackHoleMass_ID']))/len(gals_bulges)
        BHgrowth_Radio[ii]=(np.sum(gals_bulges['BlackHoleMass_Radio'])-np.sum(gals_bulges_past['BlackHoleMass_Radio']))/len(gals_bulges)
    axes.plot(np.array(list(redshift.values())),BHgrowth_ID/BHgrowth,lines[jj],label='Instability-Driven Quasar-Mode',color=colors[4],lw=linewidth[jj])
    axes.plot(np.array(list(redshift.values()))[BHgrowth_MD[:,0]>0],BHgrowth_MD[BHgrowth_MD[:,0]>0]/BHgrowth[BHgrowth_MD[:,0]>0],lines[jj],label='Merger-Driven Quasar-Mode',color=colors[3],lw=linewidth[jj])
    axes.plot(np.array(list(redshift.values()))[BHgrowth_Coalescence[:,0]>0],BHgrowth_Coalescence[BHgrowth_Coalescence[:,0]>0]/BHgrowth[BHgrowth_Coalescence[:,0]>0],lines[jj],label='BH-BH Coalescence',color=colors[2],lw=linewidth[jj])
    axes.plot(np.array(list(redshift.values())),BHgrowth_Radio/BHgrowth,lines[jj],label='Radio-Mode',color=colors[1],lw=linewidth[jj])
    if jj==0:
        axes.plot(0,0,linestyle='-',lw=2.5,color=colors[4],label='Tiamat')
        axes.plot(0,0,linestyle='--',lw=1.5,color=colors[4],label='Tiamat-125-HR')
        lgd=plt.legend(loc='upper center', bbox_to_anchor=(0.45, -0.2))
    jj+=1
  
  
  axes.set_xlabel('Redshift')
  axes.set_yscale('log')
  axes.set_xlim([1,7])
  axes.set_ylabel(r'$\log(\Delta M_{\rm{BH,~ growth~ mode}}/ \Delta M_{\rm{BH,~ total}})$')
  plt.tight_layout()
  plt.savefig('/home/mmarshal/results/plots/Paper2/BHGrowthModes_redshift_T125.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

