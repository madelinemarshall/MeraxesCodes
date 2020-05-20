import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
from _load_data import load_data
from astropy.cosmology import Planck15

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.6,4)
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4

data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'
cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  } 

if __name__=="__main__":
  #snapshots=range(63,250,5) #every 5 snapshots
  #snapshots=range(63,158,15) #every 5 snapshots
  snapshots=np.arange(50,158,20) #every 8 snapshots
  redshift={}
  #props=['BlackHoleMass','GhostFlag','ID']
  props=['ID','BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag','BlackHoleMass_Coalescence','BlackHoleMass_Radio']
  fig, axes = plt.subplots(1, 1)
  filenames=['paper2','paper2_high_eta']#'paper1_T125_both']#'draft2_reion_T125']#'tuned_reion_T125','tuned_T125_MR','tuned_T125_NR','tuned_reion']
  volume=100**3#(125/cosmo['h'])**3
  for filename in filenames:
    ii=-1
    BHgrowth=np.zeros((len(snapshots),1))
    BHgrowth_up=np.zeros((len(snapshots),1))
    BHgrowth_low=np.zeros((len(snapshots),1))
    BHgrowth_def=np.zeros((len(snapshots),1))
    time_diff=np.zeros((len(snapshots),1))
    sim_props = meraxes.io.read_input_params(data_folder+filename+meraxes_loc,h=cosmo['h'],quiet=True)
    EddingtonRatio = sim_props['EddingtonRatio']
    if 'RadioAccretionEff' in sim_props.keys():
      ETA=sim_props['RadioAccretionEff']
    else:
      ETA = 0.06
    
    for snapshot in snapshots:
      ii+=1
      redshift[snapshot] = meraxes.io.grab_redshift(data_folder+filename+meraxes_loc,snapshot)
      redshift_next = meraxes.io.grab_redshift(data_folder+filename+meraxes_loc,snapshot-1)
       
      time_diff[ii]=abs(Planck15.age(redshift_next).value-Planck15.age(redshift[snapshot]).value)*1e9 #yr
      EddFactor=(np.exp(EddingtonRatio*time_diff[ii]/(ETA*450*1e6))-1)

      gals_bulges=load_data(filename,snapshot,props)
      gals_bulges_past=load_data(filename,snapshot-1,props)
      BHlimit=6
      gals_bulges=gals_bulges[gals_bulges['BlackHoleMass']*1e10>10**BHlimit]
      gals_bulges_past=gals_bulges_past[gals_bulges_past['BlackHoleMass']*1e10>10**BHlimit]
      
      IDs=np.zeros(len(gals_bulges),dtype=bool) 
      for jj in np.arange(0,len(gals_bulges['ID'])): IDs[jj]=gals_bulges['ID'][jj] in gals_bulges_past['ID']
      gals_bulges=gals_bulges[IDs]      
      
      IDs=np.zeros(len(gals_bulges_past),dtype=bool) 
      for jj in np.arange(0,len(gals_bulges_past['ID'])): IDs[jj]=gals_bulges_past['ID'][jj] in gals_bulges['ID']
      gals_bulges_past=gals_bulges_past[IDs]      

      individualBHgrowth=((gals_bulges['BlackHoleMass']-gals_bulges['BlackHoleMass_Coalescence'])-\
                          (gals_bulges_past['BlackHoleMass']-gals_bulges_past['BlackHoleMass_Coalescence']))*1e10
      individualEddLim=gals_bulges_past['BlackHoleMass']*1e10*EddFactor

      #BHgrowth[ii]=(np.sum(gals_bulges['BlackHoleMass'])-np.sum(gals_bulges_past['BlackHoleMass']))*1e10/time_diff[ii]/volume
      #Eddington_Lim[ii]=np.sum(gals_bulges_past['BlackHoleMass'])*1e10*EddFactor/time_diff[ii]/volume
      BHgrowth[ii,0]=np.median(individualBHgrowth/individualEddLim)
      BHgrowth_up[ii,0]=np.percentile(individualBHgrowth/individualEddLim,84)
      BHgrowth_low[ii,0]=np.percentile(individualBHgrowth/individualEddLim,16)
    lines=['-','--',':']
    axes.plot((np.array(list(redshift.values()))),BHgrowth[:,0],lines[0],label=filename,lw=2.5)
    axes.fill_between(np.array(list(redshift.values())),BHgrowth_up[:,0],BHgrowth_low[:,0],alpha=0.15)

 
  lgd=plt.legend()
  axes.set_xlabel('z')
  axes.set_yscale('log')
  #axes.set_xlim([2,7])
  #axes.set_ylabel(r'$\log(\Delta M_{\rm{BH,~ growth~ mode}})$')
  axes.set_ylabel('Median Eddington Ratio')
  #plt.tight_layout()
#  axes.set_ylim(10**-5.9,10**-3.1)

 # ax2 = axes.twiny()
  #ax2.set_xlabel("z")
  #ax2.set_xticks(np.log10(1+np.linspace(0,6,7)))
  #ax2.set_xticklabels(['0','1','2','3','4','5','6'])
  
  #ax3 = axes.twinx()
  #ax3.set_ylabel(r'$3300 \log(\Delta M_{\rm{BH,~ growth~ mode}})$')
  #ax3.set_yticks(np.array([-5.9,-5.5,-5.1,-4.7,-4.3,-3.9]))
  #ax3.set_yticklabels(['-2.4','-2','-1.6','-1.2','-0.8','-0.4'])
  #ax3.set_yticks(np.array([0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1]))
  #ax3.set_yticklabels(['-2.4','-2','-1.6','-1.2','-0.8','-0.4','0','0.4'])
  
  #plt.savefig('/home/mmarshal/results/plots/BHGrowthWithRedshift.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

