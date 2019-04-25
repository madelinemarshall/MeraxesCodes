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
plt.rc('text', usetex=True)
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
  snapshots=np.arange(30,158,20) #every 5 snapshots
  redshift={}
  props=['ID','BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag','BlackHoleMass_Coalescence','BlackHoleMass_Radio']
  props_default=['ID','BlackHoleMass','GhostFlag']
  fig, axes = plt.subplots(1, 1)
  filename_def='dragons10'
  filenames=['dragons10','paper2','paper2_low_edd']#'paper1_T125_both']#'draft2_reion_T125']#'tuned_reion_T125','tuned_T125_MR','tuned_T125_NR','tuned_reion']
  volume=100**3#(125/cosmo['h'])**3
  ll=0
  for filename in filenames:
    ii=-1
    BHgrowth=np.zeros((len(snapshots),3))
    BHgrowth_def=np.zeros((len(snapshots),3))
    BHgrowth_MD=np.zeros((len(snapshots),3))
    BHgrowth_ID=np.zeros((len(snapshots),3))
    #BHgrowth_Coalescence=np.zeros((len(snapshots),3))
    #BHgrowth_Radio=np.zeros((len(snapshots),3))
    time_diff=np.zeros((len(snapshots),1))
    for snapshot in snapshots:
      ii+=1
      redshift[snapshot] = meraxes.io.grab_redshift(data_folder+filename+meraxes_loc,snapshot)
      redshift_next = meraxes.io.grab_redshift(data_folder+filename+meraxes_loc,snapshot-1)
       
      time_diff[ii]=abs(Planck15.age(redshift_next).value-Planck15.age(redshift[snapshot]).value)*1e9 #yr
      if filename!='dragons10':
        gals_bulges=load_data(filename,snapshot,props)
        gals_bulges_past=load_data(filename,snapshot-1,props)
      else:
        gals_bulges=load_data(filename,snapshot,props_default)
        gals_bulges_past=load_data(filename,snapshot-1,props_default)
      #gals_def=load_data(filename_def,snapshot,['BlackHoleMass'])
      #gals_def_past=load_data(filename_def,snapshot-1,['BlackHoleMass'])
      jj=0
      for BHlimit in [6]:
        gals_bulges=gals_bulges[gals_bulges['BlackHoleMass']*1e10>10**BHlimit]
        gals_bulges_past=gals_bulges_past[gals_bulges_past['BlackHoleMass']*1e10>10**BHlimit]
      
        IDs=np.zeros(len(gals_bulges),dtype=bool) 
        for kk in np.arange(0,len(gals_bulges['ID'])): IDs[kk]=gals_bulges['ID'][kk] in gals_bulges_past['ID']
        gals_bulges=gals_bulges[IDs]      
      
        IDs=np.zeros(len(gals_bulges_past),dtype=bool) 
        for kk in np.arange(0,len(gals_bulges_past['ID'])): IDs[kk]=gals_bulges_past['ID'][kk] in gals_bulges['ID']
        gals_bulges_past=gals_bulges_past[IDs]      
 

        if filename!='dragons10':
          individualBHgrowth=((gals_bulges['BlackHoleMass']-gals_bulges['BlackHoleMass_Coalescence'])-\
                          (gals_bulges_past['BlackHoleMass']-gals_bulges_past['BlackHoleMass_Coalescence']))*1e10
          #individualEddLim=gals_bulges_past['BlackHoleMass']*1e10*EddFactor
          #BHgrowth[ii,jj]=(np.sum(gals_bulges['BlackHoleMass'])-np.sum(gals_bulges_past['BlackHoleMass']))*1e10/time_diff[ii]/volume
          BHgrowth[ii,jj]=np.sum(individualBHgrowth)/time_diff[ii]/volume
          BHgrowth_MD[ii,jj]=np.sum(gals_bulges['BlackHoleMass_MD']-gals_bulges_past['BlackHoleMass_MD'])*1e10/time_diff[ii]/volume
          #BHgrowth_Coalescence[ii,jj]=(np.sum(gals_bulges['BlackHoleMass_Coalescence'])-np.sum(gals_bulges_past['BlackHoleMass_Coalescence']))*1e10/time_diff[ii]/volume
          BHgrowth_ID[ii,jj]=(np.sum(gals_bulges['BlackHoleMass_ID']-gals_bulges_past['BlackHoleMass_ID']))*1e10/time_diff[ii]/volume
          #BHgrowth_Radio[ii,jj]=(np.sum(gals_bulges['BlackHoleMass_Radio'])-np.sum(gals_bulges_past['BlackHoleMass_Radio']))*1e10/time_diff[ii]/volume
        else:
          individualBHgrowth=((gals_bulges['BlackHoleMass'])-\
                          (gals_bulges_past['BlackHoleMass']))*1e10
          BHgrowth[ii,jj]=np.sum(individualBHgrowth)/time_diff[ii]/volume


        #gals_def=gals_def[gals_def['BlackHoleMass']*1e10>10**BHlimit]
        #gals_def_past=gals_def_past[gals_def_past['BlackHoleMass']*1e10>10**BHlimit]
        #BHgrowth_def[ii,jj]=(np.sum(gals_def['BlackHoleMass'])-np.sum(gals_def_past['BlackHoleMass']))*1e10/time_diff[ii]/volume
        jj+=1
    jj=0
    lines=['-','--',':']
    for BHlimit in [6]:
      axes.plot(np.log10(1+np.array(list(redshift.values()))),BHgrowth[:,jj],lines[0],label=filename,lw=2.5,color=colors[ll])
      #axes.plot(np.log10(1+np.array(list(redshift.values()))),BHgrowth_def[:,jj],lines[jj],label='Q17',color=[0.5,0.5,0.5],lw=2.5)

      if filename!='dragons10':
        axes.plot(np.log10(1+np.array(list(redshift.values()))),BHgrowth_ID[:,jj],lines[1],label='Instability-Driven Quasar-Mode',color=colors[ll],lw=2.5)
        axes.plot(np.log10(1+np.array(list(redshift.values()))),BHgrowth_MD[:,jj],lines[2],label='Merger-Driven Quasar-Mode',color=colors[ll],lw=2.5)
      #axes.plot(np.array(list(redshift.values())),BHgrowth_Coalescence[:,jj],lines[jj],label='BH-BH Coalescence',color=colors[2],lw=2.5)
      #axes.plot(np.array(list(redshift.values())),BHgrowth_Radio[:,jj],lines[jj],label='Radio-Mode',color=colors[1],lw=2.5)
      #lgd=plt.legend(loc='upper center', bbox_to_anchor=(0.45, -0.2),ncol=2,usetex=False)
      jj+=1
    ll+=1
  
  axes.set_xlabel('log(1+z)')
  axes.set_yscale('log')
  #axes.set_xlim([2,7])
  axes.set_ylabel(r'$\log(\Delta M_{\rm{BH,~ growth~ mode}})$')
  #plt.tight_layout()
  axes.set_ylim(10**-5.9,10**-3.1)

  ax2 = axes.twiny()
  ax2.set_xlabel("z")
  #ax2.set_xticks(np.log10(1+np.linspace(0,6,7)))
  #ax2.set_xticklabels(['0','1','2','3','4','5','6'])
  z=np.array(list(redshift.values()))
  ax2.set_xticks(np.log10(1+np.array([min(z),max(z)])))
  ax2.set_xticklabels([str(min(z)),str(max(z))])
  
  ax3 = axes.twinx()
  ax3.set_ylabel(r'$3300 \log(\Delta M_{\rm{BH,~ growth~ mode}})$')
  #ax3.set_yticks(np.array([-5.9,-5.5,-5.1,-4.7,-4.3,-3.9]))
  #ax3.set_yticklabels(['-2.4','-2','-1.6','-1.2','-0.8','-0.4'])
  ax3.set_yticks(np.array([0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1]))
  ax3.set_yticklabels(['-2.4','-2','-1.6','-1.2','-0.8','-0.4','0','0.4'])
  
  #plt.savefig('/home/mmarshal/results/plots/BHGrowthWithRedshift.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

