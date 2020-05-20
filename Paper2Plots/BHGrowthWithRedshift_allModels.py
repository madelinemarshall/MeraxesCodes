import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
from astropy.cosmology import Planck15

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.6,4)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4

data_folder='/home/mmarshal/data_dragons/'
cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  } 

def load_data(filename,snapshot,prop):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals

def load_data_def(filename_def,snapshot,prop):
  gals=meraxes.io.read_gals(data_folder+filename_def,\
      snapshot=snapshot,\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals

def plot_obs(axes):
  #Delvecchio
  X=np.array([0.21,0.53,0.95,1.46,2.02,2.75])
  Y=10**np.array([-5.25,-4.75,-4.55,-4.39,-4.17,-4.50])
  Y_low=Y-10**np.array([-5.54,-4.78,-4.6,-4.59,-4.40,-4.76])
  Y_up=10**np.array([-5.19,-4.65,-4.50,-4.21,-3.99,-4.27])-Y
  X_low=X-np.array([0.1,0.3,0.7,1.2,1.8,2.5])
  X_up=np.array([0.3,0.7,1.2,1.8,2.5,3.8])-X
  axes.errorbar(np.log10(1+X),Y,yerr=[Y_low,Y_up],xerr=[np.log10(1+X_low),np.log10(1+X_up)],fmt='o',color=colors[1],label='Delvecchio et al. (2014)')
  ##Vita+18
  #X=np.array([3.25,3.75,5.0,5.85])
  #Y=np.array([1.93e-5,7.76e-6,1.21e-6,7.68e-7])
  #X_low=np.array([3.00,3.50,4.50,5.50])
  #X_up=np.array([3.50,4.50,5.50,6.0])
  #Y_low=np.array([2.8e-5,1.1e-5,1.7e-6,1.15e-6])
  #Y_up=np.array([1.25e-5,5.05e-6,7.7e-7,4.7e-7])
  #axes.errorbar(np.log10(1+X),Y,yerr=[Y_low,Y_up],xerr=[np.log10(1+X_low),np.log10(1+X_up)],fmt='o',color=colors[3],label='Vito et al. (2018)')


if __name__=="__main__":
  #snapshots=range(63,250,5) #every 5 snapshots
  #snapshots=range(63,158,15) #every 5 snapshots
  snapshots=np.array([52,63,78,100,116,134,158])
  redshift={}
  props=['ID','BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag','BlackHoleMass_Coalescence','BlackHoleMass_Radio']
  props_default=['ID','BlackHoleMass','GhostFlag']
  fig, axes = plt.subplots(1, 1)
  #filename_def='dragons10'
  #filenames=['dragons10','paper2','paper2_low_edd']#'paper1_T125_both']#'draft2_reion_T125']#'tuned_reion_T125','tuned_T125_MR','tuned_T125_NR','tuned_reion']
  filename='tuning_Tiamat/output/'
  meraxes_locs=['meraxes_00'+str(f)+'.hdf5' for f in [0,2,4,6]]
  lines=['-','--','-.',':']
  label=['$k_c=0.005$','$k_c=0.01$','$k_c=0.03$','$k_c=0.09$']
  volume=100**3#(125/cosmo['h'])**3

  filename_def='paper1/output/meraxes.hdf5'
  #plot_SMF(gals_default,prop,vol_def,ax,**{'linestyle':'-','label':'Q17 Meraxes\n(Tiamat)','linewidth':2,'color':[0.35,0.35,0.35],'zorder':999})

  ll=0
    
  ii=-1
  BHgrowth_def=np.zeros((len(snapshots),1))
  time_diff_def=np.zeros((len(snapshots),1))
  for snapshot in snapshots:
      ii+=1
      redshift[snapshot] = meraxes.io.grab_redshift(data_folder+filename_def,snapshot)
      redshift_next = meraxes.io.grab_redshift(data_folder+filename_def,snapshot-1)
       
      time_diff_def[ii]=abs(Planck15.age(redshift_next).value-Planck15.age(redshift[snapshot]).value)*1e9 #yr
      gals_bulges=load_data_def(filename_def,snapshot,props_default)
      gals_bulges_past=load_data_def(filename_def,snapshot-1,props_default)
      for BHlimit in [6]:
        gals_bulges=gals_bulges[gals_bulges['BlackHoleMass']*1e10>10**BHlimit]
        gals_bulges_past=gals_bulges_past[gals_bulges_past['BlackHoleMass']*1e10>10**BHlimit]
      
        IDs=np.zeros(len(gals_bulges),dtype=bool) 
        for kk in np.arange(0,len(gals_bulges['ID'])): IDs[kk]=gals_bulges['ID'][kk] in gals_bulges_past['ID']
        gals_bulges=gals_bulges[IDs]      
      
        IDs=np.zeros(len(gals_bulges_past),dtype=bool) 
        for kk in np.arange(0,len(gals_bulges_past['ID'])): IDs[kk]=gals_bulges_past['ID'][kk] in gals_bulges['ID']
        gals_bulges_past=gals_bulges_past[IDs]      
 

        individualBHgrowth_def=((gals_bulges['BlackHoleMass'])-\
                          (gals_bulges_past['BlackHoleMass']))*1e10
        BHgrowth_def[ii]=np.sum(individualBHgrowth_def)/time_diff_def[ii]/volume
  axes.plot(np.log10(1+np.array(list(redshift.values()))),BHgrowth_def,'-',label='M19 Meraxes',color=colors[2],lw=2)
  #plot_SMF(gals_default,prop,vol_def,ax,**{'linestyle':'-','label':'Q17 Meraxes\n(Tiamat)','linewidth':2,'color':[0.35,0.35,0.35],'zorder':999})


  BHgrowth=np.zeros((len(snapshots),3))
  time_diff=np.zeros((len(snapshots),1)) 
  for meraxes_loc in meraxes_locs:
    ii=-1
    BHgrowth=np.zeros((len(snapshots),3))
    #BHgrowth_def=np.zeros((len(snapshots),3))
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
        else:
          individualBHgrowth=((gals_bulges['BlackHoleMass'])-\
                          (gals_bulges_past['BlackHoleMass']))*1e10
          BHgrowth[ii,jj]=np.sum(individualBHgrowth)/time_diff[ii]/volume


        #gals_def=gals_def[gals_def['BlackHoleMass']*1e10>10**BHlimit]
        #gals_def_past=gals_def_past[gals_def_past['BlackHoleMass']*1e10>10**BHlimit]
        #BHgrowth_def[ii,jj]=(np.sum(gals_def['BlackHoleMass'])-np.sum(gals_def_past['BlackHoleMass']))*1e10/time_diff[ii]/volume
        jj+=1
    jj=0
    for BHlimit in [6]:
      axes.plot(np.log10(1+np.array(list(redshift.values()))),BHgrowth[:,jj],lines[ll],lw=2.5,color='k',label=label[ll])

      jj+=1
    ll+=1
  plot_obs(axes)
  plt.legend(loc='lower left',fontsize='small')
  
  axes.set_xlabel('log(1+z)')
  axes.set_yscale('log')
  axes.set_xlim([0.35,1])
  axes.set_ylabel(r'$\Psi_{BHAR}/M_\odot\textrm{yr}^{-1} \textrm{Mpc}^{-3}$')
  #plt.tight_layout()
  axes.set_ylim(10**-5.9,10**-3.3)

  ax2 = axes.twiny()
  ax2.set_xlabel("z")
  z=np.array(list(redshift.values()))
  ax2.set_xticks(np.log10(1+np.array([1,2,3,4,5,6,7,8])))
  ax2.set_xticklabels(['1','2','3','4','5','6','7','8'])
  ax2.set_xlim([0.35,1])
  #ax2.set_xticks(np.log10(1+np.linspace(0,6,7)))
  #ax2.set_xticklabels(['0','1','2','3','4','5','6'])
 
  plt.tight_layout() 
  plt.savefig('/home/mmarshal/results/plots/Paper2/BHGrowthWithRedshift.pdf',format='pdf')
  plt.show()

