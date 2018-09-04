##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import ContourPlot as cp
import pandas as pd
import magcalc as mc

import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,5)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_gals(filename,snapshot):
  #Setup
  cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  }
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,props=['GhostFlag','StellarMass','StellarDiskScaleLength','BulgeStellarMass','GasDiskScaleLength'],\
                                          h=cosmo['h'],quiet=True)
  return gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)]


def load_mags(filename,snapshot):
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def plot_obs(axes):
  #Lower L bin:
  #Oesch+10
  z=[5,6,6.8,8]
  r=np.array([0.71097815,0.7700958,0.48959744,0.38306573])
  r_low=r-np.array([0.6490305,0.6285012,0.33914837,0.22696671])
  r_up=np.array([0.7693811,0.9099205,0.638267,0.528518])-r
  plt.errorbar(z,r,np.array([r_low,r_up]),marker='o',color='k',label='Oesch et al. (2010)',capsize=6,linestyle='None')
  #Ono+13
  z=[6.743,7.945]
  r=np.array([0.2843476,0.31214276])
  r_low=r-np.array([0.19584617,0.2296599])
  r_up=np.array([0.37461418,0.39462563])-r
  axes[1].errorbar(z,r,np.array([r_low,r_up]),marker='s',color='k',label='Ono et al. (2013)',capsize=6,linestyle='None')
  #Holwerda+15
  z=[9.4]
  z_up=[0.5]
  z_low=[0.5]
  r=np.array([0.379])
  r_low=r-np.array([0])
  r_up=np.array([0.997])-r
  axes[1].errorbar(z,r,np.array([r_low,r_up]),np.array([z_low,z_up]),marker='d',color='k',label='Holwerda et al. (2015)',capsize=6,linestyle='None')
 
  #Upper L bin:
  #Oesch+10
  z=[5,6,6.8]
  r=np.array([1.092249,0.86277,0.75069886])
  r_low=r-np.array([0.951715,0.750699,0.631512])
  r_up=np.array([1.2327827,0.989072,0.8716645])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='o',color='k',label='Oesch et al. (2010)',capsize=6,linestyle='None')
  #Bouwens+04
  z=[4.9,6]
  r=np.array([1.16190,0.8220581])
  r_low=r-np.array([1.0658449,0.6550512])
  r_up=np.array([1.2452867,0.96952015])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='*',color='k',label='Bouwens et al. (2004)',capsize=6,linestyle='None')
  #Kawamata+15
  z=[6.5,7.94]
  r=np.array([0.54116875,0.5406563])
  r_low=r-np.array([0.43990123,0.4464954])
  r_up=np.array([0.6584278,0.96952015])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='^',color='k',label='Kawamata et al. (2015)',capsize=6,linestyle='None')
  #Ono+13
  z=[6.89,7.77]
  r=np.array([0.6440796,0.65620136])
  r_low=r-np.array([0.3242807,0.37016043])
  r_up=np.array([0.9585454,0.93513566])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='s',color='k',label='Ono et al. (2013)',capsize=6,linestyle='None')
  #Shibuya+15
  z=[10.4]
  r=np.array([0.37633])
  r_up=np.array([0.648160])-r
  r_low=r-np.array([0.100951])
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='p',color='k',label='Shibuya et al. (2015)',capsize=6,linestyle='None')
  #Holwerda+15 
  z=[9.68]
  z_up=[0.38]
  z_low=[0.36]
  r=np.array([0.536])
  r_low=r-np.array([0.3946])
  r_up=np.array([0.689])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),np.array([z_low,z_up]),marker='d',color='k',label='Holwerda et al. (2015)',capsize=6,linestyle='None')


def plot_original_model(axes):
  #Lower L bin
  z=[5,6,7,8,9,10]
  r=np.array([0.8879714,0.65150553,0.54250836,0.40690947,0.29432926,0.2510376])
  axes[1].plot(z,r,'b--',label='DRAGONS VII')
  #Upper L bin
  z=[5.020564,6.000951,6.9679995,7.998389,9.0062475,10.005273]
  r=np.array([1.0764694,0.80429316,0.6298375,0.4766792,0.35550866,0.32513568])
  axes[0].plot(z,r,'b--',label='DRAGONS VII')



if __name__=="__main__":
  filename='tuned_reion_T125'
  redshift={37:10,43:9,52:8,63:7,78:6,100:5}
  snapshots=np.flip([37,43,52,63,78,100],0)
  mean_r_upper=np.zeros(6)
  mean_r_lower=np.zeros(6)
  ii=0
  for snapshot in snapshots:
    gals=load_gals(filename,snapshot)
    mag=load_mags(filename,snapshot)
    BT=gals['BulgeStellarMass']/gals['StellarMass']
    upper=gals[(mag<-19.7)&(mag>-21)&(gals['StellarDiskScaleLength']>0)]['StellarDiskScaleLength']*1000*1.678
    lower=gals[(mag>-19.7)&(mag<-18.7)&(gals['StellarDiskScaleLength']>0)]['StellarDiskScaleLength']*1000*1.678
    #plt.plot(mag[(mag<-19.7)&(mag>-21)&(BT<0.3)&(gals['StellarDiskScaleLength']>0)],upper,'.')
    #plt.plot([-24,-10],[np.nanmean(upper),np.nanmean(upper)],'-')
    #plt.hist(gals[(gals['StellarDiskScaleLength']>0)&(BT<0.3)]['StellarDiskScaleLength']*1000,range=[0,1])
    #plt.xlim([-24,-10])
    #plt.ylim([-1.5,1])
    #plt.show()
    mean_r_upper[ii]=np.nanmean(upper)
    mean_r_lower[ii]=np.nanmean(lower)
    ii+=1
  fig,axes=plt.subplots(2,1,gridspec_kw = {'hspace':0})
  axes[0].plot(np.flip(np.array(list(redshift.values())),0),mean_r_upper,'b',label='New Model')
  axes[1].plot(np.flip(np.array(list(redshift.values())),0),mean_r_lower,'b',label='New Model')
  plot_original_model(axes)
  plot_obs(axes)
  axes[0].set_ylabel('$R_e$ (kpc)')
  axes[1].set_xlabel('Redshift')
  axes[1].set_ylabel('$R_e$ (kpc)')
  axes[0].set_xticklabels([])
  axes[0].legend(ncol=2,loc=1)
  axes[1].legend(ncol=2,loc=1)
  axes[0].set_ylim(0,1.55)
  axes[1].set_ylim(0,1.55)
  axes[0].set_xlim(4.6,10.3)
  axes[1].set_xlim(4.6,10.3)
  axes[0].text(5, 0.1, r'$(0.3-1)L^\ast_{z=3}$',weight='bold',size='x-large')
  axes[1].text(5, 0.1, r'$(0.12-0.3)L^\ast_{z=3}$',weight='bold',size='x-large')
  plt.show()
