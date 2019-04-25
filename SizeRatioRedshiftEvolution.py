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
from scipy.optimize import curve_fit
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (4.5,8)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4


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
                                          snapshot=snapshot,\
#props=['GhostFlag','StellarMass','StellarDiskScaleLength','BulgeStellarMass','GasDiskScaleLength'],\
                                          h=cosmo['h'],quiet=True)
  return gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-3)]


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
  axes[1].errorbar(z,r,np.array([r_low,r_up]),marker='o',color=[0.5,0.5,0.5],label='Oesch et al. (2010)',capsize=4,linestyle='None',markersize=4)
  #Ono+13
  z=[6.743,7.945]
  r=np.array([0.2843476,0.31214276])
  r_low=r-np.array([0.19584617,0.2296599])
  r_up=np.array([0.37461418,0.39462563])-r
  axes[1].errorbar(z,r,np.array([r_low,r_up]),marker='s',color=[0.5,0.5,0.5],label='Ono et al. (2013)',capsize=4,linestyle='None',markersize=4)
  #Holwerda+15 (Taken from their paper, plot digitizer)
  z=[9.4]
  z_up=[0.5]
  z_low=[0.5]
  r=np.array([0.388])
  r_low=r-np.array([0])
  r_up=np.array([1.026])-r
  axes[1].errorbar(z,r,np.array([r_low,r_up]),np.array([z_low,z_up]),marker='d',color=[0.5,0.5,0.5],label='Holwerda et al. (2015)',capsize=4,linestyle='None',markersize=4)
  #Shibuya+15 (Taken from their paper, plot digitizer)
  z=[3.8,4.9,5.9,6.8,7.9]
  r=np.array([0.6021,0.5318,0.5672,0.4528,0.3312])
  r_up=np.array([1.1438,1.0103,1.1620,1.0260,1.7654])-r
  r_low=r-np.array([0.3333,0.2835,0.2381,0.2414,0.1977])
  axes[1].errorbar(z,r,np.array([r_low,r_up]),marker='p',color=[0.5,0.5,0.5],label='Shibuya et al. (2015)',capsize=4,linestyle='None',markersize=4)
 
   #Upper L bin:
  #Oesch+10
  z=[5,6,6.8]
  r=np.array([1.092249,0.86277,0.75069886])
  r_low=r-np.array([0.951715,0.750699,0.631512])
  r_up=np.array([1.2327827,0.989072,0.8716645])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='o',color=[0.5,0.5,0.5],label='Oesch et al. (2010)',capsize=4,linestyle='None',markersize=4)
  #Bouwens+04
  z=[4.9,6]
  r=np.array([1.16190,0.8220581])
  r_low=r-np.array([1.0658449,0.6550512])
  r_up=np.array([1.2452867,0.96952015])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='*',color=[0.5,0.5,0.5],label='Bouwens et al. (2004)',capsize=4,linestyle='None',markersize=4)
  #Kawamata+15
  z=[6.5,7.94]
  r=np.array([0.54116875,0.5406563])
  r_low=r-np.array([0.43990123,0.4464954])
  r_up=np.array([0.6584278,0.96952015])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='^',color=[0.5,0.5,0.5],label='Kawamata et al. (2015)',capsize=4,linestyle='None',markersize=4)
  #Kawamata+18 (From their paper, plot digitizer)
  z=[6.2,7.8,8.5]
  r=np.array([0.9495,0.6946,0.5362])
  r_low=r-np.array([1.1218,0.9289,0.8083])
  r_up=np.array([0.8152,0.5500,0.4087])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='<',color=[0.5,0.5,0.5],label='Kawamata et al. (2018)',capsize=4,linestyle='None',markersize=4)
  #Ono+13
  z=[6.89,7.77]
  r=np.array([0.6440796,0.65620136])
  r_low=r-np.array([0.3242807,0.37016043])
  r_up=np.array([0.9585454,0.93513566])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='s',color=[0.5,0.5,0.5],label='Ono et al. (2013)',capsize=4,linestyle='None',markersize=4)
  #Shibuya+15 (Taken from their paper, plot digitizer)
  z=[3.8,4.9,5.9,6.8,7.9,10.4]
  r=np.array([0.8072,0.6695,0.5851,0.6047,0.4860,0.3834])
  r_up=np.array([1.5471,1.2049,0.9682,1.4599,1.4780,0.6616])-r
  r_low=r-np.array([0.4347,0.3605,0.3151,0.2427,0.1409,0.1033])
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='p',color=[0.5,0.5,0.5],label='Shibuya et al. (2015)',capsize=4,linestyle='None',markersize=4)
  #Holwerda+15 (Taken from their paper, plot digitizer)
  z=[9.67]
  z_up=[0.38]
  z_low=[0.36]
  r=np.array([0.5598])
  r_low=r-np.array([0.4100])
  r_up=np.array([0.7179])-r
  axes[0].errorbar(z,r,np.array([r_low,r_up]),np.array([z_low,z_up]),marker='d',color=[0.5,0.5,0.5],label='Holwerda et al. (2015)',capsize=4,linestyle='None',markersize=4)
  #Laporte+16 - from their text
  z=[7,8]
  r=np.array([0.8,0.45])
  r_low=np.array([0.18,0.15])
  r_up=np.array([0.18,0.15])
  axes[0].errorbar(z,r,np.array([r_low,r_up]),marker='x',color=[0.5,0.5,0.5],label='Laporte et al. (2016)',capsize=4,linestyle='None',markersize=4)


def plot_original_model(axes):
  #Lower L bin
  z=[5,6,7,8,9,10]
  r=np.array([0.8879714,0.65150553,0.54250836,0.40690947,0.29432926,0.2510376])
  axes[1].plot(z,r,'b--',label='DRAGONS VII')
  #Upper L bin
  z=[5.020564,6.000951,6.9679995,7.998389,9.0062475,10.005273]
  r=np.array([1.0764694,0.80429316,0.6298375,0.4766792,0.35550866,0.32513568])
  axes[0].plot(z,r,'b--',label='DRAGONS VII')


def func(z,a,b):
  return a*((1+z)/(1+7))**-b


def fit_equation(zz,rad,axes):
  popt,pcov = curve_fit(func,zz,rad)
  print("SLOPE: popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))
  #axes.plot(zz,func(zz,*popt),'--',color='green')


if __name__=="__main__":
  filename='draft2_reion'
  filename_def='dragons10'
  disk_length={filename:'StellarDiskScaleLength',filename_def:'DiskScaleLength'}
  label={filename:"M18",filename_def:"Q17"}
  linestyle={filename:"b-",filename_def:"b--"}

  #redshift={37:10,43:9,52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,250:0}
  redshift={37:10,43:9,52:8,63:7,78:6,100:5}
  snapshots=np.flip([37,43,52,63,78,100],0)
  fig,axes=plt.subplots(3,gridspec_kw = {'hspace':0})
  for fname in [filename,filename_def]:
    mean_r_upper=np.zeros(len(snapshots))
    mean_r_lower=np.zeros(len(snapshots))
    ii=0
    for snapshot in snapshots:
      gals=load_gals(fname,snapshot)
      mag=load_mags(fname,snapshot)
      gals_up=gals[(mag<-19.7)&(mag>-21)&(gals[disk_length[fname]]>0)]
      gals_low=gals[(mag>-19.7)&(mag<-18.7)&(gals[disk_length[fname]]>0)]
      upper=gals_up[disk_length[fname]]*1.678/gals_up['Rvir']
      lower=gals_low[disk_length[fname]]*1.678/gals_low['Rvir']
      
      mean_r_upper[ii]=np.nanmedian(upper)
      mean_r_lower[ii]=np.nanmedian(lower)
      ii+=1
    axes[0].plot(np.flip(np.array(list(redshift.values())),0),mean_r_upper,linestyle[fname],label=label[fname],color='black')
    axes[1].plot(np.flip(np.array(list(redshift.values())),0),mean_r_lower,linestyle[fname],label=label[fname],color='black')#'New Model')

    ##Fit curve to all data points (galaxies) not just their median at each z

  #plot_obs(axes)
  
  axes[2].axis('off')
  #axes[0].set_ylabel('$R_e$ (kpc)')
  axes[1].set_xlabel('Redshift')
  #axes[1].set_ylabel('$R_e$ (kpc)')
  axes[0].set_xticklabels([])
  axes[0].legend(ncol=2,loc=(-0.15,-1.9))
  #axes[0].set_ylim(0,1.6)
  #axes[1].set_ylim(0,1.6)
  axes[0].set_xlim(4.6,10.3)
  axes[1].set_xlim(4.6,10.3)
  #axes[0].text(8.1, 1.4, r'$(0.3-1)L^\ast_{z=3}$',weight='bold',size='large')
  #axes[1].text(8.1, 1.4, r'$(0.12-0.3)L^\ast_{z=3}$',weight='bold',size='large')
  #plt.savefig('/home/mmarshal/results/plots/SizeRedshift.pdf',format='pdf')
  plt.show()
