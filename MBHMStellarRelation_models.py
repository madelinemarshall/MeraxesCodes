import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF
from scipy.optimize import curve_fit


#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,4)
matplotlib.rcParams['lines.linewidth'] = 2.5
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,snapshot):
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
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','BulgeStellarMass'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['BlackHoleMass']*1e10>1e6]
  gals=gals[gals['StellarMass']*1e10>1e6]
  return gals


def plot_observations(snapshot,axes):
  #Kormendy & Ho (2013):
  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**7,10**12])/10**11)+9
  axes.errorbar([7,12],logMBH,yerr=0.28,linestyle='-',color='darkorange',label='Kormendy \& Ho (2013)',capsize=3,linewidth=2.5, zorder=101)

  logMBH=(8.43+1.16*(np.log10(np.array([1e7,1e12])/1e11)))
  HU=axes.plot([7,12],logMBH,'-',color='gold',linewidth=2.5,label="BlueTides: Huang et al. (2018)", zorder=102)
  ##Sijaki+15 Illustris sims
  logMBH=1.23*(np.log10(np.array([1e8,1e12])))-4.85
  axes.plot([8,12],logMBH,'--',color='gold',linewidth=2.5,label="Illustris: Sijaki et al. (2015)",zorder=102)
  logMBH=1.28*(np.log10(np.array([1e8,1e12])))-5.04
  axes.plot([8,12],logMBH,'--',color='gold',linewidth=2.5,zorder=102)


def func(x,a,c):
  return a*x+c


if __name__=='__main__':
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  #snapshots=[52,63,78,100,116,134,158]
  snapshots=[78]
  prop='StellarMass'
  #color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'black',213:'pink'}
  color={0:'C0',1:'C1',2:'C2',3:'C3',4:'C4',5:'black',6:'pink'}
 
  #filenames={0:'bulges_correctBHMF',1:'bulges_noreion',2:'bulges_2edd',3:'bulges_0p5edd',4:'bulges_lowseed',5:'bulges_highseed',6:'bulges_nodiskinstability'} 
  filenames={0:'tuned_reion'}
  models={0:'New model',1:'No reionization',2:'Twice Eddington',3:'Half Eddington',4:'Lighter BH seed',5:'Heavier BH seed',6:'No disk instabilities'} 
  #fig,axes = plt.subplots(len(filenames),2,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)
  fig,axes = plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)

  for snap in snapshots:
    for ii in range(0,len(filenames)):
      gals=load_data(filenames[ii],snap)
      #df=pd.DataFrame(gals)
      #grouped=df.groupby('CentralGal')
      #grouped_sum=grouped.sum()
      #Mstellar=np.array(grouped_sum['Mstellar']*1e10)
      #MBH=np.array(grouped_sum['StellarMass']*1e10)
      Mstellar=gals['StellarMass']*1e10
      Mbulge=gals['BulgeStellarMass']*1e10
      MBH=gals['BlackHoleMass']*1e10
      logMstellar=np.log10(Mstellar)
      logMbulge=np.log10(Mbulge)
      logMBH=np.log10(MBH) 
   
      #Total stellar mass
      bin_width=0.3
      min_mass=np.min(logMstellar)
      max_mass=np.max(logMstellar)
      n_bins=np.int((max_mass-min_mass)/bin_width)
      med_bh=np.zeros(n_bins)
      middle_sm=np.zeros(n_bins)
      pctl_bh=np.zeros((n_bins,2))
      for nn in range(0,n_bins):
        if np.size(MBH[(logMstellar>min_mass+(nn*bin_width))&(logMstellar<min_mass+(nn+1)*bin_width)])>5:
          med_bh[nn]=np.median(MBH[(logMstellar>min_mass+(nn*bin_width))&(logMstellar<min_mass+(nn+1)*bin_width)])
          pctl_bh[nn,:]=np.percentile(MBH[(logMstellar>min_mass+(nn*bin_width))&(logMstellar<min_mass+(nn+1)*bin_width)],[16,84])
          middle_sm[nn]=min_mass+(nn+0.5)*bin_width
        else:
          med_bh[nn]=np.nan
      #if snap==158:
      #  axes[0].plot(middle_sm,np.log10(pctl_bh[:,0]),'k:')
      #  axes[0].plot(middle_sm,np.log10(pctl_bh[:,1]),'k:')
      axes[0].plot(middle_sm,np.log10(med_bh),label='{}'.format(models[ii]),color=color[ii])
      middle_sm=middle_sm[np.logical_not(np.isnan(med_bh))]
      med_bh=med_bh[np.logical_not(np.isnan(med_bh))]
      #popt,pcov = curve_fit(func,middle_sm,np.log10(med_bh))
      #print("snap {}, popt {}, perr {}".format(snap,popt,np.sqrt(np.diag(pcov))))
      #plt.plot(middle_sm,func(middle_sm,*popt),'--')


      ##Bulge Mass
      gals=gals[gals['BulgeStellarMass']>0]
      Mbulge=gals['BulgeStellarMass']*1e10
      MBH=gals['BlackHoleMass']*1e10
      logMbulge=np.log10(Mbulge)
      logMBH=np.log10(MBH) 
      #Bulge stellar mass
      bin_width=0.3
      max_mass=np.max(logMbulge)
      n_bins=np.int((max_mass-min_mass)/bin_width)
      med_bh=np.zeros(n_bins)
      middle_sm=np.zeros(n_bins)
      for nn in range(0,n_bins):
        if np.size(MBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)])>5:
          med_bh[nn]=np.median(MBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)])
          middle_sm[nn]=min_mass+(nn+0.5)*bin_width
        else:
          med_bh[nn]=np.nan
      axes[1].plot(logMbulge,logMBH,'.k')
      axes[1].plot(middle_sm,np.log10(med_bh),label='$z={}$'.format(redshift[snap]),color=color[ii])
      middle_sm=middle_sm[np.logical_not(np.isnan(med_bh))]
      med_bh=med_bh[np.logical_not(np.isnan(med_bh))]
      popt,pcov = curve_fit(func,middle_sm,np.log10(med_bh))
      #print("snap {}, popt {}, perr {}".format(snap,popt,np.sqrt(np.diag(pcov))))
      plt.plot(middle_sm,func(middle_sm,*popt),'--')
      popt,pcov = curve_fit(func,logMbulge,logMBH)
      #print("snap {}, popt {}, perr {}".format(snap,popt,np.sqrt(np.diag(pcov))))
      plt.plot(np.array([9,12]),func(np.array([9,12]),*popt),'--')
      #axes[0].set_xlim([6,12])
      axes[0].set_ylabel(r'$\log(M_{BH})$')
      #axes[1].set_xlim([6,12])
      #axes[0].text(6, 7, ' {}'.format(str(models[ii])),usetex=False)

  axes[0].legend()
  axes[0].set_xlabel(r'$\log(M_\ast)$')
  axes[1].set_xlabel(r'$\log(M_{\textrm{bulge}})$')
  plt.show()
