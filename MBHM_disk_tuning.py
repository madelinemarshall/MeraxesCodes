import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/simulation_codes/Yuxiang/')
sys.path.append('/home/mmarshal/simulation_codes')
from _plot_obsGSMF import plot_obsGSMF
from scipy.optimize import curve_fit
import ContourPlot as cp
import matplotlib.lines as mlines


#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (6.5,3.2)
matplotlib.rcParams['lines.linewidth'] = 2.5
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def load_data(filename,snapshot):
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
      snapshot=snapshot,\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals

def plot_observations(axes,color):
  #Kormendy & Ho (2013):
  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
  obs=axes.errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)',capsize=3,linewidth=2.5, zorder=101,color=color[250])
  return obs


def plot_MBHMstellar(filename,snapshots,mass_bulge,color,axes):
  slope=np.zeros(len(snapshots))
  slope_errs=np.zeros(len(snapshots))
  inter=np.zeros(len(snapshots))
  inter_errs=np.zeros(len(snapshots))
  ii=0
  for snap in snapshots:
    gals=load_data(filename,snap)
    gals=gals[gals['BlackHoleMass']*1e10>10**5.8]
    if mass_bulge==0:
      Mstellar=gals['StellarMass']*1e10
      Mbulge=gals['BulgeStellarMass']*1e10
      BT=Mbulge/Mstellar
      Mstel=Mstellar
    else: ##Bulge Mass
      gals=gals[gals['BulgeStellarMass']>0]
      Mstellar=gals['StellarMass']*1e10
      Mbulge=gals['BulgeStellarMass']*1e10
      BT=Mbulge/Mstellar
      Mstel=Mbulge      

    MBH=gals['BlackHoleMass']*1e10
    logMstel=np.log10(Mstel)
    logMBH=np.log10(MBH) 
    #sc=axes.scatter(logMstel,logMBH,c=BT,marker='.',s=2)
    
    logMBH_disk=logMBH[BT<0.3]
    logMBH_bulge=logMBH[BT>0.7]
    logMstel_disk=logMstel[BT<0.3]
    logMstel_bulge=logMstel[BT>0.7]
    cp.contour_plot(logMstel_bulge,logMBH_bulge,xlab=None,ylab=None,xlims=[8,11.6],ylims=[6,9.2],axes=axes,colors='#e41a1c',linewidth=0.9,levels=[0.2,0.4,0.6,0.8])
    cp.contour_plot(logMstel_disk,logMBH_disk,xlab=None,ylab=None,xlims=[8,11.6],ylims=[6,9.2],axes=axes,colors='#377eb8',linewidth=0.9,levels=[0.2,0.4,0.6,0.8])#levels=np.logspace(-2,1,7)
    axes.plot(logMstel_disk,logMBH_disk,'.',color='#377eb8')
    ii+=1
  axes.set_xlim([8.05,11.6])
  axes.set_ylim([6,9.2])


def plot_MBHMstellar_line(filename,snapshots,mass_bulge,color,axes):
  slope=np.zeros(len(snapshots))
  slope_errs=np.zeros(len(snapshots))
  inter=np.zeros(len(snapshots))
  inter_errs=np.zeros(len(snapshots))
  ii=0
  for snap in snapshots:
    gals=load_data(filename,snap,['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','BulgeStellarMass'])
    Mstellar=gals['StellarMass']*1e10
    if mass_bulge==0:
      Mstel=Mstellar
    else: ##Bulge Mass
      gals=gals[gals['BulgeStellarMass']>0]
      Mbulge=gals['BulgeStellarMass']*1e10
      Mstel=Mbulge
      
    MBH=gals['BlackHoleMass']*1e10
    logMstel=np.log10(Mstel)
    logMBH=np.log10(MBH) 
    #Bulge stellar mass
    bin_width=0.5
    min_mass=np.min(logMstel)
    max_mass=np.max(logMstel)
    n_bins=np.int((max_mass-min_mass)/bin_width)
    med_bh=np.zeros(n_bins)
    middle_sm=np.zeros(n_bins)
    pctl_bh=np.zeros((n_bins,2))
    for nn in range(0,n_bins):
      if np.size(logMBH[(logMstel>min_mass+(nn*bin_width))&(logMstel<min_mass+(nn+1)*bin_width)])>10:
        med_bh[nn]=np.median(logMBH[(logMstel>min_mass+(nn*bin_width))&(logMstel<min_mass+(nn+1)*bin_width)])
        middle_sm[nn]=min_mass+(nn+0.5)*bin_width
        pctl_bh[nn,:]=np.percentile(logMBH[(logMstel>min_mass+(nn*bin_width))&(logMstel<min_mass+(nn+1)*bin_width)],[16,84])
      else:
        med_bh[nn]=np.nan 
    middle_sm=middle_sm[np.logical_not(np.isnan(med_bh))]
    med_bh=med_bh[np.logical_not(np.isnan(med_bh))]
    axes.plot(middle_sm,med_bh,label=r'All Galaxies'.format(redshift[snap]),color=color[snap])



if __name__=='__main__':
  ##SETUP
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,213:0.55,250:0}
  snapshots=[250]
  prop='StellarMass'
  color={52:'C0',63:'C1',78:'C2',100:'C3',116:'C4',134:'aqua',158:'pink',213:'k',250:'k'}
  data_folder='/home/mmarshal/data_dragons/'
  filename='tuning/'
  meraxes_loc=str(sys.argv[1])
  
  ##PLOT
  fig,ax = plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)
  plot_MBHMstellar(filename,snapshots,True,color,ax[0])
  #plot_MBHMstellar_line(filename,snapshots,True,color,ax[0])
  pink_line = mlines.Line2D([], [], color='#e41a1c', label='Bulge Dominated Galaxies')
  blue_line = mlines.Line2D([], [], color='#377eb8', label='Disc Dominated Galaxies')
  hand=[pink_line,blue_line]#,obs]
  lgd=ax[0].legend(handles=hand,fontsize='small',loc='upper left')  
  
  plot_MBHMstellar(filename,snapshots,False,color,ax[1])
  #plot_MBHMstellar_line(filename,snapshots,False,color,ax[1])
  ax[1].set_xlabel(r'$\log(\textrm{M}_\ast/M_\odot)$')
  ax[0].set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}}/M_\odot)$')

  ax[0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}}/M_\odot)$') 
  plt.tight_layout()
  plt.savefig('/home/mmarshal/results/plots/tuning_paper2/BHMF_{}.pdf'.format(int(meraxes_loc[-8:-5])),format='pdf')
  plt.show()

  ##PERFORM STATISTICS
  #find_fit()
