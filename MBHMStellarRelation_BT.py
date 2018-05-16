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
import ContourPlot as cp
import matplotlib.lines as mlines

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (8,4)
#matplotlib.rcParams['figure.figsize'] = (7.2,4)
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
  gals=gals[gals['BlackHoleMass']*1e10>1e3]
  gals=gals[gals['StellarMass']*1e10>1e8]
  return gals


def plot_observations(axes,color):
  #Kormendy & Ho (2013):
  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
  obs=axes.errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)',capsize=3,linewidth=2.5, zorder=101,color=color[213])
  return obs


def plot_MBHMstellar(filename,snapshots,mass_bulge,color,axes):
  slope=np.zeros(len(snapshots))
  slope_errs=np.zeros(len(snapshots))
  inter=np.zeros(len(snapshots))
  inter_errs=np.zeros(len(snapshots))
  ii=0
  for snap in snapshots:
    gals=load_data(filename,snap)

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
    cp.contour_plot(logMstel_bulge,logMBH_bulge,xlab=None,ylab=None,xlims=[8.15,11.6],ylims=[4,9.5],axes=axes,colors='pink',levels=np.logspace(-2,1,7),linewidth=0.9)
    cp.contour_plot(logMstel_disk,logMBH_disk,xlab=None,ylab=None,xlims=[8.15,11.6],ylims=[4,9.5],axes=axes,colors='lightblue',levels=np.logspace(-2,1,7),linewidth=0.9)
    ii+=1
  axes.set_xlim([8.15,11.6])
  axes.set_ylim([4,9.5])


if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,213:0.55}
  snapshots=[213]
  prop='StellarMass'
  color={52:'C0',63:'C1',78:'C2',100:'C3',116:'C4',134:'aqua',158:'pink',213:'k'}
  filename='tuned_t125'
  meraxes_loc='/output/meraxes.hdf5'
  #filename='tune_higherseed'
  #meraxes_loc='/output/run1/meraxes_001.hdf5'
  
  ##PLOT
  fig,ax = plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)
  #fig,axes=plt.subplots(1,1)
  plot_MBHMstellar(filename,snapshots,True,color,ax[0])
  obs=plot_observations(ax[0],color)
  pink_line = mlines.Line2D([], [], color='pink', label='Bulge Dominated Galaxies')
  blue_line = mlines.Line2D([], [], color='lightblue', label='Disk Dominated Galaxies')
  hand=[pink_line,blue_line,obs]
  lgd=ax[0].legend(handles=hand,fontsize='small',loc='upper left')  
  
  plot_MBHMstellar(filename,snapshots,False,color,ax[1])
  plot_observations(ax[1],color)
  ax[1].set_xlabel(r'$\log(\textrm{M}_\ast)$')
  ax[0].set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}})$')

  #plt.legend()
  ax[0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 
  plt.savefig('MBHMStellarRelation_BT.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

  ##PERFORM STATISTICS
  #find_fit()
