import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF
import magcalc as mc
import ContourPlot as cp

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (3.5,3)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def load_data(filename,snapshot,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['StellarMass','GhostFlag','BulgeStellarMass','Type','MergerBulgeStellarMass'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['StellarMass']*1e10>1e9]#**8.9]
  #gals=gals[gals['Type']==0]
  return gals


def plot_component_over_total(gals,axes):
  IDB=(gals['BulgeStellarMass']-gals['MergerBulgeStellarMass'])*1e10
  MDB=gals['MergerBulgeStellarMass']*1e10
  disk=(gals['StellarMass']-gals['BulgeStellarMass'])*1e10
  SM=gals['StellarMass']*1e10
  #cp.contour_plot(np.log10(SM),disk/SM,xlab=None,ylab=None,xlims=[9,12.5],ylims=[0,1],axes=axes[0])
  axes[0].plot(SM,disk/SM,'b.',label='disks',markersize=0.2)
  axes[1].plot(SM,IDB/SM,'g.',label='IDB',markersize=0.2)
  axes[2].plot(SM,MDB/SM,'r.',label='MDB',markersize=0.2)
  for ii in [0,1,2]:
    axes[ii].set_xscale('log')
    axes[ii].set_xlim([1e9,10**12.5])
    axes[ii].set_ylabel('Component / Total')
    axes[ii].set_xlabel('Stellar Mass')
    axes[ii].legend()


def plot_component_over_total_bulge(gals,axes):
  bulge=gals['BulgeStellarMass']*1e10
  disk=(gals['StellarMass']-gals['BulgeStellarMass'])*1e10
  SM=gals['StellarMass']*1e10
  #cp.contour_plot(np.log10(SM),disk/SM,xlab=None,ylab=None,xlims=[9,12.5],ylims=[0,1],axes=axes[0])
  axes[0].plot(SM,disk/SM,'b.',label='disk',markersize=0.2)
  axes[1].plot(SM,bulge/SM,'g.',label='bulge',markersize=0.2)
  for ii in [0,1]:
    axes[ii].set_xscale('log')
    axes[ii].set_xlim([1e9,10**12.5])
    axes[ii].set_ylabel('Component / Total')
    axes[ii].set_xlabel('Stellar Mass')
    axes[ii].legend()


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
  snapshot=158
  filename='bulges_IDBH_tiamat125_best'
  #filename='bulges_mod_merger_frac'
  total_bulge=0
  gals_bulges=load_data(filename,snapshot,cosmo)
   
  if total_bulge==0: 
    fig,axes=plt.subplots(3,1)
    plot_component_over_total(gals_bulges,axes)
  else:
    fig,axes=plt.subplots(2,1)
    plot_component_over_total_bulge(gals_bulges,axes)
  #plt.savefig('DiskFraction.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
