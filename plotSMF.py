import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF


def load_data(filename,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals


def plot_SMF(gals,prop,axes,label):
    maxval=np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10)) 
    minval=np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=30)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    axes.plot(Max,np.log10(hist/(bin_edges[1]-bin_edges[0])/100.**3),linewidth=2.0,label=label)
    return axes


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
  redshift={63:7,78:6,100:5,116:4,134:3,158:2}
  prop='StellarMass'
  
  filename='bulges_update1102_full'
  filename2='bulges_tiamat125_doubleSFefficiency'

  fig, axes = plt.subplots(2, 3)
  ii=-1
  j=0
  for snapshot in [63,78,100,116,134,158]:
    ii+=1
    if ii==3:
      j+=1
      ii=0
    gals_default=load_data('default',snapshot,prop,cosmo)
    gals_bulges=load_data(filename,snapshot,prop,cosmo)
    gals_2=load_data(filename2,snapshot,prop,cosmo)
  
    plot_SMF(gals_default,prop,axes[j,ii],'Default Meraxes')
    plot_SMF(gals_bulges,prop,axes[j,ii],'Bulge Model')
    plot_SMF(gals_2,prop,axes[j,ii],'Bulge Model, Larger Stellar Mass')
    if snapshot==63:
      plt.legend()
    plot_obsGSMF(axes[j,ii],redshift[snapshot],hubble_h=cosmo['h'],markersize=7,legend=True,silent=False,color=[0.5,0.5,0.5],alpha=1.0)

    axes[j,ii].set_xlabel('log(Mass)')
    axes[j,ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    axes[j,ii].set_title('Redshift {}'.format(redshift[snapshot]))
    axes[j,ii].set_xlim([7.5,11])
    axes[j,ii].set_ylim([-7,-1])
  plt.tight_layout()
  plt.show()

