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


def plot_hist2d(xdata,ydata,axes,xlims,ylims):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='BuPu',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu')
  #cb=plt.colorbar(im,ax=axes,use_gridspec=True)
  #cb.set_label('N, tot N = {}'.format(np.size(xdata)))
  axes.set_aspect('auto')
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)



if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  snap=78
  filename='tune_tiamat'
  ##PLOT
  fig,ax = plt.subplots(1,3,gridspec_kw = {'wspace':0, 'hspace':0},sharey=True)
  #fig,axes=plt.subplots(1,1)
  gals=load_data(filename,snap)


  ##StellarMass 
  x=np.log10(gals['StellarMass']*1e10)
  y=np.log10(gals['BlackHoleMass']*1e10)
  plot_hist2d(x,y,ax[1],[8,np.nanmax(x)],[4,9.5])
  
  x=np.log10(gals['Mvir']*1e10)
  plot_hist2d(x,y,ax[2],[10,np.nanmax(x)],[4,9.5])
    
  gals=gals[gals['BulgeStellarMass']>0]
  x=np.log10(gals['BulgeStellarMass']*1e10)
  y=np.log10(gals['BlackHoleMass']*1e10)
  plot_hist2d(x,y,ax[0],[4,np.nanmax(x)],[4,9.5])

  ax[2].set_xlabel(r'$\log(\textrm{M}_{\textrm{vir}})$')
  ax[1].set_xlabel(r'$\log(\textrm{M}_\ast)$')
  ax[0].set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}})$')

  plt.legend()
  lgd=plt.legend(fontsize='small',loc='upper center', bbox_to_anchor=(1.25, 0.8))  
  ax[0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 
  #plt.savefig('MBHBulge_z.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
