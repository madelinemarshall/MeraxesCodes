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
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


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
  gals=gals[gals['StellarMass']*1e10>1e7]
  return gals


def plot_hist2d(xdata,ydata,axes,xlims,ylims):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='Greys',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='Greys')
  #cb=plt.colorbar(im,ax=axes,use_gridspec=True)
  #cb.set_label('N, tot N = {}'.format(np.size(xdata)))
  axes.set_aspect('auto')
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)


def plot_colorline(xdata,ydata,z_list,axes,colorbar,i):
    points = np.array([xdata, ydata]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(z_list.min(), z_list.max())
    lc = LineCollection(segments, cmap='jet', norm=norm)
    lc.set_array(z_list)
    lc.set_linewidth(0.5)
    line = axes.add_collection(lc)
    if (i==0) & colorbar:
      fig.colorbar(line, ax=ax[0][2])
      fig.colorbar(line, ax=ax[1][2])
    

if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  snap=158
  filename='tuned_reion'
  ##PLOT
  #fig,ax = plt.subplots(6,3,sharey=True)
  fig,ax = plt.subplots(2,3,gridspec_kw = {'wspace':0, 'hspace':0},sharey=True)
  #fig,axes=plt.subplots(1,1)
  gals=load_data(filename,snap)


  ##StellarMass 
  x=np.log10(gals['StellarMass']*1e10)
  y=np.log10(gals['BlackHoleMass']*1e10)
  for i in range(0,2):
    axes=ax[i]
    plot_hist2d(x,y,axes[1],[7,9.2],[4,7])

    
  gals=gals[gals['BulgeStellarMass']>0]
  x=np.log10(gals['BulgeStellarMass']*1e10)
  y=np.log10(gals['BlackHoleMass']*1e10)
  for i in range(0,2):
    axes=ax[i]
    plot_hist2d(x,y,axes[0],[7,9.2],[4,7])
  
  x=gals['BulgeStellarMass']/gals['StellarMass']
  for i in range(0,2):
    axes=ax[i]
    plot_hist2d(x,y,axes[2],[0,1],[4,7])

  snapshots=np.arange(37,159)
  EXPANSION_FACTOR_PATH = "/fred/oz013/simulations/Tiamat/a_list.txt"
  a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
  z_list = np.array(1.0/a_list - 1.0)
  z_list = z_list[snapshots][78-37:]
  #Tracks
  #track_loc='/home/mmarshal/data_dragons/histories/tuned_reion/history_byBlackHoleMass/100/BlackHoleMass07-09/'
  track_loc='/home/mmarshal/data_dragons/histories/tuned_reion/history_byStellarMass/158/StellarMass09-09'
  ##Disks
  BH=np.log10(np.fromfile(track_loc+'_disks/'+'BlackHoleMass.bin').reshape(-1,(121+1)))+10
  bulge=np.log10(np.fromfile(track_loc+'_disks/'+'BulgeStellarMass.bin').reshape(-1,(121+1)))+10
  stellar=np.log10(np.fromfile(track_loc+'_disks/'+'StellarMass.bin').reshape(-1,(121+1)))+10
  mvir=np.log10(np.fromfile(track_loc+'_disks/'+'Mvir.bin').reshape(-1,(121+1)))+10
  colors=['red','orange','yellow','green','blue','purple']
  for i in range(0,np.shape(BH)[0]):
    plot_colorline(bulge[i][78-37:],BH[i][78-37:],z_list,ax[0][0],1,i)
    plot_colorline(stellar[i][78-37:],BH[i][78-37:],z_list,ax[0][1],0,i)
    plot_colorline(10**bulge[i][78-37:]/10**stellar[i][78-37:],BH[i][78-37:],z_list,ax[0][2],0,i)
    #axes[0].plot(bulge[i],BH[i],color=color,linewidth=0.2)
  ##Bulges
  BH=np.log10(np.fromfile(track_loc+'_bulges/'+'BlackHoleMass.bin').reshape(-1,(121+1)))+10
  bulge=np.log10(np.fromfile(track_loc+'_bulges/'+'BulgeStellarMass.bin').reshape(-1,(121+1)))+10
  stellar=np.log10(np.fromfile(track_loc+'_bulges/'+'StellarMass.bin').reshape(-1,(121+1)))+10
  mvir=np.log10(np.fromfile(track_loc+'_bulges/'+'Mvir.bin').reshape(-1,(121+1)))+10
  colors=['red','orange','yellow','green','blue','purple']
  for i in range(0,np.shape(BH)[0]):
    plot_colorline(bulge[i][78-37:],BH[i][78-37:],z_list,ax[1][0],1,i)
    plot_colorline(stellar[i][78-37:],BH[i][78-37:],z_list,ax[1][1],0,i)
    plot_colorline(10**bulge[i][78-37:]/10**stellar[i][78-37:],BH[i][78-37:],z_list,ax[1][2],0,i)
    #axes[0].plot(bulge[i],BH[i],color=color,linewidth=0.2)

  ax[0][0].set_xticks([])
  ax[0][1].set_xticks([])
  ax[0][2].set_xticks([])
  ax[0][0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 
  ax[0][1].text(7, 6.5, 'Disks',weight='bold',size='large')
  ax[1][1].text(7, 6.5, 'Bulges',weight='bold',size='large')


  ax[1][2].set_xlabel(r'$\textrm{M}_{\textrm{bulge}}/\textrm{M}_\ast$')
  ax[1][1].set_xlabel(r'$\log(\textrm{M}_\ast)$')
  ax[1][0].set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}})$')
  ax[1][0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 
  #plt.legend()
  #lgd=plt.legend(fontsize='small',loc='upper center', bbox_to_anchor=(1.25, 0.8))  
  #plt.savefig('MBHBulge_z.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.tight_layout() 
  plt.show()
