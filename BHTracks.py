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
      fig.colorbar(line, ax=ax[2][2])
      fig.colorbar(line, ax=ax[3][2])
      fig.colorbar(line, ax=ax[4][2])
      fig.colorbar(line, ax=ax[5][2])
    

if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  snap=158
  filename='tuned_reion'
  ##PLOT
  #fig,ax = plt.subplots(6,3,sharey=True)
  fig,ax = plt.subplots(6,3,gridspec_kw = {'wspace':0, 'hspace':0},sharey=True)
  #fig,axes=plt.subplots(1,1)
  gals=load_data(filename,snap)


  ##StellarMass 
  x=np.log10(gals['StellarMass']*1e10)
  y=np.log10(gals['BlackHoleMass']*1e10)
  for i in range(0,6):
    axes=ax[i]
    plot_hist2d(x,y,axes[1],[7,np.nanmax(x)],[4,9.5])

    
  gals=gals[gals['BulgeStellarMass']>0]
  x=np.log10(gals['BulgeStellarMass']*1e10)
  y=np.log10(gals['BlackHoleMass']*1e10)
  for i in range(0,6):
    axes=ax[i]
    plot_hist2d(x,y,axes[0],[7,np.nanmax(x)],[4,9.5])
  
  x=gals['BulgeStellarMass']/gals['StellarMass']
  for i in range(0,6):
    axes=ax[i]
    plot_hist2d(x,y,axes[2],[0,1],[4,9.5])

  #Tracks
  #track_loc='/home/mmarshal/data_dragons/histories/tuned_reion/history_byBlackHoleMass/100/BlackHoleMass07-09/'
  track_loc='/home/mmarshal/data_dragons/histories/tuned_reion/history_byBlackHoleMass/158/BlackHoleMass08-10/'
  snapshots=np.arange(37,159)
  EXPANSION_FACTOR_PATH = "/fred/oz013/simulations/Tiamat/a_list.txt"
  a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
  z_list = np.array(1.0/a_list - 1.0)
  z_list = z_list[snapshots][78-37:]
 
  BH=np.log10(np.fromfile(track_loc+'BlackHoleMass.bin').reshape(-1,(121+1)))+10
  bulge=np.log10(np.fromfile(track_loc+'BulgeStellarMass.bin').reshape(-1,(121+1)))+10
  stellar=np.log10(np.fromfile(track_loc+'StellarMass.bin').reshape(-1,(121+1)))+10
  mvir=np.log10(np.fromfile(track_loc+'Mvir.bin').reshape(-1,(121+1)))+10
  colors=['red','orange','yellow','green','blue','purple']
  for i in range(0,np.shape(BH)[0]):
    prop=bulge[i][-1]/stellar[i][-1]
    ax[0,1].set_title('z = 2 B/T')
    cond=[1,0.8,0.6,0.4,0.2,0]
    
    #prop=mvir[i][100-37]
    #ax[0,1].set_title('z = 5 Virial Mass')
    #cond=[12.25,11.8,11.35,10.9,10.45,10]
    
    #prop=mvir[i][-1]
    #ax[0,1].set_title('z = 2 Virial Mass')
    #cond=[13,12.4,11.8,11.2,10.6,10]
    
    #prop=stellar[i][-1]
    #ax[0,1].set_title('z = 2 Stellar Mass')
    #cond=[11.50,11.25,11,10.75,10.5,10.25,10]
    
    #prop=stellar[i][100-37]
    #ax[0,1].set_title('z = 5 Stellar Mass')
    #cond=[10.35,10.2,10.05,9.9,9.75,9.4]
    
    #prop=BH[i][-1]
    #ax[0,1].set_title('z = 2 Black Hole Mass')
    #cond=[9,8.6,8.2,7.8,7.4,7]
    
    #prop=BH[i][100-37]
    #ax[0,1].set_title('z = 5 Black Hole Mass')
    #cond=[8,7.8,7.6,7.4,7.2,7]
    if prop>cond[0]:
      color=colors[0]
      axes=ax[0]
    elif prop>cond[1]:
      color=colors[1]
      axes=ax[1]
    elif prop>cond[2]:
      color=colors[2]
      axes=ax[2]
    elif prop>cond[3]:
      color=colors[3]
      axes=ax[3]
    elif prop>cond[4]:
      color=colors[4]
      axes=ax[4]
    elif prop>cond[5]:
      color=colors[5]
      axes=ax[5]
    plot_colorline(bulge[i][78-37:],BH[i][78-37:],z_list,axes[0],1,i)
    plot_colorline(stellar[i][78-37:],BH[i][78-37:],z_list,axes[1],0,i)
    plot_colorline(10**bulge[i][78-37:]/10**stellar[i][78-37:],BH[i][78-37:],z_list,axes[2],0,i)
    #axes[0].plot(bulge[i],BH[i],color=color,linewidth=0.2)

  for i in range(0,5):
    axes=ax[i]
    axes[0].set_xticks([])
    axes[1].set_xticks([])
    axes[2].set_xticks([])
    axes[0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 
    axes[1].text(7, 8.5, r'$M>{}$'.format(cond[i]),weight='bold',size='large')
  ax[5][1].text(7, 8.5, r'$M>{}$'.format(cond[i]),weight='bold',size='large')


  ax[5][2].set_xlabel(r'$\textrm{M}_{\textrm{bulge}}/\textrm{M}_\ast$')
  ax[5][1].set_xlabel(r'$\log(\textrm{M}_\ast)$')
  ax[5][0].set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}})$')
  ax[5][0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 
  #plt.legend()
  #lgd=plt.legend(fontsize='small',loc='upper center', bbox_to_anchor=(1.25, 0.8))  
  #plt.savefig('MBHBulge_z.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.tight_layout() 
  plt.show()
