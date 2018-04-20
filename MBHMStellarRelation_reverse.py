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
#matplotlib.rcParams['figure.figsize'] = (5,4)
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
  gals=gals[gals['BlackHoleMass']*1e10>1e3]
  gals=gals[gals['StellarMass']*1e10>1e7]
  gals=gals[gals['BulgeStellarMass']*1e10>1e7]
  return gals


def plot_simulations(gals,axes,prop,xlims):
  prop_vals=np.log10(gals[prop]*1e10)
  ylims=[4,8]

  #2D Histogram
  H, xedges, yedges, img=axes.hist2d(prop_vals, np.log10(gals['BlackHoleMass']*1e10), bins=25, range=[xlims,ylims], cmin=1, cmap='BuPu',norm=matplotlib.colors.LogNorm(),cmax=2e3,vmax=2e3)#,,vmin=0,vmax=25000)
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu',norm=matplotlib.colors.LogNorm())  
  #im=axes.imshow(H,extent=extent,cmap='BuPu')
  if prop=='Mvir':
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.8, 0.11, 0.03, 0.77])
    cb=fig.colorbar(im, cax=cbar_ax,use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))
    #plt.colorbar(im,ax=axes,use_gridspec=True)
  axes.set_aspect('auto')



if __name__=='__main__':
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  #snapshots=[52,63,78,100,116,134,158]
  snap=78
  prop='StellarMass'
  #color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'black',213:'pink'}
  color={0:'C0',1:'C1',2:'C2',3:'C3',4:'C4',5:'aqua',6:'pink'}
  #filename='bulges_correctBHMF' 
  filename='bulges_noreion'

  fig,axes = plt.subplots(1,3,gridspec_kw = {'wspace':0, 'hspace':0},sharey=True)
  #fig,axes=plt.subplots(1,1)
  #fig2,axes2=plt.subplots(1,1)
  #fig3,axes3=plt.subplots(1,1)
  gals=load_data(filename,snap)
  MBH=gals['BlackHoleMass']*1e10
  Mbulge=gals['BulgeStellarMass']*1e10
  logMbulge=np.log10(Mbulge)
  logMBH=np.log10(MBH) 
  logMstellar=np.log10(gals['StellarMass']*1e10)
#    axes2.hist(np.log10(gals['BulgeStellarMass'][(gals['BlackHoleMass']*1e10>10**7)&(gals['BlackHoleMass']*1e10<10**7.5)]*1e10))
#    BT=gals['BulgeStellarMass']/gals['StellarMass']
#    axes3.hist(BT[(gals['BlackHoleMass']*1e10>10**7)&(gals['BlackHoleMass']*1e10<10**7.5)])

#      cp.contour_plot(logMbulge,logMBH,xlab=None,ylab=None,xlims=[8.15,11.6],ylims=[4,9.5],axes=axes,colors='C2',levels=np.logspace(-2,1,7),linewidth=0.9)
#      cp.contour_plot(logMstellar,logMBH,xlab=None,ylab=None,xlims=[8.15,11.6],ylims=[4,9.5],axes=axes,colors='C4',levels=np.logspace(-2,1,7),linewidth=0.9)
  plot_simulations(gals,axes[0],'BulgeStellarMass',[7,10.9])
  plot_simulations(gals,axes[1],'StellarMass',[7,10.9])
  plot_simulations(gals,axes[2],'Mvir',[9.5,12.5])
  
  axes[0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$')
  axes[1].set_xlabel(r'$\log(\textrm{M}_{\ast\textrm{, total}})$')
  axes[0].set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}})$')
  axes[2].set_xlabel(r'$\log(\textrm{M}_{\textrm{vir}})$')
  #axes.set_xlim([8.15,11.6])
  for ii in [0,1,2]:
    axes[ii].set_ylim([4,8])
  
  #plt.savefig('MBHBulge_z.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

