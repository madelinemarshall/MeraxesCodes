##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import ContourPlot as cp
import pandas as pd



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
                                          h=cosmo['h'],quiet=True)
  return gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)]



def load_mags(filename,snapshot):
  return pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')
  #return np.array(lumins[mag])



def plot_hist2d(xdata,ydata,axes,xlims,ylims):
  #xlims=[min(xdata),max(xdata)]
  #ylims=[min(ydata),max(ydata)]
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='BuPu',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu')
  cb=plt.colorbar(im,ax=axes,use_gridspec=True)
  cb.set_label('N, tot N = {}'.format(np.size(xdata)))
  axes.set_aspect('auto')
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)


if __name__=="__main__":
  filename='bulges_update1102_full'
  filename2='bulges_crotonSF'

  snapshot=78
  #snapshot=213

  ##Fig: Stellar Mass vs M1600
  default=load_gals('default',snapshot)
  gals1=load_gals(filename,snapshot)
  gals2=load_gals(filename2,snapshot)

  mag_def=load_mags('default',snapshot)['M1600']
  mag_1=load_mags(filename,snapshot)['M1600']
  mag_2=load_mags(filename2,snapshot)['M1600']

  fig, axes = plt.subplots(1, 3)
  xlims=[-24,-8]
  ylims=[6,11]

  #plot_hist2d(mag_def,np.log10(default['StellarMass']*1e10),axes[0],xlims,ylims)
  #plot_hist2d(mag_1,np.log10(gals1['StellarMass']*1e10),axes[1],xlims,ylims)
  #plot_hist2d(mag_2,np.log10(gals2['StellarMass']*1e10),axes[2],xlims,ylims)

  axes[0].set_title('Default Meraxes')
  axes[1].set_title('Bulge Model')
  axes[2].set_title('Bulge Model - Croton SF')
  axes[0].set_xlabel('$M_{1600}$')
  axes[1].set_xlabel('$M_{1600}$')
  axes[2].set_xlabel('$M_{1600}$')
  axes[0].set_ylabel('log(Stellar Mass)')
  axes[1].set_ylabel('log(Stellar Mass)')
  axes[2].set_ylabel('log(Stellar Mass)')
  plt.tight_layout()
  plt.show()

  ##Fig: B vs V
  lum_def=load_mags('default',snapshot)
  lum_1=load_mags(filename,snapshot)
  lum_2=load_mags(filename2,snapshot)

  
  #fig, axes = plt.subplots(1, 3)
  #axes[0].plot(B_def-V_def,V_def,'.')
  #xlims=[-10,10]
  #ylims=[-30,-8]
  #plot_hist2d(B_def-V_def,V_def,axes[0],xlims,ylims)
  #plt.show()

  ##Fig: M* - g-i
  fig,axes=plt.subplots(1,1)
  #plot_hist2d(np.log10(default['StellarMass']*1e10),lum_def['z850']-lum_def['J125'],axes,[6,11],[7,13])
  #plot_hist2d(lum_def['J125'],lum_def['z850']-lum_def['J125'],axes,[24,40],[0,1])
  #plot_hist2d(lum_def['z850']-lum_def['H160'],lum_def['V606']-lum_def['i775'],axes,[-0.5,1.5],[0.5,3.0])
  plot_hist2d(lum_def['i775']-lum_def['z850'],lum_def['Y105']-lum_def['H160'],axes,[-0.5,1.5],[0.5,3.0])
  plt.show()
  
