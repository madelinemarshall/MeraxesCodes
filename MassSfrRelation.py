##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import ContourPlot as cp


def load_data(filename):
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

  snapshot=116
  #snapshot=213

  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,\
                                          h=cosmo['h'],quiet=True)
  return(gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)&(gals["Sfr"]>0)])



def plot_hist2d(gals,axes):
  xlims=[6,11]
  ylims=[-8,3]

  H, xedges, yedges, img=axes.hist2d(np.log10(gals['StellarMass']*1e10), np.log10(gals['Sfr']), bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='BuPu',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu')
  cb=plt.colorbar(im,ax=axes,use_gridspec=True)
  cb.set_label('N, tot N = {}'.format(np.shape(gals)))
  axes.set_aspect('auto')
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  axes.set_xlabel('log(Stellar Mass)')
  axes.set_ylabel('log(Sfr)')



if __name__=="__main__":
  filename='bulges_update1102_full'
  filename2='bulges_crotonSF'

  default=load_data('default')
  gals1=load_data('bulges_update1102_full')
  gals2=load_data('bulges_crotonSF')

  ##Fig 1
  fig, axes = plt.subplots(1, 3)
  plot_hist2d(default,axes[0])
  plot_hist2d(gals1,axes[1])
  plot_hist2d(gals2,axes[2])
  axes[0].set_title('Default Meraxes')
  axes[1].set_title('Bulge Model')
  axes[2].set_title('Bulge Model - Croton SF')
  plt.suptitle("All Galaxies") 
  plt.tight_layout(rect=[0, 0.03, 1, 0.95])
  plt.show()
  
  ###Fig 2
  #bulges=gals1[gals1['BulgeStellarMass']/gals1['StellarMass']>0.5]
  #disks=gals1[gals1['BulgeStellarMass']/gals1['StellarMass']<0.5]
  #plt.scatter(np.log10(bulges['StellarMass']*1e10),np.log10(bulges['Sfr']),c=bulges['BulgeStellarMass']/bulges['StellarMass'],vmax=1,vmin=0)
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('B/T')
  #plt.title('B/T>0.5') 
  #plt.show()

  ###Fig 3
  #plt.scatter(np.log10(disks['StellarMass']*1e10),np.log10(disks['Sfr']),c=disks['BulgeStellarMass']/disks['StellarMass'],vmax=1,vmin=0)
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('B/T')
  #plt.title('B/T<0.5')
  #plt.show()

  ###Fig 4
  #cold=gals1[gals1['ColdGas']>np.median(gals1['ColdGas'])]
  #hot=gals1[gals1['ColdGas']<np.median(gals1['ColdGas'])]
  #plt.scatter(np.log10(cold['StellarMass']*1e10),np.log10(cold['Sfr']),c=np.log10(cold['ColdGas']*1e10))
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('log(Cold Gas Mass)')
  #plt.title('Cold Gas Mass - galaxies with lots of cold gas')
  #plt.show()
  
  ###Fig 5
  #plt.scatter(np.log10(hot['StellarMass']*1e10),np.log10(hot['Sfr']),c=np.log10(hot['ColdGas']*1e10))
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('log(Cold Gas Mass)')
  #plt.title('Cold Gas Mass - galaxies with little cold gas')
  #plt.show()
 
  ##Fig 6
  central=gals1[gals1['Type']==0]
  sats=gals1[gals1['Type']!=0]

  fig,axes=plt.subplots(1,2)
  plot_hist2d(central,axes[0])
  plot_hist2d(sats,axes[1])
  axes[0].set_title('Centrals')
  axes[1].set_title('Satellites')
  plt.show()

