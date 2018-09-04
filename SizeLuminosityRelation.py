##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import ContourPlot as cp
import pandas as pd
import magcalc as mc


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
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def plot_hist2d(xdata,ydata,axes,xlims,ylims):
  #xlims=[min(xdata),max(xdata)]
  #ylims=[min(ydata),max(ydata)]
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='BuPu',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu')
  #cb=plt.colorbar(im,ax=axes,use_gridspec=True)
  #cb.set_label('N, tot N = {}'.format(np.size(xdata)))
  axes.set_aspect('auto')
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)


if __name__=="__main__":
  filename='tuned_reion_T125'
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}

  fig, axes = plt.subplots(2, 3,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  j=0
  for snapshot in np.flip([37,43,52,63,78,100],0):
    ii+=1
    if ii==3:
      j+=1
      ii=0

    gals=load_gals(filename,snapshot)
    mag=load_mags(filename,snapshot)
    mag=mag[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
    gals=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
    xlims=[-24,-10.01]
    ylims=[-1.5,0.99]

    plot_hist2d(mag,np.log10(gals['StellarDiskScaleLength']*1000),axes[j,ii],xlims,ylims)
    #axes[j,ii].scatter(mag,np.log10(gals['StellarDiskScaleLength']*1000),c=gals['BulgeStellarMass']/gals['StellarMass'])
    #axes[j,ii].set_xlim(xlims)
    #axes[j,ii].set_ylim(ylims)
    #plt.plot(mag,np.log10(gals['StellarDiskScaleLength']*1000))
    if j==1:
      axes[j,ii].set_xlabel('$M_{UV}$')
    else: 
      axes[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(R_e/\mathrm{kpc})$')
    else:
      axes[j,ii].set_yticklabels([])
    #axes[j,ii].set_title('$z=${}'.format(redshift[snapshot]))
    axes[j,ii].text(-23, 0.6, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='x-large')
  #axes[j,ii].colorbar=plt.colorbar
  #cb=plt.colorbar(ax=axes)
  plt.show()
