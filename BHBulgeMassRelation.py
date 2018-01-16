import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import pandas as pd
import ContourPlot as cp
import pylab as p


def load_data(filename,snapshot):
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

  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,props=['StellarMass','BulgeStellarMass','BlackHoleMass'],\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']*1e10>1e7]
  gals=gals[gals['BlackHoleMass']>0]
  gals=gals[gals['BulgeStellarMass']>0]
  return gals


def plot_simulations(gals,axes):
  xlims=[7,12]
  ylims=[4,12]

  #2D Histogram
  #axes.plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
  cp_gals=gals[(gals['BulgeStellarMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10>1e4)]
  #cp.contour_plot(np.log10(cp_gals['BulgeStellarMass']*1e10),np.log10(cp_gals['BlackHoleMass']*1e10),xlims=xlims,ylims=ylims,axes=axes)
  H, xedges, yedges, img=axes.hist2d(np.log10(cp_gals['BulgeStellarMass']*1e10), np.log10(cp_gals['BlackHoleMass']*1e10), bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='BuPu',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu')
  plt.colorbar(im,ax=axes,use_gridspec=True)
  axes.set_aspect('auto')
  
  ##Line of best fit
  #x=(cp_gals['BulgeStellarMass']*1e10)
  #y=(cp_gals['BlackHoleMass']*1e10)
  #axes.plot(xlims, np.poly1d(np.polyfit(np.log10(x), np.log10(y), 1))(xlims),'k--')
  axes.set_xlim(xlims)
  axes.set_ylim(ylims)
  axes.set_xlabel('log(Bulge Mass)')
  axes.set_ylabel('log(BH Mass)')


def plot_observations(snapshot,axes):
  #Haring & Rix
  MBH=10**(8.2+1.12*np.log10(np.array([10**7,10**12])/10**11))
  #MBH_low=10**(8.1+1.06*np.log10(np.array([10**8,10**12])/10**11))
  #MBH_high=10**(8.3+1.18*np.log10(np.array([10**8,10**12])/10**11))
  axes.plot([7,12],np.log10(MBH),color='orange')
  if (snapshot==78):
    Wang_BH=np.log10([2.8e9,2.1e9,8.6e8,1.7e8])
    Wang_dyn=np.log10([9.6e10,12.5e10,7.2e10,1.3e10])
    axes.plot(Wang_dyn,Wang_BH,'ko')
  axes.set_title(f'z = {redshift[snapshot]}')


def plot_median_ratio(redshift,median_BH_bulge,eightyfourth_pctile_BH_bulge,sixteenth_pctle_BH_bulge):
  axes[1,int((len(np.array(list(redshift.values())))+1)/2)-1].plot(np.array(list(redshift.values())),np.array(list(median_BH_bulge.values())),'k')
  axes[1,int((len(np.array(list(redshift.values())))+1)/2)-1].plot(np.array(list(redshift.values())),np.array(list(eightyfourth_pctile_BH_bulge.values())),'k--')
  axes[1,int((len(np.array(list(redshift.values())))+1)/2)-1].plot(np.array(list(redshift.values())),np.array(list(sixteenth_pctle_BH_bulge.values())),'k--')
  axes[1,int((len(np.array(list(redshift.values())))+1)/2)-1].set_xlabel('Redshift')
  axes[1,int((len(np.array(list(redshift.values())))+1)/2)-1].set_ylabel('Median BH/Bulge Mass Ratio')

  axes[1,int((len(np.array(list(redshift.values())))+1)/2)-1].invert_xaxis()

  #plt.plot(np.array(list(redshift.values())),np.array(list(median_BH_bulge.values())),'k')
  #plt.plot(np.array(list(redshift.values())),np.array(list(eightyfourth_pctile_BH_bulge.values())),'k--')
  #plt.plot(np.array(list(redshift.values())),np.array(list(sixteenth_pctle_BH_bulge.values())),'k--')
  
if __name__=="__main__":
  filename='bulges_update1102_full'
  #filename='bulges_nodiskinstability'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  median_BH_bulge={52:0,63:0,78:0,100:0,116:0,134:0,158:0}
  eightyfourth_pctile_BH_bulge={52:0,63:0,78:0,100:0,116:0,134:0,158:0}
  sixteenth_pctle_BH_bulge={52:0,63:0,78:0,100:0,116:0,134:0,158:0}
  snapshots=[52,63,78,100,116,134,158]
  fig, axes = plt.subplots(2, int((len(snapshots)+1)/2))
  ii=-1
  j=0
  for snapshot in snapshots:
    ii+=1
    if ii==(len(snapshots)+1)/2:
      j+=1
      ii=0
    gals=load_data(filename,snapshot)
    plot_simulations(gals,axes[j,ii])
    plot_observations(snapshot,axes[j,ii])
    median_BH_bulge[snapshot]=np.log10(np.median(gals['BlackHoleMass']/gals['BulgeStellarMass']))
    eightyfourth_pctile_BH_bulge[snapshot]=np.log10(np.percentile(gals['BlackHoleMass']/gals['BulgeStellarMass'],84))
    sixteenth_pctle_BH_bulge[snapshot]=np.log10(np.percentile(gals['BlackHoleMass']/gals['BulgeStellarMass'],16))
  plot_median_ratio(redshift,median_BH_bulge,eightyfourth_pctile_BH_bulge,sixteenth_pctle_BH_bulge)
  plt.tight_layout()
  plt.show()
