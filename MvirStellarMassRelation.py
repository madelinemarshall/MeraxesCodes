import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import pandas as pd
import ContourPlot as cp
import pylab as p
import scipy.stats as stats

matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,5)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
  snapshot=snapshot,props=['StellarMass','BulgeStellarMass','BlackHoleMass','Mvir'],\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']*1e10>1e7]
  gals=gals[gals['BlackHoleMass']>0]
  gals=gals[gals['BulgeStellarMass']>0]
  #gals=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7] ##BULGES ONLY
  return gals


def plot_simulations(gals,axes):
  xlims=[10,12.9]
  ylims=[7,10.9]

  #2D Histogram
  #axes.plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
  cp_gals=gals[(gals['BulgeStellarMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10>1e4)]
  #cp.contour_plot(np.log10(cp_gals['BulgeStellarMass']*1e10),np.log10(cp_gals['BlackHoleMass']*1e10),xlims=xlims,ylims=ylims,axes=axes)
  H, xedges, yedges, img=axes.hist2d(np.log10(cp_gals['Mvir']*1e10), np.log10(cp_gals['StellarMass']*1e10), bins=20, range=[xlims,ylims], cmin=1, cmap='BuPu',vmax=25000,norm=matplotlib.colors.LogNorm())#,,vmin=0,vmax=25000)
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu',vmin=0,vmax=25000)
  #plt.colorbar(im,ax=axes,use_gridspec=True)
  axes.set_aspect('auto')


def plot_LOBF(gals,axes):
  xlims=[10,12.9]
  ylims=[7,10.9]

  #2D Histogram
  #axes.plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
  cp_gals=gals[(gals['BulgeStellarMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10>1e4)]
  
  ###Line of best fit
  #x=np.log10(cp_gals['BulgeStellarMass']*1e10)
  #y=np.log10(cp_gals['BlackHoleMass']*1e10)
  #axes.plot(xlims, np.poly1d(np.polyfit(x, y, 1),)(xlims),'k--')
  #print(np.poly1d(np.poly1d(np.polyfit(x, y, 1))))
 
  #slope,inter,rval,pval,err=stats.linregress(x,y)
  #axes.plot(xlims,inter+slope*np.array(xlims),'r--')

  #axes.set_xlim(xlims)
  #axes.set_ylim(ylims)
  binwidth=0.1
  med=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  mvir=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  hh=xlims[0]
  y=np.log10(cp_gals['Mvir']*1e10)
  for idx in range(0,len(med)):
    samp=cp_gals[(y>=hh)&(y<hh+binwidth)]['StellarMass']
    med[idx]=np.nanmean(samp)
    mvir[idx]=hh+binwidth/2
    hh=hh+binwidth
  mvir=mvir[np.invert(np.isnan(med))]
  med=med[np.invert(np.isnan(med))]

  #axes.plot(bulge,np.log10(med)+10,'k-')
  slope,inter,rval,pval,err=stats.linregress(mvir, np.log10(med)+10)
  #print("slope = {}, inter = {}".format(slope,inter))
  axes.plot(xlims,inter+slope*np.array(xlims),'k',linewidth=2.5,label='Best Fit', zorder=106)
  #axes.text(10, 6, r'$\beta={:.2f}$'.format(slope),weight='bold',size='large')
  return [slope,inter]


  
if __name__=="__main__":
  filename='tuned_reion_T125'
  #filename='bulges_nodiskinstability'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  snapshots=[52,63,78,100,116,134,158]
  fig, axes = plt.subplots(2, int((len(snapshots)+1)/2),gridspec_kw = {'wspace':0, 'hspace':0})

  slope={}#np.zeros(len(snapshots))
  inter={}#np.zeros(len(snapshots))
  ii=-1
  j=0
  for snapshot in snapshots:
    ii+=1
    if ii==(len(snapshots)+1)/2:
      j+=1
      ii=0
    gals=load_data(filename,snapshot)
    plot_simulations(gals,axes[j,ii])
    #slope[snapshot],inter[snapshot]=plot_LOBF(gals,axes[j,ii])
    if j==1:
      axes[j,ii].set_xlabel(r'$\log(M_{\mathrm{vir}}/M_\odot$)')
    else:
      axes[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(M_{\mathrm{\ast, total}}/M_\odot)$')
    else:
      axes[j,ii].set_yticklabels([])
    #axes[j,ii].set_xlim([7,12.1])
    axes[j,ii].text(11, 7.3, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='x-large')

  snapshot=100
  j=0
  ii=3
  lgd=axes[j,ii].legend(loc=(0.02,-0.95))
  #axes[j,ii].legend(loc=(0.02,-0.72),handles=['Meraxes - Best Fit','Marconi \& Hunt (2003)','Haring \& Rix (2004)',\
  #"Huang et al. (2018)","$M_{BH}/M_{Bulge}=10^{-2},10^{-3},10^{-4}$"])

  axes[1,int((len(snapshots)+1)/2)-1].axis('off')
  
  #plt.tight_layout()
  plt.show() 
  #fig.savefig('BHBulge.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  
  #fig,axes=plt.subplots(1,1)
  #plot_median_ratio(redshift,median_BH_bulge,eightyfourth_pctile_BH_bulge,sixteenth_pctle_BH_bulge,axes)
  
  ##plot_inter(inter,redshift)
  #plt.show()
