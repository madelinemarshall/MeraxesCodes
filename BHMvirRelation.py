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
  snapshot=snapshot,props=['StellarMass','BulgeStellarMass','BlackHoleMass','Mvir','Type'],\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']*1e10>1e7]
  gals=gals[gals['BlackHoleMass']>0]
  gals=gals[gals['BulgeStellarMass']>0]
  gals=gals[gals['Type']==0] ####ONLY CONSIDERS CENTRALS
  return gals


def plot_simulations(gals,axes):
  xlims=[10,13]
  ylims=[4,9.9]

  #2D Histogram
  #axes.plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
  cp_gals=gals[(gals['Mvir']*1e10>1e8)&(gals['BlackHoleMass']*1e10>1e4)]
  #cp.contour_plot(np.log10(cp_gals['BulgeStellarMass']*1e10),np.log10(cp_gals['BlackHoleMass']*1e10),xlims=xlims,ylims=ylims,axes=axes)
  H, xedges, yedges, img=axes.hist2d(np.log10(cp_gals['Mvir']*1e10), np.log10(cp_gals['BlackHoleMass']*1e10), bins=20, range=[xlims,ylims], cmin=1, cmap='BuPu',vmax=25000,norm=matplotlib.colors.LogNorm())#,,vmin=0,vmax=25000)
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu',vmin=0,vmax=25000)
  #plt.colorbar(im,ax=axes,use_gridspec=True)
  axes.set_aspect('auto')


def plot_LOBF(gals,axes):
  xlims=[10,13]
  ylims=[4,9.9]

  #2D Histogram
  #axes.plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
  cp_gals=gals[(gals['Mvir']*1e10>1e8)&(gals['BlackHoleMass']*1e10>1e4)]
  
  binwidth=0.1
  med=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  tot=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  hh=xlims[0]
  mvir=np.log10(cp_gals['Mvir']*1e10)
  for idx in range(0,len(med)):
    samp=cp_gals[(mvir>=hh)&(mvir<hh+binwidth)]['BlackHoleMass']
    med[idx]=np.nanmean(samp)
    tot[idx]=hh+binwidth/2
    hh=hh+binwidth
  tot=tot[np.invert(np.isnan(med))]
  med=med[np.invert(np.isnan(med))]
  slope,inter,rval,pval,err=stats.linregress(tot, np.log10(med)+10)
  axes.plot(xlims,inter+slope*np.array(xlims),'k',linewidth=2.5,label='Best Fit', zorder=106)


def plot_observations(snapshot,axes):
  xlims=[10,13]
  #Multu-Pakdil+17:
  #log(MBH /M⊙) = (1.55 ± 0.02) ∗ log(MDM /M⊙) − (11.26 ± 0.20) (R2 = 0.88, RMSE= 0.25)
  logMBH=1.55*np.array([11,13])-11.26
  axes.errorbar([11,13],logMBH,yerr=0.25,linestyle='-',color='darkorange',label='Mutlu-Pakdil et al. (2017)',capsize=3,linewidth=2.5, zorder=100)

  #Bogdan+17:
  #log(MBH/10^9)=-0.473+1.17log(M200/10^13)
  logMBH=-0.24+1.17*(np.array([12.5,13])-13)+9

  axes.errorbar([12.5,13],logMBH,yerr=0.32,linestyle='-',color='gold',label='Bogdan et al. (2017)',capsize=3,linewidth=2.5, zorder=101)


def plot_const_ratio(axes):
  xlims=np.array([10,13])
  axes.plot(xlims,-5+xlims,'k:',linewidth=2.5,label=r"$\log\frac{M_{\mathrm{BH}}}{M_{\mathrm{vir}}}=-5,-6,-7$", zorder=103)
  axes.plot(xlims,-6+xlims,'k:',linewidth=2.5,label="_nolegend_", zorder=104)
  axes.plot(xlims,-7+xlims,'k:',linewidth=2.5,label="_nolegend_", zorder=105)


  
if __name__=="__main__":
  filename='bulges_correctBHMF'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  snapshots=[52,63,78,100,116,134,158]
  fig, axes = plt.subplots(2, int((len(snapshots)+1)/2),gridspec_kw = {'wspace':0, 'hspace':0})

  ii=-1
  j=0
  for snapshot in snapshots:
    ii+=1
    if ii==(len(snapshots)+1)/2:
      j+=1
      ii=0
    gals=load_data(filename,snapshot)
    plot_simulations(gals,axes[j,ii])
    plot_LOBF(gals,axes[j,ii])
    plot_const_ratio(axes[j,ii])
    plot_observations(snapshot,axes[j,ii])

    if j==1:
      axes[j,ii].set_xlabel(r'$\log(M_{\mathrm{vir}}/M_\odot$)')
    else:
      axes[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(M_{\mathrm{BH}}/M_\odot)$')
    else:
      axes[j,ii].set_yticklabels([])
    axes[j,ii].set_xlim([10,12.9])
    axes[j,ii].text(10.9, 9.25, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='x-large')
    if snapshot==100:
      hand=(axes[j,ii].get_legend_handles_labels()[0])
      lab=(axes[j,ii].get_legend_handles_labels()[1])
      #hand.append(W[0])
      #lab.append("Wang et al. (2013)")
      #hand.append(HU[0])
      #ab.append("Huang et al. (2018)")
      lgd=axes[j,ii].legend(hand,lab,loc=(0.02,-0.72),fontsize='small')
      #axes[j,ii].legend(loc=(0.02,-0.72),handles=['Meraxes - Best Fit','Marconi \& Hunt (2003)','Haring \& Rix (2004)',\
      #"Huang et al. (2018)","$M_{BH}/M_{Bulge}=10^{-2},10^{-3},10^{-4}$"])
  axes[1,int((len(snapshots)+1)/2)-1].axis('off')
  
  fig.savefig('BHMvir.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show() 
  
