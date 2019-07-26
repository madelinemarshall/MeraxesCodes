##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import pandas as pd
import magcalc as mc
from _load_data import load_data
from scipy.optimize import curve_fit
#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.6,4)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4


def plot_hist2d(xdata,ydata,axes,xlims,ylims,cmax=None,cbar=False):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, vmin=1, vmax=cmax, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  if cbar:
    cb=plt.colorbar(img, cax=cbar,use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))


def plot_avg(xdata,ydata,axes,xlims,bin_width,**kwargs):
  min_bin=np.floor(xlims[0])
  max_bin=np.floor(xlims[1])
  n_bins=np.int((max_bin-min_bin)/bin_width)
  avg_r=np.zeros(n_bins)
  pct_r_16=np.zeros(n_bins)
  pct_r_84=np.zeros(n_bins)
  bin_centre=np.zeros(n_bins)

  bin_start=min_bin
  bin_end=min_bin+1
  bin_num=0
  while bin_num<n_bins:
    y=ydata[(xdata<bin_end)&(xdata>=bin_start)]
    if np.size(y)>=25:
      avg_r[bin_num]=np.median(y)
      pct_r_16[bin_num]=np.percentile(y,16)
      pct_r_84[bin_num]=np.percentile(y,84)
    else:
      avg_r[bin_num]=np.nan
      pct_r_16[bin_num]=np.nan
      pct_r_84[bin_num]=np.nan
    bin_centre[bin_num]=bin_start+bin_width/2
    bin_start+=bin_width
    bin_end+=bin_width
    bin_num+=1
  
  axes.errorbar(bin_centre,avg_r,yerr=np.array([avg_r-pct_r_16,pct_r_84-avg_r]),**kwargs)#,color='k',marker='s',markersize=4,label='M19 - Median')


if __name__=="__main__":
  filename='paper1'
  filename_default='dragonsIII'#'dragons10'
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  ylims=[-2.8,1]
  xlims=[0.5,1.7]
  
  zz=[5,6,7,8,9,10]

  fig, axes = plt.subplots(3,4,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1],'width_ratios':[4,4,4,0.4]})

  ii=-1
  j=0
  for snapshot in np.flip([37,43,52,63,78,100],0):
    ii+=1
    if ii==3:
      j+=1
      ii=0
    
    if (ii==2):
      cbar=axes[j,3]
    else:
      cbar=False
    
    gals=load_data(filename,snapshot,'All',centrals=False)
    gals_default=load_data(filename_default,snapshot,'All',centrals=False)
    selection=(gals['Type']==0)&(gals['StellarDiskScaleLength']>0)
    selection_default=(gals_default['Type']==0)&(gals_default['DiskScaleLength']>0)
    gals=gals[selection]  
    gals_default=gals_default[selection_default]
  
    plot_hist2d(np.log10(gals['Rvir']*1000),np.log10(gals['StellarDiskScaleLength']*1000),axes[j,ii],xlims,ylims,cbar=cbar,cmax=2e3)
    #axes[j,ii].plot(np.log10(gals['Rvir']*1000),np.log10(gals['StellarDiskScaleLength']*1000),'.')
    #plot_hist2d(np.log10(gals_default['Rvir']),np.log10(gals_default['DiskScaleLength']*1000*1.67835),axes[j,ii],masslims,ylims,cbar=cbar,cmax=2e3)
    #plot_avg(np.log10(gals_no_cut_default['StellarMass']*1e10),np.log10(gals_no_cut_default['DiskScaleLength']*1000*1.67835),axes[j,ii],masslims,bin_width=0.5,**{'color':[0.5,0.5,0.5],'marker':'o','markersize':4,'label':'L16 - Median'})


    if j==1:
      axes[j,ii].set_xlabel(r'$\log(R_{vir}/\mathrm{kpc})$')
    else: 
      axes[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(R_e/\mathrm{kpc})$')
    else:
      axes[j,ii].set_yticklabels([])
    axes[j,ii].text(1.4, -2.4, r'$z={}$'.format(redshift[snapshot]))
  axes[2,0].axis('off')
  axes[2,1].axis('off')
  axes[2,2].axis('off')
  axes[2,3].axis('off')

  fig.savefig('/home/mmarshal/results/plots/Paper1/HaloDiscSize.pdf',format='pdf')
  
  plt.show()
