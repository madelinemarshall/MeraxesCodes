import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import pandas as pd
sys.path.append('/home/mmarshal/simulation_codes')
import ContourPlot as cp
import pylab as p
import scipy.stats as stats
from scipy.optimize import curve_fit
import itertools
from _load_data import load_data

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.2,3.2)#(3.7,5.5)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def plot_median_ratio(redshift,med,pctile84,pctile16,axes,color='b',**kwargs):
  axes.fill_between(np.array(list(redshift.values())),np.array(list(pctile84.values())),\
    np.array(list(pctile16.values())),alpha=0.15,color=color)
  axes.plot(np.array(list(redshift.values())),np.array(list(med.values())),'-',linewidth=2.5,color=color,**kwargs)


def find_med(gals,stellarmass):
  gals=gals[gals[stellarmass]>0]
  if np.size(gals)==0:
    return [np.nan,np.nan,np.nan]
  else:
    pct=np.log10(np.percentile(gals['BlackHoleMass']/gals[stellarmass],[50,84,16]))
    return pct
  

def func(x,a,b):
  return b**np.log10(1+x)+np.log10(a)


def func_plot(x,a,b):
  return a *(1+x)**b


def plot_schulze(axes):
  obs=[0.025,0.34,0.52,1.29]
  intrinsic=[0.025,0.15,0.15,-0.02]
  z=[0.56,1.5,4.2,6]
  z_err_1=[-0.175,0.175]
  z_err_2=[-0.5,0.5]
  z_err_3=[-0.1,0.1]
  z_err_4=[-0.25,0.25]
  z_err=[[0.175,0.5,0.1,0.25],[0.175,0.5,0.1,0.25]]
  in_err_1=[-0.09,0.09]
  in_err_2=[-0.08,0.27]
  in_err_3=[-1.49,0.98]
  in_err_4=[-2.76,1.36]
  in_err=[[0.09,0.08,1.49,2.76],[0.09,0.27,0.98,1.36]]
  axes.errorbar(z,np.array(intrinsic)-2.85,xerr=np.array(z_err),yerr=np.array(in_err),marker='o',capsize=2,linestyle='None',color='gray',label='Schulze \& Wisotzki (2014)')
  axes.plot(z,np.array(obs)-2.85,'o',markerfacecolor='white',markeredgecolor='gray')


if __name__=="__main__":
  filename='paper2'
  filename_HE='paper2_high_eta'
  #snapshots=np.linspace(52,158,10) #30
  snapshots=np.linspace(37,158,10) #30
  #snapshots=np.linspace(52,158,8)
  #snapshots2=np.linspace(158,250,3)
  bulge_split=1 #plot BH bulge relation split into BH mass bins

  colors=['#e41a1c','#377eb8','#4daf4a']

  med={}
  pctile84={}
  pctile16={}
  med_HE={}
  pctile84_HE={}
  pctile16_HE={}

  redshift={}
  for snapshot in snapshots:
    redshift[snapshot]=meraxes.io.grab_redshift('/home/mmarshal/data_dragons/'+filename+'/output/meraxes.hdf5', int(snapshot))
    #paper 2 model
    gals=load_data(filename,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    gals=gals[(gals['BlackHoleMass']*1e10>1e6)]
    med[snapshot],pctile84[snapshot],pctile16[snapshot]=find_med(gals,'StellarMass')

    #High eta comparison
    gals_HE=load_data(filename_HE,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    gals_HE=gals_HE[(gals_HE['BlackHoleMass']*1e10>1e6)]
    med_HE[snapshot],pctile84_HE[snapshot],pctile16_HE[snapshot]=find_med(gals_HE,'StellarMass')

  fig,axes=plt.subplots(1,1,gridspec_kw = {'hspace':0})
  ##Bulge masses
  axes.set_xlabel('Redshift')

##Total stellar masses
  plot_median_ratio(redshift,med,pctile84,pctile16,axes,color='orchid',**{'label':r'$\eta=0.06$'})
  plot_median_ratio(redshift,med_HE,pctile84_HE,pctile16_HE,axes,color='cornflowerblue',**{'label':r'$\eta=0.2$','linestyle':'--'})

  lgd=axes.legend()
  
  axes.set_ylabel(r'$\log(M_{\mathrm{BH}}/M_{\ast})$')
  axes.invert_xaxis()
  axes.set_xlim([8,2])  
  axes.set_ylim([-3.8,-2.5])  

  plt.tight_layout()
  fig.savefig('/home/mmarshal/results/plots/Paper2/MeanBHBulge_etaComparison.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
