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
from _load_data import load_data

matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.5,3.5)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors=['#e41a1c','#377eb8','#4daf4a']


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
  axes.errorbar(z,np.array(intrinsic)-2.85,xerr=np.array(z_err),yerr=np.array(in_err),marker='o',capsize=2,linestyle='None',color='gray')
  axes.plot(z,np.array(obs)-2.85,'o',markerfacecolor='white',markeredgecolor='gray')

if __name__=="__main__":
  filename='paper2'
  filename2='paper2_T125'
  filenameMR='paper2_T125MR'
  filenameNR='draft2_T125NR'
  #snapshots=np.linspace(52,158,30)
  snapshots=np.arange(30,80,10)
  med={}
  pctile84={}
  pctile16={}
  med_b={}
  pctile84_b={}
  pctile16_b={}
  med_125_b={}
  pctile84_125_b={}
  pctile16_125_b={}
  med_125={}
  pctile84_125={}
  pctile16_125={}
  med_MR={}
  pctile84_MR={}
  pctile16_MR={}
  med_NR={}
  pctile84_NR={}
  pctile16_NR={}

  redshift={}
  for snapshot in snapshots:
    redshift[snapshot]=meraxes.io.grab_redshift('/home/mmarshal/data_dragons/'+filename+'/output/meraxes.hdf5', int(snapshot))
    gals=load_data(filename,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    gals=gals[(gals['BlackHoleMass']*1e10>1e6)]

    med[snapshot],pctile84[snapshot],pctile16[snapshot]=find_med(gals,'StellarMass')
 
    gals_125=load_data(filename2,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    gals_125=gals_125[(gals_125['BlackHoleMass']*1e10>1e6)]
    med_125[snapshot],pctile84_125[snapshot],pctile16_125[snapshot]=find_med(gals_125,'StellarMass')

    gals_MR=load_data(filenameMR,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    gals_MR=gals_MR[(gals_MR['BlackHoleMass']*1e10>1e6)]
    med_MR[snapshot],pctile84_MR[snapshot],pctile16_MR[snapshot]=find_med(gals_MR,'StellarMass')
     
    #gals_NR=load_data(filenameNR,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    #med_NR[snapshot],pctile84_NR[snapshot],pctile16_NR[snapshot]=find_med(gals_NR,'StellarMass')
    #med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125,'BulgeStellarMass')
  

  fig,axes=plt.subplots(1,1)
  axes.set_xlabel('Redshift')
  
  plot_median_ratio(redshift,med,pctile84,pctile16,axes,color=colors[0],**{'label':'Tiamat'})
  plot_median_ratio(redshift,med_125,pctile84_125,pctile16_125,axes,color=colors[1],**{'label':'Tiamat-125-HR','linestyle':'-','lw':1})
  plot_median_ratio(redshift,med_MR,pctile84_MR,pctile16_MR,axes,color=colors[2],**{'label':'Tiamat-125-MR','linestyle':'-','lw':1})

  lgd=axes.legend(fontsize='small',ncol=3,loc='upper center', bbox_to_anchor=(0.42, -0.2))

  axes.set_ylabel(r'$\log(M_{\mathrm{BH}}/M_{\ast})$')
  axes.invert_xaxis()
  axes.set_xlim([8,2])  
  axes.set_ylim([-4.1,-1.9])  
  axes.set_yticks(np.arange(-4.0,-1.99, 0.5))

  plt.tight_layout()
  fig.savefig('/home/mmarshal/results/plots/Paper2/MeanBHBulge_resolutionTest.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
