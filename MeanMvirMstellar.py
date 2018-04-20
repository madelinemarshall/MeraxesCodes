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
matplotlib.rcParams['figure.figsize'] = (7.2,3.4)
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
  #gals=gals[gals['BulgeStellarMass']>0]
  #gals=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7] ##BULGES ONLY
  return gals


def plot_median_ratio(redshift,med,pctile84,pctile16,axes,color='b',**kwargs):
  axes.fill_between(np.array(list(redshift.values())),np.array(list(pctile84.values())),\
    np.array(list(pctile16.values())),alpha=0.15,color=color)
  axes.plot(np.array(list(redshift.values())),np.array(list(med.values())),'-',linewidth=2.5,color=color,**kwargs)
  #axes.plot(np.array(list(redshift.values())),np.array(list(pctile84.values())),'--',**kwargs)
  #axes.plot(np.array(list(redshift.values())),np.array(list(pctile16.values())),'--',**kwargs)
  #plt.plot(np.array(list(redshift.values())),np.array(list(median_BH_bulge.values())),'k')
  #plt.plot(np.array(list(redshift.values())),np.array(list(eightyfourth_pctile_BH_bulge.values())),'k--')
  #plt.plot(np.array(list(redshift.values())),np.array(list(sixteenth_pctle_BH_bulge.values())),'k--')


def find_med(gals,stellarmass):
  gals=gals[gals[stellarmass]>0]
  med=np.log10(np.median(gals[stellarmass]/gals['Mvir']))
  eightyfourth_pctile=np.log10(np.percentile(gals[stellarmass]/gals['Mvir'],84))
  sixteenth_pctile=np.log10(np.percentile(gals[stellarmass]/gals['Mvir'],16))
  return [med,eightyfourth_pctile,sixteenth_pctile]
  


if __name__=="__main__":
  filename='bulges_correctBHMF'
  filename2='bulges_correctBHMF_tiamat125'
  #filename='bulges_nodiskinstability'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  snapshots=[52,63,78,100,116,134,158]
  redshift2={173:1.5,192:1,213:0.55}
  snapshots2=[173,192,213]
  
  med={}
  pctile84={}
  pctile16={}
  med_b={}
  pctile84_b={}
  pctile16_b={}
  med_125_b={}
  pctile84_125_b={}
  pctile16_125_b={}

  for snapshot in snapshots:
    gals=load_data(filename,snapshot)
    gals_125=load_data(filename2,snapshot)
    gals_125_b=gals_125[gals_125['BulgeStellarMass']*1e10>1e7]

  
    med[snapshot],pctile84[snapshot],pctile16[snapshot]=find_med(gals,'StellarMass')
    med_b[snapshot],pctile84_b[snapshot],pctile16_b[snapshot]=find_med(gals,'BulgeStellarMass')
    med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125_b,'BulgeStellarMass')

  for snapshot in snapshots2:
    gals_125=load_data(filename2,snapshot)
    med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125_b,'BulgeStellarMass')

  redshift_tot=redshift.copy()
  redshift_tot.update(redshift2)


  fig,axes=plt.subplots(1,1)

  plot_median_ratio(redshift,med,pctile84,pctile16,axes,color='orchid',**{'label':r'Total Stellar Mass'})
  plot_median_ratio(redshift,med_b,pctile84_b,pctile16_b,axes,color='aquamarine',**{'label':r'Bulge Stellar Mass'})

  axes.plot(np.array(list(redshift.values())),np.log10(10**med[158]*(1+2)**(0.3)/(1+np.array(list(redshift.values())))**(0.3)),'k:',label=r'$(1+z)^{-0.3}$')
  axes.plot(np.array(list(redshift.values())),np.log10(10**med_b[158]*(1+2)**(0.6)/(1+np.array(list(redshift.values())))**(0.6)),'k--',label=r'$(1+z)^{-0.6}$')

  #box = axes[0].get_position()
  axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15))

  axes.set_xlabel('Redshift')
  axes.set_ylabel(r'Median $M_{\mathrm{\ast}}/M_{\mathrm{vir}}$')
  axes.invert_xaxis()
  axes.set_xlim([8,2])  
  axes.set_ylim([-4.5,-1])  

  plt.tight_layout()
  plt.show()
  #fig.savefig('MeanBHBulge.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
