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
matplotlib.rcParams['figure.figsize'] = (3.5,5.5)
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
  snapshot=snapshot,props=['Mvir','StellarMass','BlackHoleMass'],\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['BlackHoleMass']*1e10>1e4]
  #gals=gals[(gals['StellarMass']*1e10>10**7.9)&(gals['StellarMass']*1e10<10**8)]
  #gals=gals[gals['StellarMass']*1e10>1e7]
  #gals=gals[(gals['Mvir']*1e10>10**9.9)&(gals['Mvir']*1e10<10**10)]
  return gals


def plot_median_ratio(redshift,med,axes,**kwargs):
  #axes.fill_between(np.array(list(redshift.values())),np.array(list(pctile84.values())),\
  #  np.array(list(pctile16.values())),alpha=0.15,color=color)
  axes.plot(np.array(list(redshift.values())),np.array(list(med.values())),'-',linewidth=2.5,**kwargs)


def find_med(gals,BH=True):
  med=np.log10(np.median(gals['BlackHoleMass']*1e10))
  #eightyfourth_pctile=np.log10(np.percentile(gals['BlackHoleMass']*1e10,84))
  #sixteenth_pctile=np.log10(np.percentile(gals['BlackHoleMass']*1e10,16))
  return med #[med,eightyfourth_pctile,sixteenth_pctile]



if __name__=="__main__":
  filename='bulges_correctBHMF'
  filename2='bulges_noreion'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  snapshots=[52,63,78,100,116,134,158]
  color={8:'C4',8.5:'pink',9:'C0',9.5:'C1',10:'C2',10.5:'C3',11:'C4',11.5:'pink',12:'black'}
  fig,axes=plt.subplots(1,1,gridspec_kw = {'hspace':0})
  gals_dict={}
  gals_def_dict={}
  
  for snapshot in snapshots:
    gals_dict[snapshot]=load_data(filename,snapshot)
    gals_def_dict[snapshot]=load_data(filename2,snapshot)

  for up in [8,8.5,9,9.5,10]:
    med={}
    #pctile84={}
    #pctile16={}
    med_def={}
    #pctile84_def={}
    #pctile16_def={}

    for snapshot in snapshots:
      gals=gals_dict[snapshot]
      gals_def=gals_def_dict[snapshot]

      med[snapshot]=find_med(gals[(gals['StellarMass']*1e10>10**(up-0.1))&(gals['StellarMass']*1e10<10**up)],True)
      med_def[snapshot]=find_med(gals_def[(gals_def['StellarMass']*1e10>10**(up-0.1))&(gals_def['StellarMass']*1e10<10**up)],True)

    plot_median_ratio(redshift,med,axes,**{'label':r'$\log(M_v)={}$'.format(up),'color':color[up]})
    plot_median_ratio(redshift,med_def,axes,**{'label':'__nolegend__','linestyle':'--','color':color[up]})
    index=1
    axes.plot(np.array(list(redshift.values())),np.log10(((1+np.array(list(redshift.values()),float))/(1+2))**(index)*10**med[158]),'k:',label='__nolegend__')


  lgd=axes.legend(loc='lower right',fontsize='small')
  axes.set_ylabel(r'$\log(M_\ast)$')
  axes.invert_xaxis()
  axes.set_xlabel('Redshift')
  axes.set_xlim([8,2])  

  plt.tight_layout()
  #fig.savefig('MeanBHBulge.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
