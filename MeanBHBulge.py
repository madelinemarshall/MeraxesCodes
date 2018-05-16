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
from scipy.optimize import curve_fit

matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (3.7,5.5)
#matplotlib.rcParams['figure.figsize'] = (7,3.6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def load_data(filename,meraxes_loc,snapshot):
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

  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,props=['StellarMass','BulgeStellarMass','BlackHoleMass'],\
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
  med=np.log10(np.median(gals['BlackHoleMass']/gals[stellarmass]))
  eightyfourth_pctile=np.log10(np.percentile(gals['BlackHoleMass']/gals[stellarmass],84))
  sixteenth_pctile=np.log10(np.percentile(gals['BlackHoleMass']/gals[stellarmass],16))
  return [med,eightyfourth_pctile,sixteenth_pctile]
  
def func(x,a,b):
  return b*np.log10(1+x)+np.log10(a)

def func_plot(x,a,b):
  return a *(1+x)**b


if __name__=="__main__":
  filename='tune_higherseed'
  filename2='tune_higherseed'
  meraxes_loc='/output/run1/meraxes_001.hdf5'
  #redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  redshift={63:7,78:6,100:5,116:4,134:3,158:2}
  #snapshots=[52,63,78,100,116,134,158]
  snapshots=[63,78,100,116,134,158]
  #redshift2={158:2,173:1.5,192:1,213:0.55}
  redshift2={158:2,213:0.55}
  snapshots2=[158,213]
  
  med={}
  pctile84={}
  pctile16={}
  med_b={}
  pctile84_b={}
  pctile16_b={}
  med_125_b={}
  pctile84_125_b={}
  pctile16_125_b={}
  med_1e5bulges={}
  med_125={}
  pctile84_125={}
  pctile16_125={}
  pctile84_1e5bulges={}
  pctile16_1e5bulges={}
  med_1e7bulges={}
  pctile84_1e7bulges={}
  pctile16_1e7bulges={}
  med_1e9bulges={}
  pctile84_1e9bulges={}
  pctile16_1e9bulges={}
  med_onlybulges={}
  pctile84_onlybulges={}
  pctile16_onlybulges={}
  med_onlydisks={}
  pctile84_onlydisks={}
  pctile16_onlydisks={}

  for snapshot in snapshots:
    gals=load_data(filename,meraxes_loc,snapshot)
    #gals_125=load_data(filename2,snapshot)
    #gals_125_b=gals_125[gals_125['BulgeStellarMass']*1e10>1e7]
    gals_1e5bulges=gals[(gals['BulgeStellarMass']*1e10>1e5)&(gals['BulgeStellarMass']*1e10<1e7)]
    gals_1e7bulges=gals[(gals['BulgeStellarMass']*1e10>1e7)&(gals['BulgeStellarMass']*1e10<1e9)]
    gals_1e9bulges=gals[(gals['BulgeStellarMass']*1e10>1e9)]
    gals_onlybulges=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
    gals_onlydisks=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]

    med_1e5bulges[snapshot],pctile84_1e5bulges[snapshot],pctile16_1e5bulges[snapshot]=find_med(gals_1e5bulges,'BulgeStellarMass')
    med_1e7bulges[snapshot],pctile84_1e7bulges[snapshot],pctile16_1e7bulges[snapshot]=find_med(gals_1e7bulges,'BulgeStellarMass')
    med_1e9bulges[snapshot],pctile84_1e9bulges[snapshot],pctile16_1e9bulges[snapshot]=find_med(gals_1e9bulges,'BulgeStellarMass')
    med_onlybulges[snapshot],pctile84_onlybulges[snapshot],pctile16_onlybulges[snapshot]=find_med(gals_onlybulges,'StellarMass')
    med_onlydisks[snapshot],pctile84_onlydisks[snapshot],pctile16_onlydisks[snapshot]=find_med(gals_onlydisks,'StellarMass')
  
    med[snapshot],pctile84[snapshot],pctile16[snapshot]=find_med(gals,'StellarMass')
    med_b[snapshot],pctile84_b[snapshot],pctile16_b[snapshot]=find_med(gals,'BulgeStellarMass')
    #med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125_b,'BulgeStellarMass')

  for snapshot in snapshots2:
    gals_125=load_data(filename2,meraxes_loc,snapshot)
    med_125[snapshot],pctile84_125[snapshot],pctile16_125[snapshot]=find_med(gals_125,'StellarMass')
    med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125,'BulgeStellarMass')

  #redshift_tot=redshift.copy()
  #redshift_tot.update(redshift2)


  fig,axes=plt.subplots(2,1,gridspec_kw = {'hspace':0})
  ##Bulge masses
  plot_median_ratio(redshift,med_1e5bulges,pctile84_1e5bulges,pctile16_1e5bulges,axes[1],color='turquoise',**{'label':r'$10^5M_\odot<M_{\mathrm{bulge}}<10^7M_\odot$'})
  plot_median_ratio(redshift,med_1e7bulges,pctile84_1e7bulges,pctile16_1e7bulges,axes[1],color='purple',**{'label':r'$10^7M_\odot<M_{\mathrm{bulge}}<10^9M_\odot$'})
  plot_median_ratio(redshift,med_1e9bulges,pctile84_1e9bulges,pctile16_1e9bulges,axes[1],color='r',**{'label':r'$M_{\mathrm{bulge}}>10^9M_\odot$'})

##Total stellar masses
  plot_median_ratio(redshift,med,pctile84,pctile16,axes[0],color='orchid',**{'label':r'$M_{\ast\textrm{, total}}$'})
  plot_median_ratio(redshift,med_b,pctile84_b,pctile16_b,axes[0],color='aquamarine',**{'label':r'$M_{\textrm{bulge}}$'})
  plot_median_ratio(redshift2,med_125,pctile84_125,pctile16_125,axes[0],color='orchid',**{'label':r'$M_{\ast\textrm{, total}}$'})
  plot_median_ratio(redshift2,med_125_b,pctile84_125_b,pctile16_125_b,axes[0],color='aquamarine',**{'label':r'$M_{\textrm{bulge}}$'})
  #lab='N={}'.format(np.size(gals_onlybulges))
  #plot_median_ratio(redshift,med_onlybulges,pctile84_onlybulges,pctile16_onlybulges,axes[1],color='aquamarine',**{'label':r'$B/T>0.7,$ '+lab})
  #lab='N={}'.format(np.size(gals_onlydisks))
  #plot_median_ratio(redshift,med_onlydisks,pctile84_onlydisks,pctile16_onlydisks,axes[1],color='orange',**{'label':r'$B/T<0.3,$ '+lab})
      
  #TSM
  zz=np.array(list(redshift.values()))
  sig_l=np.array(list(med.values()))-np.array(list(pctile16.values()))
  sig_u=np.array(list(med.values()))+np.array(list(pctile84.values()))
  sig=(sig_l+sig_u)/2
  popt,pcov = curve_fit(func,zz,np.array(list(med.values())))
  print("popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))  
  lab="$(1+z)^{"+f"{popt[1]:.2f}"+"}$"
  axes[0].plot(zz,np.log10(func_plot(zz,*popt)),'k:',label=lab)
  #axes[0].plot(zz,func(zz,*popt),'k:')

  #BSM
  sig_l=np.array(list(med_b.values()))-np.array(list(pctile16_b.values()))
  sig_u=-np.array(list(med_b.values()))+np.array(list(pctile84_b.values()))
  sig=(sig_l+sig_u)/2
  popt,pcov = curve_fit(func,zz,np.array(list(med_b.values())))
  print("popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))  
  lab="$(1+z)^{"+f"{popt[1]:.2f}"+"}$"
  axes[0].plot(zz,np.log10(func_plot(zz,*popt)),'k--',label=lab)

  #axes[0].plot(np.array(list(redshift.values())),np.log10(10**med[158]*(1+2)**(1.25)/(1+np.array(list(redshift.values())))**(1.25)),'k:',label=r'$(1+z)^{-1.25}$')
  #axes[0].plot(np.array(list(redshift.values())),np.log10(10**med_b[158]*(1+2)**(0.8)/(1+np.array(list(redshift.values())))**(0.8)),'k--',label=r'$(1+z)^{-0.8}$')

  #box = axes[0].get_position()
  lgd=axes[0].legend(bbox_to_anchor=(1.06,0),loc='lower right',fontsize='small',ncol=2)
  axes[1].legend(loc='lower right',fontsize='small')

  #axes[0].set_xlabel('Redshift')
  axes[0].set_ylabel(r'$\log(M_{\mathrm{BH}}/M_{\ast})$')
  axes[0].invert_xaxis()
  axes[0].set_xticklabels([])
  #axes[0].set_xlabel('Redshift')
  axes[1].set_xlabel('Redshift')
  axes[1].set_ylabel(r'$\log(M_{\mathrm{BH}}/M_{\mathrm{Bulge}})$')
  axes[1].invert_xaxis()
  axes[0].set_xlim([8,0.5])  
  axes[1].set_xlim([8,0.5])  
  axes[0].set_ylim([-4.79,-1.9])  
  axes[1].set_ylim([-5.49,-1.51])  

  from matplotlib.ticker import FormatStrFormatter
  axes[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

  plt.tight_layout()
  fig.savefig('MeanBHBulge.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
