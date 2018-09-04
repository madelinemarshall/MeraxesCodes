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
matplotlib.rcParams['figure.figsize'] = (4.4,5.5)#(3.7,5.5)
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
  gals=gals[gals['BlackHoleMass']*1e10>1e6]
  #gals=gals[gals['BulgeStellarMass']>0]
  #gals=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7] ##BULGES ONLY
  return gals


def plot_median_ratio(redshift,med,pctile84,pctile16,axes,color='b',med_0={250:1},**kwargs):
  axes.fill_between(np.log10(np.array(list(redshift.values()))+1),np.array(list(pctile84.values()))-med_0[250],\
    np.array(list(pctile16.values()))-med_0[250],alpha=0.15,color=color)
  axes.plot(np.log10(np.array(list(redshift.values()))+1),np.array(list(med.values()))-med_0[250],'-',linewidth=2.5,color=color,**kwargs)
  #axes.plot(np.array(list(redshift.values())),np.array(list(pctile84.values())),'--',**kwargs)
  #axes.plot(np.array(list(redshift.values())),np.array(list(pctile16.values())),'--',**kwargs)
  #plt.plot(np.array(list(redshift.values())),np.array(list(median_BH_bulge.values())),'k')
  #plt.plot(np.array(list(redshift.values())),np.array(list(eightyfourth_pctile_BH_bulge.values())),'k--')
  #plt.plot(np.array(list(redshift.values())),np.array(list(sixteenth_pctle_BH_bulge.values())),'k--')


def find_med(gals,stellarmass):
  gals=gals[gals[stellarmass]>0]
  if np.size(gals)==0:
    return [np.nan,np.nan,np.nan]
  else:
    pct=np.log10(np.percentile(gals['BlackHoleMass']/gals[stellarmass],[50,84,16]))
    #eightyfourth_pctile=np.log10(np.percentile(gals['BlackHoleMass']/gals[stellarmass],84))
    #sixteenth_pctile=np.log10(np.percentile(gals['BlackHoleMass']/gals[stellarmass],16))
    return pct#[med,eightyfourth_pctile,sixteenth_pctile]

  
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
  filename='tuned_reion'
  filename2='tuned_reion_T125'
  meraxes_loc='/output/meraxes.hdf5'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  snapshots=[52,63,78,100,116,134,158]
  #redshift={63:7,78:6,100:5,116:4,134:3,158:2}
  #snapshots=[63,78,100,116,134,158]
  #redshift2={158:2,173:1.5,192:1,213:0.55}
  redshift2={158:2,192:1,250:0}
  snapshots2=[158,192,250]
  bulge_split=1

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
  if bulge_split:
    med_1e5bulges={}
    pctile84_1e5bulges={}
    pctile16_1e5bulges={}
    med_1e6bulges={}
    pctile84_1e6bulges={}
    pctile16_1e6bulges={}
    med_1e7bulges={}
    pctile84_1e7bulges={}
    pctile16_1e7bulges={}
    med_1e8bulges={}
    pctile84_1e8bulges={}
    pctile16_1e8bulges={}
    med_1e9bulges={}
    pctile84_1e9bulges={}
    pctile16_1e9bulges={}
    pctile84_1e5bulges_125={}
    med_1e5bulges_125={}
    pctile16_1e5bulges_125={}
    pctile84_1e6bulges_125={}
    med_1e6bulges_125={}
    pctile16_1e6bulges_125={}
    med_1e7bulges_125={}
    pctile84_1e7bulges_125={}
    pctile16_1e7bulges_125={}
    med_1e8bulges_125={}
    pctile84_1e8bulges_125={}
    pctile16_1e8bulges_125={}
    med_1e9bulges_125={}
    pctile84_1e9bulges_125={}
    pctile16_1e9bulges_125={}
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
    if bulge_split:
      gals_1e5bulges=gals[(gals['BlackHoleMass']*1e10>1e6)]
      gals_1e6bulges=gals[(gals['BlackHoleMass']*1e10>10**6.5)]
      gals_1e7bulges=gals[(gals['BlackHoleMass']*1e10>1e7)]
      gals_1e8bulges=gals[(gals['BlackHoleMass']*1e10>10**7.5)]
      gals_1e9bulges=gals[(gals['BlackHoleMass']*1e10>1e8)]
    gals_onlybulges=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
    gals_onlydisks=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]

    if bulge_split:
      med_1e5bulges[snapshot],pctile84_1e5bulges[snapshot],pctile16_1e5bulges[snapshot]=find_med(gals_1e5bulges,'StellarMass')
      med_1e6bulges[snapshot],pctile84_1e6bulges[snapshot],pctile16_1e6bulges[snapshot]=find_med(gals_1e6bulges,'StellarMass')
      med_1e7bulges[snapshot],pctile84_1e7bulges[snapshot],pctile16_1e7bulges[snapshot]=find_med(gals_1e7bulges,'StellarMass')
      med_1e8bulges[snapshot],pctile84_1e8bulges[snapshot],pctile16_1e8bulges[snapshot]=find_med(gals_1e8bulges,'StellarMass')
      med_1e9bulges[snapshot],pctile84_1e9bulges[snapshot],pctile16_1e9bulges[snapshot]=find_med(gals_1e9bulges,'StellarMass')
    med_onlybulges[snapshot],pctile84_onlybulges[snapshot],pctile16_onlybulges[snapshot]=find_med(gals_onlybulges,'StellarMass')
    med_onlydisks[snapshot],pctile84_onlydisks[snapshot],pctile16_onlydisks[snapshot]=find_med(gals_onlydisks,'StellarMass')
  
    med[snapshot],pctile84[snapshot],pctile16[snapshot]=find_med(gals,'StellarMass')
    med_b[snapshot],pctile84_b[snapshot],pctile16_b[snapshot]=find_med(gals,'BulgeStellarMass')
    #med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125_b,'BulgeStellarMass')

  for snapshot in snapshots2:
    gals_125=load_data(filename2,meraxes_loc,snapshot)
    med_125[snapshot],pctile84_125[snapshot],pctile16_125[snapshot]=find_med(gals_125,'StellarMass')
    med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125,'BulgeStellarMass')
    if bulge_split:
      gals_1e5bulges=gals_125[(gals_125['BlackHoleMass']*1e10>1e6)]
      gals_1e6bulges=gals_125[(gals_125['BlackHoleMass']*1e10>10**6.5)]
      gals_1e7bulges=gals_125[(gals_125['BlackHoleMass']*1e10>1e7)]
      gals_1e8bulges=gals_125[(gals_125['BlackHoleMass']*1e10>10**7.5)]
      gals_1e9bulges=gals_125[(gals_125['BlackHoleMass']*1e10>1e8)]
      med_1e5bulges_125[snapshot],pctile84_1e5bulges_125[snapshot],pctile16_1e5bulges_125[snapshot]=find_med(gals_1e5bulges,'StellarMass')
      med_1e6bulges_125[snapshot],pctile84_1e6bulges_125[snapshot],pctile16_1e6bulges_125[snapshot]=find_med(gals_1e6bulges,'StellarMass')
      med_1e7bulges_125[snapshot],pctile84_1e7bulges_125[snapshot],pctile16_1e7bulges_125[snapshot]=find_med(gals_1e7bulges,'StellarMass')
      med_1e8bulges_125[snapshot],pctile84_1e8bulges_125[snapshot],pctile16_1e8bulges_125[snapshot]=find_med(gals_1e8bulges,'StellarMass')
      med_1e9bulges_125[snapshot],pctile84_1e9bulges_125[snapshot],pctile16_1e9bulges_125[snapshot]=find_med(gals_1e9bulges,'StellarMass')



  redshift_tot=redshift.copy()
  redshift_tot.update(redshift2)


  fig,axes=plt.subplots(2,1,gridspec_kw = {'hspace':0})
  ##Bulge masses
  if bulge_split:
    plot_median_ratio(redshift,med_1e5bulges,pctile84_1e5bulges,pctile16_1e5bulges,axes[1],color='turquoise',med_0=med_1e5bulges_125,**{'label':r'$>10^{6}$'})
    plot_median_ratio(redshift,med_1e6bulges,pctile84_1e6bulges,pctile16_1e6bulges,axes[1],color='green',med_0=med_1e6bulges_125,**{'label':r'$>10^{6.5}$'})
    plot_median_ratio(redshift,med_1e7bulges,pctile84_1e7bulges,pctile16_1e7bulges,axes[1],color='purple',med_0=med_1e7bulges_125,**{'label':r'$>10^{7}$'})
    plot_median_ratio(redshift,med_1e8bulges,pctile84_1e8bulges,pctile16_1e8bulges,axes[1],color='orange',med_0=med_1e8bulges_125,**{'label':r'$>10^{7.5}$'})
    plot_median_ratio(redshift,med_1e9bulges,pctile84_1e9bulges,pctile16_1e9bulges,axes[1],color='r',med_0=med_1e9bulges_125,**{'label':r'$>10^{8}$'})
    plot_median_ratio(redshift2,med_1e5bulges_125,pctile84_1e5bulges_125,pctile16_1e5bulges_125,axes[1],color='turquoise',med_0=med_1e5bulges_125,**{'label':'__nolabel__','linestyle':'--'})
    plot_median_ratio(redshift2,med_1e6bulges_125,pctile84_1e6bulges_125,pctile16_1e6bulges_125,axes[1],color='green',med_0=med_1e6bulges_125,**{'label':'__nolabel__','linestyle':'--'})
    plot_median_ratio(redshift2,med_1e7bulges_125,pctile84_1e7bulges_125,pctile16_1e7bulges_125,axes[1],color='purple',med_0=med_1e7bulges_125,**{'label':'__nolabel__','linestyle':'--'})
    plot_median_ratio(redshift2,med_1e8bulges_125,pctile84_1e8bulges_125,pctile16_1e8bulges_125,axes[1],color='orange',med_0=med_1e8bulges_125,**{'label':'__nolabel__','linestyle':'--'})
    plot_median_ratio(redshift2,med_1e9bulges_125,pctile84_1e9bulges_125,pctile16_1e9bulges_125,axes[1],color='r',med_0=med_1e9bulges_125,**{'label':'__nolabel__','linestyle':'--'})
    axes[1].set_xlabel('Redshift')
    axes[1].set_ylabel(r'$\Delta \log(M_{\mathrm{BH}})$')
    #axes[1].invert_xaxis()
    #axes[1].set_xlim([8,-0.04])  
    #axes[1].set_ylim([-3.9,-1.9])  
    from matplotlib.ticker import FormatStrFormatter
    axes[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axes[1].legend(loc='lower right',fontsize='small')
  else:
    axes[1].axis('off')
  axes[0].axis('off')

  #lab='N={}'.format(np.size(gals_onlybulges))
  #plot_median_ratio(redshift,med_onlybulges,pctile84_onlybulges,pctile16_onlybulges,axes[1],color='aquamarine',**{'label':r'$B/T>0.7,$ '+lab})
  #lab='N={}'.format(np.size(gals_onlydisks))
  #plot_median_ratio(redshift,med_onlydisks,pctile84_onlydisks,pctile16_onlydisks,axes[1],color='orange',**{'label':r'$B/T<0.3,$ '+lab})
    

  #TSM
  zz=np.array(list(redshift.values()))
  sig_l=np.array(list(med.values()))-np.array(list(pctile16.values()))
  sig_u=np.array(list(med.values()))+np.array(list(pctile84.values()))
  sig=(sig_l+sig_u)/2
  mv=np.array(list(med.values()))
  zz=zz[np.isfinite(mv)]
  mv=mv[np.isfinite(mv)]
  popt,pcov = curve_fit(func,zz,mv)
  #print("popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))  
  #lab="$(1+z)^{"+f"{popt[1]:.2f}"+"}$"

  #BSM
  sig_l=np.array(list(med_b.values()))-np.array(list(pctile16_b.values()))
  sig_u=-np.array(list(med_b.values()))+np.array(list(pctile84_b.values()))
  sig=(sig_l+sig_u)/2
  mv=np.array(list(med_b.values()))
  zz=np.array(list(redshift.values()))
  zz=zz[np.isfinite(mv)]
  mv=mv[np.isfinite(mv)]
  #popt,pcov = curve_fit(func,zz,mv)
  #print("popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))  
  #lab="$(1+z)^{"+f"{popt[1]:.2f}"+"}$"


  plt.tight_layout()
  fig.savefig('/home/mmarshal/results/plots/MeanBHBulge.pdf', format='pdf')
  plt.show()
