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
matplotlib.rcParams['figure.figsize'] = (3.5,5.5)#(3.7,5.5)
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
  filename2='paper2_T125'
  #snapshots=np.linspace(52,158,10) #30
  snapshots=np.linspace(37,158,10) #30
  snapshots2=np.linspace(158,250,5)
  #snapshots=np.linspace(52,158,8)
  #snapshots2=np.linspace(158,250,3)
  bulge_split=0 #plot BH bulge relation split into BH mass bins

  colors=['#e41a1c','#377eb8','#4daf4a']

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
    med_1e7bulges={}
    pctile84_1e7bulges={}
    pctile16_1e7bulges={}
    med_1e9bulges={}
    pctile84_1e9bulges={}
    pctile16_1e9bulges={}
    pctile84_1e5bulges_125={}
    med_1e5bulges_125={}
    pctile16_1e5bulges_125={}
    med_1e7bulges_125={}
    pctile84_1e7bulges_125={}
    pctile16_1e7bulges_125={}
    med_1e9bulges_125={}
    pctile84_1e9bulges_125={}
    pctile16_1e9bulges_125={}
  med_onlybulges={}
  pctile84_onlybulges={}
  pctile16_onlybulges={}
  med_onlydisks={}
  pctile84_onlydisks={}
  pctile16_onlydisks={}

  redshift={}
  for snapshot in snapshots:
    redshift[snapshot]=meraxes.io.grab_redshift('/home/mmarshal/data_dragons/'+filename+'/output/meraxes.hdf5', int(snapshot))
    gals=load_data(filename,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    gals=gals[(gals['BlackHoleMass']*1e10>1e6)]
    if bulge_split:
      gals_1e5bulges=gals[(gals['BlackHoleMass']*1e10>1e6)&(gals['BlackHoleMass']*1e10<1e7)]
      gals_1e7bulges=gals[(gals['BlackHoleMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10<1e8)]
      gals_1e9bulges=gals[(gals['BlackHoleMass']*1e10>1e8)]
    gals_onlybulges=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
    gals_onlydisks=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]

    if bulge_split:
      med_1e5bulges[snapshot],pctile84_1e5bulges[snapshot],pctile16_1e5bulges[snapshot]=find_med(gals_1e5bulges,'StellarMass')
      med_1e7bulges[snapshot],pctile84_1e7bulges[snapshot],pctile16_1e7bulges[snapshot]=find_med(gals_1e7bulges,'StellarMass')
      med_1e9bulges[snapshot],pctile84_1e9bulges[snapshot],pctile16_1e9bulges[snapshot]=find_med(gals_1e9bulges,'StellarMass')
    med_onlybulges[snapshot],pctile84_onlybulges[snapshot],pctile16_onlybulges[snapshot]=find_med(gals_onlybulges,'StellarMass')
    med_onlydisks[snapshot],pctile84_onlydisks[snapshot],pctile16_onlydisks[snapshot]=find_med(gals_onlydisks,'StellarMass')
  
    med[snapshot],pctile84[snapshot],pctile16[snapshot]=find_med(gals,'StellarMass')
    med_b[snapshot],pctile84_b[snapshot],pctile16_b[snapshot]=find_med(gals,'BulgeStellarMass')
   
  redshift2={}
  for snapshot in snapshots2:
    redshift2[snapshot]=meraxes.io.grab_redshift('/home/mmarshal/data_dragons/'+filename2+'/output/meraxes.hdf5', int(snapshot))
    gals_125=load_data(filename2,snapshot,['StellarMass','BulgeStellarMass','BlackHoleMass'])
    gals_125=gals_125[(gals_125['BlackHoleMass']*1e10>1e6)]
    med_125[snapshot],pctile84_125[snapshot],pctile16_125[snapshot]=find_med(gals_125,'StellarMass')
    med_125_b[snapshot],pctile84_125_b[snapshot],pctile16_125_b[snapshot]=find_med(gals_125,'BulgeStellarMass')
    if bulge_split:
      gals_1e5bulges=gals_125[(gals_125['BlackHoleMass']*1e10>1e6)&(gals_125['BlackHoleMass']*1e10<1e7)]
      gals_1e7bulges=gals_125[(gals_125['BlackHoleMass']*1e10>1e7)&(gals_125['BlackHoleMass']*1e10<1e8)]
      gals_1e9bulges=gals_125[(gals_125['BlackHoleMass']*1e10>1e8)]
      med_1e5bulges_125[snapshot],pctile84_1e5bulges_125[snapshot],pctile16_1e5bulges_125[snapshot]=find_med(gals_1e5bulges,'StellarMass')
      med_1e7bulges_125[snapshot],pctile84_1e7bulges_125[snapshot],pctile16_1e7bulges_125[snapshot]=find_med(gals_1e7bulges,'StellarMass')
      med_1e9bulges_125[snapshot],pctile84_1e9bulges_125[snapshot],pctile16_1e9bulges_125[snapshot]=find_med(gals_1e9bulges,'StellarMass')

  redshift_tot=redshift.copy()
  redshift_tot.update(redshift2)


  fig,axes=plt.subplots(2,1,gridspec_kw = {'hspace':0})
  ##Bulge masses
  if bulge_split:
    plot_median_ratio(redshift,med_1e5bulges,pctile84_1e5bulges,pctile16_1e5bulges,axes[1],color=colors[0],**{'label':r'$10^6M_\odot<M_{\mathrm{BH}}<10^7M_\odot$'})
    plot_median_ratio(redshift,med_1e7bulges,pctile84_1e7bulges,pctile16_1e7bulges,axes[1],color=colors[1],**{'label':r'$10^7M_\odot<M_{\mathrm{BH}}<10^8M_\odot$'})
    plot_median_ratio(redshift,med_1e9bulges,pctile84_1e9bulges,pctile16_1e9bulges,axes[1],color=colors[2],**{'label':r'$M_{\mathrm{BH}}>10^8M_\odot$'})
    plot_median_ratio(redshift2,med_1e5bulges_125,pctile84_1e5bulges_125,pctile16_1e5bulges_125,axes[1],color=colors[0],**{'label':'__nolabel__','linestyle':'--'})
    plot_median_ratio(redshift2,med_1e7bulges_125,pctile84_1e7bulges_125,pctile16_1e7bulges_125,axes[1],color=colors[1],**{'label':'__nolabel__','linestyle':'--'})
    plot_median_ratio(redshift2,med_1e9bulges_125,pctile84_1e9bulges_125,pctile16_1e9bulges_125,axes[1],color=colors[2],**{'label':'__nolabel__','linestyle':'--'})
    axes[1].set_xlabel('Redshift')
    axes[1].set_ylabel(r'$\log(M_{\mathrm{BH}}/M_\ast)$')
    axes[1].invert_xaxis()
    axes[1].set_xlim([8,-0.04])  
    axes[1].set_ylim([-3.9,-1.9])  
    from matplotlib.ticker import FormatStrFormatter
    axes[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axes[1].legend(loc='lower right',fontsize='small')
    axes[0].set_xticklabels([])
  else:
    axes[1].axis('off')
    axes[0].set_xlabel('Redshift')

##Total stellar masses
  plot_median_ratio(redshift,med,pctile84,pctile16,axes[0],color='orchid',**{'label':r'$M_{\ast\textrm{, total}}$'})
  plot_median_ratio(redshift,med_b,pctile84_b,pctile16_b,axes[0],color='aquamarine',**{'label':r'$M_{\ast\textrm{, bulge}}$'})
  plot_median_ratio(redshift2,med_125,pctile84_125,pctile16_125,axes[0],color='orchid',**{'label':'__nolabel__','linestyle':'--','lw':1})
  plot_median_ratio(redshift2,med_125_b,pctile84_125_b,pctile16_125_b,axes[0],color='aquamarine',**{'label':'__nolabel__','linestyle':'--','lw':1})
  #lab='N={}'.format(np.size(gals_onlybulges))
  #plot_median_ratio(redshift,med_onlybulges,pctile84_onlybulges,pctile16_onlybulges,axes[1],color='aquamarine',**{'label':r'$B/T>0.7,$ '+lab})
  #lab='N={}'.format(np.size(gals_onlydisks))
  #plot_median_ratio(redshift,med_onlydisks,pctile84_onlydisks,pctile16_onlydisks,axes[1],color='orange',**{'label':r'$B/T<0.3,$ '+lab})
    
  #axes[0].errorbar(0,-2.85,yerr=0.15,color='black',fmt='o',capsize=2,label=r'Kormendy \& Ho (2013)')
  plot_schulze(axes[0])

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
  #axes[0].plot(zz,np.log10(func_plot(zz,*popt)),'k:',label=lab)
  #axes[0].plot(zz,func(zz,*popt),'k:')

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
  #axes[0].plot(zz,np.log10(func_plot(zz,*popt)),'k--',label=lab)

  #axes[0].plot(np.array(list(redshift.values())),np.log10(10**med_b[158]*(1+2)**(0.8)/(1+np.array(list(redshift.values())))**(0.8)),'k--',label=r'$(1+z)^{-0.8}$')

  #box = axes[0].get_position()
  axes[0].plot(0,0,linestyle='-',lw=2.5,color='orchid',label='Tiamat')
  axes[0].plot(0,0,linestyle='--',lw=2.5,color='orchid',label='Tiamat-125-HR')
  if bulge_split==0:
    handles, labels = axes[0].get_legend_handles_labels()
    lgd=axes[0].legend(flip(handles, 2), flip(labels, 2),fontsize='small',ncol=2,loc='upper center', bbox_to_anchor=(0.46, -0.2))
  else:
    axes[0].plot(0,0,linestyle='-',lw=2.5,color=colors[0],label=r'$10^6M_\odot<M_{\mathrm{BH}}<10^7M_\odot$')
    axes[0].plot(0,0,linestyle='-',lw=2.5,color=colors[1],label=r'$10^7M_\odot<M_{\mathrm{BH}}<10^8M_\odot$')
    axes[0].plot(0,0,linestyle='-',lw=2.5,color=colors[2],label=r'$M_{\mathrm{BH}}>10^8M_\odot$')
    handles, labels = axes[0].get_legend_handles_labels()
    lgd=axes[1].legend(flip(handles, 2), flip(labels, 2),fontsize='small',ncol=2,loc='upper center', bbox_to_anchor=(0.46, -0.2))
    #axes[1].axvspan(8,6, alpha=0.2, color='gray')

  #axes[0].set_xlabel('Redshift')
  axes[0].set_ylabel(r'$\log(M_{\mathrm{BH}}/M_{\ast})$')
  axes[0].invert_xaxis()
  #axes[0].set_xlabel('Redshift')
  axes[0].set_xlim([8,-0.04])  
  axes[0].set_ylim([-4.1,-1.4])  
  if bulge_split:
    axes[1].set_xlim([8,-0.04])  
    axes[1].set_ylim([-4.1,-1.4])  
  axes[0].set_yticks(np.arange(-4.0,-1.49, 0.5))
  axes[1].set_yticks(np.arange(-4.0,-1.49, 0.5))

  #axes[0].axvspan(8,6, alpha=0.2, color='gray')



  plt.tight_layout()
  fig.savefig('/home/mmarshal/results/plots/Paper2/MeanBHBulge.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
