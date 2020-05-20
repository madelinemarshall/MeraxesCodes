import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _calculateQLF import calculateQLF
from KulkarniLF import plot_Kulkarni
import random
random.seed(2)

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.6,2.5)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
import itertools

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

from astrodatapy.number_density import number_density

colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D']*4
linestyles     = ['-','--','-.',':']*4


def plot_obsQLF_X(axes,z):
    L=np.logspace(41,47)
    erg_conv=1/(3.839*1e33)
    ##Ebrero+09
    if z<=3:
      Ebrero={'A':4.78*(cosmo['h']/0.7)**3,'L0':10**(43.91)*(cosmo['h']/0.7)**-2,'gamma1':0.96,'gamma2':2.35,'p1':4.07,'p2':-1.5,'zc':1.9,'La':10**(44.6)*(cosmo['h']/0.7)**-2,'alpha':0.245}
      axes.plot(np.log10(L*erg_conv),obsLF(L,z,Ebrero),label='Ebrero et al. (2009)',linestyle=':',color='#ffbf00',linewidth=2)

    # La Franca + (model 4) ??no h dependence of L0 and La mentioned in text??
    if z<=4:
      LaFranca={'A':1.21*(cosmo['h']/0.7)**3,'L0':10**(44.25)*(cosmo['h']/0.7)**-2,'gamma1':1.01,'gamma2':2.38,'p1':4.62,'p2':-1.15,'zc':2.29,'La':10**(45.74)*(cosmo['h']/0.7)**-2,'alpha':0.22}
      axes.plot(np.log10(L*erg_conv),obsLF(L,z,LaFranca),label='La Franca et al. (2005)',linestyle='--',color=colors[2],linewidth=2)

    
    Ueda={'A':5.04*(cosmo['h']/0.7)**3,'L0':10**(43.94)*(cosmo['h']/0.7)**-2,'gamma1':0.86,'gamma2':2.23,'p1':4.23,'p2':-1.5,'zc':1.9,'La':10**(44.6)*(cosmo['h']/0.7)**-2,'alpha':0.335}
    axes.plot(np.log10(L*erg_conv),obsLF(L,z,Ueda),label='Ueda et al. (2003)',linestyle='-.',color=colors[3],linewidth=2)
    
    #Yencho={'A':(1e6*10**(-6.140))*(cosmo['h']/0.7)**3,'L0':10**(44.40)*(cosmo['h']/0.7)**-2,'gamma1':0.872,'gamma2':2.36,'p1':3.61,'p2':-2.83,'zc':2.18,'La':10**(45.09)*(cosmo['h']/0.7)**-2,'alpha':0.208}
    #axes.plot(np.log10(L*erg_conv),obsLF(L,z,Yencho))
    
    #Silverman={'A':(1e6*10**(-6.163))*(cosmo['h']/0.7)**3,'L0':10**(44.33)*(cosmo['h']/0.7)**-2,'gamma1':2.15,'gamma2':1.1,'p1':4.22,'p2':-3.27,'zc':1.89,'La':10**(44.6)*(cosmo['h']/0.7)**-2,'alpha':0.333}
    #axes.plot(np.log10(L*erg_conv),obsLF(L,z,Silverman))


def obsLF(L,z,params):
    return np.log10(params['A']*1e-6*((L/params['L0'])**params['gamma1']+(L/params['L0'])**params['gamma2'])**-1*ez(L,z,params))


def ez(L,z,params):
    zc=zc_func(L,params)
    ez=np.ones_like(L)*(1+z)**params['p1']
    ez[z>=zc]=(1+zc[z>=zc])**params['p1'] * ((1+z)/(1+zc[z>=zc]))**params['p2']
    return ez

def zc_func(L,params):
    zc=np.ones_like(L)*params['zc']
    zc[L<params['La']]=params['zc']*(L[L<params['La']]/params['La'])**params['alpha']
    return zc

def make_legend(ax):
    #calculateQLF(gals_bulges,fname_1,'UV',ax,**{'linestyle':'-','label':'Bulge Model','linewidth':1.5,'color':'k','zorder':0})
    #calculateQLF(gals_default,fname_d,'UV',ax,**{'linestyle':'--','label':'Default Meraxes','linewidth':1.5,'color':'gray','zorder':1})
    #labels_dict={}
    #kk=0
    for snapshot in [78,100,116,134]:
      #kk=plot_obs(ax,redshift[snapshot],1,labels_dict,kk)
      plot_Kulkarni(redshift[snapshot]-0.2,redshift[snapshot]+0.2,ax,10)
    ax.plot([0,0.1],[0,0.1],**{'linestyle':'-','label':'M18 Meraxes','linewidth':1.5,'color':'k','zorder':0})
    ax.plot([0,0.1],[0,0.1],**{'linestyle':'--','label':'Q17 Meraxes','linewidth':1.5,'color':'gray','zorder':1})
    ax.plot([0,0.1],[0,0.1],**{'linestyle':'-','label':'Kulkarni et al. 2018','linewidth':1.5,'color':'#ffbf00','zorder':3})
    handles, labels = ax.get_legend_handles_labels()
    #ax.legend(flip(handles, 4), flip(labels, 4),fontsize='small',ncol=4,loc=(0,0))
    lgd=ax.legend(fontsize='small',loc=(0.02,-0.1))
    ax.axis('off')
    ax.set_xlim(8,8.1)
    ax.set_ylim(-6,-5.8)
    #ax.set_visible(False)
    return lgd


def load_data(full_filename,snapshot,props,centrals=False):
  cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  } 
  
  if (props=='All')|(props==None):
    gals=meraxes.io.read_gals(full_filename,\
    snapshot=snapshot,\
    h=cosmo['h'],quiet=True)
  else:
    if type(props)==str:
      props=[props]
    if 'GhostFlag' not in props:
      props.append('GhostFlag')
    if 'Type' not in props:
      props.append('Type')
    if 'StellarMass' not in props:
      props.append('StellarMass')
    #props.append('GhostFlag')
    gals=meraxes.io.read_gals(full_filename,\
    snapshot=snapshot,props=props,\
    h=cosmo['h'],quiet=True)
  gals=gals[gals['GhostFlag']==0]
  gals=gals[gals['StellarMass']*1e10>1e7]
  #if 'BlackHoleMass' in props:
  #  gals=gals[gals['BlackHoleMass']*1e10>1e5]
  if centrals==True:
    gals=gals[gals['Type']==0]
  return gals


if __name__=="__main__":
  #Setup
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
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  prop='BlackHoleMass'
  
  filename='tuning_Tiamat/'
  meraxes_loc=str(sys.argv[1])
  sim_props = meraxes.io.read_input_params(data_folder+filename+meraxes_loc,h=cosmo['h'],quiet=True)
  eta = sim_props['RadioAccretionEff']

  
##XRAY
  figX, axesX = plt.subplots(1, 5,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)
  ii=-1
  for snapshot in [100,116,134,158]:
    ii+=1
    gals_bulges=load_data(data_folder+filename+meraxes_loc,snapshot,[prop,'BlackHoleAccretedColdMass','GhostFlag','dt'])

    gals_bulges=gals_bulges[gals_bulges['BlackHoleMass']*1e10>1e6]

    calculateQLF(gals_bulges,data_folder+filename+meraxes_loc,'hardX',axesX[ii],eta=eta,**{'linestyle':'-','label':'M18 Model','linewidth':2,'color':'k','zorder':100})
    plot_obsQLF_X(axesX[ii],redshift[snapshot])
    axesX[ii].set_xlabel(r'$L_{\textrm{2-10keV}}/L_\odot$')

    if ii==0:
      axesX[ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$',usetex=True)
    #else:
    #  axesX[ii].set_yticklabels([])
    axesX[ii].set_xlim([7.5,13])
    axesX[ii].set_ylim([-8,-2])
    axesX[ii].text(10.8, -2.7, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
    #axes[j,ii].grid(color=[0.8,0.8,0.8],linestyle='--') 

  lgdX=axesX[-2].legend(loc=(1.05,0.2))
  axesX[-1].axis('off')
  plt.tight_layout()
  plt.savefig('/home/mmarshal/results/plots/tuning_paper2/Tiamat/XrayQLF_{}.pdf'.format(int(meraxes_loc[-8:-5])),format='pdf',bbox_extra_artists=(lgdX,), bbox_inches='tight')
  plt.show()

