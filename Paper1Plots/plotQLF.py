import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _calculateQLF import calculateQLF

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (6.5,2.5)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import itertools

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

from astrodatapy.number_density import number_density

colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D']*4
linestyles     = ['-','--','-.',':']*4


def load_data(filename,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'BlackHoleAccretedColdMass','GhostFlag','dt'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['BlackHoleMass']*1e10>1e6]
  return gals


def plot_obsQLF_UV(ax,z,make_legend,labels_dict,jj):
    feature = 'QLF_UV'
    xlim    = (-18, -29)
    ylim    = (-10, -3.5)

    obs    = number_density(feature=feature,z_target=z,quiet=1,h=cosmo['h'])
    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Qin2017_Tiamat") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        if (label[:label.rfind('20')]=='Gilkman')&(label[-4:]=='SDSS'):
            label = 'Gilkman et al. (2011)\nSDSS'
        elif (label[:label.rfind('20')]=='Gilkman'):
            label = 'Gilkman et al. (2011)\nNDWFS DLS'
        else:
            label      = label[:label.rfind('20')]+' et al. ('+label[label.rfind('20'):]+')'
        if label in labels_dict:
            if make_legend:
              continue
            datatype   = labels_dict[label][0]
            color   = labels_dict[label][1]
            marker   = labels_dict[label][2]
            linestyle   = labels_dict[label][3]
        else:
            datatype   = obs.target_observation['DataType'][ii]
            color      = colors[jj]#'gray'
            marker     = markers[jj]
            linestyle  = linestyles[k_func]
            labels_dict[label]=[datatype,color,marker,linestyle]
            jj += 1

        data[:,1:] = np.log10(data[:,1:])
        if datatype == 'data':
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        elif datatype == 'dataULimit':
            ax.errorbar(data[:,0],  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        else:
            ax.plot(data[:,0],data[:,1],label=label,color=color,linestyle=linestyle,lw=1)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1
    return jj


def plot_obsQLF_B(ax,z,make_legend,labels_dict,jj):
    feature = 'QLF_optical'
    xlim    = (-18, -29)
    ylim    = (-10, -3.5)

    obs    = number_density(feature=feature,z_target=z,quiet=1,h=cosmo['h'])
    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Qin2017_Tiamat") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        if (label[:label.rfind('20')]=='Gilkman')&(label[-4:]=='SDSS'):
            label = 'Gilkman et al. (2011)\nSDSS'
        elif (label[:label.rfind('20')]=='Gilkman'):
            label = 'Gilkman et al. (2011)\nNDWFS DLS'
        else:
            label      = label[:label.rfind('20')]+' et al. ('+label[label.rfind('20'):]+')'
        if label in labels_dict:
            if make_legend:
              continue
            datatype   = labels_dict[label][0]
            color   = labels_dict[label][1]
            marker   = labels_dict[label][2]
            linestyle   = labels_dict[label][3]
        else:
            datatype   = obs.target_observation['DataType'][ii]
            color      = colors[jj]#'gray'
            marker     = markers[jj]
            linestyle  = linestyles[k_func]
            labels_dict[label]=[datatype,color,marker,linestyle]
            jj += 1

        data[:,1:] = np.log10(data[:,1:])
        if datatype == 'data':
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        elif datatype == 'dataULimit':
            ax.errorbar(data[:,0],  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        else:
            ax.plot(data[:,0],data[:,1],label=label,color=color,linestyle=linestyle,lw=1)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1
    return jj


def plot_obsQLF_X(axes,z):
    L=np.logspace(41,47)
    erg_conv=1/(3.839*1e33)
    ##Ebrero+09
    if z<=3:
      Ebrero={'A':4.78*(cosmo['h']/0.7)**3,'L0':10**(43.91)*(cosmo['h']/0.7)**-2,'gamma1':0.96,'gamma2':2.35,'p1':4.07,'p2':-1.5,'zc':1.9,'La':10**(44.6)*(cosmo['h']/0.7)**-2,'alpha':0.245}
      axes.plot(np.log10(L*erg_conv),obsLF(L,z,Ebrero),label='Ebrero et al. (2009)',linestyle=':',color=colors[0],linewidth=1.5)

    # La Franca + (model 4) ??no h dependence of L0 and La mentioned in text??
    if z<=4:
      LaFranca={'A':1.21*(cosmo['h']/0.7)**3,'L0':10**(44.25)*(cosmo['h']/0.7)**-2,'gamma1':1.01,'gamma2':2.38,'p1':4.62,'p2':-1.15,'zc':2.29,'La':10**(45.74)*(cosmo['h']/0.7)**-2,'alpha':0.22}
      axes.plot(np.log10(L*erg_conv),obsLF(L,z,LaFranca),label='La Franca et al. (2005)',linestyle='--',color=colors[2],linewidth=1.5)

    
    Ueda={'A':5.04*(cosmo['h']/0.7)**3,'L0':10**(43.94)*(cosmo['h']/0.7)**-2,'gamma1':0.86,'gamma2':2.23,'p1':4.23,'p2':-1.5,'zc':1.9,'La':10**(44.6)*(cosmo['h']/0.7)**-2,'alpha':0.335}
    axes.plot(np.log10(L*erg_conv),obsLF(L,z,Ueda),label='Ueda et al. (2003)',linestyle='-.',color=colors[3],linewidth=1.5)
    
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
    labels_dict={}
    kk=0
    for snapshot in [78,100,116,134]:
      kk=plot_obsQLF_UV(ax,redshift[snapshot],1,labels_dict,kk)
    ax.plot([0,0.1],[0,0.1],**{'linestyle':'-','label':'M18 Meraxes','linewidth':1.5,'color':'k','zorder':0})
    ax.plot([0,0.1],[0,0.1],**{'linestyle':'--','label':'Q17 Meraxes','linewidth':1.5,'color':'gray','zorder':1})
    handles, labels = ax.get_legend_handles_labels()
    #ax.legend(flip(handles, 4), flip(labels, 4),fontsize='small',ncol=4,loc=(0,0))
    lgd=ax.legend(fontsize='small',loc=(0.02,-0.1))
    ax.axis('off')
    ax.set_xlim(8,8.1)
    ax.set_ylim(-6,-5.8)
    #ax.set_visible(False)
    return lgd


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
  meraxes_loc='/output/meraxes.hdf5'
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  prop='BlackHoleMass'
  
  filename='med_seed'#'M18_reion'
  filename_def='M18_reion'#'dragons10'
  fname_d=data_folder+filename_def+meraxes_loc
  fname_1=data_folder+filename+meraxes_loc
  #fname_2=data_folder+filename125+meraxes_loc
  
  fig, axes = plt.subplots(1, 5,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  j=0
  kk=0
  labels_dict={}
  for snapshot in [78,100,116,134]:
  #for snapshot in [78,116,158,213]:
    ii+=1
    if (snapshot<160):
      gals_default=load_data(filename_def,snapshot,prop,cosmo)
      gals_bulges=load_data(filename,snapshot,prop,cosmo)

    if (snapshot<160):
      calculateQLF(gals_bulges,fname_1,'UV',axes[ii],**{'linestyle':'-','label':'M18 Meraxes','linewidth':1.5,'color':'k','zorder':100})
      calculateQLF(gals_default,fname_d,'UV',axes[ii],**{'linestyle':'--','label':'Q17 Meraxes','linewidth':1.5,'color':'gray','zorder':102})
    if (snapshot==158)|(snapshot==213):
      kk=plot_obsQLF_B(axes[ii],redshift[snapshot],0,labels_dict,kk)
      axes[ii].set_xlabel(r'$M_{B}$)')
    else:
      kk=plot_obsQLF_UV(axes[ii],redshift[snapshot],0,labels_dict,kk)
      axes[ii].set_xlabel(r'$M_{1450}$')

    if ii==0:
      axes[ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$',usetex=True)
    else:
      axes[ii].set_yticklabels([])
    #axes[j,ii].set_title('$z=${}'.format(redshift[snapshot]))
    axes[ii].set_xlim([-18,-29])
    axes[ii].set_ylim([-10,-3.7])
    axes[ii].text(-23.5, -4.4, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
  
  lgd=make_legend(axes[4])
  axes[4].axis('off')
  plt.tight_layout()
  plt.savefig('/home/mmarshal/results/plots/UVQLF.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
  
##XRAY
  figX, axesX = plt.subplots(1, 4,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  for snapshot in [116,134,158]:
    ii+=1
    gals_default=load_data(filename_def,snapshot,prop,cosmo)
    gals_bulges=load_data(filename,snapshot,prop,cosmo)

    calculateQLF(gals_bulges,fname_1,'hardX',axesX[ii],**{'linestyle':'-','label':'M18 Model','linewidth':2,'color':'k','zorder':100})
    calculateQLF(gals_default,fname_d,'hardX',axesX[ii],**{'linestyle':'--','label':'Q17 Meraxes','linewidth':2,'color':'gray','zorder':102})
    plot_obsQLF_X(axesX[ii],redshift[snapshot])
    axesX[ii].set_xlabel(r'$L_{\textrm{2-10keV}}/L_\odot$')

    if ii==0:
      axesX[ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$',usetex=True)
    else:
      axesX[ii].set_yticklabels([])
    axesX[ii].set_xlim([7.5,12])
    axesX[ii].set_ylim([-8,-2])
    axesX[ii].text(10.8, -2.7, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
    #axes[j,ii].grid(color=[0.8,0.8,0.8],linestyle='--') 

  lgdX=axesX[-2].legend(loc=(1.05,0.2))
  axesX[-1].axis('off')
  plt.tight_layout()
  plt.savefig('/home/mmarshal/results/plots/XrayQLF.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

