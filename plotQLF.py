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
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,5)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

from astrodatapy.number_density import number_density

colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','--','-.',':']*4


def load_data(filename,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'BlackHoleAccretedColdMass','GhostFlag','dt'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['BlackHoleMass']*1e10>1e6]
  return gals


def plot_obsQLF(ax,z,make_legend,labels_dict):
    feature = 'QLF_UV'
    xlim    = (-15, -29)
    ylim    = (-10, -3.5)

    obs    = number_density(feature=feature,z_target=z,quiet=1,h=cosmo['h'])
    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Qin2017_Tiamat") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        if label in labels_dict:
            if make_legend:
              continue
            datatype   = labels_dict[label][0]
            color   = labels_dict[label][1]
            marker   = labels_dict[label][2]
            linestyle   = labels_dict[label][3]
        else:
            datatype   = obs.target_observation['DataType'][ii]
            color      = colors[ii]#'gray'
            marker     = markers[j_data]
            linestyle  = linestyles[k_func]
            labels_dict[label]=[datatype,color,marker,linestyle]
        data[:,1:] = np.log10(data[:,1:])
        if datatype == 'data':
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker)
            j_data +=1
        elif datatype == 'dataULimit':
            ax.errorbar(data[:,0],  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker)
            j_data +=1
        else:
            ax.plot(data[:,0],data[:,1],label=label,color=color,linestyle=linestyle,lw=3)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1


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
  
  filename='tuned_reion'#'bulges_correctBHMF'
  filename125='tuned_reion_T125'
  fname_d=data_folder+'default'+meraxes_loc
  fname_1=data_folder+filename+meraxes_loc
  fname_2=data_folder+filename125+meraxes_loc
  
  fig, axes = plt.subplots(1, 5,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  j=0
  labels_dict={}
  for snapshot in [63,78,100,116]:#,134]:#[78,116,158,213]:
    ii+=1
    #if (snapshot!=194)&(snapshot!=213):
      #gals_default=load_data('default',snapshot,prop,cosmo)
    gals_bulges=load_data(filename,snapshot,prop,cosmo)
    #gals_125=load_data(filename125,snapshot,prop,cosmo)

    #if (snapshot!=194)&(snapshot!=213):
    calculateQLF(gals_bulges,fname_1,'UV',axes[ii],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'k','zorder':100,'alpha':0.3})
      #calculateQLF(gals_default,fname_d,'UV',axes[ii],**{'linestyle':'-','label':'Default Meraxes','linewidth':2.5,'color':'C9','zorder':102})
      #calculateQLF(gals_125,fname_2,'UV',axes[ii],**{'linestyle':'-','label':'Bulge Model\n (Tiamat-125-HR)','linewidth':0.8,'color':'Purple','zorder':101})
    #else:
      #calculateQLF(gals_125,fname_2,'UV',axes[ii],**{'linestyle':'-','label':'Bulge Model\n (Tiamat-125-HR)','linewidth':0.8,'color':'Purple','zorder':100})
    #plot_obsQLF(axes[ii],redshift[snapshot],markersize=10,legend=True,hubble_h=0.678,silent=False,color='gray',alpha=1.0)
    plot_obsQLF(axes[ii],redshift[snapshot],0,labels_dict)


    axes[ii].set_xlabel(r'$\log(M_\ast (M_\odot)$)')
    if ii==0:
      axes[ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    else:
      axes[ii].set_yticklabels([])
    #axes[j,ii].set_title('$z=${}'.format(redshift[snapshot]))
    #axes[j,ii].set_xlim([7.5,12])
    #axes[j,ii].set_ylim([-5.8,-1.2])
    axes[ii].text(-18, -9.5, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='x-large')
    #axes[j,ii].grid(color=[0.8,0.8,0.8],linestyle='--') 

  lab=axes[0].get_legend_handles_labels()[1]
  hand=axes[0].get_legend_handles_labels()[0]
  lab.append(axes[1].get_legend_handles_labels()[1])
  hand.append(axes[1].get_legend_handles_labels()[0])
  lab.append(axes[2].get_legend_handles_labels()[1])
  hand.append(axes[2].get_legend_handles_labels()[0])
  #lab.append(axes[3].get_legend_handles_labels()[1])
  #hand.append(axes[3].get_legend_handles_labels()[0])
  #lab.append(axes[4].get_legend_handles_labels()[1])
  #hand.append(axes[4].get_legend_handles_labels()[0])
  plt.legend(hand,lab)
  #plt.tight_layout()
  plt.legend(loc=(1,0.5))
  plt.show()

