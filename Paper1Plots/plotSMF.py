import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
from astrodatapy.number_density import number_density
import itertools
from _load_data import load_data

def flip(items, ncol):
  return itertools.chain(*[items[i::ncol] for i in range(ncol)])

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.3,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','8']*4
linestyles     = ['-','-.']*4

cosmo = {'omega_M_0' : 0.308,
'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
'omega_b_0' : 0.04839912,
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}

def plot_SMF(gals,prop,boxwidth,axes,**kwargs):
    maxval=np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10))
    minval=np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=40)

    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    Max=Max[hist>10]
    hist=hist[hist>10]
    hist_plus=hist+np.sqrt(hist)
    hist_minus=hist-np.sqrt(hist)
    phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    phi_plus=np.log10(hist_plus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    phi_minus=np.log10(hist_minus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    axes.plot(Max,phi,**kwargs)
    axes.fill_between(Max, phi_plus, phi_minus, alpha=0.5,color=kwargs['color'],label='__nolegend__')
    return axes


def plot_obs(ax,z,make_legend,labels_dict,jj):
    obs    = number_density(feature='GSMF',z_target=z,quiet=1,h=cosmo['h'])
    xlim    = (7, 13)
    ylim    = (-7, -0.5)
    xlabel  = r"$\log_{10}[M_*/{\rm M_{\odot}}]$"
    ylabel  = r"$\log_{10}[\rm \phi/Mpc^{-3} dex^{-1}]$"

    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Qin2017_Tiamat") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
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
            color      = colors[jj]
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


def make_legend(ax):
    plot_SMF(gals_bulges,prop,vol,ax,**{'linestyle':'-','label':'M19 Meraxes\n(Tiamat)','linewidth':2,'color':'k','zorder':1000})
    plot_SMF(gals_default,prop,vol_def,ax,**{'linestyle':'-','label':'Q17 Meraxes\n(Tiamat)','linewidth':2,'color':[0.35,0.35,0.35],'zorder':999})
    plot_SMF(gals_125,prop,125/cosmo['h'],ax,**{'linestyle':'--','label':'M19 Meraxes\n(Tiamat-125-HR)','linewidth':1,'color':'k','zorder':1001})
    plot_SMF(gals_default,prop,125/cosmo['h'],ax,**{'linestyle':'--','label':'Q17 Meraxes\n(Tiamat-125-HR)','linewidth':1,'color':[0.35,0.35,0.35],'zorder':1001})
    labels_dict={}
    kk=0
    for snapshot in [52,63,78,100,116,134,158,250]:
      kk=plot_obs(ax,redshift[snapshot],1,labels_dict,kk)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(flip(handles, 4), flip(labels, 4),fontsize='small',ncol=4,loc=(-1.635,0))
    ax.axis('off')
    ax.set_xlim(8,8.1)
    ax.set_ylim(-6,-5.8)


if __name__=="__main__":
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,250:0}
  prop='StellarMass'

  filename='paper1_false_eta_0p4_2'
  vol=100
  filename125='paper1_T125_false_eta_0p4_2'
  default='dragons10'
  vol_def=100
  default_125='dragons10_T125'
  labels_dict={}
  fig, axes = plt.subplots(3, 4,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  j=0
  kk=0
  for snapshot in [52,63,78,100,116,134,158,250]:
    ii+=1
    if ii==4:
      j+=1
      ii=0
    if (snapshot!=194)&(snapshot!=250):
      gals_bulges=load_data(filename,snapshot,prop)
      gals_default=load_data(default,snapshot,prop)
    else:
      gals_default=load_data(default_125,snapshot,prop)
    gals_125=load_data(filename125,snapshot,prop)

    if (snapshot!=194)&(snapshot!=250):
      plot_SMF(gals_default,prop,vol_def,axes[j,ii],**{'linestyle':'-','label':'Q17 Meraxes','linewidth':1,'color':[0.35,0.35,0.35],'zorder':999})
      plot_SMF(gals_bulges,prop,vol,axes[j,ii],**{'linestyle':'-','label':'M19 Meraxes','linewidth':2,'color':'k','zorder':1000})
    else:
      plot_SMF(gals_default,prop,125/cosmo['h'],axes[j,ii],**{'linestyle':'--','label':'Q17 Meraxes\n (Tiamat-125-HR)','linewidth':1,'color':[0.35,0.35,0.35],'zorder':1001})
      plot_SMF(gals_125,prop,125/cosmo['h'],axes[j,ii],**{'linestyle':'--','label':'M19 Meraxes\n (Tiamat-125-HR)','linewidth':2,'color':'k','zorder':1002})
    kk=plot_obs(axes[j,ii],redshift[snapshot],0,labels_dict,kk)


    if j==1:
      axes[j,ii].set_xlabel(r'$\log(M_\ast/M_\odot$)')
    else:
      axes[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
    else:
      axes[j,ii].set_yticklabels([])
    axes[j,ii].set_xlim([7.5,12.3])
    axes[j,ii].set_ylim([-5.8,-1.5])
    axes[j,ii].text(8.1, -5.6, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
    axes[j,ii].grid(color=[0.8,0.8,0.8],linestyle='--')

  make_legend(axes[2,1])
  axes[2,0].axis('off')
  axes[2,1].axis('off')
  axes[2,2].axis('off')
  axes[2,3].axis('off')
  #plt.tight_layout()
  plt.savefig('/home/mmarshal/results/plots/Paper1/SMF.pdf',format='pdf')
  plt.show()
