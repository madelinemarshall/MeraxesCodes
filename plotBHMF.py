import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
from astrodatapy.number_density import number_density

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,3)
#matplotlib.rcParams['figure.figsize'] = (7.2,5)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','--','-.',':']*4


def load_data(filename,meraxes_loc,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'Mvir','StellarMass','GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['StellarMass']*1e10>1e7]
  #gals=gals[gals['StellarMass']*1e10<10**11]
  #gals=gals[gals['Mvir']*1e10<10**13.5]
  return gals


def plot_obs(ax,z):
    obs    = number_density(feature='BHMF',z_target=z,quiet=0,h=cosmo['h'])
    xlim    = (6, 10)
    ylim    = (-6, -0.5)
    xlabel  = r"$\log_{10}[M_{\rm BH}/{\rm M_{\odot}}]$"
    ylabel  = r"$\log_{10}[\rm \phi/Mpc^{-3} dex^{-1}]$"

    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Shen2012") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        label      = label[:label.rfind('20')]+' et al. ('+label[label.rfind('20'):]+')'
        datatype   = obs.target_observation['DataType'][ii]
        color      = obs_color[z]
        marker     = markers[j_data]
        linestyle  = linestyles[k_func]
        data[:,1:] = np.log10(data[:,1:])
        if datatype == 'data':
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker,markersize=4)
            j_data +=1
        elif datatype == 'dataULimit':
            ax.errorbar(data[:,0],  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker)
            j_data +=1
        else:
            ax.plot(data[:,0],data[:,1],label=label,color=color,linestyle=linestyle,lw=3)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1



def plot_BHMF(gals,prop,boxwidth,axes,**kwargs):
    maxval=np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10)) 
    minval=np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=30)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    Max=Max[hist>5]
    hist=hist[hist>5]
    hist_plus=hist+np.sqrt(hist)
    hist_minus=hist-np.sqrt(hist)
    phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    phi_plus=np.log10(hist_plus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    phi_minus=np.log10(hist_minus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    axes.plot(Max,phi,**kwargs)
    axes.fill_between(Max, phi_plus, phi_minus, alpha=0.5,color=kwargs['color'],label='__nolegend__')

    return axes


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
  #meraxes_loc='/output/meraxes.hdf5'
  meraxes_loc='/output/'+'meraxes'+'.hdf5'
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,192:1,250:0}
  color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'pink',192:'purple',250:'black'}
  obs_color={0.5:[0.2,0.2,0.2],0:'gray'}
  color={63:'#e41a1c',78:'#377eb8',100:'#4daf4a',116:'#984ea3',\
                  134:'#ff7f00',158:'#f781bf',192:'#a65628',250:'black'}
  prop='BlackHoleMass'
  
  
  filename_default='dragons10'
  filename_default125='dragons10_T125'
  filename='tuned_reion'
  meraxes_loc2='/output/'+'meraxes'+'.hdf5'
  boxwidth=100
  filename125='tuned_reion_T125'
  plot_z0=1

  if plot_z0:
    snapshots=[63,78,100,116,134,158,192,250]
  else:
    snapshots=[63,78,100,116,134,158]

  fig, axes = plt.subplots(1,1)
  for snapshot in snapshots:
    if (snapshot<164):
      #gals_default=load_data(filename_default,meraxes_loc,snapshot,prop,cosmo)
      gals_bulges=load_data(filename,meraxes_loc2,snapshot,prop,cosmo)
      plot_BHMF(gals_bulges,prop,boxwidth,axes,**{'linestyle':'-','label':'$z={}$'.format(redshift[snapshot]),'linewidth':2.5,'color':color[snapshot],'zorder':101})
      #gals_125=load_data(filename125,meraxes_loc2,snapshot,prop,cosmo)
      #plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':':','label':'$z={}$ (Tiamat-125-HR)'.format(redshift[snapshot]),'linewidth':2.5,'color':color[snapshot],'zorder':101})
      #plot_BHMF(gals_125,prop,125/cosmo['h'],axes[1],**{'linestyle':'-','label':'Bulge Model\n (Tiamat-125-HR)','linewidth':0.8,'color':'Purple','zorder':101})
      #plot_BHMF(gals_default,prop,100,axes,**{'linestyle':'-','label':'Default Meraxes','linewidth':2.5,'color':'C9','zorder':102})
    else:
      gals_125=load_data(filename125,meraxes_loc2,snapshot,prop,cosmo)
      plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$z={}$ (Tiamat-125-HR)'.format(redshift[snapshot]),'linewidth':2.5,'color':color[snapshot],'zorder':101})
      #gals_def=load_data(filename_default125,meraxes_loc,snapshot,prop,cosmo)
      #plot_BHMF(gals_def,prop,125/cosmo['h'],axes,**{'linestyle':':','label':'$z={}$ (Tiamat-125-HR) \n- Default Meraxes'.format(redshift[snapshot]),'linewidth':2.5,'color':'k','zorder':101})
  plot_obs(axes,0)
  #plot_obs(axes,0.5)
  axes.set_xlabel(r'$\log (\mathrm{M}_{\mathrm{BH}}/M_\odot)$')
  axes.set_ylabel(r'$\log (\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
  axes.set_xlim([6,10])
  axes.set_ylim([-5,-2])
  
  #plot_obsBHMF(axes[1],0.5,hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color='gray',alpha=1.0)
  #plot_obsBHMF(axes[1],0,hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color='lightgreen',alpha=1.0)
  fig.subplots_adjust(right=0.62,bottom=0.2,top=0.95)
  lgd=plt.legend(loc=(1.02,0.1),fontsize='small')
  #fig.subplots_adjust(hspace=0, wspace=0)
  #plt.tight_layout()
  ##plt.savefig('/home/mmarshal/PhD/plots/BHMF.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.savefig('/home/mmarshal/results/plots/BHMF.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()


