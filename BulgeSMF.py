import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
import magcalc as mc
from _load_data import load_data
from astrodatapy.number_density import number_density

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.6,4)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','--','-.',':']*4


def load_mags(filename,snapshot):
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  M6300 = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['F625W']
  M1600   = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  M6300_dust = M6300 + mc.reddening(6300., M1600, z = redshift[snapshot])
  return M6300_dust


def plot_SMF(gals,prop,boxwidth,axes,**kwargs):
  if prop=='DiskStellarMass':
    dsm=gals['StellarMass'][gals['StellarMass']>0]-gals['BulgeStellarMass'][gals['StellarMass']>0]
    dsm=dsm[dsm>0]
    maxval=np.nanmax(np.log10(dsm*1e10)) 
    minval=np.nanmin(np.log10(dsm*1e10))
    hist, bin_edges = np.histogram(np.log10(dsm*1e10),range=(minval,maxval),bins=30)
  elif prop=='InstabilityBulgeMass':
    selection=(gals['BulgeStellarMass']-gals['MergerBulgeStellarMass'])>0
    values=np.log10(gals[selection]['BulgeStellarMass']-gals[selection]['MergerBulgeStellarMass'])+10
    maxval=np.nanmax(values)
    minval=np.nanmin(values)
    hist, bin_edges = np.histogram(values,range=(minval,maxval),bins=30)
  else:
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
  axes.set_xlabel(r'$\log(M_\ast (M_\odot)$)')
  axes.set_xlim([8.2,12])
  axes.set_ylim([-5,-1.05])
  return axes



def plot_SMFs(vol,axes):
      BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
      logM=np.linspace(8,13)
      plot_SMF(gals_bulges,prop,vol,axes,**{'linestyle':'-','label':'Total','linewidth':2.5,'color':colors[0],'zorder':100})
      plot_SMF(gals_bulges[BT<0.7],prop,vol,axes,**{'linestyle':':','label':'Disk-dominated','linewidth':2.5,'color':colors[4],'zorder':100})
      plot_SMF(gals_bulges[BT>0.3],prop,vol,axes,**{'linestyle':'-.','label':'Bulge-dominated','linewidth':2.5,'color':colors[1],'zorder':100})
#      plot_SMF(gals_bulges,'BulgeStellarMass',vol,axes,**{'linestyle':'--','label':'Total Bulge','linewidth':2.5,'color':colors[1],'zorder':100})
#      plot_SMF(gals_bulges,'MergerBulgeStellarMass',vol,axes,**{'linestyle':':','label':'Merger-Driven Bulge','linewidth':2.5,'color':colors[2],'zorder':100})
#      plot_SMF(gals_bulges,'InstabilityBulgeMass',vol,axes,**{'linestyle':':','label':'Instability-Driven Bulge','linewidth':2.5,'color':colors[3],'zorder':100})
#      plot_SMF(gals_bulges,'DiskStellarMass',vol,axes,**{'linestyle':'-.','label':'Stellar Disc','linewidth':2.5,'color':colors[4],'zorder':100})



def plot_obs(ax):
  features = ['GSMF','GSMF_Disk', 'GSMF_Bulge']
  for z in [0.0,]:
    for feature, color in zip(features, [colors[0],colors[4], colors[1]]):
        obs    = number_density(feature=feature,z_target=z,quiet=1,h=cosmo['h'])
        for ii in range(obs.n_target_observation):
            data       = obs.target_observation['Data'][ii]
            label      = obs.target_observation.index[ii]
            if (feature=='GSMF')&(label!='Bell2003'):
              continue
            datatype   = obs.target_observation['DataType'][ii]
            marker     = markers[ii]
            data[:,1:] = np.log10(data[:,1:])
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker)


def plot_GAMMA(ax):
  M=[8.05,8.25,8.35,8.45,8.55,8.65,8.75,8.85,8.95,9.05,9.15,9.25,9.35,9.45,9.55,9.65,9.75,9.85,9.95,10.05,10.15,10.25,10.35,10.45,10.55,10.65,10.75,10.85,10.95,11.05,11.15,11.25,11.35,11.45,11.55]
  bulge=[1.75e-04,3.5e-04,3.5e-04,6.47e-04,1.61e-03,8.99e-04,8.99e-04,1.52e-03,2.39e-03,1.7e-03,1.55e-03,2.47e-03,2.65e-03,3e-03,3.44e-03,3.14e-03,3.98e-03,3.85e-03,4.51e-03,4.12e-03,3.17e-03,3.64e-03,2.77e-03,2.21e-03,1.37e-03,1.24e-03,9.96e-04,9.31e-04,9.96e-04,7.5e-04,7.76e-04,6.04e-04,2.46e-04,1.77e-04,3.57e-05]
  ax.plot(M,np.log10(np.array(bulge)),'*',color=colors[1])

  M=[8.05,8.15,8.25,8.35,8.45,8.55,8.65,8.75,8.85,8.95,9.05,9.15,9.25,9.35,9.45,9.65,9.85,9.95,10.05,10.15,10.25,10.35,10.45,10.55,10.65,10.75,10.85,10.95,11.05,11.15,11.25]
  disk=[2.81e-02,1.85e-02,2.98e-02,2.69e-02,1.39e-02,1.76e-02,1.83e-02,1.01e-02,1.36e-02,1.67e-02,9.34e-03,7.79e-03,7.97e-03,7.36e-03,6.14e-03,5.23e-03,4.83e-03,5.29e-03,4.12e-03,4.62e-03,4.67e-03,3.07e-03,3.03e-03,2.62e-03,2.08e-03,1.59e-03,1.13e-03,5.64e-04,2.82e-04,2.46e-04,1.06e-04]
  ax.plot(M,np.log10(np.array(disk)),'*',color=colors[4])


if __name__=="__main__":
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,242:0.1,250:0}
  snapshot=250
  prop='StellarMass'
  
  filename='draft2_reion_T125'
  vol=125/cosmo['h']

  fig,axes=plt.subplots(2,1,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,1]})
  gals_bulges=load_data(filename,snapshot,[prop,'GhostFlag','BulgeStellarMass','Type','MergerBulgeStellarMass'],centrals=True)
  BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
    
  plot_SMFs(vol,axes[0])
  plot_obs(axes[0])
  plot_GAMMA(axes[0])
 
  axes[0].set_xlabel(r'$\log(M_{\mathrm{component}} /M_\odot$)')
  axes[0].set_ylabel(r'$\log(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
  plt.tight_layout()
  axes[0].legend(loc=(-0.2,-0.4),ncol=2)
  axes[1].axis('off')
  
  plt.savefig('/home/mmarshal/results/plots/Paper1/BulgeSMF.pdf', format='pdf')
  plt.show()
