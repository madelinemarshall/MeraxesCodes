import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/PhD/simulation_codes/Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF
from _plot_obsBHMF import plot_obsBHMF

#Sets plot defaults
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,4)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,meraxes_loc,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals


def plot_SMF(gals,prop,boxwidth,axes,**kwargs):
    maxval=np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10)) 
    minval=np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=40)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    axes.plot(Max,np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3),**kwargs)
    return axes


def plot_BHMF(gals,prop,boxwidth,axes,**kwargs):
    maxval=np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10)) 
    minval=np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=30)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    Max=Max[hist>5]
    hist=hist[hist>5]
    axes.plot(Max,(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3),**kwargs)
    axes.set_yscale('log')
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
  meraxes_loc='/output/meraxes.hdf5'
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  prop='StellarMass'
  
  filename125=''
  default='default_t125'
  fig, axes = plt.subplots(1, 1,gridspec_kw = {'wspace':0, 'hspace':0})
  for snapshot in [213]:
    gals_default=load_data(default,meraxes_loc,snapshot,prop,cosmo)
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p2/meraxes.hdf5',snapshot,prop,cosmo)
    plot_SMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.2$','linewidth':2.5})
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p1/meraxes.hdf5',snapshot,prop,cosmo)
    plot_SMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.1$','linewidth':2.5})
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p1_RME0p006/meraxes.hdf5',snapshot,prop,cosmo)
    plot_SMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.1$, RME=0.006','linewidth':2.5})
    gals_125=load_data(filename125,'tune_RFE/output/eta0p06_RME0p006/meraxes.hdf5',snapshot,prop,cosmo)
    plot_SMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.06$, RME=0.006','linewidth':2.5})
    gals_125=load_data(filename125,'tune_RFE/output/eta0p1_RME0p01/meraxes.hdf5',snapshot,prop,cosmo)
    plot_SMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.1$, RME=0.01','linewidth':2.5})
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p06/meraxes.hdf5',snapshot,prop,cosmo)
    plot_SMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'Default, $\eta=0.06$','linewidth':2.5})
    #gals_125=load_data(filename125,'/output/0p5times.hdf5',snapshot,prop,cosmo)
    #plot_SMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'0.5 times'})
  
    
    plot_SMF(gals_default,prop,125/cosmo['h'],axes,**{'linestyle':':','label':'Default Meraxes','linewidth':2.5,'color':'C9','zorder':102})
    
    plot_obsGSMF(axes,redshift[snapshot],hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color=[0.5,0.5,0.5],alpha=1.0)

    axes.legend(fontsize='small')
    axes.set_xlabel(r'$\log(M_\ast/M_\odot$)')
    axes.set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    axes.set_xlim([7,12.5])
    axes.set_ylim([-5.8,-1.2])
    axes.grid(color=[0.8,0.8,0.8],linestyle='--') 

  #fig.subplots_adjust(hspace=0, wspace=0)
  plt.tight_layout()
  #plt.savefig('SMF.pdf',format='pdf')
  plt.show()
  

  prop='BlackHoleMass'
  
  filename125=''
  default='default_t125'
  fig, axes = plt.subplots(1, 1,gridspec_kw = {'wspace':0, 'hspace':0})
  for snapshot in [213]:
    gals_default=load_data(default,meraxes_loc,snapshot,prop,cosmo)
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p2/meraxes.hdf5',snapshot,prop,cosmo)
    plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.2$','linewidth':2.5})
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p1/meraxes.hdf5',snapshot,prop,cosmo)
    plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.1$','linewidth':2.5})
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p1_RME0p006/meraxes.hdf5',snapshot,prop,cosmo)
    plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.1$, RME=0.006','linewidth':2.5})
    gals_125=load_data(filename125,'tune_RFE/output/eta0p06_RME0p006/meraxes.hdf5',snapshot,prop,cosmo)
    plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.06$, RME=0.006','linewidth':2.5})
    gals_125=load_data(filename125,'tune_RFE/output/eta0p1_RME0p01/meraxes.hdf5',snapshot,prop,cosmo)
    plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$\eta=0.1$, RME=0.01','linewidth':2.5})
    
    gals_125=load_data(filename125,'tune_RFE/output/eta0p06/meraxes.hdf5',snapshot,prop,cosmo)
    plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'Default, $\eta=0.06$','linewidth':2.5})
    
    plot_BHMF(gals_default,prop,125/cosmo['h'],axes,**{'linestyle':':','label':'Default Meraxes','linewidth':2.5,'color':'C9','zorder':102})
    
    plot_obsBHMF(axes,0.5,hubble_h=cosmo['h'],markersize=3,legend=True,silent=False,color='gray',alpha=1.0)
    plot_obsBHMF(axes,0,hubble_h=cosmo['h'],markersize=3,legend=True,silent=False,color=[0.4,0.4,0.4],alpha=1.0)

    axes.legend(fontsize='small')
    axes.set_xlabel(r'$\log(M_\ast/M_\odot$)')
    axes.set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    axes.set_xlim([7,10])
    axes.set_ylim([1e-6,1e-2])
    axes.grid(color=[0.8,0.8,0.8],linestyle='--') 

  #fig.subplots_adjust(hspace=0, wspace=0)
  plt.tight_layout()
  #plt.savefig('SMF.pdf',format='pdf')
  plt.show()
