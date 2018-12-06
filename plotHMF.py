import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/simulation_codes/Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF

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
  gals=gals[gals['Mvir']*1e10>1e10]
  return gals


def plot_SMF(gals,prop,boxwidth,axes,**kwargs):
    maxval=15#np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10)) 
    minval=12#np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),bins=np.log10(dat[:,0]/0.678))
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    axes.plot(Max,np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3),**kwargs)
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
  prop='Mvir'
  
  filename='tuned_reion'
  meraxes_loc2='/output/meraxes.hdf5'
  vol=100
  filename125='tuned_reion_T125'
  default='default_reion'
  fig, axes = plt.subplots(2, 2 ,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  j=0
  for snapshot in [63,78,100,116]:#,134,158,213]:
    ii+=1
    if ii==2:
      j+=1
      ii=0
    if snapshot==213:
      dat=np.loadtxt('data/hmf_z0p5.txt')
    elif snapshot == 158:
      dat=np.loadtxt('data/hmf_z2.txt')
    elif snapshot == 116:
      dat=np.loadtxt('data/hmf_z4.txt')
    else:
      dat=np.loadtxt('data/hmf_z6.txt')
    axes[j,ii].plot(np.log10(dat[:,0]/0.678),np.log10(dat[:,7]*(0.678)**4))
    axes[j,ii].set_xlim([10,13])
    axes[j,ii].set_ylim([-6,0])
    
    if (snapshot!=194)&(snapshot!=213):
      gals_bulges=load_data(filename,meraxes_loc2,snapshot,prop,cosmo)
    gals_125=load_data(filename125,meraxes_loc2,snapshot,prop,cosmo)


    if (snapshot!=194)&(snapshot!=213):
      plot_SMF(gals_bulges,prop,vol,axes[j,ii],**{'linestyle':'-','label':'Modified Meraxes','linewidth':2.5,'color':'Purple','zorder':100})
      plot_SMF(gals_125,prop,125/cosmo['h'],axes[j,ii],**{'linestyle':'-','label':'Modified Meraxes\n (Tiamat-125-HR)','linewidth':0.8,'color':'Purple','zorder':101})
    else:
      plot_SMF(gals_125,prop,125/cosmo['h'],axes[j,ii],**{'linestyle':'-','label':'Modified Meraxes\n (Tiamat-125-HR)','linewidth':0.8,'color':'Purple','zorder':100})

    ##TIAMAT 125, boxwidth=125/cosmo['h']
    if j==1:
      axes[j,ii].set_xlabel(r'$\log(M_{\mathrm{halo}}/M_\odot$)')
    else:
      axes[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    else:
      axes[j,ii].set_yticklabels([])
    #axes[j,ii].set_xlim([7.5,12.3])
    #axes[j,ii].set_ylim([-5.8,-1.2])
    axes[j,ii].text(11, -1, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
    axes[j,ii].grid(color=[0.8,0.8,0.8],linestyle='--') 
  
  #fig.subplots_adjust(hspace=0, wspace=0)
  plt.tight_layout()
  #plt.savefig('SMF.pdf',format='pdf')
  plt.show()
