import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/PhD/simulation_codes/Yuxiang/')
from _plot_obsBHMF import plot_obsBHMF
import warnings
warnings.filterwarnings("ignore")

#Sets plot defaults
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,2.8)
#matplotlib.rcParams['figure.figsize'] = (7.2,5)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals


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
  #meraxes_loc='/output/meraxes.hdf5'
  meraxes_loc='/output/'+str(sys.argv[1])+'.hdf5'
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'pink',213:'black'}
  prop='BlackHoleMass'
  
  filename_default='default'
  filename='tune_higherseed'
  boxwidth=125/cosmo['h']
  #filename='bulges_correctBHMF'#_split'#correctBHMF'
  filename125='tune_higherseed'
  plot_z0=1

  if plot_z0:
    snapshots=[213]#[63,78,100,116,134,158,213]
  else:
    snapshots=[63,78,100,116,134,158]

  fig, axes = plt.subplots(1,1)
  for snapshot in snapshots:
    if snapshot!=213:
      #gals_default=load_data(filename_default,snapshot,prop,cosmo)
      gals_bulges=load_data(filename,snapshot,prop,cosmo)
      plot_BHMF(gals_bulges,prop,boxwidth,axes,**{'linestyle':'-','label':'$z={}$'.format(redshift[snapshot]),'linewidth':2.5,'color':color[snapshot],'zorder':101})
      #plot_BHMF(gals_125,prop,125/cosmo['h'],axes[1],**{'linestyle':'-','label':'Bulge Model\n (Tiamat-125-HR)','linewidth':0.8,'color':'Purple','zorder':101})
      #plot_BHMF(gals_default,prop,100,axes[0],**{'linestyle':'-','label':'Default Meraxes','linewidth':2.5,'color':'C9','zorder':102})
    else:
      gals_125=load_data(filename125,snapshot,prop,cosmo)
      plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$z={}$ (Tiamat-125-HR)'.format(redshift[snapshot]),'linewidth':2.5,'color':'black','zorder':101})
  axes.set_xlabel(r'$\log (\mathrm{M}_{\mathrm{BH}}/M_\odot)$')
  axes.set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
  axes.set_xlim([5,10])
  axes.set_ylim([1e-6,0.5])
  plot_obsBHMF(axes,0.5,hubble_h=cosmo['h'],markersize=3,legend=True,silent=False,color='gray',alpha=1.0)
  plot_obsBHMF(axes,0,hubble_h=cosmo['h'],markersize=3,legend=True,silent=False,color=[0.4,0.4,0.4],alpha=1.0)

  #plot_obsBHMF(axes[1],0.5,hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color='gray',alpha=1.0)
  #plot_obsBHMF(axes[1],0,hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color='lightgreen',alpha=1.0)
  fig.subplots_adjust(right=0.62,bottom=0.2,top=0.95)
  lgd=plt.legend(loc=(1.02,-0.08),fontsize='small')
  #fig.subplots_adjust(hspace=0, wspace=0)
  #plt.tight_layout()
  plt.savefig('/home/mmarshal/PhD/plots/tune_higherseed/BHMF_'+str(sys.argv[2])+'.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  #plt.savefig('BHMF.pdf',format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()


