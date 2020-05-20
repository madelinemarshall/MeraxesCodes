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


def plot_Shankar(axes):
    shan=pd.read_csv('/home/mmarshal/simulation_codes/data/Shankar_BHMF.csv',header=None)
    axes.plot(shan[0],10**shan[1],'--',color=[0.5,0.5,0.5],label=r"Shankar et al. (2009) $z=0$",markersize=3,lw=1.5)
    axes.fill_between(shan[0],10**shan[3],10**shan[2],color=[0.5,0.5,0.5],alpha=0.4)


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
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,250:0}
  color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'pink',250:'black'}
  prop='BlackHoleMass'
  
  #filename_default='paper1_T125'
  boxwidth=125/cosmo['h']
  filename125='tuning/'
  meraxes_loc=str(sys.argv[1])
  plot_z0=1

  snapshots=[250]#[63,78,100,116,134,158,213]

  fig, axes = plt.subplots(1,1)
  for snapshot in snapshots:
    gals_125=load_data(filename125,snapshot,prop,cosmo)
    plot_BHMF(gals_125,prop,125/cosmo['h'],axes,**{'linestyle':'-','label':'$z={}$ (Tiamat-125-HR)'.format(redshift[snapshot]),'linewidth':2.5,'color':'black','zorder':101})
  axes.set_xlabel(r'$\log (\mathrm{M}_{\mathrm{BH}}/M_\odot)$')
  axes.set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
  axes.set_xlim([5,10])
  axes.set_ylim([1e-6,0.5])
  #plot_obsBHMF(axes,0.0,hubble_h=cosmo['h'],markersize=3,legend=True,silent=False,color=[0.4,0.4,0.4],alpha=1.0)
  plot_Shankar(axes)
  fig.subplots_adjust(right=0.62,bottom=0.2,top=0.95)
  lgd=plt.legend(loc=(1.02,-0.08),fontsize='small')
  plt.savefig('/home/mmarshal/results/plots/tuning_paper2/BHMF_{}.pdf'.format(int(meraxes_loc[-8:-5])),format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()


