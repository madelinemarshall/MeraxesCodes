import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
#import cosmolopy
#from _plot_obsGLF import plot_obsGLF

def load_data(filename,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)]
  return gals


def load_mags(filename,snapshot):
  lumins=pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')
  return np.array(lumins['M1600'])


def plot_LF(lums,axes,lab,sty):
    maxval=-7
    minval=-23
    hist, bin_edges = np.histogram(lums,range=(minval,maxval),bins=20)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    axes.plot(Max,np.log10(hist/(bin_edges[1]-bin_edges[0])/100.**3),label=lab,linewidth=2.0,linestyle=sty)
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
  redshift={63:7,78:6,100:5,115:4,134:3,158:2}
  col={63:'turquoise',78:'r',100:'b',115:'g',134:'purple',158:'orange'} 
  filename='bulges_update1102_full'
  filename2='bulges_crotonSF'

  fig,axes=plt.subplots(1,5)

  ii=-1
  j=0
  for snapshot in [63,78,100,115,158]:
    ii+=1
  #  if ii==3:
  #    j+=1
  #    ii=0
    mags_default=load_mags('default',snapshot)
    mags_bulges=load_mags(filename,snapshot)
    mags_2=load_mags(filename2,snapshot)

    plot_LF(mags_default,axes[ii],'Default Meraxes','-')
    plot_LF(mags_bulges,axes[ii],'Bulge Model','--')
    plot_LF(mags_2,axes[ii],'Bulge Model - Croton SF',':')

 #   plot_obsGLF(axes[ii],redshift[snapshot],hubble_h=cosmo['h'],markersize=7,legend=True,silent=False,color=[0.5,0.5,0.5],alpha=1.0)

    axes[ii].set_xlabel(r'$M_{1600}$')
    axes[ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    axes[ii].set_xlim([-23,-7])
    axes[ii].set_ylim([-5,0])
    axes[ii].set_title('z={}'.format(redshift[snapshot]))

  #axes[2].set_ylim([-7,-1])
  plt.legend(['Default Meraxes','Bulge Model','Bulge Model - Croton SF'])
  plt.tight_layout()
  plt.show()

