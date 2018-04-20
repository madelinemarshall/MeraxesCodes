import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
#import cosmolopy
from _plot_obsGLF import plot_obsGLF
import magcalc as mc

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,2.8)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def load_data(filename,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)]
  return gals


def load_mags(filename,snapshot):
  redshift={52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  #if (filename=='bulges_update1102_full')&(snapshot==158):
  #  dust=pd.read_csv('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'_dust.txt',sep=' ')
  #  i775=np.array(dust)[:,0]
  #  m1600=np.array(dust)[:,1]
  #  return m1600
  #else:
  MUV=pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def plot_LF(lums,axes,**kwargs):
    maxval=-7
    minval=-23
    hist, bin_edges = np.histogram(lums,range=(minval,maxval),bins=40)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    axes.plot(Max,hist/(bin_edges[1]-bin_edges[0])/100.**3,**kwargs)
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
  redshift={52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  col={63:'turquoise',78:'r',100:'b',115:'g',134:'purple',158:'orange'} 
  filename='bulges_fullreion'

  fig,axes=plt.subplots(1,5,gridspec_kw = {'wspace':0})

  ii=-1
  j=0
  for snapshot in [52,63,78,100,115]:
    ii+=1
  #  if ii==3:
  #    j+=1
  #    ii=0
    mags_default=load_mags('default_fullreion',snapshot)
    mags_bulges=load_mags(filename,snapshot)
    #mags_2=load_mags(filename2,snapshot)

    plot_LF(mags_bulges,axes[ii],**{'linestyle':'-','label':'Modified Meraxes','linewidth':2.5,'color':'Purple'})
    plot_LF(mags_default,axes[ii],**{'linestyle':':','label':'Default Meraxes','linewidth':2.5,'color':'C9'})
    #plot_LF(mags_2,axes[ii],'Bulge Model - Croton SF',':')

    plot_obsGLF(axes[ii],redshift[snapshot],hubble_h=cosmo['h'],markersize=3,legend=True,silent=True,color=[0.5,0.5,0.5],alpha=1.0)

    if ii==0:
      axes[ii].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    else:
      axes[ii].set_yticklabels([])
    axes[ii].set_xlabel(r'$M_{1600}$')
    axes[ii].set_xlim([-23,-9])
    axes[ii].set_ylim([9e-7,1])
    #axes[ii].set_title('z={}'.format(redshift[snapshot]))
    axes[ii].text(-21.5, 0.2, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
    axes[ii].grid(color=[0.8,0.8,0.8],linestyle='--') 
  #plt.tight_layout()
  fig.subplots_adjust(right=0.9,bottom=0.2)
  lgd=plt.legend(fontsize='small',loc=[0.1,0.02])
  #axes[2].set_ylim([-7,-1])

  #axes[4].legend(['Bulge Model','Default Meraxes'])#,'Bulge Model - Croton SF'])
  plt.savefig('GLF.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

