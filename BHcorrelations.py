import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,6)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','-.']*4


def load_data(filename,meraxes_loc,snapshot,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals


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
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,250:0}
  
  filename='tuned_reion'
  meraxes_loc2='/output/meraxes.hdf5'
  vol=100

  snapshot=158
  for mass_cut in [9,10]:
    gals=load_data(filename,meraxes_loc2,snapshot,cosmo)
    gals=gals[gals['StellarMass']*1e10>10**mass_cut]
    gals=gals[gals['StellarMass']*1e10<10**(mass_cut+0.1)]
    if mass_cut == 10:
      gals=gals[gals['BlackHoleMass']*1e10>1e6]
    gals=gals[gals['Type']==0]
    
    props=['Spin','Mvir','ColdGas','GasDiskScaleLength','StellarDiskScaleLength','BulgeStellarMass','Sfr','MetalsStellarMass']
 
    fig, axes = plt.subplots(int(np.ceil(len(props)/4)), 4,gridspec_kw = {'wspace':0})

    ii=-1
    jj=0
    for prop in props:
      ii+=1
      if ii==4:
        jj+=1
        ii=0
      if prop in ['Mvir','HotGas','ColdGas','BulgeStellarMass','MergerBulgeStellarMass','MetalsStellarMass','Sfr']:
        axes[jj,ii].plot(np.log10(gals[prop][gals[prop]>0]),np.log10(gals['BlackHoleMass'][gals[prop]>0]*1e10),'.')
        axes[jj,ii].set_xlabel('log('+prop+')')
        corr=np.corrcoef(np.log10(gals[prop][gals[prop]>0]),np.log10(gals['BlackHoleMass'][gals[prop]>0]*1e10))[0,1]
        axes[jj,ii].text(0.1,0.1,'r={0:.3f}'.format(corr),weight='bold',size='large',transform=axes[jj,ii].transAxes)
      else:
        axes[jj,ii].plot(gals[prop],np.log10(gals['BlackHoleMass']*1e10),'.')
        axes[jj,ii].set_xlabel(prop)
        corr=np.corrcoef(gals[prop],np.log10(gals['BlackHoleMass']*1e10))[0,1]
        axes[jj,ii].text(0.1,0.1,'r={0:.3f}'.format(corr),weight='bold',size='large',transform=axes[jj,ii].transAxes)
      if ii==0:
        axes[jj,ii].set_ylabel('log(BlackHoleMass)')
      else:
        axes[jj,ii].set_yticks([])


    plt.tight_layout()
    #plt.savefig('/home/mmarshal/results/plots/SMF.pdf',format='pdf')
    plt.show()
