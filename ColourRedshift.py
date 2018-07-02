import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import magcalc as mc
sys.path.append('Yuxiang/')
from _function import _function

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7,3.2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
color=['#377eb8','#ff7f00','#4daf4a','#984ea3']#,\'#377eb8'
                  #134:,158:'#f781bf',194:'#a65628',213:'black'}
cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
}



def load_UVmags(filename,snapshot):
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def load_JWSTmags(filename,snapshot):
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  mags= pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'_JWSTfilters.hdf5')
  #print(list(mags.columns.values))
  mags_dust={}
  filters=['SDSSg','SDSSi','SDSSr']
  central_wavelength[4700,7500,6200]
  i=-1
  for filt in filters:
    i+=1
    mags_dust[filt] = np.array(mags[filt]) + mc.reddening(central_wavelength[i]/(1.+redshift[snapshot]), MUV, z = redshift[snapshot])
  return mags_dust


if __name__=="__main__":
  #Setup
  filename='tuned_reion_T125'
  for snapshot in range(213,250):

    JWST_mags=load_JWSTmags(filename,snapshot)
    filters=['JWST_F277W','JWST_F444W'] 
    central_wavelength=[27700,44400]
    mag=np.zeros((np.size(gals['ID']),len(filters)))
    for ii in range(0,np.size(gals['ID'])):
      if AGN_condition[ii]:
        i=-1
        for ff in filters:
          i+=1
          mag[ii,i]=JWST_mags[ff][ii]
  
