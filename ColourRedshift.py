import numpy as np
from dragons import meraxes
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import magcalc as mc
sys.path.append('Yuxiang/')
from _function import _function

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (8,3.2)
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


def load_mags(filename,snapshot,zlist):
  redshift=meraxes.io.grab_redshift('/home/mmarshal/data_dragons/'+filename+'/output/meraxes.hdf5', snapshot)
  #if np.mod(redshift,0.05)>0.015:
  #  return 0
  print(snapshot)
  mags=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')
  #print(list(mags.columns.values))
  zlist.append(redshift)
  mags_dust={}
  filters=['SDSSg','SDSSi','SDSSr']
  central_wavelength=[4700,7500,6200]
  i=-1
  for filt in filters:
    i+=1
    mags_dust[filt] = np.array(mags[filt]) + mc.reddening(central_wavelength[i]/(1.+redshift), mags['M1600-100'], z = redshift)
  return mags_dust


if __name__=="__main__":
  #Setup
  filename='tuned_reion_T125'
  g_minus_i={}
  i=0
  zlist=[]
  d={}
  snaps=[249,245,241,237,234,231,228,225,222,219,216]
  snaps=range(213,250)
  for snapshot in snaps:
    mags=load_mags(filename,snapshot,zlist)
    if mags!=0:
      g_minus_i=mags['SDSSg']-mags['SDSSi']
      r=mags['SDSSr']
      g_minus_i=g_minus_i[r<19.8]
      d[zlist[i]]=g_minus_i
      #plt.plot(1,1)
      #plt.plot(np.ones(len(g_minus_i))*zlist[i],g_minus_i,'b.',markersize=2)
      plt.violinplot(np.array(g_minus_i),[zlist[i]],widths=0.05)
      i+=1
  plt.xlim(-0.02,0.52)
  plt.ylim(-0.1,3.1)
  plt.xlabel('Redshift')
  plt.ylabel('g-i colour')
  plt.tight_layout()
  #plt.show()
  #plt.savefig('ColourRedshift_violin.png',format='png')
  df = pd.DataFrame.from_dict(d, orient='index').transpose().fillna('') 
  df.to_csv('ColourRedshift.csv', index=False)
