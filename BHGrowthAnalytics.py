import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
from astrodatapy.number_density import number_density
import itertools

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.3,6)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
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
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','-.']*4

###Best fit params
kc=0.03
ki=0.02


def load_data(filename,snapshot):
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['Vvir','StellarMass','BlackHoleMass','GhostFlag','BulgeStellarMass'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['StellarMass']*1e10>1e7]
  #gals=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
  return gals

def merger_stellar(v,gamma):
  return (1+(280/v)**2)**-1*kc*gamma/(0.57*gamma**0.7+gamma)

def major_bulge(v,gamma):
  return (1+(280/v)**2)**-1*kc*gamma/(0.57*gamma**0.7+gamma+0.4)

def disk_bulge(v):
  return ki/(2*(1+(280/v)**2)-ki)

def disk_stellar(v):
  return ki/((1+(280/v)**2)-ki)


if __name__=="__main__":
  filename='tuned_reion_T125'
  snapshot=250
  gals=load_data(filename,snapshot)

  v=gals['Vvir']

  disk_to_min_merge_S=np.zeros(len(v))
  disk_to_max_merge_S=np.zeros(len(v))
  disk_to_min_merge_B=np.zeros(len(v))
  disk_to_max_merge_B=np.zeros(len(v))
  for ii in range(0,len(v)):
    disk_to_min_merge_S[ii]=disk_stellar(v[ii])/merger_stellar(v[ii],0.3)
    disk_to_max_merge_S[ii]=disk_stellar(v[ii])/merger_stellar(v[ii],1)
    disk_to_min_merge_B[ii]=disk_bulge(v[ii])/major_bulge(v[ii],0.3)
    disk_to_max_merge_B[ii]=disk_bulge(v[ii])/major_bulge(v[ii],1)

  fig,axes=plt.subplots(3,2)
  axes[0,0].hist((disk_to_min_merge_S[disk_to_min_merge_S!=0]),log=True)
  axes[0,1].hist((disk_to_max_merge_S[disk_to_max_merge_S!=0]),log=True)
  axes[1,0].hist(disk_to_min_merge_B[disk_to_min_merge_B!=0],log=True)
  axes[1,1].hist(disk_to_max_merge_B[disk_to_max_merge_B!=0],log=True)
  
  axes[0,0].set_xlabel('Stellar: Disk to Smallest Merger')
  axes[0,1].set_xlabel('Stellar: Disk to Largest Merger')
  axes[1,0].set_xlabel('Bulge: Disk to Smallest Merger')
  axes[1,1].set_xlabel('Bulge: Disk to Largest Merger')
  
  axes[2,0].hist(v[v<150],log=True)

  plt.show()
  print(merger_stellar(100,0.1))
  print(merger_stellar(100,0.3))
  print(merger_stellar(100,1))
  print(major_bulge(100,1))
  print(disk_bulge(100))
  print(disk_stellar(100))

