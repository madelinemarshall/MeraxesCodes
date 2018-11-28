import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,4)
matplotlib.rcParams['lines.linewidth'] = 2.5
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,snapshot):
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
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','Vvir'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['StellarMass']*1e10>0.98*1e10]
  gals=gals[gals['StellarMass']*1e10<1.02*1e10]
  return gals

if __name__=='__main__':
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'

  redshift={63:7,78:6,100:5,116:4,134:3,158:2,250:0}
  snapshots=[63,78,100,116,134,158]
  prop='StellarMass'
  color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'pink',213:'black'}
 
  filename2='tuned_reion_T125'  
  filename='tuned_reion'

  max_v=np.zeros(len(snapshots)+1)
  min_v=np.zeros(len(snapshots)+1)
  med_v=np.zeros(len(snapshots)+1)
  ii=0
  for snap in snapshots:
    #New model
    gals=load_data(filename,snap)
    if np.size(gals)>0:
      max_v[ii]=np.max(gals['Vvir'])
      min_v[ii]=np.min(gals['Vvir'])
      med_v[ii]=np.median(gals['Vvir'])
    ii+=1

  gals=load_data(filename2,250)
  if np.size(gals)>0:
    max_v[ii]=np.max(gals['Vvir'])
    min_v[ii]=np.min(gals['Vvir'])
    med_v[ii]=np.median(gals['Vvir'])

  plt.plot(np.array(list(redshift.values())),max_v,'-') 
  plt.plot(np.array(list(redshift.values())),min_v,'-')    
  plt.plot(np.array(list(redshift.values())),med_v,'-')    
  plt.legend()
  plt.xlabel(r'z')
  plt.ylabel(r'$V_{vir}$')
  ##plt.ylabel(r'$M_\ast/M_{vir}$')
  #plt.yscale('log')
  plt.show()
