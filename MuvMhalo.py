import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/simulation_codes/Yuxiang/')
sys.path.append('/home/mmarshal/simulation_codes')
from _plot_obsGSMF import plot_obsGSMF
from scipy.optimize import curve_fit
import ContourPlot as cp
import matplotlib.lines as mlines
import magcalc as mc

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,3.2)
#matplotlib.rcParams['figure.figsize'] = (7.2,4)
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
      snapshot=snapshot,props=['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','BulgeStellarMass','Type','ColdGas'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[(gals['StellarMass']*1e10>1e6)]
  return gals


def load_mags(filename,snapshot):
  redshift={52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  #if (filename=='bulges_update1102_full')&(snapshot==158):
  #  dust=pd.read_csv('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'_dust.txt',sep=' ')
  #  i775=np.array(dust)[:,0]
  #  m1600=np.array(dust)[:,1]
  #  return m1600
  #else:
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust




if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,213:0.55,250:0}
  snapshot=52
  prop='StellarMass'
  color={52:'C0',63:'C1',78:'C2',100:'C3',116:'C4',134:'aqua',158:'pink',213:'k',250:'k'}
  filename='tuned_reion_T125'
  meraxes_loc='/output/meraxes.hdf5'
 
  gals=load_data(filename,snapshot)
  mags=load_mags(filename,snapshot)
  mags=mags[gals["Type"]==0]
  gals=gals[gals["Type"]==0]
  ##PLOT
  #plt.scatter(gals['Mvir']*1e10,mags,s=1,c=gals['BulgeStellarMass']/gals['StellarMass'])
  plt.scatter(gals['Mvir']*1e10,mags,s=2,c=np.log10(gals['ColdGas']/gals['StellarMass']))
  plt.colorbar()
  plt.ylim([-8,-23])
  plt.xlim([10**9.5,1e12])
  plt.xscale('log') 
  plt.show()
