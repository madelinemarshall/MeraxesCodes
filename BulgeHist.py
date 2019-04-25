import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
import magcalc as mc
sys.path.append('Paper1Plots/')
from _load_data import load_data
from astrodatapy.number_density import number_density

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.5,3)
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
                  '#ff7f00','#a65628','#f781bf','#999999']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','--','-.',':']*4


def load_mags(filename,snapshot):
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  M6300 = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['F625W']
  M1600   = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  M6300_dust = M6300 + mc.reddening(6300., M1600, z = redshift[snapshot])
  return M6300_dust


if __name__=="__main__":
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,242:0.1,250:0}
  prop='StellarMass'
  
  filename='draft2_reion'
  vol=125/cosmo['h']

  fig, axes = plt.subplots(1, 4,gridspec_kw = {'wspace':0})#,sharey=True)
  for snapshot in [100]:
    gals_bulges=load_data(filename,snapshot,[prop,'GhostFlag','BulgeStellarMass','Type','VStellarDisk','StellarDiskScaleLength','AMstars'])
    disk_mass=(gals_bulges['StellarMass']-gals_bulges['BulgeStellarMass'])*1e10
    index=(gals_bulges['StellarDiskScaleLength']>0)&(disk_mass>0)
    gals_bulges=gals_bulges[index]
    disk_mass=disk_mass[index]
    BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
    AM=np.sqrt(gals_bulges['AMstars'][:,0]**2+gals_bulges['AMstars'][:,1]**2+gals_bulges['AMstars'][:,2]**2)
    #ratio=AM/(disk_mass*gals_bulges['VStellarDisk']*gals_bulges['StellarDiskScaleLength'])
    #print(np.median(np.log10(ratio)),np.min(np.log10(ratio)),np.max(np.log10(ratio)))

    axes[0].hist(np.log10(gals_bulges[BT>0.7]['StellarDiskScaleLength']*1000*1.678),range=(-2,2),histtype='step',density=True,color='red',linestyle='--')
    axes[0].hist(np.log10(gals_bulges[BT<0.3]['StellarDiskScaleLength']*1000*1.678),range=(-2,2),histtype='step',density=True,color='blue',linestyle=':')
    axes[0].hist(np.log10(gals_bulges['StellarDiskScaleLength']*1000*1.678),range=(-2,2),histtype='step',density=True,color='black')
    axes[0].set_xlabel('log(Re)')
    axes[0].set_ylabel('Probability Density')
    axes[1].hist(np.log10(disk_mass[BT>0.7]),range=(6.5,11.5),histtype='step',density=True,color='red',linestyle='--')
    axes[1].hist(np.log10(disk_mass[BT<0.3]),range=(6.5,11.5),histtype='step',density=True,color='blue',linestyle=':')
    axes[1].hist(np.log10(disk_mass),range=(6.5,11.5),histtype='step',density=True,color='black')
    axes[1].set_xlabel('log(M disk)')
    axes[2].hist(np.log10(gals_bulges[BT>0.7]['VStellarDisk']),range=(1.4,2.8),histtype='step',density=True,color='red',linestyle='--')
    axes[2].hist(np.log10(gals_bulges[BT<0.3]['VStellarDisk']),range=(1.4,2.8),histtype='step',density=True,color='blue',linestyle=':')
    axes[2].hist(np.log10(gals_bulges['VStellarDisk']),range=(1.4,2.8),histtype='step',density=True,color='black')
    axes[2].set_xlabel('log(V)')
    axes[3].hist(np.log10(AM[BT>0.7]),range=(-6,2),histtype='step',density=True,color='red',linestyle='--')
    axes[3].hist(np.log10(AM[BT<0.3]),range=(-6,2),histtype='step',density=True,color='blue',linestyle=':')
    axes[3].hist(np.log10(AM),range=(-6,2),histtype='step',density=True,color='black')
    axes[3].set_xlabel('log(J)')
    axes[3].legend(['Bulge-dominated','Disk-dominated','All galaxies'])
    #axes[4].hist(np.log10(ratio[BT>0.7]),range=(-9.9,-9.25),histtype='step',log=True,color='red')#,density=True,color='black')
    #axes[4].hist(np.log10(ratio[BT<0.3]),range=(-9.9,-9.25),histtype='step',log=True,color='blue')#,density=True,color='black')
    #axes[4].hist(np.log10(ratio),range=(-9.9,-9.25),histtype='step',log=True)#,density=True,color='black')



  plt.tight_layout()
  
  plt.show()
