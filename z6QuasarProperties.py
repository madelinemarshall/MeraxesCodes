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
from _meraxes_util import _quasar_luminosity_boot,\
_Lbol2MB, _Lbol2M1450

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


def load_data(filename,snapshot):
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[(gals['StellarMass']*1e10>1e6)]
  return gals


def load_UVmags(filename,snapshot):
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def load_JWSTmags(filename,snapshot):
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  mags= pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_078_JWSTfilters.hdf5')
  #print(list(mags.columns.values))
  mags_dust={}
  filters=[ 'Y105','H160',"JWST_F070W","JWST_F090W","JWST_F115W","JWST_F150W","JWST_F200W","JWST_F277W"]#,"JWST_F356W","JWST_F444W"]
  central_wavelength=[10500,16000,7000,9000,11500,15000,20000,27700]#,35600,44400]
  i=-1
  for filt in filters:
    i+=1
    mags_dust[filt] = np.array(mags[filt]) + mc.reddening(central_wavelength[i]/(1.+redshift[snapshot]), MUV, z = redshift[snapshot])
  return mags_dust


def calc_Magn(gals,filename):
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  band='UV'
  eta = 0.06
  alpha_q = 1.57 
  alpha_q_op = 0.44
  sim_props = meraxes.io.read_input_params(data_folder+filename+meraxes_loc,h=cosmo['h'],quiet=True)
  volume = sim_props['Volume']
  EddingtonRatio = sim_props['EddingtonRatio']
  quasar_open_angle = sim_props['quasar_open_angle']
  observed = 1-np.cos(np.deg2rad(quasar_open_angle)/2.)
  props = ("GhostFlag","BlackHoleMass","BlackHoleAccretedColdMass",'dt')
  bins = np.linspace(-30,-15,16)
  Nboot = 1000

  bh = gals["BlackHoleMass"]*1e10
  bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
  delta_t = gals["dt"]

  Lbol  = _quasar_luminosity_boot(bh[bh_accreted_cold>0],bh_accreted_cold[bh_accreted_cold>0],delta_t[bh_accreted_cold>0],\
                                Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
  if band == 'B':
      Magn = _Lbol2MB(Lbol)
  elif band == 'UV':
      Magn = _Lbol2M1450(Lbol)
  Magn[np.isinf(Magn)]=np.nan
  return(np.nanmedian(Magn,axis=0))


if __name__=="__main__":
  #Setup
  filename='tuned_reion'
  redshift={78:6}
  snapshot=78  

  gals=load_data(filename,snapshot)
  Mgal=load_UVmags(filename,snapshot)
  JWST_mags=load_JWSTmags(filename,snapshot)
  Magn=calc_Magn(gals[gals['BlackHoleMass']*1e10>5e6],filename)
  Mgal=np.array(Mgal[gals['BlackHoleMass']*1e10>5e6])
  filters=[ 'Y105','H160',"JWST_F070W","JWST_F090W","JWST_F115W","JWST_F150W","JWST_F200W","JWST_F277W"]#,"JWST_F356W","JWST_F444W"]
  central_wavelength=[10500,16000,7000,9000,11500,15000,20000,27700]#,35600,44400]
  mag=np.zeros((np.size(gals['ID']),len(filters)))
  for ii in range(0,np.size(gals['ID'])):
    if gals[ii]['BlackHoleMass']*1e10>5e6:
      i=-1
      for ff in filters:
        i+=1
        mag[ii,i]=JWST_mags[ff][ii]
  mag=mag[gals['BlackHoleMass']*1e10>5e6]
  gals=gals[gals['BlackHoleMass']*1e10>5e6]
  Magn[np.isnan(Magn)]=0
  print(filters,np.median(mag,axis=0))
  
  #plt.errorbar(central_wavelength,np.median(mag,axis=0)) 
  #plt.gca().invert_yaxis()
  #plt.show() 
  

  fig,ax=plt.subplots(1,4, sharey=True,gridspec_kw = {'wspace':0, 'hspace':0})

  ax[0].hist(np.log10(gals['BlackHoleMass']*1e10),color=color[0],edgecolor=color[0],alpha=0.4,hatch='\\')
  ax[0].set_xlabel(r'$\log M_{BH}/M_\odot$')
  ax[1].hist(np.log10(gals['StellarMass']*1e10),range=(5.5,12),color=color[0],edgecolor=color[0],alpha=0.4,hatch='\\')
  ax[1].hist(np.log10(gals['ColdGas'][gals['ColdGas']>0]*1e10),range=(5.5,12),color=color[1],edgecolor=color[1],alpha=0.4,hatch='/')
  ax[1].set_xlabel(r'$\log M/M_\odot$')
  ax[1].legend(['Stellar','Gas'],loc='upper left',fontsize='small')
  ax[2].hist(gals['StellarDiskScaleLength']*1e3/5.843,color=color[0],edgecolor=color[0],alpha=0.4,range=(0,0.4),hatch='\\') #5.843 kpc/" at z=6
  ax[2].hist(gals['GasDiskScaleLength']*1e3/5.843,color=color[1],edgecolor=color[1],alpha=0.4,range=(0,0.4),hatch='/') #5.843 kpc/" at z=6
  ax[2].legend(['Stellar','Gas'],loc='upper left',fontsize='small')
  ax[2].set_xlabel(r'Disk Scale Length (arcsec)')
  ax[3].hist(gals['BulgeStellarMass']/gals['StellarMass'],color=color[0],edgecolor=color[0],alpha=0.4,hatch='\\')
  ax[3].set_xlabel(r'B/T')
  ax[0].set_ylabel(r'Number of Quasar Hosts')
  plt.tight_layout() 
  plt.savefig('/home/mmarshal/results/plots/z6Hosts_props.pdf',format='pdf')
  plt.show()
  
  fig,ax=plt.subplots(1,3, sharey=True,gridspec_kw = {'wspace':0, 'hspace':0})
  print(np.shape(Mgal))
  ax[0].hist(Mgal,color=color[0],edgecolor=color[0],alpha=0.4,hatch='/')
  ax[0].set_xlabel(r'$M_{UV,host}$ (includes dust extinction)')
  ax[1].hist(Magn,color=color[0],edgecolor=color[0],alpha=0.4,hatch='/')
  ax[1].set_xlabel(r'$M_{UV,quasar}$')
  ax[2].hist(Magn-Mgal,color=color[0],edgecolor=color[0],alpha=0.4,hatch='/')
  ax[2].set_xlabel(r'$M_{UV,quasar}-M_{UV,host}$')
  ax[0].set_ylabel(r'Number of Quasar Hosts')
  plt.tight_layout() 
  plt.savefig('/home/mmarshal/results/plots/z6Hosts_mags.pdf',format='pdf')
  plt.show()


