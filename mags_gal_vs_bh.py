import sys
sys.path.append('Yuxiang/')
import numpy as np
from dragons import meraxes
from _function import _function
from _meraxes_util import _quasar_luminosity_boot,\
_Lbol2MB, _Lbol2M1450
import matplotlib.pyplot as plt
import pandas as pd
import magcalc as mc

cosmo = {'omega_M_0' : 0.308, 
'omega_lambda_0' : 0.692, 
'omega_b_0' : 0.04839912, 
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}



def load_gals(filename,snap):
  gals = meraxes.io.read_gals(filename,snapshot=snap,h=cosmo['h'],quiet=True)
  gals = gals[(gals["GhostFlag"]==0)]
  gals = gals[(gals["StellarMass"]*1e10>1e6)]
  return gals


def calc_Magn(gals,fmeraxes):
  band='UV'
  eta = 0.06
  alpha_q = 1.57 
  alpha_q_op = 0.44
  sim_props = meraxes.io.read_input_params(fmeraxes,h=cosmo['h'],quiet=True)
  volume = sim_props['Volume']
  EddingtonRatio = sim_props['EddingtonRatio']
  quasar_open_angle = sim_props['quasar_open_angle']
  observed = 1-np.cos(np.deg2rad(quasar_open_angle)/2.)
  props = ("GhostFlag","BlackHoleMass","BlackHoleAccretedColdMass",'dt')
  bins = np.linspace(-30,-15,16)
  Nboot = 1

  bh = gals["BlackHoleMass"]*1e10
  bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
  delta_t = gals["dt"]

  Lbol  = _quasar_luminosity_boot(bh[bh_accreted_cold>0],bh_accreted_cold[bh_accreted_cold>0],delta_t[bh_accreted_cold>0],\
                                Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
  if band == 'B':
      Magn = _Lbol2MB(Lbol)
  elif band == 'UV':
      Magn = _Lbol2M1450(Lbol)
  return(Magn[0])


def load_gal_mags(filename,snapshot):
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust
  #return MUV


if __name__=='__main__':
  #Setup
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  filename='tuned_reion'
  fmeraxes=data_folder+filename+meraxes_loc
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,250:0}
  snap=78

  gals=load_gals(fmeraxes,snap)
  Magn=calc_Magn(gals,fmeraxes)
  Mgals=load_gal_mags(filename,snap)
  
  bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
  gals=gals[bh_accreted_cold>0]
  Mgals=Mgals[bh_accreted_cold>0]

  plt.plot(Mgals,Magn,'.')
  plt.show()

