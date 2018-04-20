import sys
import numpy as np
from dragons import meraxes
from _function import _function
from _meraxes_util import _quasar_luminosity_boot,\
_Lbol2MB, _Lbol2M1450
import matplotlib.pyplot as plt

cosmo = {'omega_M_0' : 0.308, 
'omega_lambda_0' : 0.692, 
'omega_b_0' : 0.04839912, 
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}



def load_gals(filename,redshift):
  snap = meraxes.io.check_for_redshift(filename, redshift, tol=0.1)[0]
  gals = meraxes.io.read_gals(filename,props=props,snapshot=snap,h=cosmo['h'],quiet=True)
  return gals[(gals["GhostFlag"]==0)]



def calculateQLF(gals,fmeraxes,band,axes,**kwargs):
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
  Nboot = 10000
  fy = np.zeros([Nboot,len(bins)-1])

  bh = gals["BlackHoleMass"]*1e10
  bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
  delta_t = gals["dt"]

  Lbol  = _quasar_luminosity_boot(bh[bh_accreted_cold>0],bh_accreted_cold[bh_accreted_cold>0],delta_t[bh_accreted_cold>0],\
                                Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
  if band == 'B':
      Magn = _Lbol2MB(Lbol)
  elif band == 'UV':
      Magn = _Lbol2M1450(Lbol)

  for ii in range(Nboot):
      f,ed = _function(Magn[ii],volume,bins=bins,weights=np.zeros(len(Magn[ii]))+observed,return_edges=True)
      fy[ii] = f

  fx = ed[:-1]
  fy = np.sort(fy,axis=0)
  alpha=0.95
  y_mean = fy[int((1/2.0)*len(fy))]              
  y_err =  fy[int((1-alpha)/2.0*len(fy))] 
  y_err2 = fy[int((1+alpha)/2.0*len(fy))]     
  axes.plot(fx,y_mean,**kwargs)
  axes.set_yscale('log')
  axes.set_xlim([-28.4,-15.3])
  axes.set_ylim([10**-10,2*10**-4])
  axes.invert_xaxis()



if __name__=='__main__':
  fmeraxes = '/home/mmarshal/data_dragons/bulges_update1102_full/output/meraxes.hdf5'

  redshift = float(sys.argv[1])
  band = sys.argv[2]
  gals=load_gals(fmeraxes,redshift)
  fig,axes=plt.subplot(1,1)
  calculateQLF(gals,band,axes)
  plt.show()
