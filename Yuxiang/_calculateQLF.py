import sys
import numpy as np
from dragons import meraxes
from _function import _function
from _meraxes_util import _quasar_luminosity_boot, _AGN_luminosity_boot,\
_Lbol2MB, _Lbol2M1450, _ksoftX, _khardX
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



def calculateQLF(gals,fmeraxes,band,axes,eta=None,**kwargs):
  alpha_q = 1.57 
  alpha_q_op = 0.44
  sim_props = meraxes.io.read_input_params(fmeraxes,h=cosmo['h'],quiet=True)
  volume = sim_props['Volume']
  EddingtonRatio = sim_props['EddingtonRatio']
  
  quasar_open_angle = sim_props['quasar_open_angle']
  observed = 1-np.cos(np.deg2rad(quasar_open_angle)/2.)
  if eta==None:
    if 'RadioAccretionEff' in sim_props.keys():
      eta=sim_props['RadioAccretionEff']
    else:
      eta = 0.06

  props = ("GhostFlag","BlackHoleMass","BlackHoleAccretedColdMass",'BlackHoleAccretedHotMass','dt')
  bins = np.linspace(-30,-15,16)
  Nboot = 100#0
  fy = np.zeros([Nboot,len(bins)-1,2])

  bh = gals["BlackHoleMass"]*1e10
  bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
  bh_accreted_hot = gals["BlackHoleAccretedHotMass"]*1e10
  delta_t = gals["dt"]

  #Lbol  = _quasar_luminosity_boot(bh[bh_accreted_cold>0],bh_accreted_cold[bh_accreted_cold>0],delta_t[bh_accreted_cold>0],\
  #                              Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
  Lbol  = _AGN_luminosity_boot(bh[(bh_accreted_cold+bh_accreted_hot>0)],bh_accreted_hot[(bh_accreted_cold+bh_accreted_hot>0)],bh_accreted_cold[(bh_accreted_cold+bh_accreted_hot>0)],delta_t[(bh_accreted_cold+bh_accreted_hot>0)],\
                                Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
  if band == 'B':
      Magn = _Lbol2MB(Lbol)
      axes.invert_xaxis()
  elif band == 'UV':
      Magn = _Lbol2M1450(Lbol)
      axes.invert_xaxis()
  elif band == 'softX': #0.5-2keV
      Magn = np.log10(Lbol/_ksoftX(Lbol)) ##logL, not mag
      bins = np.linspace(8,14,16)
      observed=1
  elif band == 'hardX': #2-10keV
      Magn =  np.log10(Lbol/_khardX(Lbol)) ##logL, not mag
      bins = np.linspace(8,14,16)
      observed=1
  elif band == 'bol': #Bolometric
      Magn =  np.log10(Lbol) ##logL, not mag
      bins = np.linspace(8,14,16)
      observed=1

  for ii in range(Nboot):
      f,ed = _function(Magn[ii],volume,bins=bins,return_edges=True,weights=np.zeros(len(Magn[ii]))+observed)
      fy[ii] = f

  fx = ed[:-1]
  fy = np.sort(fy,axis=0)
  alpha=0.95
  y_mean = np.log10(np.array(fy[int((1/2.0)*len(fy))])[:,1])
  y_err =  np.array(fy[int((1-alpha)/2.0*len(fy))])[:,1] 
  y_err2 =  np.array(fy[int((1+alpha)/2.0*len(fy))])[:,1] 
  axes.plot(fx,y_mean,**kwargs)
  axes.fill_between(fx, np.log10(y_err), np.log10(y_err2),color=kwargs['color'],alpha=0.3,label='__nolabel__')
  #axes.set_yscale('log')


def calculateQLF_evolvingEdd(gals,fmeraxes,redshift,band,axes,eta=None,**kwargs):
  alpha_q = 1.57 
  alpha_q_op = 0.44
  sim_props = meraxes.io.read_input_params(fmeraxes,h=cosmo['h'],quiet=True)
  volume = sim_props['Volume']
  EddingtonRatio = 0.1#1.9374e-3*(1+redshift)**3+8.063e-3
  
  quasar_open_angle = sim_props['quasar_open_angle']
  observed = 1-np.cos(np.deg2rad(quasar_open_angle)/2.)
  if eta==None:
    if 'RadioAccretionEff' in sim_props.keys():
      eta=sim_props['RadioAccretionEff']
    else:
      eta = 0.06

  props = ("GhostFlag","BlackHoleMass","BlackHoleAccretedColdMass",'BlackHoleAccretedHotMass','dt')
  bins = np.linspace(-30,-15,16)
  Nboot = 100#0
  fy = np.zeros([Nboot,len(bins)-1,2])

  bh = gals["BlackHoleMass"]*1e10
  bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
  bh_accreted_hot = gals["BlackHoleAccretedHotMass"]*1e10
  delta_t = gals["dt"]

  #Lbol  = _quasar_luminosity_boot(bh[bh_accreted_cold>0],bh_accreted_cold[bh_accreted_cold>0],delta_t[bh_accreted_cold>0],\
  #                              Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
  Lbol  = _AGN_luminosity_boot(bh[(bh_accreted_cold+bh_accreted_hot>0)],bh_accreted_hot[(bh_accreted_cold+bh_accreted_hot>0)],bh_accreted_cold[(bh_accreted_cold+bh_accreted_hot>0)],delta_t[(bh_accreted_cold+bh_accreted_hot>0)],\
                                Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
  if band == 'B':
      Magn = _Lbol2MB(Lbol)
      axes.invert_xaxis()
  elif band == 'UV':
      Magn = _Lbol2M1450(Lbol)
      axes.invert_xaxis()
  elif band == 'softX': #0.5-2keV
      Magn = np.log10(Lbol/_ksoftX(Lbol)) ##logL, not mag
      bins = np.linspace(8,14,16)
      observed=1
  elif band == 'hardX': #2-10keV
      Magn =  np.log10(Lbol/_khardX(Lbol)) ##logL, not mag
      bins = np.linspace(8,14,16)
      observed=1
  elif band == 'bol': #Bolometric
      Magn =  np.log10(Lbol) ##logL, not mag
      bins = np.linspace(8,14,16)
      observed=1

  for ii in range(Nboot):
      f,ed = _function(Magn[ii],volume,bins=bins,return_edges=True,weights=np.zeros(len(Magn[ii]))+observed)
      fy[ii] = f

  fx = ed[:-1]
  fy = np.sort(fy,axis=0)
  alpha=0.95
  y_mean = np.log10(np.array(fy[int((1/2.0)*len(fy))])[:,1])
  y_err =  np.array(fy[int((1-alpha)/2.0*len(fy))])[:,1] 
  y_err2 =  np.array(fy[int((1+alpha)/2.0*len(fy))])[:,1] 
  axes.plot(fx,y_mean,**kwargs)
  axes.fill_between(fx, np.log10(y_err), np.log10(y_err2),color=kwargs['color'],alpha=0.3,label='__nolabel__')
  #axes.set_yscale('log')

if __name__=='__main__':
  fmeraxes = '/home/mmarshal/data_dragons/bulges_update1102_full/output/meraxes.hdf5'

  redshift = float(sys.argv[1])
  band = sys.argv[2]
  gals=load_gals(fmeraxes,redshift)
  fig,axes=plt.subplot(1,1)
  calculateQLF(gals,band,axes)
  plt.show()
