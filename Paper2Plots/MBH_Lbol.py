import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/simulation_codes/Yuxiang/')
sys.path.append('/home/mmarshal/simulation_codes')
from _load_data import load_data
from _function import _function
from _calculateQLF import calculateQLF
from _meraxes_util import _quasar_luminosity_boot, _AGN_luminosity_boot,\
_Lbol2MB, _Lbol2M1450, _ksoftX, _khardX

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,3.2)
#matplotlib.rcParams['lines.linewidth'] = 2.5
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

markers        = ['s','v','^','<','>','p','*','D','8']*4
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4

cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
}


def plot_hist2d(xdata,ydata,axes,cmax=None,cbar=False):
  xlims=[np.nanmin(xdata),np.nanmax(xdata)]
  ylims=[np.nanmin(ydata),np.nanmax(ydata)]
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, vmin=1, vmax=cmax, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  if True:#cbar:
    cb=plt.colorbar(img, cax=ax[2],use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))

 
def calc_mag(fmeraxes,gals,band):
  alpha_q = 1.57 
  alpha_q_op = 0.44
  sim_props = meraxes.io.read_input_params(fmeraxes,h=cosmo['h'],quiet=True)
  volume = sim_props['Volume']
  EddingtonRatio = sim_props['EddingtonRatio']
  
  quasar_open_angle = sim_props['quasar_open_angle']
  observed = 1-np.cos(np.deg2rad(quasar_open_angle)/2.)
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
  Lbol=np.ma.masked_equal(Lbol,0).mean(axis=0) 
  #Lbol[Lbol==0]=np.nan
  #Lbol=np.nanmean(Lbol,axis=0)
  #print(test,Lbol)
  #print(np.shape(test),np.shape(Lbol))
 
  if band == 'B':
      Magn = _Lbol2MB(Lbol)
  elif band == 'UV':
      Magn = _Lbol2M1450(Lbol)
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

  return Magn,bh[(bh_accreted_cold+bh_accreted_hot>0)]

if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,213:0.55,250:0}
  color={52:'C0',63:'C1',78:'C2',100:'C3',116:'C4',134:'aqua',158:'pink',213:'k',250:'k'}
  meraxes_loc='/output/meraxes.hdf5'

  filename_125='paper2_T125'
  fmeraxes=data_folder+filename_125+meraxes_loc

  fig,ax = plt.subplots(1,3,gridspec_kw = {'wspace':0, 'hspace':0,'width_ratios':[4,4,0.4]})#,sharex=True,sharey=True)
  ii=0
  for snap in [158,250]:
    gals=load_data(filename_125,snap,'All',centrals=True)
    mag,MBH=calc_mag(fmeraxes,gals,'UV')   
    plot_hist2d(mag,np.log10(MBH),ax[ii])
    ax[ii].text(-15,9,'N accreting={}'.format(len(mag)))
    ax[ii].set_xlabel(r'$\log(M_{BH}/M_\odot)$')
    ax[ii].set_ylabel(r'$M_{UV}$')
    ii+=1
   

  ax[0].invert_xaxis()
  ax[1].invert_xaxis()

  ax[0].set_title('z=2')
  ax[1].set_title('z=3')

  ax[2].axis('off')
  #plt.tight_layout()
  lgd=ax[0].legend(fontsize='small',loc='upper center', bbox_to_anchor=(2.94, 0.8))  
  ax[0].set_ylabel(r'$\log(M_{\rm{BH}}/M_\odot)$') 
  ax[1].set_yticks([])
  plt.subplots_adjust(bottom=0.2,left=0.08,right=0.95,top=0.95)
  plt.show()
