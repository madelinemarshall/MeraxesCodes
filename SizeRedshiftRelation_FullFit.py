##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import ContourPlot as cp
import pandas as pd
import magcalc as mc
from scipy.optimize import curve_fit
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (4.5,8)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4


def load_gals(filename,snapshot):
  #Setup
  cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  }
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,\
#props=['GhostFlag','StellarMass','StellarDiskScaleLength','BulgeStellarMass','GasDiskScaleLength'],\
                                          h=cosmo['h'],quiet=True)
  return gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-3)]


def load_mags(filename,snapshot):
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def func(z,a,b):
  return a*((1+z)/(1+7))**-b
  #log r/r7=-b log(1+z/1+7)


def fit_equation(zz,rad,axes):
  popt,pcov = curve_fit(func,zz,rad)
  print("popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))
  axes.plot(zz,func(zz,*popt),'--',color='green')


if __name__=="__main__":
  filename='draft2_reion'
  filename_def='dragons10'
  disk_length={filename:'StellarDiskScaleLength',filename_def:'DiskScaleLength'}
  label={filename:"M18",filename_def:"Q17"}
  linestyle={filename:"b-",filename_def:"b--"}

  #redshift={37:10,43:9,52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,250:0}
  redshift={37:10,43:9,52:8,63:7,78:6,100:5}
  snapshots=np.flip([37,43,52,63,78,100],0)
  fig,axes=plt.subplots(3,gridspec_kw = {'hspace':0})
  for fname in [filename]:
    ii=0
    for snapshot in snapshots:
      gals=load_gals(fname,snapshot)
      mag=load_mags(fname,snapshot)
      upper=gals[(mag<-19.7)&(mag>-21)&(gals[disk_length[fname]]>0)][disk_length[fname]]*1000*1.678
      lower=gals[(mag>-19.7)&(mag<-18.7)&(gals[disk_length[fname]]>0)][disk_length[fname]]*1000*1.678
      z_up=[redshift[snapshot]]*len(upper)
      z_low=[redshift[snapshot]]*len(lower)

      if ii>0:
        upper_radii+=upper.tolist()
        lower_radii+=lower.tolist()
        upper_z+=z_up
        lower_z+=z_low
      else:
        upper_radii=upper.tolist()
        lower_radii=lower.tolist()
        upper_z=z_up
        lower_z=z_low
      ii+=1

    axes[0].plot(upper_z,upper_radii,'.',label=label[fname],color='black')
    axes[1].plot(lower_z,lower_radii,'.',label=label[fname],color='black')

    #print(fname)
    print('Upper:') 
    fit_equation(np.array(upper_z),np.array(upper_radii),axes[0])
    #print('Lower:') 
    #fit_equation(np.flip(np.array(list(redshift.values())),0),mean_r_lower,axes[1])
 
  axes[2].axis('off')
  axes[0].set_ylabel('$R_e$ (kpc)')
  axes[1].set_xlabel('Redshift')
  axes[1].set_ylabel('$R_e$ (kpc)')
  axes[0].set_xticklabels([])
  axes[0].legend(ncol=2,loc=(-0.15,-1.9))
  axes[0].set_ylim(0,1.6)
  axes[1].set_ylim(0,1.6)
  axes[0].set_xlim(4.6,10.3)
  axes[1].set_xlim(4.6,10.3)
  axes[0].text(8.1, 1.4, r'$(0.3-1)L^\ast_{z=3}$',weight='bold',size='large')
  axes[1].text(8.1, 1.4, r'$(0.12-0.3)L^\ast_{z=3}$',weight='bold',size='large')
  plt.show()
