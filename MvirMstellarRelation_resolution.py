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
  #gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals

if __name__=='__main__':
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  #data_folder='/lustre/projects/p113_astro/yqin/dragons_results/'
  #meraxes_loc='/meraxes.hdf5'

  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,250:0}
  snapshots=[78,100,116,134,158]
  prop='StellarMass'
  color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'pink',250:'black'}
 
  #filename='bulges_correctBHMF'
  filename='paper1'
  #filename='meraxes_on_tiamat_newtrees/Q_100_000_S_006_050'
  #filename2='meraxes_on_tiamat_newtrees/Q_100_000_S_006_050_Re'
  resolution_lim=np.zeros(len(snapshots))

  ii=0
  for snap in snapshots:
    #New model
    gals=load_data(filename,snap)
    df=pd.DataFrame(gals)
    grouped=df.groupby('CentralGal')
    grouped_sum=grouped.sum()
    FOFMvir=np.array(grouped_sum['Mvir']*1e10)
    FOFStellar=np.array(grouped_sum['StellarMass']*1e10)
    logFOFMvir=np.log10(FOFMvir)
    logFOFStellar=np.log10(FOFStellar) 
    FOFFrac=FOFStellar###/FOFMvir
    #logFOFMvir=np.log10(gals['Mvir']*1e10)
    #logFOFStellar=np.log10(gals['StellarMass']*1e10) 
    #FOFFrac=gals['StellarMass']*1e10###/FOFMvir
   
    bin_width=0.2
    min_mass=6#10
    max_mass=np.max(logFOFMvir)
    n_bins=np.int((max_mass-min_mass)/bin_width)
    pctile_frac=np.zeros((n_bins,3))
    middle_mvir=np.zeros(n_bins)
    for nn in range(0,n_bins):
      selection=(logFOFMvir>min_mass+(nn*bin_width))&(logFOFMvir<min_mass+(nn+1)*bin_width)
      if len(selection[selection>0])>0:
        pctile_frac[nn,:]=np.percentile(FOFFrac[selection],[5,50,95])
        middle_mvir[nn]=min_mass+(nn+0.5)*bin_width
    plt.plot(middle_mvir,pctile_frac[:,1],label='$z={}$'.format(redshift[snap]),color=color[snap])
    plt.fill_between(middle_mvir,pctile_frac[:,0],pctile_frac[:,2],alpha=0.8,color=color[snap])
    resolution_lim[ii]=pctile_frac[abs(middle_mvir-10.27875)<0.1][:,2] 

    ii+=1

  print(np.log10(resolution_lim))

  #res_poly=np.polyfit([7,6,5,4,3],np.log10(resolution_lim),1) #fit logM vs (1+z)
  #poly = np.poly1d(res_poly)
  #print('Redshift 0 = {}'.format(poly(1))) 

  plt.plot([6.5765+2,6.5765+2],[1e6,1e11],'k-')
  plt.plot([8.27875+2,8.27875+2],[1e6,1e11],'k-')
  plt.grid(which='major')
  plt.legend()
  plt.xlabel(r'$\log(M_{vir})$')
  plt.ylabel(r'$M_\ast$')
  ##plt.ylabel(r'$M_\ast/M_{vir}$')
  plt.yscale('log')
  plt.show()
