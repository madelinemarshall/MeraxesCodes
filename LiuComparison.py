import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import pandas as pd
#import ContourPlot as cp

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
data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'
scalefactor=1
filename1='draft2_reion'

fig, axes = plt.subplots(2, 2)
ii=-1
j=0
for snapshot in [100,78,52,30]: #Find out what snapshot is z=10
  ii+=1
  if ii==2:
    j+=1
    ii=0
  #Load data:
  #Default parameters
  gals_bulges=meraxes.io.read_gals(data_folder+filename1+meraxes_loc,\
                                           snapshot=snapshot,\
                                           h=cosmo['h'],quiet=True)
  gals_bulges = gals_bulges[(gals_bulges["GhostFlag"]==0)]#remove ghosts
  mass=gals_bulges['StellarMass']*1e10-gals_bulges['BulgeStellarMass']*1e10
  axes[j,ii].plot(np.log10(mass[mass>0]),np.log10(gals_bulges['StellarDiskScaleLength'][mass>0]*1e3),'.')
#cp.contour_plot(np.log10(mass[mass>0]),np.log10(gals_bulges['StellarDiskScaleLength'][mass>0]*1e3),'$log(M)$','$log(R (kpc))$',[5,11.5],[-1.5,1])
  axes[j,ii].set_xlabel('log(Disk Stellar Mass)')
  axes[j,ii].set_ylabel('log(Stellar Disk Radius (kpc))')
  axes[j,ii].set_xlim([5,11.5])
  axes[j,ii].set_ylim([-1.5,1])
  axes[j,ii].legend()
plt.show()
