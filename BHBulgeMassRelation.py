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
#filename='bulges_update1102_full'
filename='bulges_nodiskinstability'

#Haring & Rix
MBH=10**(8.2+1.12*np.log10(np.array([10**8,10**12])/10**11))
#MBH_low=10**(8.1+1.06*np.log10(np.array([10**8,10**12])/10**11))
#MBH_high=10**(8.3+1.18*np.log10(np.array([10**8,10**12])/10**11))

fig, axes = plt.subplots(2, 3)
ii=-1
j=0
for snapshot in [63,78,100,116,134,158]:
  ii+=1
  if ii==3:
    j+=1
    ii=0
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,props=['StellarMass','BulgeStellarMass','BlackHoleMass'],\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']*1e10>1e7]
  gals=gals[gals['BlackHoleMass']>0]
  gals=gals[gals['BulgeStellarMass']>0]

  axes[j,ii].plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
#  cp.contour_plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),xlims=[8,12],ylims=[5,10],axes=axes[j,ii])
  axes[j,ii].plot([8,12],np.log10(MBH),color='orange')
  #axes[j,ii].plot([8,12],np.log10(MBH_low),'--',color='orange')
  #axes[j,ii].plot([8,12],np.log10(MBH_high),'--',color='orange')
  axes[j,ii].set_xlim([8,12])
  axes[j,ii].set_ylim([5,10])
  axes[j,ii].set_xlabel('log(Bulge Mass)')
  axes[j,ii].set_ylabel('log(BH Mass)')
  if (snapshot==63):
    axes[j,ii].set_title('z = 7')
  elif (snapshot==78):
    axes[j,ii].set_title('z = 6')
    Wang_BH=np.log10([2.8e9,2.1e9,8.6e8,1.7e8])
    Wang_dyn=np.log10([9.6e10,12.5e10,7.2e10,1.3e10])
    axes[j,ii].plot(Wang_dyn,Wang_BH,'ro')
  elif snapshot==100:
    axes[j,ii].set_title('z = 5')
  elif snapshot==116:
    axes[j,ii].set_title('z = 4')
  elif snapshot==134:
    axes[j,ii].set_title('z = 3')
  elif snapshot==158:
    axes[j,ii].set_title('z = 2')

plt.show()
