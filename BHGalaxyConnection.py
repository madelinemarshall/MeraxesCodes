import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import pandas as pd
import ContourPlot as cp
 
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
filename='bulges_update1102_full'


fig1, axes1 = plt.subplots(2, 3,num=1)
fig2, axes2 = plt.subplots(2, 3,num=2)
fig3, axes3 = plt.subplots(2, 3,num=3)
ii=-1
j=0
for snapshot in [63,78,100,116,134,158]:
  ii+=1
  if ii==3:
    j+=1
    ii=0
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['Type']==0] ##Select centrals only as in MP17
  gals=gals[gals['StellarMass']*1e10>1e10]
  gals=gals[gals['StellarMass']*1e10<1e13]

  ###Fig 7 
  Mstell=np.array([10,13])
  MBH=1.53*Mstell-8.85

  cp.contour_plot(np.log10(gals['StellarMass'][(gals['StellarMass']>0)]*1e10),np.log10(gals['BlackHoleMass'][(gals['StellarMass']>0)]*1e10),'log(Stellar Mass)','log(BH Mass)',[10,13],[5.2,10.5],axes1[j,ii])
#  axes1[j,ii].plot(np.log10(gals['StellarMass'][(gals['StellarMass']>0)]*1e10),np.log10(gals['BlackHoleMass'][(gals['StellarMass']>0)]*1e10),'.')
  axes1[j,ii].plot(Mstell,MBH)
#  axes1[j,ii].set_xlabel('log(Stellar Mass)')
#  axes1[j,ii].set_ylabel('log(BH Mass)')
#  axes1[j,ii].set_xlim([10,13])
#  axes1[j,ii].set_ylim([5.2,10.5])
  if (snapshot==63):
    axes1[j,ii].set_title('z = 7')
  elif (snapshot==78):
    axes1[j,ii].set_title('z = 6')
  elif snapshot==100:
    axes1[j,ii].set_title('z = 5')
  elif snapshot==116:
    axes1[j,ii].set_title('z = 4')
  elif snapshot==134:
    axes1[j,ii].set_title('z = 3')
  elif snapshot==158:
    axes1[j,ii].set_title('z = 2')
    Jahnke_BH=[8.52,8.07,8.25,8.05,8.43,8.35,8.40]
    Jahnke_stellar=[11.18,11.00,11.44,10.90,10.88,11.64,10.93]
    axes1[j,ii].plot(Jahnke_stellar,Jahnke_BH,'ro')
  ##Fig 8
  MDM=np.array([11,14.5])
  MBH=1.55*MDM-11.26

  cp.contour_plot(np.log10(gals['Mvir'][(gals['BlackHoleMass']>0)]*1e10),np.log10(gals['BlackHoleMass'][(gals['BlackHoleMass']>0)]*1e10),'log(DM Halo Mass)','log(BH Mass)',[11,14.5],[5.2,10.5],axes2[j,ii])
#  axes2[j,ii].plot(np.log10(gals['Mvir'][(gals['BlackHoleMass']>0)]*1e10),np.log10(gals['BlackHoleMass'][(gals['BlackHoleMass']>0)]*1e10),'.')
  axes2[j,ii].plot(MDM,MBH)
#  axes2[j,ii].set_xlabel('log(DM Halo Mass)')
#  axes2[j,ii].set_ylabel('log(BH Mass)')
#  axes2[j,ii].set_xlim([11,14.5])
#  axes2[j,ii].set_ylim([5.2,10.5])
  if (snapshot==63):
    axes2[j,ii].set_title('z = 7')
  elif (snapshot==78):
    axes2[j,ii].set_title('z = 6')
  elif snapshot==100:
    axes2[j,ii].set_title('z = 5')
  elif snapshot==116:
    axes2[j,ii].set_title('z = 4')
  elif snapshot==134:
    axes2[j,ii].set_title('z = 3')
  elif snapshot==158:
    axes2[j,ii].set_title('z = 2')

  ##Fig 9
  MDM=np.array([11,14.5])
  MBH=1.62*MDM-12.07

  cp.contour_plot(np.log10(gals['Mvir'][(gals['BlackHoleMass']>0)]*1e10+gals['StellarMass'][(gals['BlackHoleMass']>0)]*1e10+\
  gals['ColdGas'][(gals['BlackHoleMass']>0)]*1e10+gals['HotGas'][(gals['BlackHoleMass']>0)]*1e10+gals['BlackHoleMass'][(gals['BlackHoleMass']>0)]*1e10),\
  np.log10(gals['BlackHoleMass'][(gals['BlackHoleMass']>0)]*1e10),'log(Total Mass)','log(BH Mass)',[11,14.5],[5.2,10.5],axes3[j,ii])

#  axes3[j,ii].plot(np.log10(gals['Mvir'][(gals['BlackHoleMass']>0)]*1e10+gals['StellarMass'][(gals['BlackHoleMass']>0)]*1e10+\
#  gals['ColdGas'][(gals['BlackHoleMass']>0)]*1e10+gals['HotGas'][(gals['BlackHoleMass']>0)]*1e10+gals['BlackHoleMass'][(gals['BlackHoleMass']>0)]*1e10),\
#  np.log10(gals['BlackHoleMass'][(gals['BlackHoleMass']>0)]*1e10),'.')
  axes3[j,ii].plot(MDM,MBH)
#  axes3[j,ii].set_xlabel('log(Total Mass)')
#  axes3[j,ii].set_ylabel('log(BH Mass)')
#  axes3[j,ii].set_xlim([11,14.5])
#  axes3[j,ii].set_ylim([5.2,10.5])
  if (snapshot==63):
    axes3[j,ii].set_title('z = 7')
  elif (snapshot==78):
    axes3[j,ii].set_title('z = 6')
  elif snapshot==100:
    axes3[j,ii].set_title('z = 5')
  elif snapshot==116:
    axes3[j,ii].set_title('z = 4')
  elif snapshot==134:
    axes3[j,ii].set_title('z = 3')
  elif snapshot==158:
    axes3[j,ii].set_title('z = 2')

plt.tight_layout()
plt.show()
