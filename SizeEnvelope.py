import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import pandas as pd

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

redshift=[8,7,6,5,4,3,2]

GasDisk=np.zeros(len(redshift))
StellarDisk=np.zeros(len(redshift))

ii=0
for snapshot in [52,63,78,100,116,134,158]:
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']*1e10>1e9]
  gals=gals[gals['ColdGas']*1e10>10**(9.5)]
  #gals=gals[gals['StellarMass']>0]
  GasDisk[ii]=np.nanmedian(gals['ColdGas'])
  StellarDisk[ii]=np.nanmedian(gals['StellarMass'])
  ii+=1


fig,(ax1)=plt.subplots(1,1)

ax1.plot(redshift,GasDisk,'b')
ax1.set_ylabel('Median Disk Scale Length')
ax1.set_yscale('log')

ax1.plot(redshift,StellarDisk,'r')
ax1.set_yscale('log')
ax1.set_xlabel('Redshift')
plt.legend(['Gas','Stellar'])
plt.gca().invert_xaxis()


filename='bulges_tiamat125_crotonSF'
redshift=[6,5,4,3,2,1.5,1,0.5]

GasDisk=np.zeros(len(redshift))
StellarDisk=np.zeros(len(redshift))

ii=0
for snapshot in [78,100,116,134,158,173,192,213]:
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']*1e10>1e9]
  gals=gals[gals['ColdGas']*1e10>10**(9.5)]
  #gals=gals[gals['StellarMass']>0]
  GasDisk[ii]=np.nanmedian(gals['ColdGas'])
  StellarDisk[ii]=np.nanmedian(gals['StellarMass'])
  ii+=1

ax1.plot(redshift,GasDisk,'b--')
ax1.plot(redshift,StellarDisk,'r--')
plt.show()
