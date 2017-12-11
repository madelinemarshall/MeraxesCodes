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

MMBH=np.zeros(len(redshift))
MMBulge=np.zeros(len(redshift))
MMBH_2=np.zeros(len(redshift))
MMBulge_2=np.zeros(len(redshift))
MMBH_3=np.zeros(len(redshift))
MMBulge_3=np.zeros(len(redshift))
MMBHBulge=np.zeros(len(redshift))
MMBHBulge_2=np.zeros(len(redshift))
MMBHBulge_3=np.zeros(len(redshift))
MedBHBulge=np.zeros(len(redshift))
MedBHBulge_2=np.zeros(len(redshift))
MedBHBulge_3=np.zeros(len(redshift))
MinBHBulge=np.zeros(len(redshift))
MinBHBulge_2=np.zeros(len(redshift))
MinBHBulge_3=np.zeros(len(redshift))

ii=0
for snapshot in [52,63,78,100,116,134,158]:
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['BulgeStellarMass']*1e10>1e4]
  gals=gals[gals['BlackHoleMass']*1e10>1e6]
  MMBH[ii]=np.nanmax(gals['BlackHoleMass']*1e10)
  MMBulge[ii]=np.nanmax(gals['BulgeStellarMass']*1e10)
  MMBHBulge[ii]=np.nanmax(gals['BlackHoleMass']/gals['BulgeStellarMass'])
  MedBHBulge[ii]=np.nanmedian(gals['BlackHoleMass']/gals['BulgeStellarMass'])
  MinBHBulge[ii]=np.min(gals['BlackHoleMass']/gals['BulgeStellarMass'])
  
  gals=gals[gals['BlackHoleMass']*1e10>1e6]
  gals=gals[gals['BulgeStellarMass']*1e10>1e8]
  if np.size(gals)>0:
    MMBH_2[ii]=np.nanmax(gals['BlackHoleMass']*1e10)
    MMBulge_2[ii]=np.nanmax(gals['BulgeStellarMass']*1e10)
    MMBHBulge_2[ii]=np.nanmax(gals['BlackHoleMass']/gals['BulgeStellarMass'])
    MedBHBulge_2[ii]=np.nanmedian(gals['BlackHoleMass']/gals['BulgeStellarMass'])
    MinBHBulge_2[ii]=np.min(gals['BlackHoleMass']/gals['BulgeStellarMass'])

  gals=gals[gals['BlackHoleMass']*1e10>1e7]
  gals=gals[gals['BulgeStellarMass']*1e10>1e8]
  if np.size(gals)>0: 
    MMBH_3[ii]=np.nanmax(gals['BlackHoleMass']*1e10)
    MMBulge_3[ii]=np.nanmax(gals['BulgeStellarMass']*1e10)
    MMBHBulge_3[ii]=np.nanmax(gals['BlackHoleMass']/gals['BulgeStellarMass'])
    MedBHBulge_3[ii]=np.nanmedian(gals['BlackHoleMass']/gals['BulgeStellarMass'])
    MinBHBulge_3[ii]=np.min(gals['BlackHoleMass']/gals['BulgeStellarMass'])

  ii+=1


fig,(ax1,ax2,ax3)=plt.subplots(3,1,sharex=True)

ax1.plot(redshift,MMBH,'b')
ax1.plot(redshift,MMBH_2,'r--')
ax1.plot(redshift,MMBH_3,'g:')
ax1.set_ylabel('Max BH Mass')
ax1.set_yscale('log')
ax1.legend(['BH > 1e4, Bulge > 1e6','BH > 1e6, Bulge > 1e6','BH > 1e7, Bulge > 1e8'])

ax2.plot(redshift,MMBulge,'b')
ax2.plot(redshift,MMBulge_2,'r--')
ax2.plot(redshift,MMBulge_3,'g:')
ax2.set_ylabel('Max Bulge Mass')
ax2.set_yscale('log')

ax3.plot(redshift,MMBHBulge,'b')
ax3.plot(redshift,MMBHBulge_2,'r--')
ax3.plot(redshift,MMBHBulge_3,'g:')
ax3.plot(redshift,MedBHBulge,'b')
ax3.plot(redshift,MedBHBulge_2,'r--')
ax3.plot(redshift,MedBHBulge_3,'g:')
ax3.plot(redshift,MinBHBulge,'b')
ax3.plot(redshift,MinBHBulge_2,'r--')
ax3.plot(redshift,MinBHBulge_3,'g:')
ax3.set_ylabel('Max, Median and Min BH/Bulge Mass')
ax3.set_yscale('log')
ax3.set_xlabel('Redshift')
plt.gca().invert_xaxis()
plt.show()
