#Does a quick comparison with the Okumara+17 j/m for stellar disks.
#Ours isn't crazy (0.5 relative to their 0.77)
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys

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

#Import simulation data
data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'
filename='bulges_update0915_ddsf'
snapshot=158

gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
	snapshot=snapshot,\
	h=cosmo['h'],quiet=True)

#Implement cuts as in paper
cuts=(gals['StellarMass']>10**(8.3)/1e10) & (gals['StellarMass']<10**(11.1)/1e10) \
& (gals['Mvir']>8/cosmo['h']) & (gals['Mvir']<2*1e2/cosmo['h']) \
& (gals['Sfr']>0)
z34cuts=gals['StellarMass']<10**(10.4)/1e10
samp=gals[cuts]

#Make estimate of SFR cuts
bestSFRminus2sig=-10+np.log10(samp['StellarMass']*1e10)*1 #from simple fit to bottom curve of sample in Figure 2
samp=samp[np.log10(samp['Sfr'])-bestSFRminus2sig>0]

AM=np.sqrt(samp['AMstars'][:,0]**2+samp['AMstars'][:,1]**2+samp['AMstars'][:,2]**2)
AMhalo=samp['Spin']*np.sqrt(2)*samp['Mvir']*samp['Vvir']*samp['Rvir']
AMrat=AM/AMhalo
mrat=(samp['StellarMass']-samp['BulgeStellarMass'])/samp['Mvir']

#Plot j/m against Mvir in mass bins (Fig 8 for 1 redshift)
Mbins=[8.3,9.0,9.7,10.4,11.1]
SM=np.log10(samp['StellarMass']*1e10)

AMonm=np.zeros(len(Mbins)-1)
Mdm=np.zeros(len(Mbins)-1)
for b in range(0,len(Mbins)-1):
  AMonm[b]=np.median(AMrat[(samp['StellarMass']>Mbins[b])&(samp['StellarMass']<Mbins[b+1])]\
  /mrat[(samp['StellarMass']>Mbins[b])&(samp['StellarMass']<Mbins[b+1])])
  Mdm[b]=np.median(samp['Mvir'][(samp['StellarMass']>Mbins[b])&(samp['StellarMass']<Mbins[b+1])])
plt.plot(np.log10(Mdm*1e10*cosmo['h']),AMonm,'o')
plt.xlabel('$log(M_{DM}/h)$')
plt.ylabel('$j_*/m_*$')
plt.xlim([10.5,12.5])
plt.ylim([0.45,1.3])
plt.show()

import ContourPlot as cp
cp.contour_plot(np.log10(mrat[(mrat>0)&(AMrat>0)]),np.log10(AMrat[(mrat>0)&(AMrat>0)]),'$log(j_*)$','$log(m_*)$')

