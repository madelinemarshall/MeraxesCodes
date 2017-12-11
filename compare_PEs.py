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

snaps=[52,63,78,100,116,134,158];
n_red=len(snaps)
n_split=5
f1_bin=np.zeros((n_split+1,n_red+1));
f2_bin=np.zeros((n_split+1,n_red+1));
f3_bin=np.zeros((n_split+1,n_red+1));
f4_bin=np.zeros((n_split+1,n_red+1));
f5_bin=np.zeros((n_split+1,n_red+1));
M_bin=np.zeros((n_split+1,n_red+1));
G=43.01

#Import simulation data
data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'
filename='bulges_update1014_proj'

ss=0
color=['k','b','r','g','m','purple','orange']
fig, axes = plt.subplots(1,1)
for snapshot in snaps:
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']>0]

  f1=G*(gals['StellarMass']-gals['BulgeStellarMass'])*gals['ColdGas']*gals['StellarDiskScaleLength']/(gals['StellarDiskScaleLength']+gals['GasDiskScaleLength'])**2
  f2=G*(gals['StellarMass']-gals['BulgeStellarMass'])**2/(4*gals['StellarDiskScaleLength'])
  f3=(gals['StellarMass']-gals['BulgeStellarMass'])*(gals['Vvir']**2)           
  f4=(gals['StellarMass']-gals['BulgeStellarMass'])*(gals['VStellarDisk']**2)/2
#  f5=G*(gals['StellarMass']-gals['BulgeStellarMass'])*gals['BulgeStellarMass']/gals['StellarDiskScaleLength']
  sigmav=10**(0.31*np.log10(gals['BulgeStellarMass']*1e10)-1.15)
  print(np.mean(sigmav))
  f5=(gals['StellarMass']-gals['BulgeStellarMass'])*((np.sqrt(3)*sigmav)**2)         
  Mvir=np.log10(gals['Mvir'])+10
  M_split=(max(Mvir)-min(Mvir))/n_split
  for ii in range(0,n_split):
    f1_bin[ii,ss]=np.nanmedian(f1[(Mvir>min(Mvir)+M_split*ii)&(Mvir<min(Mvir)+M_split*(ii+1))])
    f2_bin[ii,ss]=np.nanmedian(f2[(Mvir>min(Mvir)+M_split*ii)&(Mvir<min(Mvir)+M_split*(ii+1))])
    f3_bin[ii,ss]=np.nanmedian(f3[(Mvir>min(Mvir)+M_split*ii)&(Mvir<min(Mvir)+M_split*(ii+1))])
    f4_bin[ii,ss]=np.nanmedian(f4[(Mvir>min(Mvir)+M_split*ii)&(Mvir<min(Mvir)+M_split*(ii+1))])
    f5_bin[ii,ss]=np.nanmedian(f5[(Mvir>min(Mvir)+M_split*ii)&(Mvir<min(Mvir)+M_split*(ii+1))])
    M_bin[ii,ss]=(min(Mvir)+M_split*ii+min(Mvir)+M_split*(ii+1))/2
  axes.plot(M_bin[:,ss],np.log10(f1_bin[:,ss]/f4_bin[:,ss]),':',color=color[ss])
  axes.plot(M_bin[:,ss],np.log10(f2_bin[:,ss]/f4_bin[:,ss]),'-.',color=color[ss])
  axes.plot(M_bin[:,ss],np.log10(f3_bin[:,ss]/f4_bin[:,ss]),'-',color=color[ss])
  axes.plot(M_bin[:,ss],np.log10(f5_bin[:,ss]/f4_bin[:,ss]),'--',color=color[ss])
  ss+=1

lines=axes.get_lines()
legend1=plt.legend([lines[i] for i in [0,1,2,3]],['PE Gas Disk / PE Halo','PE Stellar Disk / PE Halo','KE / PE Halo','PE Bulge / PE Halo'],loc=2)
axes.add_artist(legend1)
plt.legend([lines[i] for i in [0,4,8,12,16,20,24]],['z=8','z=7','z=6','z=5','z=4','z=3','z=2'],loc=4)
#axes.add_artist(legend2)
plt.xlabel('log(Halo Mass)')
plt.ylabel('log(Fraction)')
plt.show()
