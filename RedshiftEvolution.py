##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import ContourPlot as cp

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
#filename='bulges_update1102_full'#'bulges_update0925_ddsf'#update0901_ddsf'
#filename='bulges_tiamat125'#'bulges_update0925_ddsf'#update0901_ddsf'
#sim='bulges_crotonSF'
#sim = 'bulges_update1102_full' #str(sys.argv[1])
sim='bulges_originalSFcriteria'
fmeraxes = '/home/mmarshal/data_dragons/'+sim+'/output/meraxes.hdf5'
snapshots=[52,78,115,158]
redshift={52:8,78:6,115:4,158:2}
gals={}
for snapshot in snapshots:
  gals[snapshot]=meraxes.io.read_gals(fmeraxes,\
                                snapshot=snapshot,\
                                h=cosmo['h'],quiet=True)
  gals[snapshot]=gals[snapshot][(gals[snapshot]["GhostFlag"]==0)&(gals[snapshot]["StellarMass"]*1e10>1e9)]
  gals[snapshot]=gals[snapshot][(gals[snapshot]["Type"]==0)&(gals[snapshot]["StellarMass"]*1e10<1e10)]

f, ((ax1, ax2,ax3), (ax4, ax5, ax6)) = plt.subplots(2,3)
for snapshot in snapshots:
  ax1.plot(redshift[snapshot],np.median(gals[snapshot]['StellarMass']*1e10),'ro')
  ax2.plot(redshift[snapshot],np.median(gals[snapshot]['ColdGas']*1e10),'ro')
  ax3.plot(redshift[snapshot],np.median(gals[snapshot]['StellarDiskScaleLength']),'ro')
  ax4.plot(redshift[snapshot],np.median(gals[snapshot]['GasDiskScaleLength']),'ro')
  ax5.plot(redshift[snapshot],np.median(gals[snapshot]['StellarDiskScaleLength']/gals[snapshot]['Rvir']),'r*')
  ax5.plot(redshift[snapshot],np.median(gals[snapshot]['GasDiskScaleLength']/gals[snapshot]['Rvir']),'ro')
  ax6.plot(redshift[snapshot],np.median(gals[snapshot]['StellarMass']/gals[snapshot]['Mvir']),'r*')
  ax6.plot(redshift[snapshot],np.median(gals[snapshot]['ColdGas']/gals[snapshot]['Mvir']),'ro')
ax1.plot(np.array(list(redshift.values())),(1/(1+np.array(list(redshift.values())))*(np.median(gals[52]['StellarMass']*1e10))*(1+8)),'r')
ax2.plot(np.array(list(redshift.values())),(1/(1+np.array(list(redshift.values())))*(np.median(gals[52]['ColdGas']*1e10))*(1+8)),'r')
ax3.plot(np.array(list(redshift.values())),((1+8)/(1+np.array(list(redshift.values()))))*((1+0.25*8)/(1+0.25*np.array(list(redshift.values()))))*(np.median(gals[52]['StellarDiskScaleLength'])),'r')
ax4.plot(np.array(list(redshift.values())),((1+8)/(1+np.array(list(redshift.values()))))*((1+0.25*8)/(1+0.25*np.array(list(redshift.values()))))*(np.median(gals[52]['GasDiskScaleLength'])),'r')

sim = 'default'
fmeraxes = '/home/mmarshal/data_dragons/'+sim+'/output/meraxes.hdf5'
gals={}
for snapshot in snapshots:
  gals[snapshot]=meraxes.io.read_gals(fmeraxes,\
                                snapshot=snapshot,\
                                h=cosmo['h'],quiet=True)
  #gals[snapshot]=gals[snapshot][(gals[snapshot]["GhostFlag"]==0)&(gals[snapshot]["StellarMass"]>0)]
  #gals[snapshot]=gals[snapshot][(gals[snapshot]["Type"]==0)&(gals[snapshot]["ColdGas"]>0)]
  gals[snapshot]=gals[snapshot][(gals[snapshot]["GhostFlag"]==0)&(gals[snapshot]["StellarMass"]*1e10>1e9)]
  gals[snapshot]=gals[snapshot][(gals[snapshot]["Type"]==0)&(gals[snapshot]["StellarMass"]*1e10<1e10)]
for snapshot in snapshots:
  ax1.plot(redshift[snapshot],np.median(gals[snapshot]['StellarMass']*1e10),'bo')
  ax2.plot(redshift[snapshot],np.median(gals[snapshot]['ColdGas']*1e10),'bo')
  ax3.plot(redshift[snapshot],np.median(gals[snapshot]['DiskScaleLength']),'bo')
  ax4.plot(redshift[snapshot],np.median(gals[snapshot]['DiskScaleLength']),'bo')
  ax5.plot(redshift[snapshot],np.median(gals[snapshot]['DiskScaleLength']/gals[snapshot]['Rvir']),'bo')
  ax6.plot(redshift[snapshot],np.median(gals[snapshot]['StellarMass']/gals[snapshot]['Mvir']),'b*')
  ax6.plot(redshift[snapshot],np.median(gals[snapshot]['ColdGas']/gals[snapshot]['Mvir']),'bo')

ax1.plot(np.array(list(redshift.values())),(1/(1+np.array(list(redshift.values())))*(np.median(gals[52]['StellarMass']*1e10))*(1+8)),'b')
ax2.plot(np.array(list(redshift.values())),(1/(1+np.array(list(redshift.values())))*(np.median(gals[52]['ColdGas']*1e10))*(1+8)),'b')
ax3.plot(np.array(list(redshift.values())),((1+8)/(1+np.array(list(redshift.values()))))*((1+0.25*8)/(1+0.25*np.array(list(redshift.values()))))*(np.median(gals[52]['DiskScaleLength'])),'b')
ax4.plot(np.array(list(redshift.values())),((1+8)/(1+np.array(list(redshift.values()))))*((1+0.25*8)/(1+0.25*np.array(list(redshift.values()))))*(np.median(gals[52]['DiskScaleLength'])),'b')

ax1.set_ylabel('StellarMass')
ax2.set_ylabel('Cold Gas Mass')
ax3.set_ylabel('Stellar R')
ax4.set_ylabel('Gas R')
ax5.set_ylabel('R/Rvir')
ax6.set_ylabel('M/Mvir')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax4.set_yscale('log')
ax5.set_yscale('log')
ax6.set_yscale('log')
ax1.invert_xaxis()
ax2.invert_xaxis()
ax3.invert_xaxis()
ax4.invert_xaxis()
ax5.invert_xaxis()
ax6.invert_xaxis()
plt.show()



