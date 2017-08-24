## Compares the distributions of various galaxy properties when altering certain
##parameters in the meraxes model.

import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt

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
data_folder='/home/mmarshal/data_dragons'
meraxes_loc='/output/meraxes.hdf5'
snapshot=81

#Load data:
#Default parameters
gals_default=meraxes.io.read_gals(data_folder+'/snap81_default'+meraxes_loc,\
                                        snapshot=81,\
                                        h=cosmo['h'],quiet=True)
gals_default = gals_default[(gals_default["GhostFlag"]==0)]#remove ghosts


gals_bulges=meraxes.io.read_gals(data_folder+'/bulges_update0823_ddsf'+meraxes_loc,\
                                         snapshot=81,\
                                         h=cosmo['h'],quiet=True)
gals_bulges = gals_bulges[(gals_bulges["GhostFlag"]==0)]#remove ghosts

gals_bulges_ddsf=meraxes.io.read_gals(data_folder+'/bulges_crotonSF'+meraxes_loc,\
                                         snapshot=81,\
                                         h=cosmo['h'],quiet=True)
gals_bulges_ddsf = gals_bulges_ddsf[(gals_bulges_ddsf["GhostFlag"]==0)]#remove ghosts


#Plot histograms of the BH mass distribution for each of the different models
maxval=np.nanmax(np.log10(gals_default['BlackHoleMass'][gals_default['BlackHoleMass']>0]*1e10))
minval=np.nanmin(np.log10(gals_default['BlackHoleMass'][gals_default['BlackHoleMass']>0]*1e10))
hist_default, bin_edges = np.histogram(np.log10(gals_default['BlackHoleMass']*1e10),range=(minval,maxval))
hist_bulges, bin_edges = np.histogram(np.log10(gals_bulges['BlackHoleMass']*1e10),range=(minval,maxval))
hist_bulges_ddsf, bin_edges = np.histogram(np.log10(gals_bulges_ddsf['BlackHoleMass']*1e10),range=(minval,maxval))

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.plot(bin_edges[:-1],hist_default/(bin_edges[1]-bin_edges[0])/100**3,label='Default')
ax1.plot(bin_edges[:-1],hist_bulges/(bin_edges[1]-bin_edges[0])/100**3,'--',label='Bulges (Correct Units SF)')
ax1.plot(bin_edges[:-1],hist_bulges_ddsf/(bin_edges[1]-bin_edges[0])/100**3,'--',label='Bulges (Density dependent SF)')
ax1.set_title('BH Mass Function')
ax1.set_xlabel('log(Mass)')
ax1.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax1.set_yscale('log')
ax1.legend()

#Plot histograms of stellar mass distribution for each of the different models
prop='StellarMass'
#prop='ColdGas'
#prop='HotGas'
#prop='Mvir'
maxval=np.nanmax(np.log10(gals_default[prop][gals_default[prop]>0]*1e10))
minval=np.nanmin(np.log10(gals_default[prop][gals_default[prop]>0]*1e10))
hist_default, bin_edges = np.histogram(np.log10(gals_default[prop]*1e10),range=(minval,maxval))
hist_bulges, bin_edges = np.histogram(np.log10(gals_bulges[prop]*1e10),range=(minval,maxval))
hist_bulges_ddsf, bin_edges = np.histogram(np.log10(gals_bulges_ddsf[prop]*1e10),range=(minval,maxval))
ax2.plot(bin_edges[:-1],hist_default/(bin_edges[1]-bin_edges[0])/100**3,label='Default')
ax2.plot(bin_edges[:-1],hist_bulges/(bin_edges[1]-bin_edges[0])/100**3,'--',label='Bulges')
ax2.plot(bin_edges[:-1],hist_bulges_ddsf/(bin_edges[1]-bin_edges[0])/100**3,'--',label='Bulges (Density dependent SF)')
ax2.set_xlabel('log(Mass)')
ax2.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax2.set_title('Stellar Mass Function')
ax2.set_yscale('log')

#Plot histograms of cold gas mass distribution for each of the different models
prop='ColdGas'
maxval=np.nanmax(np.log10(gals_default[prop][gals_default[prop]>0]*1e10))
minval=np.nanmin(np.log10(gals_default[prop][gals_default[prop]>0]*1e10))
hist_default, bin_edges = np.histogram(np.log10(gals_default[prop]*1e10),range=(minval,maxval))
hist_bulges, bin_edges = np.histogram(np.log10(gals_bulges[prop]*1e10),range=(minval,maxval))
hist_bulges_ddsf, bin_edges = np.histogram(np.log10(gals_bulges_ddsf[prop]*1e10),range=(minval,maxval))
ax3.plot(bin_edges[:-1],hist_default/(bin_edges[1]-bin_edges[0])/100**3,label='Default')
ax3.plot(bin_edges[:-1],hist_bulges/(bin_edges[1]-bin_edges[0])/100**3,'--',label='Bulges')
ax3.plot(bin_edges[:-1],hist_bulges_ddsf/(bin_edges[1]-bin_edges[0])/100**3,'--',label='Bulges (Density dependent SF)')
ax3.set_xlabel('log(Mass)')
ax3.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax3.set_title('Cold Gas Mass Function')
ax3.set_yscale('log')

#Plot histograms of hot gas mass distribution for each of the different models
prop='BulgeStellarMass'
maxval=np.nanmax(np.log10(gals_bulges[prop][gals_bulges[prop]>0]*1e10))
minval=np.nanmin(np.log10(gals_bulges[prop][gals_bulges[prop]>0]*1e10))
hist_bulges, bin_edges = np.histogram(np.log10(gals_bulges[prop]*1e10),range=(minval,maxval))
hist_disks, bin_edges = np.histogram(np.log10(gals_bulges['StellarMass']*1e10-gals_bulges[prop]*1e10),range=(minval,maxval))
hist_bulges_ddsf, bin_edges = np.histogram(np.log10(gals_bulges_ddsf[prop]*1e10),range=(minval,maxval))
hist_disks_ddsf, bin_edges = np.histogram(np.log10(gals_bulges_ddsf['StellarMass']*1e10-gals_bulges_ddsf[prop]*1e10),range=(minval,maxval))
ax4.plot(bin_edges[:-1],hist_bulges/(bin_edges[1]-bin_edges[0])/100**3,'r--',label='Bulges (Correct Units SF)')
ax4.plot(bin_edges[:-1],hist_disks/(bin_edges[1]-bin_edges[0])/100**3,'k-',label='Disks (Correct Units SF)')
ax4.plot(bin_edges[:-1],hist_bulges_ddsf/(bin_edges[1]-bin_edges[0])/100**3,'b--',label='Bulges (Density Dependent SF)')
ax4.plot(bin_edges[:-1],hist_disks_ddsf/(bin_edges[1]-bin_edges[0])/100**3,'g-',label='Disks (Density Dependent SF)')
ax4.set_xlabel('log(Mass)')
ax4.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax4.set_title('Stellar Mass Function')
ax4.set_yscale('log')
plt.legend()
plt.show()
