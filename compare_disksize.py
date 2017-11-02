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
data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'
snapshot=60#158

#Load data:
#Default parameters
gals_default=meraxes.io.read_gals(data_folder+'default'+meraxes_loc,\
                                        snapshot=snapshot,\
                                        h=cosmo['h'],quiet=True)
gals_default = gals_default[(gals_default["GhostFlag"]==0)]#remove ghosts
gals_default = gals_default[(gals_default['StellarMass']>1e-4)]#remove ghosts
#gals_default = gals_default[(gals_default['BulgeStellarMass']>1e-3)]#/gals_bulges['StellarMass']<0.1)]gals_default = gals_default[(gals_default['BulgeStellarMass']/gals_default['StellarMass']<0.1)]

filename1=str(sys.argv[1])
filename2=str(sys.argv[2])

gals_bulges=meraxes.io.read_gals(data_folder+filename1+meraxes_loc,\
                                         snapshot=snapshot,\
                                         h=cosmo['h'],quiet=True)
gals_bulges = gals_bulges[(gals_bulges["GhostFlag"]==0)]#remove ghosts
gals_bulges = gals_bulges[(gals_bulges['StellarMass']>1e-4)]#remove ghosts
#gals_bulges = gals_bulges[(gals_bulges['BulgeStellarMass']<1e-2)]#/gals_bulges['StellarMass']<0.1)]

gals_bulges_2=meraxes.io.read_gals(data_folder+filename2+meraxes_loc,\
                                         snapshot=snapshot,\
                                         h=cosmo['h'],quiet=True)
gals_bulges_2 = gals_bulges_2[(gals_bulges_2["GhostFlag"]==0)]#remove ghosts
gals_bulges_2 = gals_bulges_2[(gals_bulges_2['StellarMass']>1e-4)]#remove ghosts
#gals_bulges_2 = gals_bulges_2[(gals_bulges_2['BulgeStellarMass'])<1e-2]#/gals_bulges_2['StellarMass']<0.1)]

#Plot Disk Scale Length Function
#prop='Sfr'
prop='StellarDiskScaleLength'
maxval=np.nanmax(np.log10(gals_bulges[prop][gals_bulges[prop]>0]*1000)) 
minval=np.nanmin(np.log10(gals_bulges[prop][gals_bulges[prop]>0]*1000))
hist_default, bin_edges = np.histogram(np.log10(gals_default['DiskScaleLength'][gals_default['DiskScaleLength']>0]*1000),range=(minval,maxval),bins=30)
hist_bulges, bin_edges = np.histogram(np.log10(gals_bulges[prop][gals_bulges[prop]>0]*1000),range=(minval,maxval),bins=30)
hist_bulges_2, bin_edges = np.histogram(np.log10(gals_bulges_2[prop][gals_bulges_2[prop]>0]*1000),range=(minval,maxval),bins=30)

bin_edges=np.array(bin_edges, dtype=np.float128)

Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
plt.plot(Max,hist_default/(bin_edges[1]-bin_edges[0])/100.**3,label='Default')
plt.plot(Max,hist_bulges/(bin_edges[1]-bin_edges[0])/100.**3,'--',label=filename1)
plt.plot(Max,hist_bulges_2/(bin_edges[1]-bin_edges[0])/100.**3,'--',label=filename2)
plt.xlabel('log(Disk Radius (kpc))')
plt.ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
plt.title('Disk Radius Function')
plt.yscale('log')
#plt.xlim([7.5,10.5])
#plt.ylim([10.**-4,10.**-0.8])
plt.legend()
plt.show()
