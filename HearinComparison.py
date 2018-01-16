import numpy as np
from dragons import meraxes
#import os
import matplotlib.pyplot as plt
import sys
import ContourPlot as cp
import scipy.stats as stats

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


if len(sys.argv)==2:
  filename=str(sys.argv[1])
else:
  filename='bulges_update1102_full'

snapshot=163

gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,\
                                          h=cosmo['h'],quiet=True)
gals=gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-1)&(gals["Mvir"]*1e10>8.6e8*100/cosmo['h'])]
#gals=gals[gals["Type"]==0]
BtoT=gals["BulgeStellarMass"]/gals["StellarMass"]
gals_B=gals[BtoT>0.7]
gals_D=gals[BtoT<0.3]

halfrad=1.67835*gals['StellarDiskScaleLength']*1000
halfrad_B=1.67835*gals_B['StellarDiskScaleLength']*1000
halfrad_D=1.67835*gals_D['StellarDiskScaleLength']*1000
#Plot simulated galaxies
#plt.plot(np.log10(gals_D['StellarMass']*1e10),np.log10(halfrad_D),'r.')
#plt.plot(np.log10(gals_B['StellarMass']*1e10),np.log10(halfrad_B),'b.')
#cp.contour_plot(np.log10(gals_D['StellarMass'][halfrad_D>0]*1e10),np.log10(halfrad_D[halfrad_D>0]),colors='b')
#cp.contour_plot(np.log10(gals_B['StellarMass'][halfrad_B>0]*1e10),np.log10(halfrad_B[halfrad_B>0]),colors='r')
cp.contour_plot(np.log10(gals['StellarMass'][halfrad>0]*1e10),np.log10(halfrad[halfrad>0]))

plt.xlim([np.log10(3.1e9),np.log10(4.5e11)])
plt.ylim([np.log10(0.5),np.log10(31)])

#Overplot median
bin_means, bin_edges, binnumber = stats.binned_statistic(np.log10(gals['StellarMass']*1e10),\
  np.log10(halfrad),statistic='median', bins=10)
bin_width=bin_edges[1]-bin_edges[0]
plt.plot((bin_edges[:-1]+bin_width),bin_means,'k--')

bin_means, bin_edges, binnumber = stats.binned_statistic(np.log10(gals_D['StellarMass']*1e10),\
  np.log10(halfrad_D),statistic='median', bins=10)
bin_width=bin_edges[1]-bin_edges[0]
plt.plot((bin_edges[:-1]+bin_width),bin_means,'b--')

bin_means, bin_edges, binnumber = stats.binned_statistic(np.log10(gals_B['StellarMass']*1e10),\
  np.log10(halfrad_B),statistic='median', bins=10)
bin_width=bin_edges[1]-bin_edges[0]
plt.plot((bin_edges[:-1]+bin_width),bin_means,'r--')


plt.xlabel('Stellar Mass (solar masses)')
plt.ylabel('Half-Light Radius (kpc)')
plt.legend(['Median - all galaxies','Median - disk galaxies','Median - bulge galaxies'])
plt.show()


