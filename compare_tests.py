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
#BH growth rate double the default amount (was 0.05, 0.1 here)
gals_doubleBHGR=meraxes.io.read_gals(data_folder+'/snap81_BHGR0p1'+meraxes_loc,\
                                        snapshot=81,\
                                        h=cosmo['h'],quiet=True)
gals_doubleBHGR = gals_doubleBHGR[(gals_doubleBHGR["GhostFlag"]==0)]
#BH growth rate ten times the default amount (was 0.05, 0.5 here)
gals_tentimesBHGR=meraxes.io.read_gals(data_folder+'/snap81_BHGR0p5'+meraxes_loc,\
                                        snapshot=81,\
                                        h=cosmo['h'],quiet=True)
gals_tentimesBHGR = gals_tentimesBHGR[(gals_tentimesBHGR["GhostFlag"]==0)]
#Eddington ratio half the default amount (was 1, 0.5 here)
gals_halfER=meraxes.io.read_gals(data_folder+'/snap81_EddRat0p5'+meraxes_loc,\
                                        snapshot=81,\
                                        h=cosmo['h'],quiet=True)
gals_halfER = gals_halfER[(gals_halfER["GhostFlag"]==0)]
#BH seed mass an order of magnitude smaller than the default amount (was 1e-7, 1e-8 here)
gals_1em8BHseed=meraxes.io.read_gals(data_folder+'/snap81_BHseed1e-8'+meraxes_loc,\
                                        snapshot=81,\
                                        h=cosmo['h'],quiet=True)
gals_1em8BHseed = gals_1em8BHseed[(gals_1em8BHseed["GhostFlag"]==0)]
#BH seed mass an order of magnitude larger than the default amount (was 1e-7, 1e-6 here)
gals_1em6BHseed=meraxes.io.read_gals(data_folder+'/snap81_BHseed1e-6'+meraxes_loc,\
                                        snapshot=81,\
                                        h=cosmo['h'],quiet=True)
gals_1em6BHseed = gals_1em6BHseed[(gals_1em6BHseed["GhostFlag"]==0)]

#Plot histograms of the BH mass distribution for each of the different models
maxval=np.nanmax(np.log10(gals_default['BlackHoleMass'][np.isfinite(np.log10(gals_default['BlackHoleMass']))]*1e10))
minval=np.nanmin(np.log10(gals_default['BlackHoleMass'][np.isfinite(np.log10(gals_default['BlackHoleMass']))]*1e10))
hist_default, bin_edges = np.histogram(np.log10(gals_default['BlackHoleMass']*1e10),range=(minval,maxval))
hist_doubleBHGR, bin_edges = np.histogram(np.log10(gals_doubleBHGR['BlackHoleMass']*1e10),range=(minval,maxval))
hist_tentimesBHGR, bin_edges = np.histogram(np.log10(gals_tentimesBHGR['BlackHoleMass']*1e10),range=(minval,maxval))
hist_halfER, bin_edges = np.histogram(np.log10(gals_halfER['BlackHoleMass']*1e10),range=(minval,maxval))
hist_1em6BHseed, bin_edges = np.histogram(np.log10(gals_1em6BHseed['BlackHoleMass']*1e10),range=(minval,maxval))
hist_1em8BHseed, bin_edges = np.histogram(np.log10(gals_1em8BHseed['BlackHoleMass']*1e10),range=(minval,maxval))

plt.plot(bin_edges[:-1],np.log10(hist_default),label='Default')
plt.plot(bin_edges[:-1],np.log10(hist_doubleBHGR),label='Double BHGR')
plt.plot(bin_edges[:-1],np.log10(hist_tentimesBHGR),label='Ten times BHGR')
plt.plot(bin_edges[:-1],np.log10(hist_halfER),label='Half EddRatio')
plt.plot(bin_edges[:-1],np.log10(hist_1em6BHseed),label='10 times BHseed')
plt.plot(bin_edges[:-1],np.log10(hist_1em8BHseed),label='0.1 times BHseed')
plt.legend()
plt.title('Histogram of BH Mass')
plt.xlabel('log(BH Mass)')
plt.ylabel('log(Number of BHs)')
plt.show()

#Notes: - A larger BH growth rate leads to more massive BHs (double BHGR - ~
# factor of 2 increase, 10*BHGR - ~ order of magnitude increase)
# - A smaller Eddington ratio leads to slightly less massive BHs, but converges
# on the low-mass end
# - Having larger BH seed masses significantly increases the number of BHs
# (especially at the low-mass end). The reverse is true for smaller BH seed masses

#Plot histograms of any other mass distribution for each of the different models
prop='StellarMass'
#prop='ColdGas'
#prop='HotGas'
#prop='Mvir' 
maxval=np.nanmax(np.log10(gals_default[prop][np.isfinite(np.log10(gals_default[prop]))]*1e10))
minval=np.nanmin(np.log10(gals_default[prop][np.isfinite(np.log10(gals_default[prop]))]*1e10))
hist_default, bin_edges = np.histogram(np.log10(gals_default[prop]*1e10),range=(minval,maxval))
hist_doubleBHGR, bin_edges = np.histogram(np.log10(gals_doubleBHGR[prop]*1e10),range=(minval,maxval))
hist_tentimesBHGR, bin_edges = np.histogram(np.log10(gals_tentimesBHGR[prop]*1e10),range=(minval,maxval))
hist_halfER, bin_edges = np.histogram(np.log10(gals_halfER[prop]*1e10),range=(minval,maxval))
hist_1em6BHseed, bin_edges = np.histogram(np.log10(gals_1em6BHseed[prop]*1e10),range=(minval,maxval))
hist_1em8BHseed, bin_edges = np.histogram(np.log10(gals_1em8BHseed[prop]*1e10),range=(minval,maxval))

plt.plot(bin_edges[:-1],np.log10(hist_default),label='Default')
plt.plot(bin_edges[:-1],np.log10(hist_doubleBHGR),label='Double BHGR')
plt.plot(bin_edges[:-1],np.log10(hist_tentimesBHGR),label='Ten times BHGR')
plt.plot(bin_edges[:-1],np.log10(hist_halfER),label='Half EddRatio')
plt.plot(bin_edges[:-1],np.log10(hist_1em6BHseed),label='10 times BHseed')
plt.plot(bin_edges[:-1],np.log10(hist_1em8BHseed),label='0.1 times BHseed')
plt.legend()
plt.xlabel('log('+prop+')')
plt.ylabel('log(Number of galaxies)')
plt.title('Histogram of '+prop)
plt.show()

#Notes: - Changing BH properties doesn't alter Mvir (good)
# - Hot gas mass is less (~0.5-0.1 magnitude) when the BH seed is larger, for lower
# hot gas masses (<10^4)
# - The reverse is ~ true when the BH seed is smaller - but varies from larger to smaller
# - Cold gas mass is less for larger BH seed and more for smaller BH seed
# - Cold gas mass are similar for all models, and converge at high-mass end
# - Stellar mass is very similar for all models
