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

data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'
filename='bulges_update09091_ddsf'

gals=meraxes.io.read_gals(data_folder+filename1+meraxes_loc,\
	snapshot=158,\
	h=cosmo['h'],quiet=True)
AMsamp=gals[(np.log10(gals['StellarMass']*1e10)>8.3)&(gals['StellarMass']>0)&(gals['Sfr']>0)]
bestSFRminus2sig=-10+np.log10(AMsamp['StellarMass']*1e10)*1 #from simple fit to bottom curve of sample in Figure 2
AMsamp_SF=AMsamp[np.log10(AMsamp['Sfr'])-bestSFR>0]
AM=np.sqrt(AMsamp_SF['AMstars'][:,0]**2+AMsamp_SF['AMstars'][:,1]**2+AMsamp_SF['AMstars'][:,2]**2)
np.mean(AM/AMsamp_SF['StellarMass'])

