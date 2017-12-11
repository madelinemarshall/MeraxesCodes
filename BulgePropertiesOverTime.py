##Plots BH mass, stellar mass and bulge fraction functions, bulge fraction against BH mass, bulge fraction
##against stellar mass, and histograms of bulge fraction in bins of BH mass. Also performs KS test to see 
##if distribution of bulge fractions is the same for log(MBH)>7 and other mass samples. 
#Input: snapshot

import numpy as np
from dragons import meraxes
import os
import sys
import matplotlib.pyplot as plt
from chi_sq import chi_sq
import scipy.stats as stats

cosmo = {'omega_M_0' : 0.308,
'omega_lambda_0' : 0.692,
'omega_b_0' : 0.04839912,
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}

sim = 'bulges_update1102_full' #str(sys.argv[1])
fmeraxes = '/home/mmarshal/data_dragons/'+sim+'/output/meraxes.hdf5'
snapshots=[43,52,63,78,100,115,134,158]
redshifts=[9,8,7,6,5,4,3,2]
gals={}
bulgeFrac={}
for snapshot in snapshots:
  gals[snapshot]=meraxes.io.read_gals(fmeraxes,\
                                snapshot=snapshot,\
                                h=cosmo['h'],quiet=True)
  bulgeFrac[snapshot]=gals[snapshot]['BulgeStellarMass']/gals[snapshot]['StellarMass']

medBF_1e8=np.zeros((len(snapshots),3))
medBF_maxSM=np.zeros((len(snapshots),3))
medBF_maxBHM=np.zeros((len(snapshots),3))
medBF_1000BHM=np.zeros((len(snapshots),3))
medBF_1000SM=np.zeros((len(snapshots),3))
medBF_medBHM=np.zeros((len(snapshots),3))
medBF_medSM=np.zeros((len(snapshots),3))

ii=0
for snapshot in snapshots:
  medBF_1e8[ii]=np.nanpercentile(bulgeFrac[snapshot][gals[snapshot]['StellarMass']>1e-2],[16,50,84]) 
  medBF_maxSM[ii]=np.nanpercentile(bulgeFrac[snapshot][np.log10(gals[snapshot]['StellarMass'])>np.max(np.log10(gals[snapshot]['StellarMass']))-1],[16,50,84])
  medBF_maxBHM[ii]=np.nanpercentile(bulgeFrac[snapshot][np.log10(gals[snapshot]['BlackHoleMass'])>np.max(np.log10(gals[snapshot]['BlackHoleMass']))-1],[16,50,84])

  ind = np.argpartition(gals[snapshot]['BlackHoleMass'], -1000)[-1000:]
  medBF_1000BHM[ii]=np.nanpercentile(bulgeFrac[snapshot][ind],[16,50,84])
  ind = np.argpartition(gals[snapshot]['StellarMass'], -1000)[-1000:]
  medBF_1000SM[ii]=np.nanpercentile(bulgeFrac[snapshot][ind],[16,50,84])
 
  medBH=np.nanmedian(gals[snapshot]['BlackHoleMass'])
  stdBH=np.nanpercentile(gals[snapshot]['BlackHoleMass'],[16,84]) 
  medSM=np.nanmedian(gals[snapshot]['StellarMass'])
  stdSM=np.nanpercentile(gals[snapshot]['StellarMass'],[16,84])
  
  medBF_medBHM[ii]=np.nanpercentile(bulgeFrac[snapshot][(gals[snapshot]['BlackHoleMass']>stdBH[0])&(gals[snapshot]['BlackHoleMass']<stdBH[1])],[16,50,84])

  medBF_medSM[ii]=np.nanpercentile(bulgeFrac[snapshot][(gals[snapshot]['StellarMass']>stdSM[0])&(gals[snapshot]['StellarMass']<stdSM[1])],[16,50,84])

  ii+=1

plt.plot(redshifts,medBF_maxSM[:,1],'k-')
plt.plot(redshifts,medBF_maxSM[:,0],'k--')
plt.plot(redshifts,medBF_maxSM[:,2],'k--')
plt.plot(redshifts,medBF_maxBHM[:,1],'r-')
plt.plot(redshifts,medBF_maxBHM[:,0],'r--')
plt.plot(redshifts,medBF_maxBHM[:,2],'r--')
plt.gca().invert_xaxis()
plt.xlabel('Redshift')
plt.ylabel('Median Bulge Fraction')
plt.legend(('Median - 1000 most massive BHs','16th percentile - 1000 most massive BHs','84th percentile - 1000 most massive BHs','Median - 1000 largest stellar masses','16th percentile - 1000 largest stellar masses','84th percentile - 1000 largest stellar masses'))
plt.show()



