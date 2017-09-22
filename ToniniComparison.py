##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys

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
filename='bulges_update0921_ddsf'#update0901_ddsf'
snapshot=158
gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,\
                                          h=cosmo['h'],quiet=True)
gals=gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-1)&(gals["Mvir"]*1e10>8.6e8*100/cosmo['h'])]
gals=gals[gals["Type"]==0]
#Figure 3
#Gas disk radius is the radius where the surface density drops to 0.1*original value
diskdiam=-2*np.log(0.1)*gals['GasDiskScaleLength']*1000
plt.plot(np.log10(diskdiam),np.log10((gals['ColdGas'])*1e10),'.');
plt.xlim(0.8,2.5); plt.ylim(8.0,11.5)
plt.xlabel('Log(Gas Disk Diameter (kpc))'); plt.ylabel('Log(Gas Disk Mass)');
plt.show()

#Figure 4
diskFrac=1.0-gals['BulgeStellarMass']/gals['StellarMass']
ScMass=(gals['StellarMass'][diskFrac>=0.8]-gals['BulgeStellarMass'][diskFrac>=0.8])*1e10
ScRad=1.67835*gals['StellarDiskScaleLength'][diskFrac>=0.8]*1000
plt.plot(np.log10(ScMass),np.log10(ScRad),'.'); plt.xlim(9.5,11.7); plt.ylim(-0.2,2.2);
plt.xlabel('Log(Stellar Disk Mass)')
plt.ylabel('Log(Half-mass Radius (kpc)')
plt.show()

Mass=(gals['StellarMass']-gals['BulgeStellarMass'])*1e10
Rad=gals['StellarDiskScaleLength']*1000
plt.plot(np.log10(Mass),np.log10(Rad),'.'); plt.xlim(9.0,12.0); plt.ylim(-2,4);
plt.xlabel('Log(Stellar Disk Mass)'); plt.ylabel('Log(StellarDiskScaleLength) (kpc)');
plt.show()

#Figure 9
MaxMass=12.5
MinMass=9.0
BinWidth=0.1
nBins=int(np.ceil((MaxMass-MinMass)/BinWidth))
TotMinBin=np.zeros(nBins)
BulgeMinBin=np.zeros(nBins)
BinStart=np.zeros(nBins)
SM=gals['StellarMass']*1e10
BSM=gals['BulgeStellarMass']*1e10
for ii in range(0,nBins):
  BinStart[ii]=MinMass+ii*BinWidth
  BinEnd=MinMass+(ii+1)*BinWidth
  TotMinBin[ii]=np.nansum(SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
  BulgeMinBin[ii]=np.nansum(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
FracInBin=BulgeMinBin/TotMinBin
plt.plot(BinStart+0.5*BinWidth,1-FracInBin)
plt.xlabel('Log(Stellar Mass)')
plt.ylabel('Relative Contribution of Disk to Stellar Mass')
plt.ylim(0,1)
plt.show()
