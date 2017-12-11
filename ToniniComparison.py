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
#snapshot=163
filename='bulges_tiamat125'#'bulges_update0925_ddsf'#update0901_ddsf'
snapshot=213
#filename='bulges_originalSFcriteria'
#filename='bulges_crotonSF'

#snapshot=158
gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,\
                                          h=cosmo['h'],quiet=True)
gals=gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-1)&(gals["Mvir"]*1e10>8.6e8*100/cosmo['h'])]
gals=gals[gals["Type"]==0]
#Figure 3
#Gas disk radius is the radius where the surface density drops to 0.1*original value
diskdiam=-2*np.log(0.1)*gals['GasDiskScaleLength']*1000
cp.contour_plot(np.log10(diskdiam),np.log10(gals['ColdGas']*1e10),'Log(Gas Disk Diameter (kpc))','Log(Gas Disk Mass)',[0.8,2.5],[8.0,11.5])
MH1_low=12.88*(10**0.8*1000)**2 #Broeils & Rhee eqn
MH1_high=12.88*(10**2.5*1000)**2
MH1_low_lwr=9.03*(10**0.8*1000)**2 #Broeils & Rhee eqn
MH1_high_lwr=9.03*(10**2.5*1000)**2
MH1_low_upr=16.73*(10**0.8*1000)**2 #Broeils & Rhee eqn
MH1_high_upr=16.73*(10**2.5*1000)**2
plt.plot([0.8,2.5],np.log10([MH1_low,MH1_high]),color='orange');
plt.legend(['Meraxes Galaxies','Broeils & Rhee 1997 Relation'])
plt.plot([0.8,2.5],np.log10([MH1_low_lwr,MH1_high_lwr]),'--',color='orange');
plt.plot([0.8,2.5],np.log10([MH1_low_upr,MH1_high_upr]),'--',color='orange');

#plt.xlabel('Log(Gas Disk Diameter (kpc))'); plt.ylabel('Log(Gas Disk Mass)');
plt.show()

#Figure 4
diskFrac=1.0-gals['BulgeStellarMass']/gals['StellarMass']
ScMass=(gals['StellarMass'][diskFrac>=0.8]-gals['BulgeStellarMass'][diskFrac>=0.8])*1e10
ScRad=1.67835*gals['StellarDiskScaleLength'][diskFrac>=0.8]*1000
cp.contour_plot(np.log10(ScMass),np.log10(ScRad),[9.5,11.7],[-0.2,2.2])
Mdata=[8.69642,8.901704,9.092361,9.297737,9.493349,9.698732,9.904136,10.099755,10.300282,10.500779,10.7012205,10.901662,11.0922985,11.302496,11.497974]
Rdata=[1.3664275,1.6502042,1.8943602,2.082251,2.2395494,2.4439173,2.609627,2.7865145,2.9326444,3.177138,3.6472878,4.18701,4.9120507,5.8047624,7.1640797]
plt.plot(Mdata,np.log10(Rdata),'s')
plt.legend(['Meraxes Galaxies','Shen et al. 2003'])
plt.xlabel('Log(Stellar Disk Mass)')
plt.ylabel('Log(Half-mass Radius (kpc)')
plt.xlim(9.5,11.7); plt.ylim(-1,1.1);
plt.show()

Mass=(gals['StellarMass']-gals['BulgeStellarMass'])*1e10
Rad=gals['StellarDiskScaleLength']*1000
cp.contour_plot(np.log10(Mass),np.log10(Rad),[9.0,12.0],[-2,4]);
Mdata=[9.2,12]	
Rdata=[0.84,12.06]
Rdata_lwr=[0.25,3.68]
Rdata_upr=[2.39,35.35]
plt.plot(Mdata,np.log10(Rdata),color='orange')
plt.legend(['Meraxes Galaxies','Gadotti 2009 (From Tonini et al. 2016)'])
plt.plot(Mdata,np.log10(Rdata_lwr),'--',color='orange')
plt.plot(Mdata,np.log10(Rdata_upr),'--',color='orange')
plt.xlabel('Log(Stellar Disk Mass)'); plt.ylabel('Log(StellarDiskScaleLength) (kpc)');
plt.xlim(9.2,12); plt.ylim(-3,2);
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
