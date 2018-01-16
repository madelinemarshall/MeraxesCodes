##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import ContourPlot as cp

def plot_disk_contribution(gals):
  #Figure 9 of Tonini+16
  MaxMass=12.5
  MinMass=9.0
  BinWidth=0.2
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
  #plt.show()


def load_gals_Tonini_match(filename,snapshot):
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
  
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                            snapshot=snapshot,\
                                            h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-1)&(gals["Mvir"]*1e10>8.6e8*100/cosmo['h'])]
  gals=gals[gals["Type"]==0]
  return gals


def load_gals(filename,snapshot):
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
  
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                            snapshot=snapshot,\
                                            h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]
  return gals

if __name__=="__main__":
  ##All galaxies, no Tonini+16 matching
  #gals_tiamat125=load_gals('bulges_tiamat125',158)
  #gals_tiamat=load_gals('bulges_update1102_full',158)
  #All galaxies, Tonini+16 matching
  gals_tiamat125=load_gals_Tonini_match('bulges_tiamat125',158)
  gals_tiamat=load_gals_Tonini_match('bulges_update1102_full',158)
  
  plot_disk_contribution(gals_tiamat)
  plot_disk_contribution(gals_tiamat125)
  plt.legend(['Tiamat','Tiamat-125-HR'])
  plt.show()
