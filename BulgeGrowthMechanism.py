import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd

#Sets plot defaults
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (8.27,3.2)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,snapshot,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
    snapshot=snapshot,\
    h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['StellarMass']*1e10>10**6]#**8.9]
  #gals=gals[gals['Type']==0]
  return gals



def plot_BH_frac(gals,axes):
  #Figure 9
  MaxMass=12
  MinMass=7.25
  BinWidth=0.25
  nBins=int(np.ceil((MaxMass-MinMass)/BinWidth))
  BinStart=np.zeros(nBins)
  
  TotMinBin=np.zeros(nBins)
  TotNinBin=np.zeros(nBins)
  ID_MinBin=np.zeros(nBins)
  MD_MinBin=np.zeros(nBins)
  
  BH=gals['BulgeStellarMass']*1e10
  MD_BH=gals['MergerBulgeStellarMass']*1e10
  
  for ii in range(0,nBins):
    BinStart[ii]=MinMass+ii*BinWidth
    BinEnd=MinMass+(ii+1)*BinWidth
    TotMinBin[ii]=np.nansum(BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])
    TotNinBin[ii]=np.size(BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])
    MD_MinBin[ii]=np.nansum(MD_BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])

  MD_frac=MD_MinBin/TotMinBin
  ID_frac=1-MD_frac

  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,MD_frac[TotNinBin>5],label='Merger-Driven Bulge',linewidth=2.5)
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,ID_frac[TotNinBin>5],label='Instability-Driven Bulge',linewidth=2.5)
  axes.set_xlabel(r'$\log(M_{\mathrm{Bulge}})$')
  axes.set_xlim(MinMass,MaxMass)
  axes.set_yscale('log')
  axes.set_ylim(10**(-2.5),1.1)
  axes.text(10.3, 10**(-2.3), r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')


if __name__=="__main__":
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
  #meraxes_loc='/output/'+str(sys.argv[1])+'.hdf5'
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  
  filename='BHComponentSplit'

  fig,axes=plt.subplots(1,4,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  for snapshot in [78,116,158,213]:
    ii+=1
    gals=load_data(filename,snapshot,cosmo)
    plot_BH_frac(gals,axes[ii])
    if ii>0:
      axes[ii].set_yticks([])
    else:
      axes[ii].set_ylabel('Fraction of Bulge Mass \nfrom Growth Modes')
  plt.tight_layout()
  lgd=axes[-1].legend(fontsize='small',loc='upper center', bbox_to_anchor=(0.75, 0.53))
  
  #plt.savefig('BHGrowthModes.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
