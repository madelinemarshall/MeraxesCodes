import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,3.6)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4

def load_data(filename,snapshot,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
    snapshot=snapshot,\
    h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['BlackHoleMass']*1e10>10**6]#**8.9]
  #gals=gals[gals['Type']==0]
  return gals



def plot_BH_frac(gals,axes,linestyle='-',lw=2.5):
  #Figure 9
  MaxMass=10.5
  MinMass=6
  BinWidth=0.25
  nBins=int(np.ceil((MaxMass-MinMass)/BinWidth))
  BinStart=np.zeros(nBins)
  
  TotMinBin=np.zeros(nBins)
  TotNinBin=np.zeros(nBins)
  ID_MinBin=np.zeros(nBins)
  MD_MinBin=np.zeros(nBins)
  Coalescence_MinBin=np.zeros(nBins)
  Radio_MinBin=np.zeros(nBins)
  
  BH=gals['BlackHoleMass']*1e10
  ID_BH=gals['BlackHoleMass_ID']*1e10
  MD_BH=gals['BlackHoleMass_MD']*1e10
  Coalescence_BH=gals['BlackHoleMass_Coalescence']*1e10
  Radio_BH=gals['BlackHoleMass_Radio']*1e10
  
  for ii in range(0,nBins):
    BinStart[ii]=MinMass+ii*BinWidth
    BinEnd=MinMass+(ii+1)*BinWidth
    TotMinBin[ii]=np.nansum(BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])
    TotNinBin[ii]=np.size(BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])
    ID_MinBin[ii]=np.nansum(ID_BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])
    MD_MinBin[ii]=np.nansum(MD_BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])

  ID_frac=ID_MinBin/TotMinBin
  MD_frac=MD_MinBin/TotMinBin
  Rest=1-(ID_MinBin+MD_MinBin)/TotMinBin
  Rest_MinBin=TotMinBin-(ID_MinBin+MD_MinBin)

  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,(MD_frac[TotNinBin>5]),label='Merger-Driven Quasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[3])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,(ID_frac[TotNinBin>5]),label='Instability-Driven Quasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[4])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,(Rest[TotNinBin>5]),label='Other Modes',linewidth=lw,linestyle=linestyle,color=[0.5,0.5,0.5])
  axes.set_xlabel(r'$\log(M_{\mathrm{BH}})$')
  axes.set_xlim(MinMass-0.15,MaxMass-0.05)
  axes.set_ylim(0,1)
  #axes.set_yscale('log')
  #axes.set_ylim(10**-1.8,10**0.05)


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
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,192:1,213:0.55,250:0}
  
  filename='tuned_reion'
  filename_125='tuned_reion_T125'

  fig,axes=plt.subplots(1,3,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  for snapshot in [250]:
    ii+=1
    gals=load_data(filename_125,snapshot,cosmo)
    gals_d=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
    gals_b=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
    plot_BH_frac(gals,axes[0])
    plot_BH_frac(gals_b,axes[1])
    plot_BH_frac(gals_d,axes[2])
    
    #axes[2].plot([1,1.1],[1,1.1],linewidth=2.5,linestyle='-',color=colors[3],label='Tiamat')
    #axes[2].plot([1,1.1],[1,1.1],linewidth=1.5,linestyle='-',color=colors[3],label='Tiamat-125-HR')
    lgd=axes[1].legend(loc='lower left', bbox_to_anchor=(0, -0.4))


    axes[2].set_yticks([])
    axes[1].set_yticks([])
    axes[0].set_ylabel(r'Fraction of BH Mass from Growth Mode')
  axes[0].set_title('All\nGalaxies') 
  axes[1].set_title('Bulge-Dominated\nGalaxies') 
  axes[2].set_title('Disc-Dominated\nGalaxies') 
  
  #plt.tight_layout()
  
  plt.savefig('/home/mmarshal/results/plots/BHGrowthModes_talk.pdf', format='pdf', bbox_inches='tight')
  plt.show()
