import numpy as np
from dragons import meraxes
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,2.5)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4

def load_data(filename,snapshot):
  gals=meraxes.io.read_gals(data_folder+filename,\
      snapshot=snapshot,\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals


def plot_BH_frac(gals,axes,linestyle='-',lw=2.5):
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
    Coalescence_MinBin[ii]=np.nansum(Coalescence_BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])
    Radio_MinBin[ii]=np.nansum(Radio_BH[(np.log10(BH)>=BinStart[ii])&(np.log10(BH)<BinEnd)])

  ID_frac=ID_MinBin/TotMinBin
  MD_frac=MD_MinBin/TotMinBin
  Coalescence_frac=Coalescence_MinBin/TotMinBin
  Radio_frac=Radio_MinBin/TotMinBin
  Rest=1-(ID_MinBin+MD_MinBin+Coalescence_MinBin+Radio_MinBin)/TotMinBin
  Rest_MinBin=TotMinBin-(ID_MinBin+MD_MinBin+Coalescence_MinBin+Radio_MinBin)

  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(Rest[TotNinBin>5]),label='Seed Condition',linewidth=lw,linestyle=linestyle,color=colors[0])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(Radio_frac[TotNinBin>5]),label='Radio-Mode',linewidth=lw,linestyle=linestyle,color=colors[1])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(Coalescence_frac[TotNinBin>5]),label='BH-BH Coalescence',linewidth=lw,linestyle=linestyle,color=colors[2])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(MD_frac[TotNinBin>5]),label='Merger-Driven\nQuasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[3])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(ID_frac[TotNinBin>5]),label='Instability-Driven\nQuasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[4])
  axes.set_xlabel(r'$\log(M_{\mathrm{BH}}/M_\odot)$')
  axes.set_xlim(MinMass-0.15,MaxMass-0.05)
  axes.set_ylim(-5.3,0.1)


if __name__=="__main__":
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,250:0}
  
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
  filenames=['paper2/output/meraxes.hdf5']+['paper2_kc_models/output_0p0'+str(f)+'/meraxes.hdf5' for f in [1,3,9]] 
  filenames_T125=['paper2_T125/output/meraxes.hdf5']+['paper2_kc_models_T125/output_0p0'+str(f)+'/meraxes.hdf5' for f in [1,3,9]] 
  kc=[0.005,0.01,0.03,0.09]
  props=['BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag','BlackHoleMass_Coalescence','BlackHoleMass_Radio','BulgeStellarMass','StellarMass']

  fig,axes=plt.subplots(1,5,gridspec_kw = {'wspace':0,'hspace':0.3},sharey=True,sharex=True)
  #ii=-1
  #snapshot=158
  #for filename in filenames:
  #  ii+=1
  #  gals=load_data(filename,snapshot)
  #  plot_BH_frac(gals,axes[0,ii])
  #  axes[0,ii].set_title(r'$k_c={}$'.format(kc[ii]))
  
  ii=-1
  snapshot=250
  for filename in filenames_T125:
    ii+=1
    gals=load_data(filename,snapshot)
    plot_BH_frac(gals,axes[ii],linestyle='--',lw=1.5)
    
  axes[0].set_ylim([-4.5,0.2])
  axes[0].set_xlim([6.5,10.8])
  axes[4].axis('off')
  #axes[1,4].axis('off')
  lgd=axes[3].legend(fontsize='small',loc=(1.01,0.5))
  axes[0].set_ylabel(r'$\log\left(M_{\rm{BH,~ growth~ mode}}/M_{\rm{BH,~ total}}\right)$')
  #axes[1,0].set_ylabel(r'$\log\left(M_{\rm{BH,~ growth~ mode}}/M_{\rm{BH,~ total}}\right)$')
  #plt.tight_layout()
  plt.subplots_adjust(bottom=0.2,left=0.08,right=0.95)

  #plt.savefig('/home/mmarshal/results/plots/Paper2/BHGrowthModes.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.savefig('/home/mmarshal/results/plots/Paper2/BHGrowthModes_modelComparison_T125.pdf',format='pdf')
  plt.show()
