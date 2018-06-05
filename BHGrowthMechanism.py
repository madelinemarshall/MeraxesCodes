import numpy as np
from dragons import meraxes
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd

#Sets plot defaults
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (8,7)
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
  MaxMass=10.75
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
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(MD_frac[TotNinBin>5]),label='Merger-Driven Quasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[3])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(ID_frac[TotNinBin>5]),label='Instability-Driven Quasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[4])
  axes.set_xlabel(r'$\log(M_{\mathrm{BH}})$')
  axes.set_xlim(MinMass-0.15,MaxMass-0.05)
  axes.set_ylim(-5.3,0.1)


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
  
  filename='tuned_reion'
  filename_125='tuned_reion_T125'

  fig,axes=plt.subplots(3,5,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  for snapshot in [78,100,116,158,213]:
    ii+=1
    if snapshot!=213:
      gals=load_data(filename,snapshot,cosmo)
      gals_d=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
      gals_b=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
      plot_BH_frac(gals,axes[0,ii])
      plot_BH_frac(gals_b,axes[1,ii])
      plot_BH_frac(gals_d,axes[2,ii])
    
    if ii==0:      
      axes[2,ii].plot([1,1.1],[1,1.1],linewidth=2.5,linestyle='-',color=colors[3],label='Tiamat')
      axes[2,ii].plot([1,1.1],[1,1.1],linewidth=1.5,linestyle='--',color=colors[3],label='Tiamat-125-HR')
      lgd=axes[2,ii].legend(fontsize='small',loc='lower left')#, bbox_to_anchor=(1.75, 0.73))

    gals=load_data(filename_125,snapshot,cosmo)
    gals_d=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
    gals_b=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
    plot_BH_frac(gals,axes[0,ii],'--',1.5)
    plot_BH_frac(gals_b,axes[1,ii],'--',1.5)
    plot_BH_frac(gals_d,axes[2,ii],'--',1.5)

    axes[0,ii].set_xticks([])
    if ii>1:
      axes[1,ii].set_xticks([])
    if ii>0:
      axes[0,ii].set_yticks([])
      axes[1,ii].set_yticks([])
      axes[2,ii].set_yticks([])
    else:
      axes[1,ii].set_ylabel(r'Fraction of BH Mass from Growth Modes')
    axes[0,ii].set_title(r'$z={}$'.format(redshift[snapshot]))
  axes[0,-1].set_ylabel(r'All Galaxies') 
  axes[1,-1].set_ylabel(r'B/T$>0/7$') 
  axes[2,-1].set_ylabel(r'B/T$<0.3$') 
  axes[0,-1].yaxis.set_label_position("right") 
  axes[1,-1].yaxis.set_label_position("right") 
  axes[2,-1].yaxis.set_label_position("right") 
  axes[2,0].axis('off')
  axes[2,1].axis('off')
  
  plt.tight_layout()
  
  plt.savefig('/home/mmarshal/results/plots/BHGrowthModes_MixedTrees.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
