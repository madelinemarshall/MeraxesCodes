import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
from _load_data import load_data

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,5.6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4


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
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(MD_frac[TotNinBin>5]),label='Merger-Driven Quasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[3])
  axes.plot(BinStart[TotNinBin>5]+0.5*BinWidth,np.log10(ID_frac[TotNinBin>5]),label='Instability-Driven Quasar-Mode',linewidth=lw,linestyle=linestyle,color=colors[4])
  axes.set_xlabel(r'$\log(M_{\mathrm{BH}}/M_\odot)$')
  axes.set_xlim(MinMass-0.15,MaxMass-0.05)
  axes.set_ylim(-5.3,0.1)


if __name__=="__main__":
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,250:0}
  
  filename='paper1_matching_eta_0p06'
  filename_125='paper1_T125_matching_eta_0p06'#'med_seed_T125'#'draft2_reion_T125'
  props=['BlackHoleMass','BlackHoleMass_ID','BlackHoleMass_MD','GhostFlag','BlackHoleMass_Coalescence','BlackHoleMass_Radio','BulgeStellarMass','StellarMass']

  fig,axes=plt.subplots(3,5,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  for snapshot in [78,100,116,158,250]:
    ii+=1
    if snapshot<160:
      gals=load_data(filename,snapshot,props)
      gals_d=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
      gals_b=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
      plot_BH_frac(gals,axes[0,ii])
      plot_BH_frac(gals_b,axes[1,ii])
      plot_BH_frac(gals_d,axes[2,ii])
    else:
      gals=load_data(filename_125,snapshot,props)
      gals_d=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
      gals_b=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7]
      plot_BH_frac(gals,axes[0,ii],'--',1.5)
      plot_BH_frac(gals_b,axes[1,ii],'--',1.5)
      plot_BH_frac(gals_d,axes[2,ii],'--',1.5)
    
    if ii==0:      
      axes[2,ii].plot([1,1.1],[1,1.1],linewidth=2.5,linestyle='-',color=colors[3],label='Tiamat')
      axes[2,ii].plot([1,1.1],[1,1.1],linewidth=1.5,linestyle='--',color=colors[3],label='Tiamat-125-HR')
      lgd=axes[2,ii].legend(fontsize='small',loc='lower left')


    axes[0,ii].set_xticks([])
    if ii>1:
      axes[1,ii].set_xticks([])
    if (ii>0):
      axes[0,ii].set_yticks([])
      axes[1,ii].set_yticks([])
      if ii!=2:
        axes[2,ii].set_yticks([])
      else:
        axes[2,ii].set_yticks(np.arange(-5,-0.5,1))
    else:
      axes[1,ii].set_ylabel(r'$\log(M_{\rm{BH,~ growth~ mode}}/M_{\rm{BH,~ total}})$')
    axes[0,ii].set_title(r'$z={}$'.format(redshift[snapshot]))
    axes[1,0].set_xlabel('')
    axes[1,1].set_xlabel('')
  axes[0,-1].set_ylabel(r'All Galaxies') 
  axes[1,-1].set_ylabel(r'B/T$>0/7$') 
  axes[2,-1].set_ylabel(r'B/T$<0.3$') 
  axes[0,-1].yaxis.set_label_position("right") 
  axes[1,-1].yaxis.set_label_position("right") 
  axes[2,-1].yaxis.set_label_position("right") 
  axes[2,0].axis('off')
  axes[2,1].axis('off')
  
  plt.tight_layout()
  
  plt.savefig('/home/mmarshal/results/plots/Paper2/BHGrowthModes.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
