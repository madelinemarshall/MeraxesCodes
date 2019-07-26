import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import magcalc as mc
sys.path.append('/home/mmarshal/simulation_codes/Paper1Plots')
from _load_data import load_data

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.6,4)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
color=['#ff7f00','#4daf4a','#984ea3','#98ff98']#,\'#377eb8'
                  #134:,158:'#f781bf',194:'#a65628',213:'black'}


def plot_disk_frac(gals,axes,label,split,bulge_frac,mags):
  MaxMass=11.5
  MinMass=7.0
  BinWidth=0.2 #Same as obs
  nBins=int(np.ceil((MaxMass-MinMass)/BinWidth))
  TotMinBin=np.zeros(nBins)
  TotNinBin=np.zeros(nBins)
  BulgeMinBin=np.zeros(nBins)
  IBulgeMinBin=np.zeros(nBins)
  MBulgeMinBin=np.zeros(nBins)
  IFracInBin=np.zeros((nBins,3))
  MFracInBin=np.zeros((nBins,3))
  BFracInBin=np.zeros((nBins,3))
  
  BinStart=np.zeros(nBins)
  BinFull=np.zeros(nBins, dtype=bool)
  SM=gals['StellarMass']*1e10
  BSM=gals['BulgeStellarMass']*1e10
  if split:
    MBSM=gals['MergerBulgeStellarMass']*1e10
    IBSM=BSM-MBSM
  for ii in range(0,nBins):
    BinStart[ii]=MinMass+ii*BinWidth
    BinEnd=MinMass+(ii+1)*BinWidth
    condition=(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)
    if (np.size(BSM[condition])>5):
      BinFull[ii]=True
      BFracInBin[ii,:]=np.percentile(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]/SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)],[16,50,84])
     # print("Num in bin: {}, bin start: {}, bin end: {}".format(np.size(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]),BinStart[ii],BinEnd))
      if split:
        IFracInBin[ii,:]=np.percentile(IBSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]/SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)],[16,50,84])
        MFracInBin[ii,:]=np.percentile(MBSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]/SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)],[16,50,84])

  if split:
    if  not mags:
      axes.plot(BinStart[BinFull]+0.5*BinWidth,IFracInBin[BinFull][:,1],label=label+r'Instability-driven Bulge',lw=2,color=color[0],linestyle='--')
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,IFracInBin[BinFull][:,0],IFracInBin[BinFull][:,2],alpha=0.2,color=color[0])
      var, = axes.plot(BinStart[BinFull]+0.5*BinWidth,MFracInBin[BinFull][:,1],label=label+r'Merger-driven Bulge',lw=2,color=color[1],linestyle='-.')
      var.set_dashes([2,2,5,2])
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,MFracInBin[BinFull][:,0],MFracInBin[BinFull][:,2],alpha=0.2,color=color[1])
      axes.plot(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,1],label=label+r'Total Bulge',lw=2,color=color[2])
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,0],BFracInBin[BinFull][:,2],alpha=0.2,color=color[2])
    else:
      #axes.plot(BinStart+0.5*BinWidth,IFracInBin[:,1],label=label+' - Instability-driven Bulge',lw=2.5,color='firebrick',linestyle='--')
      #axes.fill_between(BinStart+0.5*BinWidth,IFracInBin[:,0],IFracInBin[:,2],alpha=0.2,color='firebrick')
      #axes.plot(BinStart+0.5*BinWidth,MFracInBin[:,1],label=label+' - Merger-driven Bulge',lw=2.5,color='gold',linestyle='--')
      #axes.fill_between(BinStart+0.5*BinWidth,MFracInBin[:,0],MFracInBin[:,2],alpha=0.2,color='gold')
      axes.plot(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,1],label=label+r'Total Bulge',lw=2,color=color[2],linestyle='--')
      #axes.fill_between(BinStart+0.5*BinWidth,BFracInBin[:,0],BFracInBin[:,2],alpha=0.2,color='royalblue')
      #axes.errorbar(BinStart+0.5*BinWidth,IFracInBin[:,1],[IFracInBin[:,1]-IFracInBin[:,0],IFracInBin[:,2]-IFracInBin[:,1]],label=label+' - Instability-driven Bulge',color='firebrick',linestyle='--')
      #axes.errorbar(BinStart+0.5*BinWidth,MFracInBin[:,1],[MFracInBin[:,1]-MFracInBin[:,0],MFracInBin[:,2]-MFracInBin[:,1]],label=label+' - Merger-driven Bulge',color='gold',linestyle='--')
      #axes.errorbar(BinStart+0.5*BinWidth,BFracInBin[:,1],[BFracInBin[:,1]-BFracInBin[:,0],BFracInBin[:,2]-BFracInBin[:,1]],label=label+' - Total Bulge',color='royalblue',linestyle='--')

  else:
    axes.plot(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,1],label=label+r'Total Bulge',lw=2,color=color[1],linestyle='-')
    axes.fill_between(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,0],BFracInBin[BinFull][:,2],alpha=0.1,color=color[1])
  #axes.set_ylabel('Stellar Mass Fraction of Components')
  #print(err)
  #print(1-FracInBin)
  #axes.plot([9,12],[0.5,0.5],'k:',label='_nolegend_')
  #axes.plot(BinStart+0.5*BinWidth,1-FracInBin,label=label)
  axes.set_ylim(0.01,1)
  axes.set_xlim(MinMass+0.5*BinWidth,MaxMass-0.5*BinWidth)



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
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  prop='StellarMass'
  snapshots=[63,78,100,115,134,158]
  filename='paper1'
  split=True
  bulge_frac=True


  fig_frac,axes=plt.subplots(2,4,gridspec_kw = {'wspace':0, 'hspace':0,'width_ratios':[4,4,4,2.4]},sharey=True,sharex=True)
  ii=-1
  jj=0
  for snapshot in snapshots:
    ii+=1
    if ii==3:
      jj+=1
      ii=0
    gals_bulges=load_data(filename,snapshot,[prop,'GhostFlag','BulgeStellarMass','Type','MergerBulgeStellarMass'],centrals=True)

    BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
    plot_disk_frac(gals_bulges,axes[jj,ii],'',split,bulge_frac,0)
    if jj==1:
      axes[jj,ii].set_xlabel(r'$\log(M_\ast/M_\odot)$')
    #else:
    #  axes[jj,ii].set_xticklabels([])
    if ii==0:
      axes[jj,ii].set_ylabel(r'$M_{\mathrm{Bulge}}/M_\ast$')
    #else:
    #  axes[jj,ii].set_yticklabels([])
    axes[jj,ii].text(10, 0.1, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
  
  axes[0,3].axis('off')
  axes[1,3].axis('off')
  #plt.tight_layout()
  lgd=axes[0,2].legend(fontsize='small',loc=(1.05,-0.15))#,loc='upper center', bbox_to_anchor=(0.5, -0.2))
  
  plt.savefig('/home/mmarshal/results/plots/Paper1/BulgeFractionEvolution.pdf', format='pdf')
  plt.show()
