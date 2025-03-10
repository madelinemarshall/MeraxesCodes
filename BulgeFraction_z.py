import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF
import magcalc as mc

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (3.5,3)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,snapshot,prop,cosmo,split):
  if split:
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag','BulgeStellarMass','Type','MergerBulgeStellarMass'],\
      h=cosmo['h'],quiet=True)
  else:
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag','BulgeStellarMass','Type'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['StellarMass']*1e10>10**6]#**8.9]
  #gals=gals[gals['Type']==0]
  return gals


def load_mags(filename,snapshot):
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  #MUV=pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  #AUV=mc.reddening(1600,MUV,redshift[snapshot])
  #MUV_dust=MUV+AUV
  MUV=pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['Y098']
  AUV=mc.reddening(6300,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust
  #return MUV


def plot_disk_frac(gals,axes,label,split,bulge_frac):
  #Figure 9
  MaxMass=12
  MinMass=9.0
  BinWidth=0.2
  nBins=int(np.ceil((MaxMass-MinMass)/BinWidth))
  TotMinBin=np.zeros(nBins)
  TotNinBin=np.zeros(nBins)
  BulgeMinBin=np.zeros(nBins)
  IBulgeMinBin=np.zeros(nBins)
  MBulgeMinBin=np.zeros(nBins)
  BinStart=np.zeros(nBins)
  SM=gals['StellarMass']*1e10
  BSM=gals['BulgeStellarMass']*1e10
  if split:
    MBSM=gals['MergerBulgeStellarMass']*1e10
    IBSM=BSM-MBSM
  for ii in range(0,nBins):
    BinStart[ii]=MinMass+ii*BinWidth
    BinEnd=MinMass+(ii+1)*BinWidth
    TotMinBin[ii]=np.nansum(SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
    TotNinBin[ii]=np.size(SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
    BulgeMinBin[ii]=np.nansum(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
    if split:
      IBulgeMinBin[ii]=np.nansum(IBSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
      MBulgeMinBin[ii]=np.nansum(MBSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])

  if split:
    IFracInBin=IBulgeMinBin/TotMinBin
    MFracInBin=MBulgeMinBin/TotMinBin
    FracInBin=BulgeMinBin/TotMinBin
    Ierr=np.sqrt(IFracInBin*(1-IFracInBin)/TotNinBin)
    Merr=np.sqrt(MFracInBin*(1-MFracInBin)/TotNinBin)
    err=np.sqrt(FracInBin*(1-FracInBin)/TotNinBin)
    axes.errorbar(BinStart+0.5*BinWidth,IFracInBin,Ierr,ls='None',label=label+' - Instability-driven Bulge',marker='*')
    axes.errorbar(BinStart+0.5*BinWidth,MFracInBin,Merr,ls='None',label=label+' - Merger-driven Bulge',marker='s')
    axes.errorbar(BinStart+0.5*BinWidth,FracInBin,err,ls='None',label=label+' - Total Bulge',marker='^')
    axes.set_ylabel('Stellar Mass Fraction of Components')
  else:
    FracInBin=BulgeMinBin/TotMinBin
    err=np.sqrt(FracInBin*(1-FracInBin)/TotNinBin)
    if bulge_frac:
      axes.errorbar(BinStart+0.5*BinWidth,FracInBin,err,label=label)
      axes.set_ylabel('Bulge Stellar Mass Fraction')
    else:
      axes.errorbar(BinStart+0.5*BinWidth,1-FracInBin,err,label=label)
      axes.set_ylabel('Disk Stellar Mass Fraction')
  axes.plot([9,12.5],[0.5,0.5],'k--',label='_nolegend_')
  #axes.plot(BinStart+0.5*BinWidth,1-FracInBin,label=label)
  axes.set_xlabel(r'$\log(M_\ast)$')
  axes.set_ylim(0,1)
  axes.set_xlim(MinMass,MaxMass)


def plot_disk_frac_SDSS(axes,bulge_frac):
  M=[8.90326,9.102003,9.300743,9.501402,9.703978,9.902666,10.101332,10.301882,10.502413,10.699095,10.901584,11.100241,11.300841,11.499548,11.702157,11.898974]
  F=[0.8897018,0.88715476,0.88215226,0.8673269,0.8414509,0.7910218,0.7209487,0.6103591,0.48380876,0.3658544,0.26385814,0.18519081,0.120027825,0.085559405,0.08914936,0.09028649]
  if bulge_frac:
    axes.plot(M,1-np.array(F),'--',label='Thanjavur et al. (2016)')
  else:
    axes.plot(M,F,'--',label='Thanjavur et al. (2016)')


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
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  prop='StellarMass'
  
  filename='bulges_IDBH_tune_t125'
  split=0
  bulge_frac=True
  plot_mags=False

  fig_frac,axes_frac=plt.subplots(1,1)
  for snapshot in [63,78,100,116,134,158,213]:
    gals_bulges=load_data(filename,snapshot,prop,cosmo,split)
    plot_disk_frac(gals_bulges,axes_frac,'$z={}$'.format(redshift[snapshot]),split,bulge_frac)
    
    if plot_mags:
      mag_def=load_mags(filename,snapshot)#['i775']
      plot_disk_frac(gals_bulges[(mag_def<21.9)],axes_frac,'$z=0.55$ Model Galaxies, $m_{Y}<-21.9$',split,bulge_frac)

  plot_disk_frac_SDSS(axes_frac,bulge_frac)
  plt.tight_layout()
  plt.legend()
  #lgd=axes_frac.legend(fontsize='small',loc='upper center', bbox_to_anchor=(0.5, -0.2))
  
  #plt.savefig('DiskFraction.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
