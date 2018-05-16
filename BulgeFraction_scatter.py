import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/PhD/simulation_codes/Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF
import magcalc as mc

#Sets plot defaults
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (3.5,3)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
#plt.rc('text', usetex=True)
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

def plot_SMF(gals,prop,boxwidth,axes,**kwargs):
  if prop=='DiskStellarMass':
    dsm=gals['StellarMass'][gals['StellarMass']>0]-gals['BulgeStellarMass'][gals['StellarMass']>0]
    dsm=dsm[dsm>0]
    maxval=np.nanmax(np.log10(dsm*1e10)) 
    minval=np.nanmin(np.log10(dsm*1e10))
    hist, bin_edges = np.histogram(np.log10(dsm*1e10),range=(minval,maxval),bins=30)
   
  else:
    maxval=np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10)) 
    minval=np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=30)
  bin_edges=np.array(bin_edges, dtype=np.float128)
  Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
  axes.set_xlabel(r'$\log(M_\ast (M_\odot)$)')
  axes.set_xlim([8.9,13])
  axes.set_ylim([-6,-1.2])
  axes.plot(Max,np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3),**kwargs)
  return axes



def plot_schechter(phi_star,alpha,logM,M_star,axes,label,color):
  phi= np.log(10)*phi_star*(10**((alpha+1)*(logM-M_star)))*(np.exp(-10**(logM-M_star)))*1#dM  
  axes.plot(logM,np.log10(phi),'--',label=label,color=color)



def plot_disk_frac(gals,axes,label,split,bulge_frac):
  #Figure 9
  MaxMass=13
  MinMass=9.0
  BinWidth=0.2
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
  SM=gals['StellarMass']*1e10
  BSM=gals['BulgeStellarMass']*1e10
  if split:
    MBSM=gals['MergerBulgeStellarMass']*1e10
    IBSM=BSM-MBSM
  for ii in range(0,nBins):
    BinStart[ii]=MinMass+ii*BinWidth
    BinEnd=MinMass+(ii+1)*BinWidth
    if (np.size(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])>0):
      BFracInBin[ii,:]=np.percentile(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]/SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)],[16,50,84])
      if split:
        IFracInBin[ii,:]=np.percentile(IBSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]/SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)],[16,50,84])
        MFracInBin[ii,:]=np.percentile(MBSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]/SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)],[16,50,84])

  if split:
    axes.errorbar(BinStart+0.5*BinWidth,IFracInBin[:,1],[IFracInBin[:,1]-IFracInBin[:,0],IFracInBin[:,2]-IFracInBin[:,1]],label=label+' - Instability-driven Bulge')
    axes.errorbar(BinStart+0.5*BinWidth,MFracInBin[:,1],[MFracInBin[:,1]-MFracInBin[:,0],MFracInBin[:,2]-MFracInBin[:,1]],label=label+' - Merger-driven Bulge')
    axes.errorbar(BinStart+0.5*BinWidth,BFracInBin[:,1],[BFracInBin[:,1]-BFracInBin[:,0],BFracInBin[:,2]-BFracInBin[:,1]],label=label+' - Total Bulge')
    axes.set_ylabel('Stellar Mass Fraction of Components')
  else:
    FracInBin=BulgeMinBin/TotMinBin
    err=np.sqrt(FracInBin*(1-FracInBin)/TotNinBin)
    if bulge_frac:
      axes.errorbar(BinStart+0.5*BinWidth,FracInBin,err,ls='None',label=label,marker='o')
      axes.set_ylabel('Bulge Stellar Mass Fraction')
    else:
      axes.errorbar(BinStart+0.5*BinWidth,1-FracInBin,err,ls='None',label=label,marker='o')
      axes.set_ylabel('Disk Stellar Mass Fraction')
  #print(err)
  #print(1-FracInBin)
  axes.plot([9,14],[0.5,0.5],'k--',label='_nolegend_')
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
  #meraxes_loc='/output/'+str(sys.argv[1])+'.hdf5'
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  prop='StellarMass'
  
  filename='tuned_quasarfeedback'#str(sys.argv[1])#'bulges_IDBHparam_tune'
  split=1
  bulge_frac=True
  plot_SMFs=False
  plot_mags=False

  if plot_SMFs:
    fig, axes = plt.subplots(1, 3,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  for snapshot in [213]:
    ii+=1
    #  gals_default=load_data('default',snapshot,prop,cosmo)
    gals_bulges=load_data(filename,snapshot,prop,cosmo,split)
    #gals_bulges=gals_bulges[(gals_bulges['StellarMass']*1e10>10**(8.9))&(mag_def<25)] 
    BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
    if plot_SMFs:
      BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
      #if snapshot!=116:
     # plot_obsGSMF(axes[j,ii],redshift[snapshot],hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color=[0.5,0.5,0.5],alpha=1.0)
      #axes[j,ii].legend()
      logM=np.linspace(8,13)
      plot_schechter(0.332*1e-3,-1.524,logM,10.740,axes[ii],label=r'$B/T<0.2$',color='C0')
      plot_schechter(0.489*1e-3,-1.270,logM,10.988,axes[ii],label=r'$0.2<B/T<0.5$',color='C1')
      plot_schechter(1.016*1e-3,-0.764,logM,10.926,axes[1],label=r'$0.5<B/T<0.8$',color='C2')
      plot_schechter(1.203*1e-3,-0.580,logM,11.023,axes[1],label=r'$0.8<B/T$',color='C3')

      plot_schechter(1.911*1e-3,-1.145,logM,11.116,axes[2],label=r'SDSS - Entire galaxy',color='C0')
      plot_schechter(1.424*1e-3,-1.073,logM,11.085,axes[2],label=r'SDSS - Bulge component',color='C1')
      plot_schechter(1.582*1e-3,-1.277,logM,10.707,axes[2],label=r'SDSS - Disk component',color='C2')

      plot_SMF(gals_bulges[BT<0.2],prop,100,axes[ii],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C0','zorder':100})
      plot_SMF(gals_bulges[(BT>0.2)&(BT<0.5)],prop,100,axes[ii],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C1','zorder':100})
      plot_SMF(gals_bulges[(BT>0.5)&(BT<0.8)],prop,100,axes[1],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C2','zorder':100})
      plot_SMF(gals_bulges[BT>0.8],prop,100,axes[1],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C3','zorder':100})
      plot_SMF(gals_bulges,prop,100,axes[2],**{'linestyle':'-','label':'All','linewidth':2.5,'color':'C0','zorder':100})
      plot_SMF(gals_bulges,'BulgeStellarMass',100,axes[2],**{'linestyle':'-','label':'Bulge components','linewidth':2.5,'color':'C1','zorder':100})
      plot_SMF(gals_bulges,'DiskStellarMass',100,axes[2],**{'linestyle':'-','label':'Disk components','linewidth':2.5,'color':'C2','zorder':100})

      #if snapshot==158:
      axes[0].legend()
      axes[1].legend()
      axes[2].legend()
      axes[0].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')

      axes[1].set_yticklabels([])
      axes[2].set_yticklabels([])
      #axes[j,ii].set_title('$z=${}'.format(redshift[snapshot]))
      axes[0].set_title('$B/T<0.5$')
      axes[1].set_title('$B/T>0.5$')
      axes[2].set_xlabel(r'$\log(M_{component} (M_\odot)$)')
    
    fig_frac,axes_frac=plt.subplots(1,1)
    plot_disk_frac(gals_bulges,axes_frac,'$z=0.55$ Model Galaxies',split,bulge_frac)
    if plot_mags:
      mag_def=load_mags(filename,snapshot)#['i775']
      plot_disk_frac(gals_bulges[(mag_def<21.9)],axes_frac,'$z=0.55$ Model Galaxies, $m_{Y}<-21.9$',split,bulge_frac)

    ##TIAMAT 125, boxwidth=125/cosmo['h']
  #plt.tight_layout()
  #plt.show()

  #print(mag_def)
  plot_disk_frac_SDSS(axes_frac,bulge_frac)
  #plt.title(r'{}'.format(sys.argv[1]))
  plt.tight_layout()
  lgd=axes_frac.legend(fontsize='small',loc='upper center', bbox_to_anchor=(0.5, -0.2))
  
  #plt.savefig('DiskFraction.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  #plt.savefig('BulgeFrac.pdf',format='pdf')
  plt.show()
