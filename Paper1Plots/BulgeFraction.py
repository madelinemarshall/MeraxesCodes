import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
import magcalc as mc
from _load_data import load_data
from astrodatapy.number_density import number_density

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.5,4)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#999999']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','--','-.',':']*4


def load_mags(filename,snapshot):
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  M6300 = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['F625W']
  M1600   = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  M6300_dust = M6300 + mc.reddening(6300., M1600, z = redshift[snapshot])
  return M6300_dust


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
  Max=Max[hist>5]
  hist=hist[hist>5]
  hist_plus=hist+np.sqrt(hist)
  hist_minus=hist-np.sqrt(hist)
  phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3)
  phi_plus=np.log10(hist_plus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
  phi_minus=np.log10(hist_minus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
  axes.plot(Max,phi,**kwargs)
  axes.fill_between(Max, phi_plus, phi_minus, alpha=0.5,color=kwargs['color'],label='__nolegend__')
  axes.set_xlabel(r'$\log(M_\ast (M_\odot)$)')
  axes.set_xlim([8,12.5])
  axes.set_ylim([-6,-1])
  return axes


def plot_obs_SMF(prop,ax):
  xlim     = (8, 12.5)
  ylim     = (-6, -1)
  xlabel   = r"$\log_{10}[M_*/{\rm M_{\odot}}]$"
  ylabel   = r"$\log_{10}[\rm \phi/Mpc^{-3} dex^{-1}]$"
  zs       = [2.0, 0.0]

  for z in [0.0,]:
    obs    = number_density(feature=prop,z_target=z,quiet=1,h=cosmo['h'])
    for ii in range(obs.n_target_observation):
      data       = obs.target_observation['Data'][ii]
      label      = obs.target_observation.index[ii]
      datatype   = obs.target_observation['DataType'][ii]
      marker     = markers[ii]
      data[:,1:] = np.log10(data[:,1:])
      ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=[0.5,0.5,0.5],fmt=marker)

    ax.text(0.95,0.95, "z=%.2f"%z,horizontalalignment='right',\
          verticalalignment='top',transform=ax.transAxes)
    leg = ax.legend(loc='lower left')
    leg.get_frame().set_alpha(0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def plot_Moffett(ax):
  Moffett_total=pd.read_csv('/home/mmarshal/simulation_codes/data/Moffett_total.csv',header=None)
  Moffett_disks=pd.read_csv('/home/mmarshal/simulation_codes/data/Moffett_spheroids.csv',header=None)
  Moffett_spheroids=pd.read_csv('/home/mmarshal/simulation_codes/data/Moffett_spheroids.csv',header=None)
  Moffett_spheroids_L=pd.read_csv('/home/mmarshal/simulation_codes/data/Moffett_bulges_lower.csv',header=None)
  Moffett_spheroids_U=pd.read_csv('/home/mmarshal/simulation_codes/data/Moffett_bulges_upper.csv',header=None)
  sph_err_U=np.log10(Moffett_spheroids_U[1])-np.log10(Moffett_spheroids[1])
  sph_err_L=np.log10(Moffett_spheroids[1])-np.log10(Moffett_spheroids_L[1])
  #ax[0].plot(Moffett_disks[0],np.log10(Moffett_disks[1]),'^',label='Moffett et al. (2016)',color=[0.5,0.5,0.5])
  #ax[1].errorbar(Moffett_spheroids[0],np.log10(Moffett_spheroids[1]),yerr=[sph_err_L,sph_err_U],label='Moffett et al. (2016)',fmt='^',color=[0.5,0.5,0.5])
  ax.plot(Moffett_total[0],np.log10(Moffett_total[1]),'^',label='Moffett et al. (2016)',color=[0.5,0.5,0.5])



def plot_schechter(phi_star,alpha,logM,M_star,axes,color):
  phi= np.log(10)*phi_star*(10**((alpha+1)*(logM-M_star)))*(np.exp(-10**(logM-M_star)))*1#dM  
  axes.plot(logM,np.log10(phi),'--',color=color)



def plot_disk_frac(gals,axes,label,split,bulge_frac,mags=False):
  #Figure 9
  MaxMass=11.8
  MinMass=9.4
  BinWidth=0.2 #Same as obs
  nBins=int(np.ceil((MaxMass-MinMass)/BinWidth))
  BMassInBin=np.zeros(nBins)
  IMassInBin=np.zeros(nBins)
  MMassInBin=np.zeros(nBins)
  TotMassInBin=np.zeros(nBins)
  TotNInBin=np.zeros(nBins)

  BinStart=np.zeros(nBins)
  SM=gals['StellarMass']*1e10
  BSM=gals['BulgeStellarMass']*1e10
  if split:
    MBSM=gals['MergerBulgeStellarMass']*1e10
    IBSM=BSM-MBSM
  for ii in range(0,nBins):
    BinStart[ii]=MinMass+ii*BinWidth
    BinEnd=MinMass+(ii+1)*BinWidth
    condition=(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)
    if (np.size(BSM[condition])>25):
      BMassInBin[ii]=np.sum(BSM[condition])
      TotMassInBin[ii]=np.sum(SM[condition])
      TotNInBin[ii]=np.size(SM[condition])
      if split:
        IMassInBin[ii]=np.sum(IBSM[condition])
        MMassInBin[ii]=np.sum(MBSM[condition])
  BFracInBin=BMassInBin/TotMassInBin  
  IFracInBin=IMassInBin/TotMassInBin  
  MFracInBin=MMassInBin/TotMassInBin  
  Berr=np.sqrt(BFracInBin*(1-BFracInBin)/TotNInBin) #Binomial errors
  Ierr=np.sqrt(IFracInBin*(1-IFracInBin)/TotNInBin)
  Merr=np.sqrt(MFracInBin*(1-MFracInBin)/TotNInBin)

  if bulge_frac:
    if mags:
      axes.errorbar(BinStart+0.5*BinWidth,BFracInBin,Berr,ls='None',label=label,marker='o',color='#984ea3')
    else:
      axes.plot(BinStart+0.5*BinWidth,BFracInBin,ls='None',label=label,marker='o',color='#984ea3',fillstyle='none')
    axes.set_ylabel(r'$M_{\rm{bulge}}/M_\ast$')
    if split:
      axes.errorbar(BinStart+0.5*BinWidth,IFracInBin,Ierr,ls='None',label=label+' - Instability-driven Bulge',marker='*')
      axes.errorbar(BinStart+0.5*BinWidth,MFracInBin,Merr,ls='None',label=label+' - Merger-driven Bulge',marker='s')
      axes.set_ylabel('Stellar Mass Fraction of Components')
  else:
    axes.errorbar(BinStart+0.5*BinWidth,1-BFracInBin,Berr,ls='None',label=label,marker='o',color='#984ea3')
    axes.set_ylabel('Disk Stellar Mass Fraction')
  axes.set_xlabel(r'$\log(M_\ast/M_\odot)$')
  axes.set_ylim(0,1)
  axes.set_xlim(MinMass,MaxMass)


def plot_disk_frac_SDSS(axes,bulge_frac):
  M=[8.90326,9.102003,9.300743,9.501402,9.703978,9.902666,10.101332,10.301882,10.502413,10.699095,10.901584,11.100241,11.300841,11.499548,11.702157,11.898974]
  F=[0.8897018,0.88715476,0.88215226,0.8673269,0.8414509,0.7910218,0.7209487,0.6103591,0.48380876,0.3658544,0.26385814,0.18519081,0.120027825,0.085559405,0.08914936,0.09028649]
  if bulge_frac:
    axes.plot(M-2*np.log10(cosmo['h']/0.7),1-np.array(F),'--',label='Thanjavur et al. (2016)',color='k',zorder=100)
  else:
    axes.plot(M-2*np.log10(cosmo['h']/0.7),F,'--',label='Thanjavur et al. (2016)',color='k',zorder=100)


def plot_SMFs(vol):
      BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
      logM=np.linspace(8,13)
      #plot_schechter(0.332*1e-3,-1.524,logM,10.740,axes[ii],label=r'$B/T<0.2$',color='C0')
      #plot_schechter(0.489*1e-3,-1.270,logM,10.988,axes[ii],label=r'$0.2<B/T<0.5$',color='C1')
      #plot_schechter(1.016*1e-3,-0.764,logM,10.926,axes[1],label=r'$0.5<B/T<0.8$',color='C2')
      #plot_schechter(1.203*1e-3,-0.580,logM,11.023,axes[1],label=r'$0.8<B/T$',color='C3')

      plot_schechter(1.911*1e-3,-1.145,logM,11.116,axes,color='C0')#label=r'Thanjavur et al. (2016) - Entire galaxy',color='C0')
      plot_schechter(1.424*1e-3,-1.073,logM,11.085,axes,color='C1')#label=r'Thanjavur et al. (2016) - Bulge component',color='C1')
      plot_schechter(1.582*1e-3,-1.277,logM,10.707,axes,color='C2')#label=r'Thanjavur et al. (2016) - Disk component',color='C2')

      #plot_SMF(gals_bulges[BT<0.2],prop,vol,axes[ii],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C0','zorder':100})
      #plot_SMF(gals_bulges[(BT>0.2)&(BT<0.5)],prop,vol,axes[ii],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C1','zorder':100})
      #plot_obs_SMF('GSMF_Disk',axes[ii])
      #plot_SMF(gals_bulges[(BT>0.5)&(BT<0.8)],prop,vol,axes[1],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C2','zorder':100})
      #plot_SMF(gals_bulges[BT>0.8],prop,vol,axes[1],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C3','zorder':100})
      #plot_obs_SMF('GSMF_Bulge',axes[1])
      plot_SMF(gals_bulges,prop,vol,axes,**{'linestyle':'-','label':'Entire galaxy','linewidth':2.5,'color':'C0','zorder':100})
      plot_SMF(gals_bulges,'BulgeStellarMass',vol,axes,**{'linestyle':'-','label':'Bulge component','linewidth':2.5,'color':'C1','zorder':100})
      plot_SMF(gals_bulges,'DiskStellarMass',vol,axes,**{'linestyle':'-','label':'Disk component','linewidth':2.5,'color':'C2','zorder':100})
      #plot_Moffett(axes)

      #axes[0].legend()
      #axes[1].legend()
      axes.set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')

      #axes[1].set_yticklabels([])
      #axes[2].set_yticklabels([])
      #axes[0].set_title('$B/T<0.5$')
      #axes[1].set_title('$B/T>0.5$')
      axes.set_xlabel(r'$\log(M_{\mathrm{component}} (M_\odot)$)')

      axes.plot(0,0,color='C0',linestyle='--',label='Thanjavur et al. (2016)',linewidth=2.5)
      axes.plot(0,0,color='C0',linestyle='-',label='M19',linewidth=2.5)
      axes.legend()
      
      ##legend


if __name__=="__main__":
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,242:0.1,250:0}
  prop='StellarMass'
  
  filename='paper1_T125'
  split=0
  bulge_frac=True
  plot_SMF_bool=True
  vol=125/cosmo['h']
  plot_mags=True

  if plot_SMF_bool:
    fig, axes = plt.subplots(1, 1,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  for snapshot in [242]:
    ii+=1
    gals_bulges=load_data(filename,snapshot,[prop,'GhostFlag','BulgeStellarMass','Type','MergerBulgeStellarMass'])
    selection=(gals_bulges['Type']==0)
    gals_bulges=gals_bulges[selection]
    BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
    if plot_SMF_bool:
      plot_SMFs(vol)

    fig_frac,axes_frac=plt.subplots(1,1)
    plot_disk_frac(gals_bulges,axes_frac,r'All $M_\ast>10^6M_\odot$ galaxies',split,bulge_frac)
    if plot_mags:
      mag_def=load_mags(filename,snapshot)
      mag_def=mag_def[selection]
      plot_disk_frac(gals_bulges[(mag_def<17.77)],axes_frac,r'$r<17.77$ galaxies',split,bulge_frac,True)
  
  plot_disk_frac_SDSS(axes_frac,bulge_frac)
  plt.tight_layout()
  lgd=axes_frac.legend(fontsize='small',loc='upper center', bbox_to_anchor=(0.5, -0.2))
  
  fig_frac.savefig('/home/mmarshal/results/plots/Paper1/BulgeFrac_z={}_points.pdf'.format(redshift[snapshot]),format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  fig.savefig('/home/mmarshal/results/plots/Paper1/SMF_morphology.pdf',format='pdf')
  plt.show()
