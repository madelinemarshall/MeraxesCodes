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
from _load_data import load_data

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,3.2)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
color=['#ff7f00','#4daf4a','#984ea3','#377eb8','#e41a1c','#e41a1c','#377eb8']#,\'#377eb8'
                  #134:,158:'#f781bf',194:'#a65628',213:'black'}

def load_mags(filename,snapshot):
  redshift={63:7,78:6,100:5,115:4,134:3,158:2,192:1,213:0.55,242:0.1,250:0}
  #MUV=pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  #AUV=mc.reddening(1600,MUV,redshift[snapshot])
  #MUV_dust=MUV+AUV
  M6300 = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['F625W']
  M1600   = pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  M6300_dust = M6300 + mc.reddening(6300., M1600, z = redshift[snapshot])
  return M6300_dust
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
  axes.set_xlim([7.5,12])#[8.9,12])
  axes.set_ylim([-6,-1.2])
  axes.plot(Max,np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3),**kwargs)
  return axes



def plot_schechter(phi_star,alpha,logM,M_star,axes,label,color):
  phi= np.log(10)*phi_star*(10**((alpha+1)*(logM-M_star)))*(np.exp(-10**(logM-M_star)))*1#dM  
  axes.plot(logM,np.log10(phi),'--',label=label,color=color)



def plot_disk_frac(gals,axes,label,split,bulge_frac,mags):
  MaxMass=11.8
  MinMass=7.2#8.8
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
  NumInBin=np.zeros(nBins)
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
    NumInBin[ii]=np.size(BSM[condition])
    if NumInBin[ii]>25:
      BinFull[ii]=True
      BFracInBin[ii,:]=np.percentile(BSM[condition]/SM[condition],[16,50,84])
      #print("Num in bin: {}, bin start: {}, bin end: {}".format(np.size(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)]),BinStart[ii],BinEnd))
      if split:
        IFracInBin[ii,:]=np.percentile(IBSM[condition]/SM[condition],[16,50,84])
        MFracInBin[ii,:]=np.percentile(MBSM[condition]/SM[condition],[16,50,84])


  if split:
    if  not mags:
      axes.plot(BinStart[BinFull]+0.5*BinWidth,IFracInBin[BinFull][:,1],label=label+r'Instability-driven Bulge',lw=2,color=color[0],linestyle='--')
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,IFracInBin[BinFull][:,0],IFracInBin[BinFull][:,2],alpha=0.2,color=color[0])
      var, =axes.plot(BinStart[BinFull]+0.5*BinWidth,MFracInBin[BinFull][:,1],label=label+r'Merger-driven Bulge',lw=2,color=color[1],linestyle='-.')
      var.set_dashes([2,2,5,2])
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,MFracInBin[BinFull][:,0],MFracInBin[BinFull][:,2],alpha=0.2,color=color[1])
      axes.plot(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,1],label=label+r'Total Bulge',lw=2,color=color[2],linestyle='-')
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,0],BFracInBin[BinFull][:,2],alpha=0.2,color=color[2])
    else:
      #axes.plot(BinStart+0.5*BinWidth,IFracInBin[:,1],label=label+' - Instability-driven Bulge',lw=2.5,color='firebrick',linestyle='--')
      #axes.fill_between(BinStart+0.5*BinWidth,IFracInBin[:,0],IFracInBin[:,2],alpha=0.2,color='firebrick')
      #axes.plot(BinStart+0.5*BinWidth,MFracInBin[:,1],label=label+' - Merger-driven Bulge',lw=2.5,color='gold',linestyle='--')
      #axes.fill_between(BinStart+0.5*BinWidth,MFracInBin[:,0],MFracInBin[:,2],alpha=0.2,color='gold')
      axes.plot(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,1],'-',label=label+r'Total Bulge',lw=2,color=color[2],zorder=100,markersize=3)
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,0],BFracInBin[BinFull][:,2],alpha=0.2,color=color[2])
      #axes.errorbar(BinStart+0.5*BinWidth,IFracInBin[:,1],[IFracInBin[:,1]-IFracInBin[:,0],IFracInBin[:,2]-IFracInBin[:,1]],label=label+' - Instability-driven Bulge',color='firebrick',linestyle='--')
      #axes.errorbar(BinStart+0.5*BinWidth,MFracInBin[:,1],[MFracInBin[:,1]-MFracInBin[:,0],MFracInBin[:,2]-MFracInBin[:,1]],label=label+' - Merger-driven Bulge',color='gold',linestyle='--')
      #axes.errorbar(BinStart+0.5*BinWidth,BFracInBin[:,1],[BFracInBin[:,1]-BFracInBin[:,0],BFracInBin[:,2]-BFracInBin[:,1]],label=label+' - Total Bulge',color='royalblue',linestyle='--')


  else:
    if  mags:
      axes.plot(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,1],label=label+r'Total Bulge',lw=2,color='k',linestyle='--',zorder=100)
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,0],BFracInBin[BinFull][:,2],alpha=0.1,color=color[2])
    else:
      axes.plot(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,1],label=label+r'Total Bulge',lw=2,color=color[2],linestyle='-')
      axes.fill_between(BinStart[BinFull]+0.5*BinWidth,BFracInBin[BinFull][:,0],BFracInBin[BinFull][:,2],alpha=0.05,color=color[2])
  #print(err)
  #print(1-FracInBin)
  #axes.plot([9,12],[0.5,0.5],'k:',label='_nolegend_')
  #axes.plot(BinStart+0.5*BinWidth,1-FracInBin,label=label)
  axes.set_xlabel(r'$\log(M_\ast)$')
  axes.set_ylim(0,1)
  axes.set_xlim(MinMass+0.5*BinWidth,MaxMass-0.5*BinWidth)


def plot_disk_frac_SDSS(axes,bulge_frac):
  M=[8.90326,9.102003,9.300743,9.501402,9.703978,9.902666,10.101332,10.301882,10.502413,10.699095,10.901584,11.100241,11.300841,11.499548,11.702157,11.898974]
  F=[0.8897018,0.88715476,0.88215226,0.8673269,0.8414509,0.7910218,0.7209487,0.6103591,0.48380876,0.3658544,0.26385814,0.18519081,0.120027825,0.085559405,0.08914936,0.09028649]
  if bulge_frac:
    #var, =axes.plot(M,1-np.array(F),'-.',label='Thanjavur et al. (2016)',color=color[3],lw=1.5)
    #var.set_dashes([2,2,5,2])
    axes.plot(M,1-np.array(F),'-o',label='Thanjavur et al. (2016)',color=[0.4,0.4,0.4],markersize=3,lw=2)
  else:
    axes.plot(M,F,'-.',label='Thanjavur et al. (2016)',color=color[3],lw=2)
    

def plot_bulge_frac_GAMMA(axes,bulge_frac):
  M=np.array([8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1,11.3,11.5])
  F=np.array([0.00,0.01,0.05,0.05,0.10,0.14,0.24,0.34,0.40,0.45,0.46,0.45,0.39,0.38,0.53,0.74,0.89,1.00])
  F_up=np.array([0.00,0.01,0.05,0.05,0.10,0.15,0.26,0.37,0.43,0.49,0.50,0.49,0.43,0.43,0.62,0.89,1.00,1.00])
  F_low=np.array([0.00,0.01,0.05,0.05,0.10,0.12,0.21,0.30,0.36,0.42,0.42,0.41,0.34,0.32,0.44,0.59,0.64,0.42])


  if bulge_frac:
    axes.errorbar(M,F,yerr=[F-F_low,F_up-F],marker='o',label='Moffett et al. (2016)',color=[0.4,0.4,0.4],linestyle='',markersize=3)
    #axes.plot(M,F,label='Moffett et al. (2016)',color=[0.5,0.5,0.5],linestyle='-.',linewidth=2)


def plot_bulge_frac_EAGLE(axes,bulge_frac):
  M=np.array([9.12,9.38,9.62,9.88,10.13,10.38,10.63,10.87,11.13,11.37,11.63])
  F=np.array([0.83,0.73,0.65,0.6,0.53,0.47,0.5,0.59,0.75,0.79,0.89])
  F_low=np.array([0.48,0.41,0.36,0.34,0.29,0.27,0.22,0.25,0.43,0.61,0.69])
  F_up=np.array([0.97,0.96,0.94,0.92,0.88,0.88,0.94,0.96,0.95,0.98,0.99])

  if bulge_frac:
    #axes.errorbar(M,F,yerr=[F-F_low,F_up-F],marker='*',label='Clauwens et al. (2018)',color='k')
    #axes.fill_between(M,F_low,F_up,color=[0.6,0.6,0.6],alpha=0.1)
    #axes.plot(M,F,label='Clauwens et al. (2018)',color='#00ced1',lw=2,linestyle='-.')
    axes.plot(M,F,label='Clauwens et al. (2018)',color=[0.4,0.4,0.4],lw=2.5,linestyle=(0, (1, 1)))


def plot_SMFs():
      BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
      #if snapshot!=116:
     # plot_obsGSMF(axes[j,ii],redshift[snapshot],hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color=[0.5,0.5,0.5],alpha=1.0)
      #axes[j,ii].legend()
      logM=np.linspace(8,12)
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
      axes[0].legend(handlelength=3)
      axes[1].legend(handlelength=3)
      axes[2].legend(handlelength=3)
      axes[0].set_ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')

      axes[1].set_yticklabels([])
      axes[2].set_yticklabels([])
      #axes[j,ii].set_title('$z=${}'.format(redshift[snapshot]))
      axes[0].set_title('$B/T<0.5$')
      axes[1].set_title('$B/T>0.5$')
      axes[2].set_xlabel(r'$\log(M_{component} (M_\odot)$)')


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
  
  filename='paper1_T125'#str(sys.argv[1])#'bulges_IDBHparam_tune'
  split=True
  bulge_frac=True
  plot_SMF=False
  plot_mags=1#True

  if plot_SMF:
    fig, axes = plt.subplots(1, 3,gridspec_kw = {'wspace':0, 'hspace':0})

  ii=-1
  for snapshot in [242]:
    ii+=1
    gals_bulges=load_data(filename,snapshot,[prop,'GhostFlag','BulgeStellarMass','Type','MergerBulgeStellarMass'])
    selection=(gals_bulges['Type']==0)
    gals_bulges=gals_bulges[selection]

    BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
    if plot_SMF:
      plot_SMFs()
    
    fig_frac,axes_frac=plt.subplots(2,3,gridspec_kw={'width_ratios':[4,4,4],'height_ratios':[4,2],'hspace':0,'wspace':0},sharey=True)
    #plot_disk_frac(gals_bulges,axes_frac[0],'',0,bulge_frac,0)
    plot_disk_frac(gals_bulges,axes_frac[0,0],'',1,bulge_frac,0)

    if plot_mags:
      mag_def=load_mags(filename,snapshot)#['i775']
      mag_def=mag_def[selection]
      plot_disk_frac(gals_bulges[(mag_def<17.77)],axes_frac[0,1],r'$r<17.77$ - ',split,bulge_frac,1)
      plot_disk_frac(gals_bulges[(mag_def<19.8)],axes_frac[0,2],r'$r<19.8$ - ',split,bulge_frac,1)
  
  plot_disk_frac_SDSS(axes_frac[0,1],bulge_frac)
  plot_bulge_frac_GAMMA(axes_frac[0,2],bulge_frac)
  plot_bulge_frac_EAGLE(axes_frac[0,0],bulge_frac)

  axes_frac[0,0].plot([8.65,8.65],[-3,2],'--',color=[0.5,0.5,0.5],label='Convergence Limit',zorder=200)
  axes_frac[0,1].plot([8.65,8.65],[-3,2],'--',color=[0.5,0.5,0.5],label='Convergence Limit',zorder=200)
  axes_frac[0,2].plot([8.65,8.65],[-3,2],'--',color=[0.5,0.5,0.5],label='Convergence Limit',zorder=200)
  axes_frac[0,0].set_ylabel(r'$M_{\textrm{bulge}}/M_\ast$')

  axes_frac[0,1].set_xlim(7.9,11.7)
  axes_frac[0,2].set_xlim(7.9,11.7)

  plt.tight_layout()
  axes_frac[1,0].axis('off')
  axes_frac[1,1].axis('off')
  axes_frac[1,2].axis('off')
  lgd=axes_frac[0,0].legend(fontsize='small',loc=(0.14, -0.7))
  lgd2=axes_frac[0,1].legend(fontsize='small',loc=(0.14, -0.6))
  lgd3=axes_frac[0,2].legend(fontsize='small',loc=(0.14, -0.6))
  
  if split:
    plt.savefig('/home/mmarshal/results/plots/Paper1/BulgeFrac_z={}.pdf'.format(redshift[snapshot]), format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  else:
    plt.savefig('/home/mmarshal/results/plots/Paper1/BulgeFrac_z={}_nosplit.pdf'.format(redshift[snapshot]), format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  #plt.savefig('/home/mmarshal/results/plots/BulgeFrac.pdf',format='pdf')
  plt.show()
