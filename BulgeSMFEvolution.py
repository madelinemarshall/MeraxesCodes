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
matplotlib.rcParams['figure.figsize'] = (7.2,4)
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
  elif prop=='InstabilityBulgeMass':
    selection=(gals['BulgeStellarMass']-gals['MergerBulgeStellarMass'])>0
    values=np.log10(gals[selection]['BulgeStellarMass']-gals[selection]['MergerBulgeStellarMass'])+10
    maxval=np.nanmax(values)
    minval=np.nanmin(values)
    hist, bin_edges = np.histogram(values,range=(minval,maxval),bins=30)
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
  axes.set_xlim([8.2,12])
  axes.set_ylim([-6,-1.05])
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


def plot_SMFs(vol,axes):
      BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
      logM=np.linspace(8,13)
      #plot_schechter(0.332*1e-3,-1.524,logM,10.740,axes[ii],label=r'$B/T<0.2$',color='C0')
      #plot_schechter(0.489*1e-3,-1.270,logM,10.988,axes[ii],label=r'$0.2<B/T<0.5$',color='C1')
      #plot_schechter(1.016*1e-3,-0.764,logM,10.926,axes[1],label=r'$0.5<B/T<0.8$',color='C2')
      #plot_schechter(1.203*1e-3,-0.580,logM,11.023,axes[1],label=r'$0.8<B/T$',color='C3')

      #plot_schechter(1.911*1e-3,-1.145,logM,11.116,axes,color='C0')#label=r'Thanjavur et al. (2016) - Entire galaxy',color='C0')
      #plot_schechter(1.424*1e-3,-1.073,logM,11.085,axes,color='C1')#label=r'Thanjavur et al. (2016) - Bulge component',color='C1')
      #plot_schechter(1.582*1e-3,-1.277,logM,10.707,axes,color='C2')#label=r'Thanjavur et al. (2016) - Disk component',color='C2')

      #plot_SMF(gals_bulges[BT<0.2],prop,vol,axes[ii],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C0','zorder':100})
      #plot_SMF(gals_bulges[(BT>0.2)&(BT<0.5)],prop,vol,axes[ii],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C1','zorder':100})
      #plot_obs_SMF('GSMF_Disk',axes[ii])
      #plot_SMF(gals_bulges[(BT>0.5)&(BT<0.8)],prop,vol,axes[1],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C2','zorder':100})
      #plot_SMF(gals_bulges[BT>0.8],prop,vol,axes[1],**{'linestyle':'-','label':'Bulge Model','linewidth':2.5,'color':'C3','zorder':100})
      #plot_obs_SMF('GSMF_Bulge',axes[1])
      plot_SMF(gals_bulges,prop,vol,axes,**{'linestyle':'-','label':'Total','linewidth':2.5,'color':colors[0],'zorder':100})
      plot_SMF(gals_bulges,'BulgeStellarMass',vol,axes,**{'linestyle':'--','label':'Total Bulge','linewidth':2.5,'color':colors[1],'zorder':100})
      plot_SMF(gals_bulges,'MergerBulgeStellarMass',vol,axes,**{'linestyle':':','label':'Merger-Driven Bulge','linewidth':2.5,'color':colors[2],'zorder':100})
      plot_SMF(gals_bulges,'InstabilityBulgeMass',vol,axes,**{'linestyle':':','label':'Instability-Driven Bulge','linewidth':2.5,'color':colors[3],'zorder':100})
      plot_SMF(gals_bulges,'DiskStellarMass',vol,axes,**{'linestyle':'-.','label':'Stellar Disc','linewidth':2.5,'color':colors[4],'zorder':100})
      #plot_Moffett(axes)

      #axes[0].legend()
      #axes[1].legend()

      #axes[1].set_yticklabels([])
      #axes[2].set_yticklabels([])
      #axes[0].set_title('$B/T<0.5$')
      #axes[1].set_title('$B/T>0.5$')

      #axes.legend()
      
      ##legend


if __name__=="__main__":
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55,242:0.1,250:0}
  snapshots=[63,78,100,116,134,158]
  prop='StellarMass'
  
  filename='draft2_reion'
  vol=125/cosmo['h']

  fig,axes=plt.subplots(3,3,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1.2]})
  ii=-1
  jj=0
  for snapshot in snapshots:
    ii+=1
    if ii==3:
      jj+=1
      ii=0
    gals_bulges=load_data(filename,snapshot,[prop,'GhostFlag','BulgeStellarMass','Type','MergerBulgeStellarMass'],centrals=True)
    BT=gals_bulges['BulgeStellarMass']/gals_bulges['StellarMass']
    
    plot_SMFs(vol,axes[jj,ii])
    if jj==1:
      axes[jj,ii].set_xlabel(r'$\log(M_{\mathrm{component}} /M_\odot$)')
    else:
      axes[jj,ii].set_xticklabels([])
    if ii==0:
      #axes[jj,ii].set_ylabel(r'$M_{\mathrm{Bulge}}/M_\ast$')
      axes[jj,ii].set_ylabel(r'$\log(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
    else:
      axes[jj,ii].set_yticklabels([])
    axes[jj,ii].text(11, -1.8, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
  for ii in [0,1,2]:
    axes[2,ii].axis('off')
 
  plt.tight_layout()
  axes[1,1].legend(loc=(-0.75,-0.6),ncol=3)
  plt.savefig('/home/mmarshal/results/plots/Paper1/BulgeSMFEvolution.pdf', format='pdf')
  plt.show()
