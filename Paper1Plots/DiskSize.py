import matplotlib
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/mmarshal/simulation_codes')
import ContourPlot as cp
import pandas as pd
from _load_data import load_data

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.5,3.5)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
color         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
# red,blue,green,purple,orange,brown,pink,mint

def plot_hist2d(xdata,ydata,axes,cmax=None,cbar=False):
  xlims=[7,12]
  ylims=[-2.1,1.9]#[np.nanmin(y),np.nanmax(y)]
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, vmin=1, vmax=cmax, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  if cbar:
    cb=plt.colorbar(img, cax=cbar,use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))
 

def load_mags(filename,snapshot): 
  return pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600']


def plot_gadotti(axes):
  dat=pd.read_csv('/home/mmarshal/simulation_codes/data/gadotti_data.csv',header=None)
  rad=np.log10(dat[0]) #scalelength
  #half_mass_rad=rad*1.67835
  diskmass=np.log10(dat[4])
  axes.plot(diskmass[dat[2]<0.3],rad[dat[2]<0.3],'o',color=color[4],label='Gadotti et al. (2009)',markersize=2)


def plot_obs(axes):
  #Dutton+10 - half-light radius
  logM=np.linspace(9,12)
  MonM0=10**(logM-10.44)
  R=10**0.72*MonM0**0.18*(0.5+0.5*MonM0**1.8)**((0.52-0.18)/1.8)
  scat=0.27+(0.47-0.27)/(1+MonM0**2.2) 
  axes.plot(logM,np.log10(R),color=color[0],label='Dutton et al. (2010)',linewidth=2.5)
  axes.fill_between(logM,np.log10(R)-scat,np.log10(R)+scat,color=color[0],label='__nolabel__',alpha=0.3)
  #axes.plot(logM,np.log10(R)-np.log10(1.67835),color=color[0],label='Dutton et al. (2010)',linewidth=2.5)
  #axes.fill_between(logM,np.log10(R)-np.log10(1.67835)-scat,np.log10(R)-np.log10(1.67835)+scat,color=color[0],label='__nolabel__',alpha=0.3)

  ##Lange+16. R=a(M/10^10)^b - half-light radius
  #For disks, a=5.141 kpc, b=0.274
  M=10**np.array([9.15,9.45,9.75,10.05,10.35,10.65,10.95])
  R=5.141*(M/1e10)**0.274
  scatter=np.array([0.185,0.176,0.17,0.151,0.133,0.15,0.202])
  #axes.plot(np.log10(M),np.log10(R)-np.log10(1.67835),color='k',linewidth=2.5,label="Lange et al. (2016)",zorder=100)
  #axes.fill_between(np.log10(M),np.log10(R)-np.log10(1.67835)-scatter,np.log10(R)-np.log10(1.67835)+scatter,alpha=0.3,color='k',label="__nolabel__",zorder=101)
  axes.plot(np.log10(M),np.log10(R),color='k',linewidth=2.5,label="Lange et al. (2016)",zorder=100)
  axes.fill_between(np.log10(M),np.log10(R)-scatter,np.log10(R)+scatter,alpha=0.3,color='k',label="__nolabel__",zorder=101)
  
  #Wu+17 - scalelength
  axes.plot([7.25,11.25],0.321*(np.array([7.25,11.25])-10)+0.343+np.log10(1.67835),color=color[2],label='Wu (2017)',linewidth=2.5)
  axes.fill_between([7.25,11.25],0.321*(np.array([7.25,11.25])-10)+0.343-0.36+np.log10(1.67835),0.321*(np.array([7.25,11.25])-10)+0.343+0.36+np.log10(1.67835),color=color[2],label='__nolabel__',alpha=0.4)

  #Lapi+18 - half-mass
  logM=np.linspace(9,11.5)
  logRe=0.7519 + 0.2333*(logM-10.5) +0.0494*(logM-10.5)**2+0.0267*(logM-10.5)**3
  #axes.plot(logM,logRe-np.log10(1.67835),color=color[3],label='Lapi et al. (2018)',linewidth=2.5)
  axes.plot(logM,logRe,color=color[3],label='Lapi et al. (2018)',linewidth=2.5)
 

if __name__=='__main__':
  filename='paper1_T125'
  snapshot=250
  gals=load_data(filename,snapshot,['StellarMass',\
    'BulgeStellarMass','StellarDiskScaleLength','GhostFlag','Type'],centrals=True)

  disks=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
  disk_mass=np.log10((disks['StellarMass']-disks['BulgeStellarMass'])*1e10)
  scale_rad=np.log10(disks['StellarDiskScaleLength']*1000*1.67835)
  ### scale radius = half_mass_rad/1.67835 -> log(Rscale)=log(halfmass)-log(1.67835)
  fig,axes=plt.subplots(1,2,gridspec_kw ={'wspace':0,'width_ratios':[4,0.4]})
  plot_hist2d(np.log10((disks['StellarMass'])*1e10),scale_rad,axes[0],cbar=axes[1])
  plot_gadotti(axes[0])
  plot_obs(axes[0])
  axes[0].plot([8,8],[-3,2],'--',color=[0.5,0.5,0.5],label='Convergence Limit')
  axes[0].set_xlabel(r'$\log (M_{\ast}/M_\odot)$')
  axes[0].set_ylabel(r'$\log(R_e/\mathrm{kpc})$')

  #plt.xlim([7.1,11.7])
  #plt.ylim([-1.55,1.8])
  lgd=axes[0].legend(fontsize='small',loc='upper center', bbox_to_anchor=(0.52, -0.2),ncol=2)
  plt.tight_layout()
  plt.savefig('/home/mmarshal/results/plots/Paper1/DiskSize.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
  

##z=2 magnitude cut
#  filename='bulges_correctBHMF'
#  snapshot=158
#  gals=load_data(filename,snapshot,['StellarMass',\
#    'BulgeStellarMass','StellarDiskScaleLength','GhostFlag','Type'])
#  mags=np.array(load_mags(filename,snapshot))
#  mags=mags[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
#  disks=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
#  disk_mass=np.log10((disks['StellarMass']-disks['BulgeStellarMass'])*1e10)
#  half_mass_rad=np.log10(1.67835*disks['StellarDiskScaleLength']*1000)
#  rad_kpc=(disks['StellarDiskScaleLength']*1000)
#
#  mu=mags+5*np.log10(rad_kpc)+1.344
#
#  fig,axes=plt.subplots(1,2)
#  sc=axes[0].scatter(np.log10((disks['StellarMass'])*1e10),np.log10(disks['StellarDiskScaleLength']*1000),c=mu,s=0.1)
#  plt.colorbar(sc)
#  #plot_hist(np.log10((disks['StellarMass'])*1e10),np.log10(disks['StellarDiskScaleLength']*1000),axes[0])
#  axes[0].set_xlabel(r'$\log M_{\ast\mathrm{\ast}}$')
#  axes[0].set_ylabel(r'$\log R_{\mathrm{scale}}$ (kpc)')
#  plot_obs(axes[0])
#  axes[0].legend(loc='lower right')

#  mags=np.array(load_mags(filename,snapshot))
#  gals=gals[mags<-18]
#  mags=mags[mags<-18]
#  disks=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
#  mags=mags[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
#  disk_mass=np.log10((disks['StellarMass']-disks['BulgeStellarMass'])*1e10)
#  plot_hist(np.log10((disks['StellarMass'])*1e10),np.log10(disks['StellarDiskScaleLength']*1000),axes[1])
#  plot_obs(axes[1])
#  axes[1].set_xlabel(r'$\log M_{\ast\mathrm{\ast}}$')
#  axes[1].set_ylabel(r'$\log R_{\mathrm{scale}}$ (kpc)')
#  axes[1].legend(loc='lower right')
  
#  plt.tight_layout()
#  plt.show()
  #fig.savefig('DiskSize.pdf', format='pdf')
