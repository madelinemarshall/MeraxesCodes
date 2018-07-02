#Sets plot defaults
import matplotlib
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import ContourPlot as cp
import pandas as pd
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.5,3.5)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
color         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
# red,blue,green,purple,orange,brown,pink,mint


def plot_hist(x,y,axes):
  xlims=[7,12]
  ylims=[-2.1,1.9]#[np.nanmin(y),np.nanmax(y)]
  H, xedges, yedges, img=axes.hist2d(x, y, bins=20, range=[xlims,ylims], cmin=1, cmap='BuPu',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu')
  #plt.colorbar(im,ax=axes,use_gridspec=True)
  axes.set_aspect('auto')


def load_data(filename,snapshot):
  #Setup
  cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  }
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,props=['StellarMass',\
                                          'BulgeStellarMass','StellarDiskScaleLength',\
                                          'GhostFlag','Type'],\
                                          h=cosmo['h'],quiet=True)
  gals=gals[gals["GhostFlag"]==0]
  #gals=gals[gals["Type"]==0]
  gals=gals[gals['StellarMass']*1e10>1e6]
  return gals
 

def load_mags(filename,snapshot): 
  #no_dust=pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')
  #dust=pd.read_csv('/home/mmarshal/PhD/results/mags_output/'+'bulges_update1102_full'+'/mags_6_'+format(snapshot,'03d')+'_dust.txt',sep=' ')
  #i775=np.array(dust)[:,0]
  #m1600=np.array(dust)[:,1]
  #print(np.median(no_dust['i775']-i775))
  return pd.read_hdf('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600']
  #return m1600


def plot_gadotti(axes):
  dat=pd.read_csv('gadotti_data.csv',header=None)
  rad=np.log10(dat[0])
  half_mass_rad=rad*1.67835
  diskmass=np.log10(dat[4])
  axes.plot(diskmass[dat[2]<0.3],half_mass_rad[dat[2]<0.3],'o',color=color[6],label='Gadotti et al. (2009)',markersize=2)


def plot_obs(axes):
#  logM=[8.751655,8.957122,9.151583,9.353381,9.551511,9.753309,9.951439,10.153237,10.351367,10.549497,10.754964,10.953094,11.154892,11.35669,11.551151,11.75295]
#  R=[1.3738238,1.653449,1.8873918,2.0433598,2.2122164,2.3713737,2.5419817,2.724864,2.8634958,3.1622777,3.597779,4.1753187,4.8455696,5.793367,6.995642,8.363996]
#  axes.plot(logM,np.log10(R)-np.log10(1.67835),'k.',label='Shen et al. (2003)')

  #Dutton+10
  logM=np.linspace(9,12)
  MonM0=10**(logM-10.44)
  R=10**0.72*MonM0**0.18*(0.5+0.5*MonM0**1.8)**((0.52-0.18)/1.8)
  scat=0.27+(0.47-0.27)/(1+MonM0**2.2) 
  axes.plot(logM,np.log10(R),color=color[0],label='Dutton et al. (2010)',linewidth=2.5)
  axes.plot(logM,np.log10(R)+scat,':',color=color[0],label='__nolabel__',linewidth=2.5)
  axes.plot(logM,np.log10(R)-scat,':',color=color[0],label='__nolabel',linewidth=2.5)

  ##Lange+16. R=a(M/10^10)^b
  #For disks, a=5.141 kpc, b=0.274
  #M=np.logspace(7,12,base=10)
  M=10**np.array([9.15,9.45,9.75,10.05,10.35,10.65,10.95])
  R=5.141*(M/1e10)**0.274
  scatter=np.array([0.185,0.176,0.17,0.151,0.133,0.15,0.202])
  axes.plot(np.log10(M),np.log10(R)-np.log10(1.67835),color=color[4],linewidth=2.5,label="Lange et al. (2016)",zorder=100)
  axes.plot(np.log10(M),np.log10(R)-np.log10(1.67835)+scatter,':',color=color[4],linewidth=2.5,label="__nolabel__",zorder=101)
  axes.plot(np.log10(M),np.log10(R)-np.log10(1.67835)-scatter,':',color=color[4],linewidth=2.5,label="__nolabel__",zorder=102)

  #Wu+17
  axes.plot([7.25,11.25],0.321*(np.array([7.25,11.25])-10)+0.343,color=color[2],label='Wu (2017)',linewidth=2.5)
  axes.plot([7.25,11.25],0.321*(np.array([7.25,11.25])-10)+0.343+0.36,':',color=color[2],label='__nolabel__',linewidth=2.5)
  axes.plot([7.25,11.25],0.321*(np.array([7.25,11.25])-10)+0.343-0.36,':',color=color[2],label='__nolabel__',linewidth=2.5)
  #axes.plot([7.25,11.25],0.321*(np.array([7.25,11.25])-10)+0.343+0.36,':k',label='Wu (2017) 2$\sigma$ scatter')
  #axes.plot([7.25,11.25],0.321*(np.array([7.25,11.25])-10)+0.343-0.36,':k',label='__nolegend__')

  
 

if __name__=='__main__':
  filename='tuned_reion_T125'
  snapshot=250
  gals=load_data(filename,snapshot)

  disks=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
  disk_mass=np.log10((disks['StellarMass']-disks['BulgeStellarMass'])*1e10)
  scale_rad=np.log10(disks['StellarDiskScaleLength']*1000)
  ### scale radius = half_mass_rad/1.67835 -> log(Rscale)=log(halfmass)-log(1.67835)

  #cp.contour_plot(disk_mass,half_mass_rad,'log(Stellar Disk Mass)','log(Half-Mass Radius/kpc)',[7.5,11.5],[-0.7,1.5])
  #plt.show()
  #fig,axes=plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0})
  fig,axes=plt.subplots(1,1)
  plot_hist(np.log10((disks['StellarMass'])*1e10),scale_rad,axes)
  #sc=axes[0].scatter(np.log10((disks['StellarMass'])*1e10),half_mass_rad,c=np.log10(disks['BulgeStellarMass']*1e10),s=0.1)
  #plt.colorbar(sc)
  #axes[0].plot(np.log10((disks['StellarMass'])*1e10),half_mass_rad,'.')
  plot_gadotti(axes)
  plot_obs(axes)
  axes.set_xlabel(r'$\log M_{\ast\mathrm{, total}}$')
  axes.set_ylabel(r'$\log R_{\mathrm{scale}}$ (kpc)')

#  plot_hist(disk_mass,scale_rad,axes[0])
#  ##Gadotti+09, disks
#  #axes[1].plot([7.04367,11.9312],np.log10([0.554176,6.20772]),'k-',linewidth=2.5,label="SDSS disks, best fit -\nGadotti et al. (2009)")
#  axes[1].set_xlabel(r'$\log M_{\ast\mathrm{, disk}}$')
#  axes[1].set_yticklabels([])
#  axes[1].legend(loc='lower right')
  plt.xlim([7.1,11.7])
  plt.ylim([-1.55,1.8])
  lgd=axes.legend(fontsize='small',loc='upper center', bbox_to_anchor=(0.5, -0.2))
  plt.savefig('/home/mmarshal/results/plots/DiskSize.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
  

##z=2 magnitude cut
#  filename='bulges_correctBHMF'
#  snapshot=158
#  gals=load_data(filename,snapshot)
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
