##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import ContourPlot as cp
#Sets plot defaults
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,2.8)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename):
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

  snapshot=173

  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
                                          snapshot=snapshot,\
                                          h=cosmo['h'],quiet=True)
  return(gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)])#&(gals["Sfr"]>0)])



def plot_hist2d(gals,axes,cbar=True,move_cbar=False):
  xlims=[7,11]
  ylims=[-6,3]
  H, xedges, yedges, img=axes.hist2d(np.log10(gals['StellarMass']*1e10), np.log10(gals['Sfr']), bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=5e5, data=None,cmap='BuPu',norm=matplotlib.colors.LogNorm(),vmax=5e5)
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu',norm=matplotlib.colors.LogNorm())
  if move_cbar:
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.511, 0.225, 0.025, 0.693])
    cb=fig.colorbar(im, cax=cbar_ax,use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))
  elif cbar:
    cb=plt.colorbar(im,ax=axes,use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))
  axes.set_aspect('auto')
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  #cp.contour_plot(np.log10(gals['StellarMass']*1e10), np.log10(gals['Sfr']),xlab=None,ylab=None,xlims=xlims,ylims=ylims,axes=axes,colors=None,levels=None,linewidth=2.5)


def plot_hist3d(gals):
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  xlims=[6,11]
  ylims=[-7,3]
  hist, xedges, yedges = np.histogram2d(np.log10(gals['StellarMass']*1e10), np.log10(gals['Sfr']), bins=100, range=[xlims,ylims])
  
  xpos, ypos = np.meshgrid(xedges[:-1] + (11-6)/41, yedges[:-1] + (3+8)/40)
  xpos = xpos.flatten('F')
  ypos = ypos.flatten('F')
  zpos = np.zeros_like(xpos)

  # Construct arrays with the dimensions for the 16 bars.
  dx = 0.5 * np.ones_like(zpos)
  dy = dx.copy()
  dz = hist.flatten()

  ax.bar3d(xpos, ypos, zpos, dx, dy, np.log10(dz), color='purple', zsort='average')
  for angle in range(0, 360):
    ax.view_init(30, angle)
    plt.draw()
    plt.pause(.0001)
  plt.show()



def plot_obs(axes):
  M=np.array([9.356677,9.659609,9.959284,10.2557,10.555374,10.858306,11.157981])
  M_low=np.array([9.2,9.5,9.8,10.1,10.4,10.7,11])
  M_high=np.array([9.5,9.8,10.1,10.4,10.7,11,11.3])
  SFR=np.log10(np.array([4.7073727,7.914956,14.433758,23.877958,31.4691,36.42118,38.865368]))
  SFR_low=np.log10(np.array([1.8960778,3.0364738,6.203913,9.775229,13.308174,15.154303,17.82616]))
  SFR_high=np.log10(np.array([11.498701,19.972073,34.689453,57.387226,76.86959,88.96604,87.53305]))
  axes.errorbar(M,SFR,yerr=(np.array([SFR-SFR_low,SFR_high-SFR])),xerr=(np.array([M-M_low,M_high-M])),color='k',label='Lee et al. (2015)\n $z\leq 1.3$',zorder=106)

  axes.errorbar([7,11],0.867*np.array([7,11])-7.484,yerr=[0.354,0.354],color='limegreen',linestyle='-',label='Kurczynski et al. (2016)\n $1.5<z\leq 2.0$',linewidth=2.5,zorder=102)

  axes.plot([7.94287,11.0],[-1.24521,0.409962],':',color='gold',label='Bisigello et al. (2017) \nQuenched galaxies, $1<z<2$',linewidth=2.5,zorder=103)
  axes.plot([7.94287,11.0],[-0.59386975,2.0421455],'--',color='gold',label='Bisigello et al. (2017) \nMain sequence, $1<z<2$',linewidth=2.5,zorder=104)

  axes.plot([7.6,11.0],1.04*(np.array([7.6,11.0])-9.7)+1.04,'darkorange',label='Santini et al. (2017) \n$1.3\leq z<2.0$',linewidth=2.5,zorder=105)


if __name__=="__main__":
  filename='tuned_reion_T125'

  default=load_data('dragons10_T125')
  gals1=load_data(filename)

  
  ##Fig 1
  fig, axes = plt.subplots(1, 2,gridspec_kw = {'wspace':0, 'hspace':0})
  plot_hist2d(default,axes[0],False)
  plot_hist2d(gals1,axes[1],True,True)
  #plot_hist2d(gals2,axes[2])
  plot_obs(axes[0])
  plot_obs(axes[1])
  #plot_obs(axes[2])
  #axes[0].set_title('Default Meraxes')
  #axes[1].set_title('Modified Meraxes')
  axes[0].text(8, -5.5, r'$  z=1.5$'+'\n '+r'Q17 Meraxes',weight='bold',size='large')
  axes[1].text(8, -5.5,  r'$  z=1.5$'+'\n '+r'M18 Meraxes',weight='bold',size='large')
  axes[0].set_xlabel(r'$\log(M_\ast/M_\odot)$')
  axes[0].set_ylabel(r'$\log(\textrm{SFR})$')
  axes[1].set_xlabel(r'$\log(M_\ast/M_\odot)$')
  #axes[2].set_title('Bulge Model - Croton SF')
  #plt.suptitle("All Galaxies") 
  axes[1].set_yticklabels([])
  handles, labels = axes[1].get_legend_handles_labels()
  order = [3,4,2,1,0]
  lgd=axes[1].legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc=[1.55,0],fontsize='small')
  plt.tight_layout(rect=[0, 0.03, 0.98, 0.98])
  
  plt.savefig('/home/mmarshal/results/plots/MassSFR.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
  
  ###Fig 2
  #bulges=gals1[gals1['BulgeStellarMass']/gals1['StellarMass']>0.5]
  #disks=gals1[gals1['BulgeStellarMass']/gals1['StellarMass']<0.5]
  #plt.scatter(np.log10(bulges['StellarMass']*1e10),np.log10(bulges['Sfr']),c=bulges['BulgeStellarMass']/bulges['StellarMass'],vmax=1,vmin=0)
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('B/T')
  #plt.title('B/T>0.5') 
  #plt.show()

  ###Fig 3
  #plt.scatter(np.log10(disks['StellarMass']*1e10),np.log10(disks['Sfr']),c=disks['BulgeStellarMass']/disks['StellarMass'],vmax=1,vmin=0)
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('B/T')
  #plt.title('B/T<0.5')
  #plt.show()

  ###Fig 4
  #cold=gals1[gals1['ColdGas']>np.median(gals1['ColdGas'])]
  #hot=gals1[gals1['ColdGas']<np.median(gals1['ColdGas'])]
  #plt.scatter(np.log10(cold['StellarMass']*1e10),np.log10(cold['Sfr']),c=np.log10(cold['ColdGas']*1e10))
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('log(Cold Gas Mass)')
  #plt.title('Cold Gas Mass - galaxies with lots of cold gas')
  #plt.show()
  
  ###Fig 5
  #plt.scatter(np.log10(hot['StellarMass']*1e10),np.log10(hot['Sfr']),c=np.log10(hot['ColdGas']*1e10))
  #plt.xlabel('log(Stellar Mass)')
  #plt.ylabel('log(Sfr)')
  #cbar=plt.colorbar()
  #cbar.set_label('log(Cold Gas Mass)')
  #plt.title('Cold Gas Mass - galaxies with little cold gas')
  #plt.show()
 
  ##Fig 6
  central=gals1[gals1['Type']==0]
  sats=gals1[gals1['Type']!=0]

  fig,axes=plt.subplots(1,2)
  plot_hist2d(central,axes[0])
  plot_hist2d(sats,axes[1])
  axes[0].set_title('Centrals')
  axes[1].set_title('Satellites')
  plt.show()

