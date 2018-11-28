import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
from scipy.optimize import curve_fit
import ContourPlot as cp


#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (8,4)
#matplotlib.rcParams['figure.figsize'] = (7.2,4)
matplotlib.rcParams['lines.linewidth'] = 2.5
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def load_data(filename,snapshot,split):
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
  if split:
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','BulgeStellarMass','MergerBulgeStellarMass'],\
      h=cosmo['h'],quiet=True)
  else:
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','BulgeStellarMass'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['BlackHoleMass']*1e10>1e3]
  gals=gals[gals['StellarMass']*1e10>1e8]
  return gals


def plot_observations(axes,color):
  #Kormendy & Ho (2013):
  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
  axes.errorbar([8.17,12],logMBH-[8.17,12],yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)',capsize=3,linewidth=2.5, zorder=101,color=color[213])

  ##BLUETIDES
  #logMBH=(8.43+1.16*(np.log10(np.array([1e7,1e12])/1e11)))
  #HU=axes.plot([7,12],logMBH,':',linewidth=2.5,label="BlueTides: Huang et al. (2018)", zorder=102,color=color[0])
  
  ##Sijaki+15 Illustris sims
  logMBH=1.23*(np.log10(np.array([1e8,1e12])))-4.85
  axes.plot([8,12],logMBH-[8,12],':',linewidth=2.5,label="Illustris, $z=4$",zorder=102,color=color[116])
  logMBH=1.28*(np.log10(np.array([1e8,1e12])))-5.04
  axes.plot([8,12],logMBH-[8,12],':',linewidth=2.5,zorder=102,color=color[158],label="Illustris, $z=2$")


def func(x,a,b):
  return a*x+b


def quad(x,a,b,c):
  return a*(x+b)**2+c


def find_fit():
  #Find best fit to slope and intercept 
  zz=np.array(list(redshift.values()))
  plt.errorbar(zz,slope,slope_errs)
  popt,pcov = curve_fit(quad,zz,slope,sigma=slope_errs)
  print("SLOPE: popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))
  plt.plot(zz,quad(zz,*popt),'--')
  popt,pcov = curve_fit(quad,zz,slope)
  #print("SLOPE: popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))
  plt.plot(zz,quad(zz,*popt),':')
  plt.show()

  plt.errorbar(zz,inter,inter_errs)
  popt,pcov = curve_fit(quad,zz,inter,sigma=inter_errs)
  print("INTERCEPT: popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))
  plt.plot(zz,quad(zz,*popt),'--')
  popt,pcov = curve_fit(quad,zz,inter)
  #print("SLOPE: popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))
  plt.plot(zz,quad(zz,*popt),':')
  plt.show()


def best_fit(z,logMstellar):
  ##binwidth 0.3
  #slope=-0.013*(z-3.7)**2+0.87
  #inter=0.10*(z-4.3)**2-2.3
  ##binwidth 0.5
  slope=-0.017*(z-4.1)**2+0.91
  inter=0.14*(z-4.6)**2-2.8
  return slope*logMstellar+inter


def plot_MBHMstellar(filename,snapshots,mass_bulge,split,bulge_type,contours,color,axes):
  slope=np.zeros(len(snapshots))
  slope_errs=np.zeros(len(snapshots))
  inter=np.zeros(len(snapshots))
  inter_errs=np.zeros(len(snapshots))
  ii=0
  for snap in snapshots:
    gals=load_data(filename,snap,split)
    Mstellar=gals['StellarMass']*1e10
    if mass_bulge==0:
      Mstel=Mstellar
    else: ##Bulge Mass
      if split:
        #Total bulge mass
        if bulge_type==0:
          gals=gals[gals['BulgeStellarMass']>0]
          Mbulge=gals['BulgeStellarMass']*1e10
        #Merger bulge mass
        if bulge_type==2:
          gals=gals[gals['MergerBulgeStellarMass']>0]
          Mbulge=gals['MergerBulgeStellarMass']*1e10
        #ID bulge mass
        if bulge_type==1:
          gals=gals[gals['BulgeStellarMass']-gals['MergerBulgeStellarMass']>0]
          Mbulge=gals['BulgeStellarMass']*1e10-gals['MergerBulgeStellarMass']*1e10
      else:
        gals=gals[gals['BulgeStellarMass']>0]
        Mbulge=gals['BulgeStellarMass']*1e10
      Mstel=Mbulge
      

    MBH=gals['BlackHoleMass']*1e10
    logMstel=np.log10(Mstel)
    logMBH=np.log10(MBH) 
    #Bulge stellar mass
    bin_width=0.5
    min_mass=np.min(logMstel)
    max_mass=np.max(logMstel)
    n_bins=np.int((max_mass-min_mass)/bin_width)
    med_bh=np.zeros(n_bins)
    middle_sm=np.zeros(n_bins)
    pctl_bh=np.zeros((n_bins,2))
    for nn in range(0,n_bins):
      if np.size(logMBH[(logMstel>min_mass+(nn*bin_width))&(logMstel<min_mass+(nn+1)*bin_width)])>10:
        med_bh[nn]=np.median(logMBH[(logMstel>min_mass+(nn*bin_width))&(logMstel<min_mass+(nn+1)*bin_width)])
        middle_sm[nn]=min_mass+(nn+0.5)*bin_width
        pctl_bh[nn,:]=np.percentile(logMBH[(logMstel>min_mass+(nn*bin_width))&(logMstel<min_mass+(nn+1)*bin_width)],[16,84])
      else:
        med_bh[nn]=np.nan 
    middle_sm=middle_sm[np.logical_not(np.isnan(med_bh))]
    med_bh=med_bh[np.logical_not(np.isnan(med_bh))]
    axes.plot(middle_sm,med_bh-middle_sm,label='$z={}$'.format(redshift[snap]),color=color[snap])

    if (snap==158):
      print('Plotting Contour')
      cp.contour_plot(logMstel,logMBH-logMstel,xlab=None,ylab=None,xlims=[8.15,11.6],ylims=[-4,-2],axes=axes,colors=color[snap],levels=np.logspace(-2,1,7),linewidth=0.9)
    ii+=1
  axes.set_xlim([8.15,11.6])
  axes.set_ylim([-4,-2])


def plot_MBHMstellar_hist(filename,snapshots,mass_bulge,split,bulge_type,contours,color,axes):
  snapshots=[158]
  for snap in snapshots:
    gals=load_data(filename,snap,split)
    Mstellar=gals['StellarMass']*1e10
    if mass_bulge==0:
      Mstel=Mstellar
    else: ##Bulge Mass
      if split:
        #Total bulge mass
        if bulge_type==0:
          gals=gals[gals['BulgeStellarMass']>0]
          Mbulge=gals['BulgeStellarMass']*1e10
        #Merger bulge mass
        if bulge_type==2:
          gals=gals[gals['MergerBulgeStellarMass']>0]
          Mbulge=gals['MergerBulgeStellarMass']*1e10
        #ID bulge mass
        if bulge_type==1:
          gals=gals[gals['BulgeStellarMass']-gals['MergerBulgeStellarMass']>0]
          Mbulge=gals['BulgeStellarMass']*1e10-gals['MergerBulgeStellarMass']*1e10
      else:
        gals=gals[gals['BulgeStellarMass']>0]
        Mbulge=gals['BulgeStellarMass']*1e10
      Mstel=Mbulge
      

    MBH=gals['BlackHoleMass']*1e10
    logMstel=np.log10(Mstel)
    logMBH=np.log10(MBH) 
   
    plot_hist2d(logMstel,logMBH,axes,[8,max(logMstel)],[6,max(logMBH)])

def plot_hist2d(xdata,ydata,axes,xlims,ylims):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=20, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='BuPu',norm=matplotlib.colors.LogNorm())
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu')
  #cb=plt.colorbar(im,ax=axes,use_gridspec=True)
  #cb.set_label('N, tot N = {}'.format(np.size(xdata)))
  axes.set_aspect('auto')
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)




if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,213:0.55}
  snapshots=[63,78,100,116,134,158]
  #snapshots=[116,134,158]
  prop='StellarMass'
  color={52:'C0',63:'C1',78:'C2',100:'C3',116:'C4',134:'aqua',158:'pink',213:'k'}
  color={52:'#a65628',63:'#e41a1c',78:'#377eb8',100:'#4daf4a',116:'#984ea3',\
                  134:'#ff7f00',158:'#f781bf',194:'#a65628',213:'black'}
  filename='tuned_reion'
  #meraxes_loc='/output/run1/meraxes_001.hdf5'
  #filename='tuned_best_t125'
  meraxes_loc='/output/meraxes.hdf5'

  ##OPTIONS	
  plot_type=2 #Plot MBH vs Mstellar and Mbulge (1) or vs Mbulge for different bulge types (2)
  z0=1 #Plot z=0.55 relation?
  #filename_125='tuned_best_t125'
  filename_125='tuned_reion_T125'
  mass_bulge=0 #Plot the total stellar mass (0) or bulge stellar mass (1)?
  split=True #Is bulge split into MDBH and IDBH?
  bulge_type=0 #If wanting only one bulge component, which type? 1:IDBH, 2:MDBH
  contour=1 #Plot contours of z=2 distribution?
  
  ##PLOT
  fig,ax = plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)
  #fig,axes=plt.subplots(1,1)
  plot_type
  if plot_type==1: 
    plot_MBHMstellar(filename,snapshots,True,split,bulge_type,contour,color,ax[0])
    plot_MBHMstellar(filename,snapshots,False,split,bulge_type,contour,color,ax[1])
    if z0:
      plot_MBHMstellar(filename_125,[213],True,split,bulge_type,contour,color,ax[0])
      plot_MBHMstellar(filename_125,[213],False,split,bulge_type,contour,color,ax[1])
    plot_observations(ax[0],color)
    plot_observations(ax[1],color)
    ax[1].set_xlabel(r'$\log(\textrm{M}_\ast)$')
    ax[0].set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}})$')
  elif plot_type==2: 
    plot_MBHMstellar_hist(filename,snapshots,True,split,1,contour,color,ax[0])
    plot_MBHMstellar_hist(filename,snapshots,True,split,2,contour,color,ax[1])
    #if z0:
    #  plot_MBHMstellar(filename_125,[213],True,split,1,contour,color,ax[0])
    #  plot_MBHMstellar(filename_125,[213],True,split,2,contour,color,ax[1])
    plot_observations(ax[0],color)
    plot_observations(ax[1],color)
    ax[1].set_xlabel(r'$\log(\textrm{M}_{\textrm{merger-driven bulge}})$')
    ax[0].set_xlabel(r'$\log(\textrm{M}_{\textrm{instability-driven bulge}})$')


  plt.legend()
  lgd=plt.legend(fontsize='small',loc='upper center', bbox_to_anchor=(1.3, 0.82))  
  ax[0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 
  if plot_type==1:
    plt.savefig('/home/mmarshal/results/plots/MBHMStellarRelation.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  else:
    plt.savefig('/home/mmarshal/results/plots/MBHMStellarRelation_BulgeType.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

  ##PERFORM STATISTICS
  #find_fit()
