import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF
from scipy.optimize import curve_fit
import ContourPlot as cp


#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (5,4)
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
  axes.errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)',capsize=3,linewidth=2.5, zorder=101,color='black')

  ##BLUETIDES
  #logMBH=(8.43+1.16*(np.log10(np.array([1e7,1e12])/1e11)))
  #HU=axes.plot([7,12],logMBH,':',linewidth=2.5,label="BlueTides: Huang et al. (2018)", zorder=102,color=color[0])
  
  ##Sijaki+15 Illustris sims
  logMBH=1.23*(np.log10(np.array([1e8,1e12])))-4.85
  axes.plot([8,12],logMBH,':',linewidth=2.5,label="Illustris, $z=4$",zorder=102,color=color[4])
  logMBH=1.28*(np.log10(np.array([1e8,1e12])))-5.04
  axes.plot([8,12],logMBH,':',linewidth=2.5,zorder=102,color=color[6],label="Illustris, $z=2$")


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


if __name__=='__main__':
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  #redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  #snapshots=[52,63,78,100,116,134,158]
  redshift={63:7,63:7,78:6,100:5,116:4,134:3,158:2}
  snapshots=[63,63,78,100,116,134,158]
  #snapshots=[158]
  prop='StellarMass'
  #color={63:'C0',78:'C1',100:'C2',116:'C3',134:'C4',158:'black',213:'pink'}
  color={0:'C0',1:'C1',2:'C2',3:'C3',4:'C4',5:'aqua',6:'pink'}
  #filename='bulges_correctBHMF' 
  #filename='bulges_noreion'
  #filename='bulges_fullreion_nodiskinstability'
  #filename='bulges_split_d9a8304'#fullreion'
  filename='bulges_IDBH_tune_t125'
  z0=1
  split=True
  bulge=0

  slope=np.zeros(len(snapshots))
  slope_errs=np.zeros(len(snapshots))
  inter=np.zeros(len(snapshots))
  inter_errs=np.zeros(len(snapshots))
  #fig,axes = plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)
  fig,axes=plt.subplots(1,1)
  ii=0
  for snap in snapshots:
    gals=load_data(filename,snap,split)
    Mstellar=gals['StellarMass']*1e10
    ##Bulge Mass
    if split:
      #Total bulge mass
      if bulge==0:
        gals=gals[gals['BulgeStellarMass']>0]
        Mbulge=gals['BulgeStellarMass']*1e10
      #Merger bulge mass
      if bulge==2:
        gals=gals[gals['MergerBulgeStellarMass']>0]
        Mbulge=gals['MergerBulgeStellarMass']*1e10
      #ID bulge mass
      if bulge==1:
        gals=gals[gals['BulgeStellarMass']-gals['MergerBulgeStellarMass']>0]
        Mbulge=gals['BulgeStellarMass']*1e10-gals['MergerBulgeStellarMass']*1e10
    else:
      gals=gals[gals['BulgeStellarMass']>0]
      Mbulge=gals['BulgeStellarMass']*1e10
    MBH=gals['BlackHoleMass']*1e10
    logMbulge=np.log10(Mbulge)
    logMBH=np.log10(MBH) 
    #Bulge stellar mass
    bin_width=0.5
    min_mass=np.min(logMbulge)
    max_mass=np.max(logMbulge)
    n_bins=np.int((max_mass-min_mass)/bin_width)
    med_bh=np.zeros(n_bins)
    middle_sm=np.zeros(n_bins)
    pctl_bh=np.zeros((n_bins,2))
    for nn in range(0,n_bins):
      if np.size(logMBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)])>10:
        med_bh[nn]=np.median(logMBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)])
        middle_sm[nn]=min_mass+(nn+0.5)*bin_width
        pctl_bh[nn,:]=np.percentile(logMBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)],[16,84])
      else:
        med_bh[nn]=np.nan 
    middle_sm=middle_sm[np.logical_not(np.isnan(med_bh))]
    med_bh=med_bh[np.logical_not(np.isnan(med_bh))]
    axes.plot(middle_sm,med_bh,label='$z={}$'.format(redshift[snap]),color=color[ii])
    #best_fit,V=np.polyfit(middle_sm,med_bh,1,cov=True)
    #print("snapshot = {}, best fit = {}".format(snap,best_fit))
    #poly = np.poly1d(best_fit)
    #axes.plot(middle_sm,poly(middle_sm),':',color=color[ii])
   
 
  ##fit line to the data 
  #  if snap==78:
  #    popt,pcov = curve_fit(func,middle_sm,med_bh)
  #    print("snap {}, popt {}, perr {}".format(snap,popt,np.sqrt(np.diag(pcov))))
  #  slope[ii]=popt[0]
  #  slope_errs[ii]=np.sqrt(np.diag(pcov))[0]
  #  inter[ii]=popt[1]
  #  inter_errs[ii]=np.sqrt(np.diag(pcov))[1]
  #  axes.plot(middle_sm,func(middle_sm,*popt),':',color=color[ii])
 
  ##Plot line given by MBH,Mstellar,z best fit relation
  #  axes.plot(middle_sm,best_fit(redshift[snap],middle_sm),'--',color=color[ii])   


    #axes[0].text(6, 7, ' {}'.format(str(models[ii])),usetex=False)
    if snap==158:
      cp.contour_plot(logMbulge,logMBH,xlab=None,ylab=None,xlims=[8.15,11.6],ylims=[4,9.5],axes=axes,colors='pink',levels=np.logspace(-2,1,7),linewidth=0.9)
    ii+=1

  if z0:
    filename=filename
    gals=load_data(filename,213,split)
    Mstellar=gals['StellarMass']*1e10
    ##Bulge Mass
    gals=gals[gals['BulgeStellarMass']*1e10>1e8]
    Mbulge=gals['BulgeStellarMass']*1e10
    MBH=gals['BlackHoleMass']*1e10
    logMbulge=np.log10(Mbulge)
    logMBH=np.log10(MBH) 
    #Bulge stellar mass
    bin_width=0.5
    min_mass=np.min(logMbulge)
    max_mass=np.max(logMbulge)
    n_bins=np.int((max_mass-min_mass)/bin_width)
    med_bh=np.zeros(n_bins)
    middle_sm=np.zeros(n_bins)
    pctl_bh=np.zeros((n_bins,2))
    for nn in range(0,n_bins):
      if np.size(logMBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)])>10:
        med_bh[nn]=np.median(logMBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)])
        middle_sm[nn]=min_mass+(nn+0.5)*bin_width
        pctl_bh[nn,:]=np.percentile(logMBH[(logMbulge>min_mass+(nn*bin_width))&(logMbulge<min_mass+(nn+1)*bin_width)],[16,84])
      else:
        med_bh[nn]=np.nan 
    middle_sm=middle_sm[np.logical_not(np.isnan(med_bh))]
    med_bh=med_bh[np.logical_not(np.isnan(med_bh))]
    axes.plot(middle_sm,med_bh,label='$z=0.55$',color='k')
       
    
  plot_observations(axes,color)
  axes.set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$')
  axes.set_xlabel(r'$\log(\textrm{M}_{\textrm{bulge}})$')
  axes.set_xlim([8.15,11.6])
  axes.set_ylim([4,9.5])
  axes.legend()
  lgd=axes.legend(fontsize='small',loc='upper center', bbox_to_anchor=(1.25, 0.8))  
  plt.savefig('MBHBulge_z.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

  #find_fit()
