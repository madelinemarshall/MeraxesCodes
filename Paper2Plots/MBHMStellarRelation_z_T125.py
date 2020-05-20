import numpy as np
from dragons import meraxes
import os
import matplotlib
matplotlib.use('Agg')
import sys
sys.path.append('/home/mmarshal/simulation_codes/Yuxiang/')
sys.path.append('/home/mmarshal/simulation_codes')
from _function import _polyfit_bootstrap
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import ContourPlot as cp
from _load_data import load_data
from collections import OrderedDict


#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (6.5,3.2)
matplotlib.rcParams['lines.linewidth'] = 2.5
plt.rc('text', usetex=True)
plt.rc('font', family='serif') 
ls = OrderedDict(
    [(250,               (0, ())),
     (134,              (0, (1, 2))),
     (78,      (0, (1, 1))),

     (158,              (0, (5, 3))),
     (116,      (0, (5, 1))),

     ('dashdotted',          (0, (3, 3, 1, 3))),
     (194,  (0, (3, 1, 1, 1))),

     (100,         (0, (3, 2, 1, 2, 1, 2)))])



def plot_hist2d(xdata,ydata,axes,cmax=None,cbar=False):
  xlims=[9.5,11.6]
  ylims=[6,9.5]#[np.nanmin(y),np.nanmax(y)]
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, vmin=1, vmax=cmax, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  if True:#cbar:
    cb=plt.colorbar(img, cax=ax[2],use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))
 

def plot_observations(axes,color):
  #Kormendy & Ho (2013):
#  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
#  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
#  axes.errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)',capsize=3,linewidth=2.5, zorder=101,color=color[250])

  ##BLUETIDES
  #logMBH=(8.43+1.16*(np.log10(np.array([1e7,1e12])/1e11)))
  #HU=axes.plot([7,12],logMBH,':',linewidth=2.5,label="BlueTides: Huang et al. (2018)", zorder=102,color=color[0])
  
  ##Sijaki+15 Illustris sims
  logMBH=1.23*(np.log10(np.array([1e8,1e12])))-4.85
  axes.plot([8,12],logMBH,':',linewidth=2.5,label="Illustris, $z=4$",zorder=102,color=color[116])
  logMBH=1.28*(np.log10(np.array([1e8,1e12])))-5.04
  axes.plot([8,12],logMBH,':',linewidth=2.5,zorder=102,color=color[158],label="Illustris, $z=2$")


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
  plt.plot(zz,quad(zz,*popt),':')
  plt.show()

  plt.errorbar(zz,inter,inter_errs)
  popt,pcov = curve_fit(quad,zz,inter,sigma=inter_errs)
  print("INTERCEPT: popt {}, perr {}".format(popt,np.sqrt(np.diag(pcov))))
  plt.plot(zz,quad(zz,*popt),'--')
  popt,pcov = curve_fit(quad,zz,inter)
  plt.plot(zz,quad(zz,*popt),':')
  plt.show()


def best_fit(z,logMstellar):
  slope=-0.017*(z-4.1)**2+0.91
  inter=0.14*(z-4.6)**2-2.8
  return slope*logMstellar+inter


def fit_parameters(x,y,mass_bulge,snap):
    #elements=(y>6)&(y<7)
    elements=(x>9.5)
    x=x[elements]
    y=y[elements]
    fit,err,cov_0,cov_1=_polyfit_bootstrap(x,y,1,sigma=[], num_samples = 10000, alpha=0.68, asymmetric=False,cov=True)
    print(cov_0,cov_1)
    resid=y-func(x,fit[0],fit[1])
    std=np.sqrt(np.sum(resid**2)/(len(resid)-1))

    fit_params.loc[redshift[snap],'N']=len(y)

    if mass_bulge:
      fit_params.loc[redshift[snap],'a_B']=fit[0] 
      fit_params.loc[redshift[snap],'a_err_B']=err[0] 
      fit_params.loc[redshift[snap],'b_B']=fit[1] 
      fit_params.loc[redshift[snap],'b_err_B']=err[1] 
      fit_params.loc[redshift[snap],'s_B']=std
    else:
      fit_params.loc[redshift[snap],'a']=fit[0] 
      fit_params.loc[redshift[snap],'a_err']=err[0] 
      fit_params.loc[redshift[snap],'b']=fit[1] 
      fit_params.loc[redshift[snap],'b_err']=err[1] 
      fit_params.loc[redshift[snap],'s']=std
    test_scatter=True
    if test_scatter:
      fit_highm,err_highm=_polyfit_bootstrap(x[x>10],y[x>10],1,sigma=[], num_samples = 1000, alpha=0.68, asymmetric=False)
      resid=y[x>10]-func(x[x>10],fit_highm[0],fit_highm[1])
      std_highm1=np.sqrt(np.sum(resid**2)/(len(resid)-1))
      fit_highm,err_highm=_polyfit_bootstrap(x[x>10.5],y[x>10.5],1,sigma=[], num_samples = 1000, alpha=0.68, asymmetric=False)
      resid=y[x>10.5]-func(x[x>10.5],fit_highm[0],fit_highm[1])
      std_highm2=np.sqrt(np.sum(resid**2)/(len(resid)-1))
      print(std,std_highm1,std_highm2)
    return fit



def plot_MBHMstellar(filename,snapshots,mass_bulge,split,bulge_type,contours,color,axes,cbar=False):
  slope=np.zeros(len(snapshots))
  slope_errs=np.zeros(len(snapshots))
  inter=np.zeros(len(snapshots))
  inter_errs=np.zeros(len(snapshots))
  ii=0
  for snap in snapshots:
    print(f"Snapshot {snap}")
    gals=load_data(filename,snap,['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','BulgeStellarMass','MergerBulgeStellarMass','BlackHoleMass_ID','BlackHoleMass_MD']) 
    Mstellar=gals['StellarMass']*1e10
    if mass_bulge==0:
      Mstel=Mstellar
      MBH=gals['BlackHoleMass']*1e10
    else: ##Bulge Mass
      if split:
        #Total bulge mass
        if bulge_type==0:
          gals=gals[gals['BulgeStellarMass']>0]
          Mbulge=gals['BulgeStellarMass']*1e10
          MBH=gals['BlackHoleMass']*1e10
        #Merger bulge mass
        if bulge_type==2:
          gals=gals[gals['MergerBulgeStellarMass']>0]
          Mbulge=gals['MergerBulgeStellarMass']*1e10
          MBH=gals['BlackHoleMass']*1e10
        #ID bulge mass
        if bulge_type==1:
          gals=gals[gals['BulgeStellarMass']-gals['MergerBulgeStellarMass']>0]
          Mbulge=gals['BulgeStellarMass']*1e10-gals['MergerBulgeStellarMass']*1e10
          MBH=gals['BlackHoleMass']*1e10
      else:
        gals=gals[gals['BulgeStellarMass']>0]
        Mbulge=gals['BulgeStellarMass']*1e10
        MBH=gals['BlackHoleMass']*1e10
      Mstel=Mbulge
      

    logMstel=np.log10(Mstel)
    logMBH=np.log10(MBH) 
    #Bulge stellar mass
    bin_width=0.5
    min_mass=np.min(logMstel)
    max_mass=np.max(logMstel)

    fit=fit_parameters(logMstel,logMBH,mass_bulge,snap)
 
    axes.plot(np.array([min_mass,max_mass]),func(np.array([min_mass,max_mass]),fit[0],fit[1]),color=color[snap],label='$z={}$'.format(redshift[snap]),linestyle=ls[snap])

    if (snap==250)&contours:
      print('Plotting Contour')
      plot_hist2d(logMstel,logMBH,axes,cmax=800,cbar=cbar)
    ii+=1
  axes.set_xlim([9.51,11.6])
  axes.set_ylim([6,9.5])



if __name__=='__main__':
  ##SETUP
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:1,250:0}
  snapshots=[78,100,116,134,158]
  prop='StellarMass'
  color={63:'pink',78:'#ff7f00',100:'#a65628',116:'#e41a1c',134:'yellow',100:'#4daf4a',158:'#f781bf',194:'#a65628',250:'black'}
  filename='paper2_T125'

  ##OPTIONS	
  plot_type=1 #Plot MBH vs Mstellar and Mbulge (1) or vs Mbulge for different bulge types (2)
  z0=1 #Plot z=0.55 relation?
  filename_125='paper2_T125'
  split=True #Is bulge split into MDBH and IDBH?
  bulge_type=0 #If wanting only one bulge component, which type? 1:IDBH, 2:MDBH
  contour=1 #Plot contours of z=2 distribution?
  
  zz=[0,1,2,3,4,5,6]
  fit_params=pd.DataFrame({'z':zz}, columns=['z','N','a','a_err','b','b_err','s','a_B','a_err_B','b_B','b_err_B','s_B']) 
  fit_params.set_index("z", inplace=True)
 
  ##PLOT
  fig,ax = plt.subplots(1,3,gridspec_kw = {'wspace':0, 'hspace':0,'width_ratios':[4,4,0.4]})
  plot_type
  if plot_type==1: 
    plot_MBHMstellar(filename,snapshots,True,split,bulge_type,contour,color,ax[0])
    plot_MBHMstellar(filename,snapshots,False,split,bulge_type,contour,color,ax[1])
    if z0:
      plot_MBHMstellar(filename_125,[194,250],True,split,bulge_type,contour,color,ax[0])
      plot_MBHMstellar(filename_125,[194,250],False,split,bulge_type,contour,color,ax[1])
    ax[1].set_xlabel(r'$\log(M_\ast/M_\odot)$')
    ax[0].set_xlabel(r'$\log(M_{\rm{bulge}}/M_\odot)$')
    ax[0].set_ylabel(r'$\log(M_{\rm{BH}}/M_\odot)$') 
  elif plot_type==2: 
    plot_MBHMstellar(filename,snapshots,True,split,1,contour,color,ax[0])
    plot_MBHMstellar(filename,snapshots,True,split,2,contour,color,ax[1])
    if z0:
      plot_MBHMstellar(filename_125,[194,250],True,split,1,contour,color,ax[0])
      plot_MBHMstellar(filename_125,[194,250],True,split,2,contour,color,ax[1])
    ax[1].set_xlabel(r'$\log(\textrm{M}_{\textrm{merger-driven bulge}})$')
    ax[0].set_xlabel(r'$\log(\textrm{M}_{\textrm{instability-driven bulge}})$')
    ax[0].set_ylabel(r'$\log(\textrm{M}_{\textrm{BH}})$') 

  ax[1].set_yticks([])
  plt.tight_layout()
  lgd=ax[0].legend(fontsize='small',loc='upper left')#, bbox_to_anchor=(1.13, 0.7))  
  if plot_type==1:
    plt.savefig('/home/mmarshal/results/plots/Paper2/MBHMStellarRelation_T125.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  else:
    plt.savefig('/home/mmarshal/results/plots/Paper2/MBHMStellarRelation_BulgeType_T125.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()
  fit_params=fit_params.astype(float)
  fit_params=fit_params.round({'a':3,'a_err': 3,'a_B':3,'a_err_B':3,'b':2,'b_err': 2,'b_B':2,'b_err_B': 2,'s': 2,'s_B':2,'N':0})
  with open('fit_params.tex', 'w') as tf:
     tf.write(fit_params.to_latex())
  print(fit_params)

  ##PERFORM STATISTICS
  #find_fit()
