import numpy as np
from dragons import meraxes
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('/home/mmarshal/simulation_codes/Yuxiang/')
sys.path.append('/home/mmarshal/simulation_codes')
from _plot_obsGSMF import plot_obsGSMF
from scipy.optimize import curve_fit
import ContourPlot as cp


#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (6.6,3.2)
#matplotlib.rcParams['figure.figsize'] = (7.2,4)
matplotlib.rcParams['lines.linewidth'] = 2.5
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4


def load_data(filename,snapshot):
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
  gals=meraxes.io.read_gals(filename,\
      snapshot=snapshot,props=['GhostFlag','Mvir','StellarMass','BlackHoleMass'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[gals['BlackHoleMass']*1e10>10**(6)]
  gals=gals[gals['StellarMass']*1e10>1e8]
  #gals=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.8]
  return gals


def plot_observations(axes,color,mass):
  #Kormendy & Ho (2013):
  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
  axes.errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)\n- Fit',capsize=3,linewidth=2.5, color=colors[1],zorder=0)
  
  kh=pd.read_csv('/home/mmarshal/simulation_codes/data/KormendyHo_ellipticals_MBHMbulge.csv')
  logMBH=np.log10(kh['log(M_BH) [Msun]'])
  logMbulge=kh['log(M) [Msun]']
  axes.plot(logMbulge,logMBH,'o',color=colors[1],label="Kormendy \& Ho (2013)\n- Ellipticals",markersize=3)
  
  kh=pd.read_csv('/home/mmarshal/simulation_codes/data/KormendyHo_classicalbulges_MBHMbulge.csv')
  logMBH=np.log10(kh['log(M_BH) [Msun]'])
  logMbulge=kh['log(M) [Msun]']
  axes.plot(logMbulge,logMBH,'^',color=colors[1],label="Kormendy \& Ho (2013)\n- Classical Bulges",markersize=3)


def plot_MBHMstellar(filename,snapshots,contours,color,axes):
  ii=0
  for snap in snapshots:
    gals=load_data(filename,snap)
    Mstellar=gals['StellarMass']*1e10
    Mstel=Mstellar

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
    #axes.plot(middle_sm,med_bh,label='$z={}$ (Tiamat-125-HR)'.format(redshift[snap]),color=color[snap])
    
    cp.contour_plot(logMstel,logMBH,xlab=None,ylab=None,xlims=[9.1,11.6],ylims=[6,9.5],axes=axes,colors=color,linewidth=0.9)#,levels=np.logspace(-2,1,7))
    ii+=1
  axes.set_xlim([8.6,11.6])
  axes.set_ylim([5.95,9.5])


def plot_SMF(filename,snapshots,prop,boxwidth,axes,**kwargs):
  for snap in snapshots:
    gals=load_data(filename,snap)
    maxval=np.nanmax(np.log10(gals[prop][gals[prop]>0]*1e10)) 
    minval=np.nanmin(np.log10(gals[prop][gals[prop]>0]*1e10))
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=40)
    #hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]),range=(minval,maxval),bins=40)
    
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    Max=Max[hist>5]
    hist=hist[hist>5]
    #axes.plot(Max[hist>0],np.log10(hist[hist>0]/(bin_edges[1]-bin_edges[0])/boxwidth**3),**kwargs)
    hist_plus=hist+np.sqrt(hist)
    hist_minus=hist-np.sqrt(hist)
    phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    phi_plus=np.log10(hist_plus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    phi_minus=np.log10(hist_minus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    axes.plot(Max,phi,**kwargs)
    axes.fill_between(Max, phi_plus, phi_minus, alpha=0.5,color=kwargs['color'],label='__nolegend__')
    return axes


if __name__=='__main__':
  ##SETUP
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,213:0.55,250:0}
  color={52:'C0',63:'C1',78:'C2',100:'C3',116:'C4',134:'aqua',158:'pink',213:'k',250:'k'}

  ##OPTIONS	
  z0=1 #Plot z=0.55 relation?
  filename_1='/fred/oz013/yqiu/projects/dust/models/IMF_Sal/meraxes.hdf5'
  filename_2='/fred/oz013/yqiu/projects/dust/models/IMF_Kro/meraxes.hdf5'
  contour=1 #Plot contours of z=2 distribution?
  
  ##PLOT
  fig,axes=plt.subplots(1,1)
  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
  axes.errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)\n- Fit',capsize=3,linewidth=2.5, color=colors[1],zorder=0)
  plot_observations(axes,color,"stellar")
  plot_MBHMstellar(filename_1,[249],contour,'red',axes)
  plot_MBHMstellar(filename_2,[249],contour,'blue',axes)
  axes.set_xlabel(r'$\log(M_\ast/M_\odot)$')
  axes.set_ylabel(r'$\log(M_{\rm{BH}}/M_\odot)$') 
  plt.show()


  fig,axes=plt.subplots(1,1)
  plot_SMF(filename_1,[249],'StellarMass',125/0.678,axes,color='red')
  plot_SMF(filename_2,[249],'StellarMass',125/0.678,axes,color='blue')
  plt.show()
