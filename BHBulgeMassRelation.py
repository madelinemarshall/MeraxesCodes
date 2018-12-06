import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import pandas as pd
import ContourPlot as cp
import pylab as p
import scipy.stats as stats

matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,5)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def load_data(filename,snapshot):
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

  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
  snapshot=snapshot,props=['StellarMass','BulgeStellarMass','BlackHoleMass'],\
  h=cosmo['h'],quiet=True)
  gals=gals[gals['StellarMass']*1e10>1e7]
  gals=gals[gals['BlackHoleMass']>0]
  gals=gals[gals['BulgeStellarMass']>0]
  #gals=gals[gals['BulgeStellarMass']/gals['StellarMass']>0.7] ##BULGES ONLY
  return gals


def plot_simulations(gals,axes):
  xlims=[7,12]
  ylims=[6,10.5]

  #2D Histogram
  #axes.plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
  cp_gals=gals[(gals['BulgeStellarMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10>1e4)]
  #cp.contour_plot(np.log10(cp_gals['BulgeStellarMass']*1e10),np.log10(cp_gals['BlackHoleMass']*1e10),xlims=xlims,ylims=ylims,axes=axes)
  H, xedges, yedges, img=axes.hist2d(np.log10(cp_gals['BulgeStellarMass']*1e10), np.log10(cp_gals['BlackHoleMass']*1e10), bins=20, range=[xlims,ylims], cmin=1, cmap='BuPu',vmax=25000,norm=matplotlib.colors.LogNorm())#,,vmin=0,vmax=25000)
  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
  im=axes.imshow(H,extent=extent,cmap='BuPu',vmin=0,vmax=25000)
  #plt.colorbar(im,ax=axes,use_gridspec=True)
  axes.set_aspect('auto')


def plot_LOBF(gals,axes):
  xlims=[7,12]
  ylims=[6,10.5]

  #2D Histogram
  #axes.plot(np.log10(gals['BulgeStellarMass']*1e10),np.log10(gals['BlackHoleMass']*1e10),'.')
  cp_gals=gals[(gals['BulgeStellarMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10>1e4)]
  
  ###Line of best fit
  #x=np.log10(cp_gals['BulgeStellarMass']*1e10)
  #y=np.log10(cp_gals['BlackHoleMass']*1e10)
  #axes.plot(xlims, np.poly1d(np.polyfit(x, y, 1),)(xlims),'k--')
  #print(np.poly1d(np.poly1d(np.polyfit(x, y, 1))))
 
  #slope,inter,rval,pval,err=stats.linregress(x,y)
  #axes.plot(xlims,inter+slope*np.array(xlims),'r--')

  #axes.set_xlim(xlims)
  #axes.set_ylim(ylims)
  binwidth=0.1
  med=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  bulge=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  hh=xlims[0]
  bulge_mass=np.log10(cp_gals['BulgeStellarMass']*1e10)
  for idx in range(0,len(med)):
    samp=np.log10(cp_gals[(bulge_mass>=hh)&(bulge_mass<hh+binwidth)]['BlackHoleMass']*1e10)
    med[idx]=np.nanmedian(samp)
    bulge[idx]=hh+binwidth/2
    hh=hh+binwidth
  bulge=bulge[np.invert(np.isnan(med))]
  med=med[np.invert(np.isnan(med))]

  #axes.plot(bulge,np.log10(med)+10,'k-')
  slope,inter,rval,pval,err=stats.linregress(bulge, med)
  #print("slope = {}, inter = {}".format(slope,inter))
  axes.plot(xlims,inter+slope*np.array(xlims),'k',linewidth=2.5,label='Best Fit', zorder=106)
  #axes.text(10, 6, r'$\beta={:.2f}$'.format(slope),weight='bold',size='large')
  ##slope,inter,rval,pval,err=stats.linregress(np.log10(cp_gals['BulgeStellarMass']*1e10), np.log10(cp_gals['BlackHoleMass']*1e10))
  ##axes.plot(xlims,inter+slope*np.array(xlims),'r',linewidth=2.5,label='Best Fit', zorder=106)
  return [slope,inter]


def plot_LOBF_totalmass(gals,axes):
  xlims=[7,12]
  ylims=[6,10.5]
  cp_gals=gals[(gals['StellarMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10>1e4)]
  
  binwidth=0.1
  med=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  tot=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  hh=xlims[0]
  tot_mass=np.log10(cp_gals['StellarMass']*1e10)
  for idx in range(0,len(med)):
    samp=np.log10(cp_gals[(tot_mass>=hh)&(tot_mass<hh+binwidth)]['BlackHoleMass']*1e10)
    med[idx]=np.nanmean(samp)
    tot[idx]=hh+binwidth/2
    hh=hh+binwidth
  tot=tot[np.invert(np.isnan(med))]
  med=med[np.invert(np.isnan(med))]

  slope,inter,rval,pval,err=stats.linregress(tot,med)
  axes.plot(xlims,inter+slope*np.array(xlims),'--',color='darkslateblue',linewidth=2.5,label=r'Best Fit -- $M_{\ast,\mathrm{total}}$', zorder=107)



def plot_LOBF_bulges(gals,axes):
  xlims=[7,12]
  ylims=[6,10.5]
  cp_gals=gals[(gals['BulgeStellarMass']*1e10>1e7)&(gals['BlackHoleMass']*1e10>1e4)&(gals['BulgeStellarMass']/gals['StellarMass']>0.7)]
  
  binwidth=0.1
  med=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  tot=np.zeros(np.int((xlims[1]-xlims[0])/binwidth))
  hh=xlims[0]
  tot_mass=np.log10(cp_gals['BulgeStellarMass']*1e10)
  for idx in range(0,len(med)):
    samp=np.log10(cp_gals[(tot_mass>=hh)&(tot_mass<hh+binwidth)]['BlackHoleMass']*1e10)
    med[idx]=np.nanmean(samp)
    tot[idx]=hh+binwidth/2
    hh=hh+binwidth
  tot=tot[np.invert(np.isnan(med))]
  med=med[np.invert(np.isnan(med))]

  slope,inter,rval,pval,err=stats.linregress(tot, med)
  axes.plot(xlims,inter+slope*np.array(xlims),'-.',color='red',linewidth=2.5,label=r'Best Fit -- $B/T>0.7$', zorder=108)



def plot_observations(snapshot,axes):
  #Marconi & Hunt (2003)
  #27 objects (entire sample of known galaxies) with secure BH mass estimates from direct gas kinematical
  #or stellar dynamical determination, and accurate Lbulge.
  #Uses virial bulge mass (kR\sigma^2/G). Same sample as Merritt & Ferrarese. 0.21 dex dispersion in relation
  #log(MBH)=(8.28\pm0.06)+(0.96\pm0.07)(log(Mbulge)-10.9), rms 0.25 dex
  #or for all galaxies in their sample:
  #log(MBH)=(8.12\pm0.09)+(1.06\pm0.10)(log(Mbulge)-10.9), rms 0.49 dex
  
#  #logMBH=(8.28+0.96*(np.log10(np.array([10**7,10**12]))-10.9))
#  logMBH=(8.12+1.06*(np.log10(np.array([10**7,10**12]))-10.9))
#  #logMBH_low=(8.12+1.06*(np.log10(np.array([10**7,10**12]))-10.9)-0.49)
#  #logMBH_high=(8.12+1.06*(np.log10(np.array([10**7,10**12]))-10.9)+0.49)
#  
#  axes.errorbar([7,12],logMBH,yerr=0.333,linestyle='-',color='darkorange',label='Marconi \& Hunt (2003)',capsize=3,linewidth=2.5, zorder=100)
#  #axes.plot([7,12],logMBH_low,'--',color='orange',label='_nolegend_')
#  #axes.plot([7,12],logMBH_high,'--',color='orange',label='_nolegend_')
 
  #Magorrian:
  #log(MBH)=-1.79\pm1.35+(0.96\pm0.12)log(Mbulge)
  #logMBH=-1.79+0.96*(np.log10(np.array([10**7,10**12])))
  #axes.plot([7,12],logMBH,'C2-',label="Magorrian")
 
  #Haring & Rix
  #logMBH=(8.2+1.12*np.log10(np.array([10**7,10**12])/10**11))
  #logMBH_high=(8.2+1.12*np.log10(np.array([10**7,10**12])/10**11)+0.3)
  #logMBH_low=(8.2+1.12*np.log10(np.array([10**7,10**12])/10**11)-0.3)
  #axes.errorbar([7,12],logMBH,yerr=0.333,linestyle='-',color='gold',label='Haring \& Rix (2004)',capsize=3,linewidth=2.5, zorder=100)
  #axes.plot([7,12],logMBH_low,'--',color='g',label='_nolegend_')
  #axes.plot([7,12],logMBH_high,'--',color='g',label='_nolegend_')


  #Kormendy & Ho (2013):
  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
  logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**7,10**12])/10**11)+9
  axes.errorbar([7,12],logMBH,yerr=0.28,linestyle='-',color='darkorange',label='Kormendy \& Ho (2013)',capsize=3,linewidth=2.5, zorder=101)

  #if (snapshot==78):
  #  Wang_BH=np.log10([2.8e9,2.1e9,8.6e8,1.7e8])
  #  Wang_dyn=np.log10([9.6e10,12.5e10,7.2e10,1.3e10])
  #  W=axes.plot(Wang_dyn,Wang_BH,'ko')
  #  return W
  if (snapshot==52):
    ##Huang+18 (BlueTides)
    ##log(MBH)=(8.43\pm0.06)+(1.16\pm0.01)log(Mbulge) (where MBH is in Mo and Mbulge is in 10^11 Mo)
    logMBH=(8.43+1.16*(np.log10(np.array([1e7,1e12])/1e11)))
    HU=axes.plot([7,12],logMBH,'-',color='gold',linewidth=2.5,label="BlueTides: Huang et al. (2018)", zorder=102)
    return HU
    #olivedrab
  ##Sijaki+15 Illustris sims
  #Here the total stellar mass within the
  #stellar half-mass radius has been adopted as a proxy for the bulge
  #mass. Note that we do not morphologically distinguish between the
  #real bulges and pseudo-bulges but we do split galaxies into different
  #categories based on their colours.4 Furthermore, from now on we
  #take into account all galaxies hosting supermassive black holes
  #with stellar half-mass greater than 108 M and we refer to the
  #MBH–Mbulge relation of the whole population, even though many of
  #these galaxies might not contain a real bulge or might be effectively
  #bulgeless.
  if (snapshot==158):
    logMBH=1.23*(np.log10(np.array([1e8,1e12])))-4.85
    SJ=axes.plot([8,12],logMBH,'--',color='gold',linewidth=2.5,label="Illustris: Sijaki et al. (2015)",zorder=102)
    return SJ
  elif (snapshot==116):
    logMBH=1.28*(np.log10(np.array([1e8,1e12])))-5.04
    axes.plot([8,12],logMBH,'--',color='gold',linewidth=2.5,zorder=102)
  elif (snapshot==78):
    #Willott+17
    M_dyn=[10.097822,10.111968,10.5474205,10.5834055,10.589627,10.614172,10.641459,10.749401,10.770727,10.851768,10.825474,10.854801,10.975973,11.03756,11.122337,11.161653,11.15687,11.30606,10.6420,11.0849,11.3478]
    M_BH=[7.918992,8.239988,8.985337,9.380784,10.100969,8.394597,8.093332,9.287499,9.6985445,9.3348675,8.9825115,8.943487,9.456685,9.053792,9.441602,9.249973,9.187332,8.917856,7.71368,8.41601,8.71450]
    WL=axes.plot(M_dyn,M_BH,'k.',label='Willott et al. (2017)')
    return WL


def plot_inter(inter,redshift):
  ##log(MBH) \propto log(Mbulge)*slope 
  plt.plot(np.array(list(redshift.values())),10**np.array(list(inter.values())))
  plt.gca().invert_xaxis()
  plt.xlabel('Redshift')
  plt.yscale('log')


def plot_const_ratio(axes):
  xlims=np.array([7,12])
  axes.plot(xlims,-2+xlims,'k:',linewidth=2.5,label=r"$\log\frac{M_{\mathrm{BH}}}{M_{\mathrm{Bulge}}}=-2,-3,-4$", zorder=103)
  axes.plot(xlims,-3+xlims,'k:',linewidth=2.5,label="_nolegend_", zorder=104)
  axes.plot(xlims,-4+xlims,'k:',linewidth=2.5,label="_nolegend_", zorder=105)


def plot_median_ratio(redshift,median_BH_bulge,eightyfourth_pctile_BH_bulge,sixteenth_pctle_BH_bulge,axes):
  axes.plot(np.array(list(redshift.values())),np.array(list(median_BH_bulge.values())),'k')
  axes.plot(np.array(list(redshift.values())),np.array(list(eightyfourth_pctile_BH_bulge.values())),'k--')
  axes.plot(np.array(list(redshift.values())),np.array(list(sixteenth_pctle_BH_bulge.values())),'k--')
  axes.set_xlabel('Redshift')
  axes.set_ylabel('Median BH/Bulge Mass Ratio')

  axes.invert_xaxis()

  #plt.plot(np.array(list(redshift.values())),np.array(list(median_BH_bulge.values())),'k')
  #plt.plot(np.array(list(redshift.values())),np.array(list(eightyfourth_pctile_BH_bulge.values())),'k--')
  #plt.plot(np.array(list(redshift.values())),np.array(list(sixteenth_pctle_BH_bulge.values())),'k--')

  
if __name__=="__main__":
  filename='tuned_reion'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2}
  #redshift={63:7,63:7,78:6,100:5,116:4,134:3,158:2}
  median_BH_bulge={52:0,63:0,78:0,100:0,116:0,134:0,158:0}
  eightyfourth_pctile_BH_bulge={52:0,63:0,78:0,100:0,116:0,134:0,158:0}
  sixteenth_pctle_BH_bulge={52:0,63:0,78:0,100:0,116:0,134:0,158:0}
  snapshots=[52,63,78,100,116,134,158]
  #snapshots=[63,63,78,100,116,134,158]
  fig, axes = plt.subplots(2, int((len(snapshots)+1)/2),gridspec_kw = {'wspace':0, 'hspace':0})

  slope={}#np.zeros(len(snapshots))
  inter={}#np.zeros(len(snapshots))
  ii=-1
  j=0
  for snapshot in snapshots:
    ii+=1
    if ii==(len(snapshots)+1)/2:
      j+=1
      ii=0
    gals=load_data(filename,snapshot)
    plot_simulations(gals,axes[j,ii])
    #slope[snapshot],inter[snapshot]=plot_LOBF(gals,axes[j,ii])
    #plot_LOBF_totalmass(gals,axes[j,ii])
    ##plot_LOBF_bulges(gals,axes[j,ii])
    plot_const_ratio(axes[j,ii])
    if snapshot==52:
      HU=plot_observations(snapshot,axes[j,ii])
    elif snapshot==78:
      WL=plot_observations(snapshot,axes[j,ii])
    elif snapshot==158:
      SJ=plot_observations(snapshot,axes[j,ii])
    #elif snapshot==78:
    #  W=plot_observations(snapshot,axes[j,ii])
    else:
      plot_observations(snapshot,axes[j,ii])
    median_BH_bulge[snapshot]=np.log10(np.median(gals['BlackHoleMass']/gals['BulgeStellarMass']))
    eightyfourth_pctile_BH_bulge[snapshot]=np.log10(np.percentile(gals['BlackHoleMass']/gals['BulgeStellarMass'],84))
    sixteenth_pctle_BH_bulge[snapshot]=np.log10(np.percentile(gals['BlackHoleMass']/gals['BulgeStellarMass'],16))
    if j==1:
      axes[j,ii].set_xlabel(r'$\log(M_{\mathrm{Bulge}}/M_\odot$)')
    else:
      axes[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(M_{\mathrm{BH}}/M_\odot)$')
    else:
      axes[j,ii].set_yticklabels([])
    axes[j,ii].set_xlim([7,12.1])
    axes[j,ii].text(8.7, 9.8, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='x-large')

  snapshot=100
  j=0
  ii=3
  hand=(axes[j,ii].get_legend_handles_labels()[0])
  lab=(axes[j,ii].get_legend_handles_labels()[1])
  #hand.append(W[0])
  #lab.append("Wang et al. (2013)")
  hand.append(HU[0])
  lab.append("BlueTides: Huang et al. (2018)")
  hand.append(SJ[0])
  lab.append("Illustris: Sijacki et al. (2015)")
  hand.append(WL[0])
  lab.append('Willott et al. (2017)')
  lgd=axes[j,ii].legend(hand,lab,loc=(0.02,-0.75),fontsize='small')
  #axes[j,ii].legend(loc=(0.02,-0.72),handles=['Meraxes - Best Fit','Marconi \& Hunt (2003)','Haring \& Rix (2004)',\
  #"Huang et al. (2018)","$M_{BH}/M_{Bulge}=10^{-2},10^{-3},10^{-4}$"])

  axes[1,int((len(snapshots)+1)/2)-1].axis('off')
  
  #plt.tight_layout()
  fig.savefig('BHBulge.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show() 
  
  #fig,axes=plt.subplots(1,1)
  #plot_median_ratio(redshift,median_BH_bulge,eightyfourth_pctile_BH_bulge,sixteenth_pctle_BH_bulge,axes)
  
  ##plot_inter(inter,redshift)
  #plt.show()
