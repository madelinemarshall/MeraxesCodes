##Creates plots like those of figure  3, 4 and 9 in Tonini+16, for comparison with that work
import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import matplotlib
import sys
import pandas as pd
import magcalc as mc
from _load_data import load_data
from scipy.optimize import curve_fit
#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.6,4)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4



def load_mags(filename,snapshot):
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def plot_hist2d(xdata,ydata,axes,xlims,ylims,cmax=None,cbar=False):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, vmin=1, vmax=cmax, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  if cbar:
    cb=plt.colorbar(img, cax=cbar,use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))


def plot_avg(xdata,ydata,axes,xlims,bin_width,**kwargs):
  min_bin=np.floor(xlims[0])
  max_bin=np.floor(xlims[1])
  n_bins=np.int((max_bin-min_bin)/bin_width)
  avg_r=np.zeros(n_bins)
  pct_r_16=np.zeros(n_bins)
  pct_r_84=np.zeros(n_bins)
  bin_centre=np.zeros(n_bins)

  bin_start=min_bin
  bin_end=min_bin+1
  bin_num=0
  while bin_num<n_bins:
    y=ydata[(xdata<bin_end)&(xdata>=bin_start)]
    if np.size(y)>=25:
      avg_r[bin_num]=np.median(y)
      pct_r_16[bin_num]=np.percentile(y,16)
      pct_r_84[bin_num]=np.percentile(y,84)
    else:
      avg_r[bin_num]=np.nan
      pct_r_16[bin_num]=np.nan
      pct_r_84[bin_num]=np.nan
    bin_centre[bin_num]=bin_start+bin_width/2
    bin_start+=bin_width
    bin_end+=bin_width
    bin_num+=1
  
  axes.errorbar(bin_centre,avg_r,yerr=np.array([avg_r-pct_r_16,pct_r_84-avg_r]),**kwargs)#,color='k',marker='s',markersize=4,label='M19 - Median')


def plot_obs_lum(axes,z):
  #Kawamata+18 (from fit)
  if z==9:
    kaw_mass=np.array([-21.5,-18])
    axes.plot(kaw_mass,-0.4*0.56*(kaw_mass+21)+np.log10(1.2),'-',color=colors[2],label='Kawamata et al. (2018)')
  if z==6 or z==7:
    kaw_mass=np.array([-21,-15])
    axes.plot(kaw_mass,-0.4*0.46*(kaw_mass+21)+np.log10(0.94),'-',color=colors[2],label='Kawamata et al. (2018)')
  if z==8:
    kaw_mass=np.array([-21.5,-16.5])
    axes.plot(kaw_mass,-0.4*0.38*(kaw_mass+21)+np.log10(0.81),'-',color=colors[2],label='Kawamata et al. (2018)')
  

  mag=np.array([-18,-21.6])
  if z==4:
    axes.plot(mag,mag_to_R(1.34,0.22,mag),'-',label='Huang et al. (2013)',color=colors[3],zorder=100) #verified same eqn from original paper
  if z==5:
    axes.plot(mag,mag_to_R(1.19,0.25,mag),'-',label='Huang et al. (2013)',color=colors[3],zorder=100)
  #if z==7:
    #axes.plot(mag,mag_to_R(1.55,1.05,mag),'--',label='Ono et al. (2013)',color=colors[0])
    #axes.plot(mag,mag_to_R(0.86,0.24,mag),'--',label='Grazian et al. (2012)',color=colors[3])
  #if z==8:
  #  axes.plot(mag,mag_to_R(1.44,1.07,mag),'--',label='Ono et al. (2013)',color=colors[0])
  #if z==9 or z==10:
  #  axes.plot(mag,mag_to_R(0.57,0.12,mag),'--',label='Holwerda et al. (2015)',color=colors[3])

  #Holwerda+15
  if z==9 or z==10:
    mag=np.array([-21.33,-20.94,-20.87,-20.84,-20.81,-20.69,-18.07])
    R=np.log10(np.array([0.601,0.639,0.790,0.382,0.459,0.494,0.313]))
    R_down=np.log10(np.array([0.331,0.550,0.621,0.306,0.372,0.352,0.173]))
    R_up=np.log10(np.array([0.869,0.731,0.953,0.459,0.545,0.637,0.454]))
    axes.errorbar(mag,R,yerr=np.array([R-R_down,R_up-R]),color=colors[3],marker='h',markersize=5,label='Holwerda et al. (2015)',linestyle='',elinewidth=0.5,capthick=0.5)


  #Ono+13 ##Stacked objects in luminosity bins
  if z==7:
    mag=np.array([-20.17,-18.94,-18.19])
    R=np.log10(np.array([0.642,0.284,0.352]))
    R_up=np.log10(np.array([0.958,0.372,0.497]))
    R_down=np.log10(np.array([0.325,0.196,0.207]))
    mag_up=np.array([-19.85,-18.80,-17.93])
    mag_down=np.array([-20.50,-19.05,-18.44])
    axes.errorbar(mag,R,yerr=np.array([R-R_down,R_up-R]),xerr=np.array([mag-mag_down,mag_up-mag]),color=colors[0],marker='o',markersize=5,label='Ono et al. (2013)',linestyle='',elinewidth=0.5,capthick=0.5)
  if z==8 or z==9:
    mag=np.array([-20.14,-19.15,-18.22])
    R=np.log10(np.array([0.654,0.310,0.362]))
    R_up=np.log10(np.array([0.937,0.392,0.494]))
    R_down=np.log10(np.array([0.371,0.226,0.231]))
    mag_up=np.array([-20.00,-19.04,-17.94])
    mag_down=np.array([-20.28,-19.26,-18.53])
    axes.errorbar(mag,R,yerr=np.array([R-R_down,R_up-R]),xerr=np.array([mag-mag_down,mag_up-mag]),color=colors[0],marker='o',markersize=5,label='Ono et al. (2013)',linestyle='',elinewidth=0.5,capthick=0.5)

  #Grazian+12
  if z==7:
    mag=np.array([-20.5,-19.5,-18.5])
    R=np.log10(10**np.array([-0.785,-1.066,-1.097])*5.227) #=log(Rh) in arcsec, 1 arcsec= 5.227kpc at z=7
    R_up=np.log10(10**np.array([-0.660,-0.965,-0.979])*5.227)
    R_down=np.log10(10**np.array([-0.850,-1.124,-1.4])*5.227)
    mag_up=np.array([-21,-20,-19])
    mag_down=np.array([-20,-19,-18])
    axes.errorbar(mag,R,yerr=np.array([R-R_down,R_up-R]),xerr=np.array([mag-mag_down,mag_up-mag]),color=colors[4],marker='p',markersize=5,label='Grazian et al. (2012)',linestyle='',elinewidth=0.5,capthick=0.5)



def lum_func(mag,r0,b):
  return -0.4*b*(mag+21)+np.log10(r0) #returns log(r)


def mass_func(mass,r0,b):
  return b*(mass-9)+np.log10(r0) #returns log(r)


def fit_equation_mass(mass,rad,axes,snapshot):
  popt,pcov = curve_fit(mass_func,mass,rad)

  fit_params.loc[redshift[snapshot],'M_r']=popt[0] 
  fit_params.loc[redshift[snapshot],'M_b']=popt[1] 
  fit_params.loc[redshift[snapshot],'M_r_err']=np.sqrt(np.diag(pcov))[0] 
  fit_params.loc[redshift[snapshot],'M_b_err']=np.sqrt(np.diag(pcov))[1] 

  #axes.plot(mass,mass_func(mass,*popt),'--',color='k',lw=0.5,label='M19 - Best Fit')


def fit_equation_lum(mag,rad,axes,snapshot):
  popt,pcov = curve_fit(lum_func,mag,rad)

  fit_params.loc[redshift[snapshot],'L_r']=popt[0] 
  fit_params.loc[redshift[snapshot],'L_b']=popt[1] 
  fit_params.loc[redshift[snapshot],'L_r_err']=np.sqrt(np.diag(pcov))[0] 
  fit_params.loc[redshift[snapshot],'L_b_err']=np.sqrt(np.diag(pcov))[1] 

  #axes.plot(mag,lum_func(mag,*popt),'--',color='k',lw=0.5,label='M19 - Best Fit')


def fit_equation_lum_cut(mag,rad,axes,snapshot):
  popt,pcov = curve_fit(lum_func,mag,rad)

  fit_params.loc[redshift[snapshot],'L_r_MC']=popt[0] 
  fit_params.loc[redshift[snapshot],'L_b_MC']=popt[1] 
  fit_params.loc[redshift[snapshot],'L_r_err_MC']=np.sqrt(np.diag(pcov))[0] 
  fit_params.loc[redshift[snapshot],'L_b_err_MC']=np.sqrt(np.diag(pcov))[1] 

  #axes.plot(mag,lum_func(mag,*popt),'--',color='k',lw=0.5,label=r'M19 - Best Fit, $M_{UV}<-14$')


def plot_obs_mass(axes,z):
  #See Holwerda+15 who list/derive these values from other studies
  #M=np.array([1e7,10**9.5])
  #if z==6:
    #axes.plot(np.log10(M),mass_to_R(0.75,0.14,M),'--',label='Mosleh et al. (2014)',color=colors[0])
  #if z==7: #Not given in original papers, masses converted by Holwerda+15
    #axes.plot(np.log10(M),mass_to_R(0.27,1.35,M),'--',label='Ono et al. (2013)',color=colors[0])
    #axes.plot(np.log10(M),mass_to_R(0.64,0.24,M),'--',label='Grazian et al. (2012)',color=colors[0])
  #if z==9 or z==10:
  #  axes.plot(np.log10(M),mass_to_R(0.57,0.12,M),'--',label='Holwerda et al. (2015)',color=colors[0])

  #Mosleh+14
  if z==4:
    axes.errorbar([8.76,9.35],[-0.112,-0.043],yerr=[np.array([0.082,0.082]),np.array([0.082,0.082])],\
    marker='D',label='Mosleh et al. (2014)',color=colors[4],markersize=5,linestyle='',elinewidth=0.5,capthick=0.5)
  if z==5:
    axes.errorbar([8.78,9.42],[-0.224,-0.121],yerr=[np.array([0.306-0.224,0.203-0.121]),np.array([0.224-0.138,0.121-0.034])],\
    marker='D',label='Mosleh et al. (2014)',color=colors[4],markersize=5,linestyle='',elinewidth=0.5,capthick=0.5)
  if z==6:
    axes.errorbar([8.76,9.35],[-0.515,-0.334],yerr=[np.array([0.602-0.515,-0.334+0.541]),np.array([+0.515-0.32,-0.015+0.334])],\
    marker='D',label='Mosleh et al. (2014)',color=colors[4],markersize=5,linestyle='',elinewidth=0.5,capthick=0.5)
  if z==7:
    axes.errorbar([7.63],[-0.554],yerr=[np.array([0.640-0.554]),np.array([0.554-0.287])],\
    marker='D',label='Mosleh et al. (2014)',color=colors[4],markersize=5,linestyle='',elinewidth=0.5,capthick=0.5)

  #Holwerda+15
  if z==10 or z==9:
    mag=np.array([7.63,7.90,8.05,8.50,8.70,9.20,9.40,9.50])
    R=np.log10(np.array([0.466,0.458,0.310,0.492,0.596,0.380,0.639,0.786]))
    R_up=np.log10(np.array([1.347,0.553,0.458,0.639,0.872,0.466,0.734,0.950]))
    R_down=np.log10(np.array([0,0.363,0.164,0.354,0.337,0.302,0.544,0.639]))
    mag_up=np.array([7.63,7.90,8.05,8.90,9.00,9.50,9.70,9.80])
    mag_down=np.array([7.63,7.90,8.05,8.10,8.40,8.90,9.10,9.20])
    axes.errorbar(mag,R,yerr=np.array([R-R_down,R_up-R]),xerr=np.array([mag-mag_down,mag_up-mag]),color=colors[3],marker='h',markersize=5,label='Holwerda et al. (2015)',linestyle='',elinewidth=0.5,capthick=0.5)


def mag_to_R(Rstar,b,mag):  #See Holwerda+15
  return np.log10(Rstar*10**(b/2.5*(-21.0-mag)))


def mass_to_R(Rstar,beta,M):
  return np.log10(Rstar*(M/1e9)**(beta)) #See Holwerda+15


def make_legend_mass(ax):
    #ax.plot([0,0],[0,0],'-',color=[0.5,0.5,0.5],lw=0.5,label='M19 - Best Fit')
    ax.errorbar([0],[0],yerr=[1],color='k',marker='s',markersize=4,label='M19 - Median')
    ax.errorbar([0],[0],yerr=[1],color=[0.5,0.5,0.5],marker='o',markersize=4,label='L16 - Median')
    ax.errorbar([0],[0],yerr=[np.array([0.640-0.554]),np.array([0.554-0.287])],\
    marker='D',label='Mosleh et al. (2014)',color=colors[4],markersize=5,linestyle='',elinewidth=0.5,capthick=0.5)
    ax.errorbar(0,0,yerr=1,xerr=1,color=colors[3],marker='h',markersize=5,label='Holwerda et al. (2015)',linestyle='',elinewidth=0.5,capthick=0.5)
    
    ax.legend(fontsize='small',ncol=4,loc=(-0.9,-0.75))
    ax.axis('off')
    ax.set_xlim(1,1.1)
    ax.set_ylim(1,1.1)


def make_legend_lum(ax):
    #ax.plot([0,0],[0,0],'-',color=[0.5,0.5,0.5],lw=0.5,label='M19 - Best Fit')
    ax.errorbar([0],[0],yerr=[1],color='k',marker='s',markersize=4,label='M19 - Median')
    ax.errorbar([0],[0],yerr=[1],color=[0.5,0.5,0.5],marker='o',markersize=4,label='L16 - Median')
    
    ax.plot([0,0],[0,0],'-',color=colors[2],label='Kawamata et al. (2018)')
    ax.plot([0,0],[0,0],'-',label='Huang et al. (2013)',color=colors[3])
    ax.errorbar(0,0,yerr=1,color=colors[3],marker='h',markersize=5,label='Holwerda et al. (2015)',linestyle='',elinewidth=0.5,capthick=0.5)
    ax.errorbar(0,0,yerr=1,xerr=1,color=colors[0],marker='o',markersize=5,label='Ono et al. (2013)',linestyle='',elinewidth=0.5,capthick=0.5)
    ax.errorbar(0,0,yerr=1,xerr=1,color=colors[4],marker='p',markersize=5,label='Grazian et al. (2012)',linestyle='',elinewidth=0.5,capthick=0.5)

    ax.legend(fontsize='small',ncol=4,loc=(-1.1,-0.75))
    ax.axis('off')
    ax.set_xlim(1,1.1)
    ax.set_ylim(1,1.1)


if __name__=="__main__":
  filename='paper1'
  filename_default='dragonsIII'#'dragons10'
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  maglims=[-22,-12.5]
  ylims=[-1.6,0.9]
  masslims=[7,10.6]
  
  zz=[5,6,7,8,9,10]
  fit_params=pd.DataFrame({'z':zz}, columns=['z','M_r','M_r_err','M_b','M_b_err','L_r_MC','L_r_err_MC','L_b_MC','L_b_err_MC']) 
  fit_params.set_index("z", inplace=True)

  fig, axes = plt.subplots(3,4,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1],'width_ratios':[4,4,4,0.4]})
  fig2, axes2 = plt.subplots(3,4,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1],'width_ratios':[4,4,4,0.4]})

  ii=-1
  j=0
  for snapshot in np.flip([37,43,52,63,78,100],0):
    ii+=1
    if ii==3:
      j+=1
      ii=0
    
    if (ii==2):
      cbar=axes2[j,3]
    else:
      cbar=False

    gals=load_data(filename,snapshot,'All',centrals=False)
    gals_default=load_data(filename_default,snapshot,'All',centrals=False)
    selection=(gals['Type']==0)&(gals['StellarDiskScaleLength']>0)
    selection_default=(gals_default['Type']==0)&(gals_default['DiskScaleLength']>0)
    gals=gals[selection]  
    gals_default=gals_default[selection_default]
  
    gals_no_cut=gals[gals['StellarDiskScaleLength']>0]#[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
    gals_no_cut_default=gals_default[gals_default['DiskScaleLength']>0]
    plot_hist2d(np.log10(gals_no_cut['StellarMass']*1e10),np.log10(gals_no_cut['StellarDiskScaleLength']*1000*1.67835),axes2[j,ii],masslims,ylims,cbar=cbar,cmax=2e3)
    plot_avg(np.log10(gals_no_cut_default['StellarMass']*1e10),np.log10(gals_no_cut_default['DiskScaleLength']*1000*1.67835),axes2[j,ii],masslims,bin_width=0.5,**{'color':[0.5,0.5,0.5],'marker':'o','markersize':4,'label':'L16 - Median'})
    plot_avg(np.log10(gals_no_cut['StellarMass']*1e10),np.log10(gals_no_cut['StellarDiskScaleLength']*1000*1.67835),axes2[j,ii],masslims,bin_width=0.5,**{'color':'k','marker':'s','markersize':4,'label':'M19 - Median'})
    plot_obs_mass(axes2[j,ii],redshift[snapshot])
    fit_equation_mass(np.log10(gals_no_cut['StellarMass']*1e10),np.log10(gals_no_cut['StellarDiskScaleLength']*1000*1.67835),axes2[j,ii],snapshot)

    #axes2[j,ii].legend(fontsize='small')
    
    mag=load_mags(filename,snapshot)
    mag=mag[selection]
    mag_default=load_mags(filename_default,snapshot)
    mag_default=mag_default[selection_default]
    #mag=mag[(gals['BulgeStellarMass']/gals['StellarMass']<0.3)&(gals['StellarDiskScaleLength']>0)]
    #gals=gals[(gals['BulgeStellarMass']/gals['StellarMass']<0.3)&(gals['StellarDiskScaleLength']>0)]
    #gals=gals[selection]
    #gals_default=gals_default[selection_default]

    if (ii==2):
      cbar=axes[j,3]
    else:
      cbar=False
   
    plot_hist2d(mag,np.log10(gals['StellarDiskScaleLength']*1000*1.67835),axes[j,ii],maglims,ylims,cbar=cbar,cmax=6e2)
    plot_avg(mag_default,np.log10(gals_default['DiskScaleLength']*1000*1.67835),axes[j,ii],maglims,bin_width=1,**{'color':[0.5,0.5,0.5],'marker':'o','markersize':4,'label':'L16 - Median'})
    plot_avg(mag,np.log10(gals['StellarDiskScaleLength']*1000*1.67835),axes[j,ii],maglims,bin_width=1,**{'color':'k','marker':'s','markersize':4,'label':'M19 - Median'})
    plot_obs_lum(axes[j,ii],redshift[snapshot])
    #fit_equation_lum(mag,np.log10(gals['StellarDiskScaleLength']*1000*1.67835),axes[j,ii],snapshot)
    fit_equation_lum_cut(mag[mag<-18],np.log10(gals['StellarDiskScaleLength'][mag<-18]*1000*1.67835),axes[j,ii],snapshot)


    #axes[j,ii].legend(fontsize='small')
    #axes[j,ii].scatter(mag,np.log10(gals['StellarDiskScaleLength']*1000),c=gals['BulgeStellarMass']/gals['StellarMass'])
    #axes[j,ii].set_xlim(xlims)
    #axes[j,ii].set_ylim(ylims)
    #plt.plot(mag,np.log10(gals['StellarDiskScaleLength']*1000))
    if j==1:
      axes[j,ii].set_xlabel('$M_{UV}$')
      axes2[j,ii].set_xlabel(r'$\log(M_{\ast \textrm{total} }/M_\odot)$')
    else: 
      axes[j,ii].set_xticklabels([])
      axes2[j,ii].set_xticklabels([])
    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(R_e/\mathrm{kpc})$')
      axes2[j,ii].set_ylabel(r'$\log(R_e/\mathrm{kpc})$')
    else:
      axes[j,ii].set_yticklabels([])
      axes2[j,ii].set_yticklabels([])
    #axes[j,ii].set_title('$z=${}'.format(redshift[snapshot]))
    axes[j,ii].text(-21.5, -1.4, r'$z={}$'.format(redshift[snapshot]))
    axes2[j,ii].text(9.5, -1.4, r'$z={}$'.format(redshift[snapshot]))
  #axes[j,ii].colorbar=plt.colorbar
  #cb=plt.colorbar(ax=axes)
  make_legend_lum(axes[2,1])
  make_legend_mass(axes2[2,1])
  axes2[2,0].axis('off')
  axes2[2,1].axis('off')
  axes2[2,2].axis('off')
  axes2[2,3].axis('off')
  axes[2,0].axis('off')
  axes[2,1].axis('off')
  axes[2,2].axis('off')
  axes[2,3].axis('off')

  fig.savefig('/home/mmarshal/results/plots/Paper1/SizeLuminosity.pdf',format='pdf')
  fig2.savefig('/home/mmarshal/results/plots/Paper1/SizeMass.pdf',format='pdf')
  
  fit_params=fit_params.astype(float)
  fit_params=fit_params.round({'M_r':3,'M_r_err': 3,'M_b':3,'M_b_err':3,'L_r':3,'L_r_err': 3,'L_b':3,'L_b_err': 3,'L_r_MC':3,'L_r_err_MC': 3,'L_b_MC':3,'L_b_err_MC': 3})
  with open('fit_params_SizeLuminosity.tex', 'w') as tf:
     tf.write(fit_params.to_latex())
  print(fit_params)
  
  plt.show()
