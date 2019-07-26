#Does a quick comparison with the Okumara+17 j/m for stellar disks.
#Ours isn't crazy (0.5 relative to their 0.77)
import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
from _load_data import load_data
import sys
from collections import OrderedDict
sys.path.append('/home/mmarshal/simulation_codes')

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (3.5,3.2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','8']*4
ls = OrderedDict(
    [('solid',               (0, ())),
     ('dotted',              (0, (1, 1.5))),
     ('densely dotted',      (0, (1, 1))),

     ('dashed',              (0, (5, 3))),
     ('densely dashed',      (0, (5, 1))),

     ('dashdotted',          (0, (3, 3, 1, 3))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 2, 1, 2, 1, 2)))])

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

def plot_hist2d(xdata,ydata,axes,xlims,ylims,cmax=None,cbar=False):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, vmin=1, vmax=cmax, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  if cbar:
    cb=plt.colorbar(img, cax=cbar,use_gridspec=True)
    cb.set_label('Number of Galaxies')#; Total N = {:.0e}'.format(np.size(gals)))


def plot_avg(xdata,ydata,axes,xlims,bin_width):
  min_bin=np.floor(xlims[0])
  max_bin=np.ceil(xlims[1])
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
    if np.size(y)>=20:
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
  
  axes.errorbar(bin_centre,avg_r,yerr=np.array([avg_r-pct_r_16,pct_r_84-avg_r]),color='k',marker='s',markersize=4,label='M19 - Median',zorder=1000)


def plot_jm_vs_m(AMratio,massRatio,samp):
  #Plot j/m against Mvir in mass bins (Okamura Fig 8 for 1 redshift)
  Mbins=[8.3,9.0,9.7,10.4,11.1]
  SM=np.log10(samp['StellarMass']*1e10)

  AMonm=np.zeros(len(Mbins)-1)
  Mdm=np.zeros(len(Mbins)-1)
  for b in range(0,len(Mbins)-1):
    AMonm[b]=np.nanmedian(AMratio[(SM>Mbins[b])&(SM<Mbins[b+1])]\
      /massRatio[(SM>Mbins[b])&(SM<Mbins[b+1])])
    Mdm[b]=np.nanmedian(np.log10(samp['Mvir'][(SM>Mbins[b])&(SM<Mbins[b+1])]*1e10))
  plt.plot(Mdm,AMonm,'s',color=colors[1])
  plt.xlabel(r'$\log(M_{vir})$') ##Need to -log(cosmo['h']) to compare with specific Okamura values
  plt.ylabel(r'$(J_{\ast}/M_{\ast})/(J_{\rm{H}}/M_{\rm{vir}}$')
  plt.plot([10.5,12.5],[np.nanmedian(AMratio/massRatio),np.nanmedian(AMratio/massRatio)],'k')
  plt.plot([10.5,12.5],[0.77]*2,color=[0.7,0.7,0.7])
  plt.legend(['M19','M19 - Median','Okamura et al. 2017'])
  plt.plot([10.5,12.5],[0.71]*2,'--',color=[0.7,0.7,0.7])
  plt.plot([10.5,12.5],[0.83]*2,'--',color=[0.7,0.7,0.7])
  plt.xlim([10.5,12.5])
  plt.ylim([0.45,1.3])
  plt.show()


def plot_y_vs_x_hist(y,x,axes,type='MassRatio',range=None,cbar=False):
  #type can be 'MassRatio' 'StellarMass' or 'VirialMass'
  if not range:
    range=([np.log10(min(x)),np.log10(max(x))],[np.log10(min(y)),np.log10(max(y))])
  xlims=range[0]
  ylims=range[1]
  axes.plot(ylims,ylims,'--',color=[0.5,0.5,0.5])
  
  if type=='MassRatio':
    l2 = np.array((-1.5, -1.7))
    # Rotate angle
    angle = 40
    trans_angle = axes.transData.transform_angles(np.array((45,)),
                                                   l2.reshape((1,2)))[0]
    th2 = axes.text(l2[0], l2[1], r'$j_\ast=m_\ast$',
               rotation=trans_angle, rotation_mode='anchor')

    Okamura_dots_m=[0.0016340608,0.0024865912,0.0043960297,0.005342158,0.009234321,0.010294778,0.009659196,0.016571935,0.017795298]
    Okamura_dots_j=[5.462987E-4,0.0022105814,0.0033079053,0.004241901,0.008070309,0.0069755013,0.011275648,0.013385341,0.013616899]
    Okamura_stars_m=[0.0026040417,0.0030831706,0.0039646695,0.0061275386,0.006909487,0.008304536,0.013176631,0.014044758,0.016943919,0.022368314]
    Okamura_stars_j=[0.0021452166,0.0026581592,0.0044850684,0.0048035714,0.0056536347,0.007831677,0.009872188,0.010304701,0.013616899,0.014646558]
    axes.plot(np.log10(Okamura_dots_m),np.log10(Okamura_dots_j),'o',label='Okamura et al. (2017)\n(clustering analysis)',color=colors[0])
    axes.plot(np.log10(Okamura_stars_m),np.log10(Okamura_stars_j),'*',label='Okamura et al. (2017)\n(abundance matching)',color=colors[0])  
    
    selection=(y>0)&(x>0)
    plot_hist2d(np.log10(x[selection]),np.log10(y[selection]),axes,xlims,ylims,cbar=cbar)
    plot_avg(np.log10(x[selection]),np.log10(y[selection]),axes,xlims,0.5)
  
  elif type=='StellarMass':
    #Lapi+18 -  km s−1 kpc #Disk-dominated galaxies, seem to use total stellar mass
    logM=np.linspace(9,11.5)
    jRatio=0.7212 -0.1963*(logM-10.5) -0.0884*(logM-10.5)**2+0.0492*(logM-10.5)**3
    axes.plot(logM,jRatio,color=colors[0],label='Lapi et al. (2018)',linewidth=2,linestyle=ls['densely dotted'])
    lower_err=jRatio-np.ones_like(jRatio)*0.2
    upper_err=jRatio+np.ones_like(jRatio)*0.2
    axes.fill_between(logM,(lower_err),(upper_err),color=colors[0],alpha=0.4)

    selection=(y>0)&(x>0)
    plot_hist2d(np.log10(x[selection]),y[selection],axes,xlims,ylims,cmax=2500,cbar=cbar)
    plot_avg(np.log10(x[selection]),y[selection],axes,xlims,0.5)
 
  elif type=='VirialMass':
    #Lapi+18 -  km s−1 kpc
    logM=np.linspace(11.2,12.1)
    jRatio=0.6865 -0.2666*(logM-12) -0.07258*(logM-12)**2+0.1277*(logM-12)**3
    lower_err=jRatio-np.ones_like(jRatio)*0.2
    upper_err=jRatio+np.ones_like(jRatio)*0.2
    axes.plot(logM,jRatio,color=colors[0],label='Lapi et al. (2018)',linewidth=2,linestyle=ls['densely dotted'])
    axes.fill_between(logM,lower_err,upper_err,color=colors[0],alpha=0.4)

    #Dutton & van den Bosch 2012
    logM=np.array([11.3,12.7])
    jRatio=np.array([0.61,0.61])
    lower_err=jRatio-np.array([0.11,0.11])
    upper_err=jRatio+np.array([0.13,0.13])
    axes.plot(logM,jRatio,color=colors[2],label='Dutton \& van den Bosch (2012)',linewidth=2,linestyle=ls['densely dashed'])
    axes.fill_between(logM,lower_err,upper_err,color=colors[2],alpha=0.4)
  
    ##Okamura+17
    #axes.plot([10.5,12.5],[0.83]*2,color=colors[3],label='Okamura et al. (2017)',linewidth=2.5)
    #axes.fill_between([10.5,12.5],[0.70]*2,[0.96]*2,color=colors[3],alpha=0.4)

    selection=(y>0)&(x>0)
    plot_hist2d(np.log10(x[selection]),y[selection],axes,xlims,ylims,cmax=2500,cbar=cbar)
    plot_avg(np.log10(x[selection]),y[selection],axes,xlims,0.5)


def plot_jstar_vs_mstar(AM,x,axes,range=None,cbar=False):
  #AM*=1e3 #Mpc->kpc
  y=(AM*1e3)/(x*1e-10) #Specific AM in km/s kpc /Mo
  if not range:
    range=([np.log10(min(x)),np.log10(max(x))],[np.log10(min(y)),np.log10(max(y))])
  xlims=range[0]
  ylims=range[1]
  selection=(y>0)&(x>0)
  plot_hist2d(np.log10(x[selection]),np.log10(y[selection]),axes,xlims,ylims,cbar=cbar)
  plot_avg(np.log10(x[selection]),np.log10(y[selection]),axes,xlims,0.5)

  #Romanowsky&Fall+13
  #for disks
  obs_mass=np.array([9,11.5])
  axes.plot(obs_mass,-3.04+0.60*(obs_mass),label='Fall \& Romanowsky (2013)\n(discs)',color=colors[0],linewidth=1.5,linestyle=ls['densely dashdotted'])
  
  #Harrison+2017
  logM=np.array([9,11])
  logj=2.59+0.6*(logM-10.10)
  lower_err=logj-0.2
  upper_err=logj+0.2
  axes.plot(logM,logj,color=colors[7],label='Harrison et al. (2017)\n($z\simeq 0.9$)',linewidth=2,linestyle=ls['dashdotdotted'])
  axes.fill_between(logM,lower_err,upper_err,color=colors[7],alpha=0.4)
 
  
  #Lapi+18 -  km s−1 kpc
  logM=np.linspace(9,11.5)
  logj=2.9834 + 0.4951*(logM-10.5) +0.0572*(logM-10.5)**2+0.0140*(logM-10.5)**3
  axes.plot(logM,logj,color=colors[3],label='Lapi et al. (2018)',linewidth=1.5,linestyle=ls['densely dotted'])
 
  #Sweet+18, same mass range as Lapi
  logj=0.56*logM-2.76
  axes.plot(logM,logj,color=colors[4],label='Sweet et al. (2018)',linewidth=2,linestyle=ls['dotted'])

  #Posti+18
  ##Disks only (disk mass, j)
  #logM=np.linspace(7,11.5)
  #logj=0.58*(logM-11)+3.43
  #err=np.ones_like(logM)*0.15
  ##Spiral galaxies (total mass, j)
  logM=np.linspace(7,11.5)
  logj=0.55*(logM-11)+3.34
  err=np.ones_like(logM)*0.17
  axes.plot(logM,logj,color=colors[2],label='Posti et al. (2018)',linewidth=1.5,linestyle=ls['dashed'])
  axes.fill_between(logM,logj-err,logj+err,color=colors[2],alpha=0.4)
  #axes.errorbar(logM,logj,err,color=colors[5],label='Posti et al. (2018)',linewidth=2.5)
 
  #Burkert+16
  obs_mass=np.array([9.3,11.8])
  axes.plot(obs_mass,3.33+(2/3)*(obs_mass-11),label='Burkert et al. (2016)\n($z\simeq 0.8$-2.6, deredshifted)',color=colors[6],linewidth=1.5,linestyle=ls['dashdotted'])
   
  #Fall&Romanowsky+18
  # log j*/j0 = a log M*/M0,  logM0=10.5 (j in kpc km/s, not Mpc!!)
  #		*  a =0.58 +- 0.10 and log j0 =3.07 +- 0.03 for pure disks
  #		*  a = 0.83+- 0.16 and log j0 = 2.20 +- 0.1 for pure bulges
  obs_mass=np.array([9.5,11.5])
  axes.plot(obs_mass,3.07+0.58*(obs_mass-10.5),label='Fall \& Romanowsky (2018)\n(discs)',color=colors[5],linewidth=1.5,linestyle=ls['densely dashed'])
  #axes.plot(obs_mass,2.20+0.83*(obs_mass-10.5),'k--',label=r'Fall \& Romanowsky (2018) - bulges')



def plot_morphology_obs(axes):
  #Rom & Fall 2018
  obs_mass=np.array([9.5,11.5])
  axes.plot(obs_mass,3.07+0.58*(obs_mass-10.5),'--',label='Fall \& Romanowsky (2018) (discs)',color='k',linewidth=2)
  axes.plot(obs_mass,2.20+0.83*(obs_mass-10.5),':',label='Fall \& Romanowsky (2018) (bulges)',color='k',linewidth=2)


def plot_gas_obs(axes):
  #D. Obreschkow and K. Glazebrook
  m_ast=np.array([10.1,9.9,9.7,10.8,9.1,10.3,10.1,10.4,10.7,10.6,10.3,10.8,10.6,10.5,10.9,9.5])
  m_gas=np.array([9.8,9.8,9.5,10.1,8.4,9.7,10.1,9.4,10.1,9.4,9.0,10.2,9.8,10.0,10.2,9.2])
  m_bary=np.array([10.27,10.16,9.91,10.88,9.18,10.41,10.41,10.44,10.81,10.62,10.32,10.91,10.66,10.62,10.99,9.69])
  j_ast=np.array([2.98,2.94,2.62,3.40,2.03,3.03,2.97,2.91,3.06,2.84,2.34,3.18,3.18,3.02,3.35,2.43]) #log kpc km /s
  #axes.scatter(m_ast,j_ast,c=10**m_gas/(10**m_gas+10**m_ast),s=20,cmap=cm,vmin=0,vmax=1,edgecolors='k')
  axes.scatter(m_ast,j_ast,c=(10**m_bary-10**m_ast)/(10**m_bary),s=30,cmap=cm,vmin=0,vmax=1,edgecolors='k',label='Obreschkow and Glazebrook (2014)')
  

if __name__=='__main__':
  filename='paper1_T125'#str(sys.argv[1])#'bulges_update0915_ddsf'
  snapshot=250
  
  gals=load_data(filename,snapshot,props='All',centrals=True)
  gals=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
  gals=gals[gals['StellarDiskScaleLength']>0] 
  #gals=gals[(gals['Sfr']>0)]
 
  #These are J not j: (in km/s Mpc) 
  AM=np.sqrt(gals['AMstars'][:,0]**2+gals['AMstars'][:,1]**2+gals['AMstars'][:,2]**2)
  AMhalo=gals['Spin']*np.sqrt(2)*gals['Mvir']*gals['Vvir']*gals['Rvir']
  
  #Convergence selection
  #selection=(AMhalo>1e-1)
  #gals=gals[selection]
  #AM=AM[selection]
  #AMhalo=AMhalo[selection]

  AMratio=AM/AMhalo
  diskMass=(gals['StellarMass']-gals['BulgeStellarMass'])*1e10
  stellarMass=gals['StellarMass']*1e10
  diskMassRatio=diskMass/(gals['Mvir']*1e10)
  totalMassRatio=gals['StellarMass']/gals['Mvir']

  #plot_jm_vs_m(AMratio,massRatio,gals)

  #Plot J*/Jh vs M*d/Mh
  #fig,axes=plt.subplots(1,1)
  #plot_y_vs_x_hist(AMratio,massRatio,axes,type='MassRatio')
  #axes.set_xlabel(r'$\log(M_{\ast~\rm{disc}}/M_{\rm{DM}})$')
  #axes.set_ylabel(r'$J_\ast/J_{\rm{DM}}$')
  #axes.legend()
  ##plt.savefig('/home/mmarshal/results/plots/AngularMomentumMass.pdf', format='pdf')
  #plt.show()
  
  #Plot j*/jh vs M*d
  fig1,axes1=plt.subplots(1,4,gridspec_kw = {'wspace':0,'width_ratios':[4,4,0.4,4.8]},figsize=(7.2,3)) 
  plot_y_vs_x_hist(AMratio/totalMassRatio,stellarMass,axes1[0],type='StellarMass',range=([7,11.9],[0,1.4]))
  axes1[0].plot([8,8],[0,1.5],'--',color=[0.5,0.5,0.5])
  axes1[0].set_xlabel(r'$\log(M_{\ast})$')
  axes1[0].set_ylabel(r'$j_\ast/j_{\rm{H}}$')
  ##plt.savefig('/home/mmarshal/results/plots/AngularMomentumMass.pdf', format='pdf')
  #plt.show()
  
  #Plot j*/jh vs Mh
  plot_y_vs_x_hist(AMratio/totalMassRatio,gals['Mvir']*1e10,axes1[1],type='VirialMass',range=([9.8,13],[0,1.4]),cbar=axes1[2])
  print("Maximum j/jh = {}".format(max(AMratio/totalMassRatio)))
  print("Minimum j/jh = {}".format(min(AMratio/totalMassRatio)))
  axes1[1].plot([10.3,10.3],[0,1.5],'--',color=[0.5,0.5,0.5],label='Convergence Limits')
  axes1[1].set_xlabel(r'$\log(M_{\rm{vir}})$')
  axes1[1].set_yticklabels([])
  axes1[3].axis('off')
  #axes[1].set_yticklabels([])#set_ylabel(r'$j_\ast/j_{\rm{DM}}$')
  axes1[1].legend(fontsize='small',loc=(1.35,0.3))
  plt.tight_layout()
  fig1.savefig('/home/mmarshal/results/plots/Paper1/AMRatioMvir.pdf', format='pdf')

  #Plot j* vs M*
  fig2,axes2=plt.subplots(2,2,figsize=(3.75,4.1),gridspec_kw ={'height_ratios':[3.5,1],'width_ratios':[4,0.3],'wspace':0})
  #plot_jstar_vs_mstar(AM,diskMass,axes2[0,0],range=([7,12],[0,4.2]),cbar=axes2[0,1])
  plot_jstar_vs_mstar(AM,stellarMass,axes2[0,0],range=([7,12],[0,4.2]),cbar=axes2[0,1])
  axes2[1,0].axis('off')
  axes2[1,1].axis('off')
  axes2[0,0].plot([8.65,8.65],[0,5],'--',color=[0.5,0.5,0.5])
  #axes2[0,0].plot([7,12],[1.6,1.6],'--',color=[0.5,0.5,0.5],label='Convergence Limits')
  axes2[0,0].legend(fontsize='small',loc=(-0.18,-0.59),ncol=2)
  axes2[0,0].set_xlabel(r'$\log(M_\ast/M_\odot)$')
  axes2[0,0].set_ylabel(r'$\log(j_\ast (\rm{km~kpc/s}))$')
  fig2.subplots_adjust(left=0.125,right=0.85,top=0.97)
  #plt.tight_layout()
  fig2.savefig('/home/mmarshal/results/plots/Paper1/AngularMomentumMass.pdf', format='pdf')

  from matplotlib.colors import LinearSegmentedColormap
  cm=LinearSegmentedColormap.from_list('maddie', [(0,'#d9bcdd'),(0.3,'#b683be'),(0.7,colors[3]),(1,'#6a3672')], N=30)

  #Plot J* vs M*, colour by B/T
  fig3,axes3=plt.subplots(2,2,gridspec_kw={'width_ratios':[4,0.3],'height_ratios':[4,1],'wspace':0,'hspace':0})
  gals=load_data(filename,snapshot,props='All',centrals=True)
  AM=np.sqrt(gals['AMstars'][:,0]**2+gals['AMstars'][:,1]**2+gals['AMstars'][:,2]**2)
  AMhalo=gals['Spin']*np.sqrt(2)*gals['Mvir']*gals['Vvir']*gals['Rvir']
  stellarMass=gals['StellarMass']*1e10
  
  #Convergence selection
  #selection=(AMhalo>1e-1)
  #gals=gals[selection]
  #AM=AM[selection]
  #AMhalo=AMhalo[selection]
  specificAM=AM/((gals['StellarMass']))
  im=axes3[0,0].scatter(np.log10(stellarMass),np.log10(specificAM)+3,c=gals['BulgeStellarMass']/gals['StellarMass'],s=2,cmap=cm,vmin=0,vmax=1)
  plot_morphology_obs(axes3[0,0])
  axes3[0,0].plot([8.65,8.65],[-4,6],'--',color=[0.5,0.5,0.5],label='Convergence Limit')
  axes3[0,0].legend(fontsize='small',loc=(0.05,-0.41))

  cb=plt.colorbar(im,cax=axes3[0,1],use_gridspec=True)
  cb.set_label(r'B/T')
  axes3[0,0].set_xlabel(r'$\log(M_\ast/M_\odot)$')
  axes3[0,0].set_ylabel(r'$\log(j_\ast (\rm{km~kpc/s}))$')
  axes3[1,0].axis('off')
  axes3[1,1].axis('off')
  axes3[0,0].set_ylim(-4,6)
  axes3[0,0].set_xlim(7,12)
  plt.tight_layout()
  fig3.subplots_adjust(left=0.125,right=0.85,top=0.97)
  plt.savefig('/home/mmarshal/results/plots/Paper1/AngularMomentumMass_BT.png', format='png',dpi=500)
  
  fig4,axes4=plt.subplots(2,2,gridspec_kw={'width_ratios':[4,0.3],'height_ratios':[4,1],'wspace':0,'hspace':0})
  im=axes4[0,0].scatter(np.log10(stellarMass),np.log10(specificAM)+3,c=gals['ColdGas']/(gals['StellarMass']+gals['ColdGas']),s=2,cmap=cm)
  plot_gas_obs(axes4[0,0])
  axes4[0,0].plot([8.65,8.65],[-4,6],'--',color=[0.5,0.5,0.5],label='Convergence Limit')
  axes4[0,0].legend(fontsize='small',loc=(0.05,-0.41))

  cb=plt.colorbar(im,cax=axes4[0,1],use_gridspec=True)
  cb.set_label(r'$f_{\rm{gas}}$')
  axes4[0,0].set_xlabel(r'$\log(M_\ast/M_\odot)$')
  axes4[0,0].set_ylabel(r'$\log(j_\ast (\rm{km~kpc/s}))$')
  axes4[1,0].axis('off')
  axes4[1,1].axis('off')
  axes4[0,0].set_ylim(-4,6)
  axes4[0,0].set_xlim(7,12)
  plt.tight_layout()
  fig4.subplots_adjust(left=0.125,right=0.85,top=0.97)
  plt.savefig('/home/mmarshal/results/plots/Paper1/AngularMomentumMass_fgas.png', format='png',dpi=500)
  plt.show()
