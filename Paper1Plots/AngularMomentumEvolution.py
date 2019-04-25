#Does a quick comparison with the Okumara+17 j/m for stellar disks.
#Ours isn't crazy (0.5 relative to their 0.77)
import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
from _load_data import load_data
import sys
sys.path.append('/home/mmarshal/simulation_codes')

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.6,3)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','8']*4

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

def plot_hist2d(xdata,ydata,axes,xlims,ylims):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  plt.tight_layout()


def plot_avg(xdata,ydata,axes,xlims,bin_width,redshiftPlot=False):
  min_bin=np.floor(xlims[0])-bin_width
  max_bin=np.floor(xlims[1])+bin_width
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
      avg_r[bin_num]=np.nanmedian(y)
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
  
  if redshiftPlot:
    if snapshot<158:
      axes.plot(bin_centre,avg_r,color=color[snapshot],label=r'$z={}$'.format(redshift[snapshot]))
    else:
      axes.plot(bin_centre,avg_r,color=color[snapshot],label=r'$z={}$'.format(redshift[snapshot]))
      axes.fill_between(bin_centre,pct_r_16,pct_r_84,color=color[snapshot],alpha=0.2)
  else:
    axes.errorbar(bin_centre,avg_r,yerr=np.array([avg_r-pct_r_16,pct_r_84-avg_r]),color='k',marker='s',markersize=4,label='M19 - Median')
  #return bin_centre,avg_r,np.array([avg_r-pct_r_16,pct_r_84-avg_r])


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
  plt.legend(['M19','M19 - Median','Okamura et al. (2018)'])
  plt.plot([10.5,12.5],[0.71]*2,'--',color=[0.7,0.7,0.7])
  plt.plot([10.5,12.5],[0.83]*2,'--',color=[0.7,0.7,0.7])
  plt.xlim([10.5,12.5])
  plt.ylim([0.45,1.3])
  plt.show()


def plot_y_vs_x_hist(y,x,axes,type='MassRatio',range=None):
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
    axes.plot(np.log10(Okamura_dots_m),np.log10(Okamura_dots_j),'o',label='Okamura et al. (2018)\n(clustering analysis)',color=colors[0])
    axes.plot(np.log10(Okamura_stars_m),np.log10(Okamura_stars_j),'*',label='Okamura et al. (2018)\n(abundance matching)',color=colors[0])  
    
    selection=(y>0)&(x>0)
    plot_hist2d(np.log10(x[selection]),np.log10(y[selection]),axes,xlims,ylims)


def plot_y_vs_x(y,x,axes,MassRatio=1,Okamura_Fig8=0,range=None,redshiftAxes=None):
  if not range:
    range=([np.log10(min(x)),np.log10(max(x))],[(min(y)),(max(y))])
  xlims=range[0]
  ylims=range[1]
  axes.plot(ylims,ylims,'--',color=[0.5,0.5,0.5])
  #axes.set_xlim([-3.5,-1.2])
  #axes.set_ylim([-3.5,-1.2])
  
  if MassRatio:
    l2 = np.array((-1.5, -1.7))
    #plt.tight_layout()
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
    axes.plot(np.log10(Okamura_dots_m),np.log10(Okamura_dots_j),'o',label='Okamura et al. (2018)\n(clustering analysis)',color=colors[0])
    axes.plot(np.log10(Okamura_stars_m),np.log10(Okamura_stars_j),'*',label='Okamura et al. (2018)\n(abundance matching)',color=colors[0])  
 
  if Okamura_Fig8:
    #abundance matching
    if redshift[snapshot]==4:
      jRatio=[1.13, 0.96, 0.80] 
      logm=[11.22, 11.45, 11.78]
    elif redshift[snapshot]==3:
      jRatio=[0.86, 0.82, 0.75] 
      logm=[11.29,11.53, 11.81]
    elif redshift[snapshot]==2:
      jRatio=[0.81, 0.79, 0.73, 0.66]
      logm=[11.30,11.51,11.79,12.23]
    if redshift[snapshot] in [2,3,4]: 
      axes.plot(np.array(logm),jRatio,'*',label='Okamura et al. (2018)\n(Abundance Matching)',color=color[snapshot])
      if redshiftAxes:
        if redshift[snapshot]==2:
          lab='Okamura et al. (2018)\n(Abundance Matching)'
        else:
          lab='__nolabel__'
        redshiftAxes.plot(np.array(logm),jRatio,'*',label=lab,color=color[snapshot],markersize=4)

    #clustering
    if redshift[snapshot]==4:
      jRatio=[0.79,0.80,0.62] 
      logm=[11.64,11.79,12.03]
      logm_low=[11.13,11.23,11.12]
      logm_up=[11.94,12.10,12.41]
      jRatio_low=[0.67,0.68,0.51]
      jRatio_up=[1.07,1.12,1.02]
    elif redshift[snapshot]==3:
      jRatio=[1.16,0.87,0.69] 
      logm=[10.79,11.40,11.92]
      logm_low=[10.23,11.19,11.69]
      logm_up=[11.35,11.57,12.08]
      jRatio_low=[0.82,0.79,0.63]
      jRatio_up=[1.50,0.99,0.78]
    elif redshift[snapshot]==2:
      jRatio=[0.89,0.75,0.78]
      logm=[11.32,11.65,11.69]
      logm_low=[11.09,11.53,11.12]
      logm_up=[11.51,11.76,12.01]
      jRatio_low=[0.79,0.71,0.66]
      jRatio_up=[1.02,0.81,1.08]
    if redshift[snapshot] in [2,3,4]: 
      logm_err=[(np.array(logm)-np.array(logm_low)),(np.array(logm_up)-np.array(logm))]
      jRatio_err=[np.array(jRatio)-np.array(jRatio_low),np.array(jRatio_up)-np.array(jRatio)]
      axes.errorbar(np.array(logm),jRatio,jRatio_err,logm_err,'o',label='Okamura et al. (2018)\n(Clustering Analysis)',color=colors[2])
      if redshiftAxes:
        if redshift[snapshot]==2:
          lab='Okamura et al. (2018)\n(Clustering Analysis)'
        else:
          lab='__nolabel__'
        redshiftAxes.errorbar(np.array(logm),jRatio,jRatio_err,logm_err,'o',label=lab,color=color[snapshot],markersize=4)


  selection=(y>0)&(x>0)
  #print(np.log10(massRatio[selection]))
  #axes.plot(np.log10(x[selection]),np.log10(y[selection]),'.')
  plot_hist2d(np.log10(x[selection]),(y[selection]),axes,xlims,ylims)
  #import ContourPlot as cp
  #cp.contour_plot(np.log10(x[selection]),np.log10(y[selection]),'$log(m_*)$','$log(j_*)$',xlims,ylims,linewidth=0.9)
  plot_avg(np.log10(x[selection]),y[selection],axes,xlims,0.25)
    
  
  if redshiftAxes:
    plot_avg(np.log10(x[selection]),(y[selection]),redshiftAxes,xlims,0.25,redshiftPlot=True)
    redshiftAxes.set_ylim(ylims)
    redshiftAxes.set_xlim(xlims)
    
  #return x,y,y_err
  #axes.set_ylabel(r'$J_\ast/J_{\rm{DM}}$')
  #plt.legend()


def plot_jstar_vs_mstar(AM,x,axes,range=None):
  #AM*=1e3 #Mpc->kpc
  y=(AM*1e3)/(x*1e-10) #Specific AM
  if not range:
    range=([np.log10(min(x)),np.log10(max(x))],[np.log10(min(y)),np.log10(max(y))])
  xlims=range[0]
  ylims=range[1]
  selection=(y>0)&(x>0)
  plot_hist2d(np.log10(x[selection]),np.log10(y[selection]),axes,xlims,ylims)
  plot_avg(np.log10(x[selection]),np.log10(y[selection]),axes,xlims,0.25)


if __name__=='__main__':
  filename='paper1'#str(sys.argv[1])#'bulges_update0915_ddsf'
  filename_T125='paper1_T125'
  
  fig1,axes1=plt.subplots(3,3,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1]})
  fig2,axes2=plt.subplots(3,3,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1]})
  fig3,axes3=plt.subplots(3,3,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1]})
  fig4,axes4=plt.subplots(1,5,gridspec_kw = {'wspace':0, 'hspace':0,'width_ratios':[3,1.03,3,3,2.1],'left':0.08})

#4)#,sharey=True,gridspec_kw={'wspace':0})
  
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,116:4,134:3,158:2,192:1,250:0}
  color={63:'#e41a1c',78:'#984ea3',100:'#377eb8',116:'#4daf4a',\
                  134:'#ff7f00',158:'#f781bf',192:'#a65628',250:'black'}
 
 
  ii=-1
  j=0
  for snapshot in np.flip([63,78,100,116,134,158],0):#,192,250],0):
    ii+=1
    if ii==3:
      j+=1
      ii=0
    if snapshot<160:
      gals=load_data(filename,snapshot,props='All',centrals=True)
    else:
      gals=load_data(filename_T125,snapshot,props='All',centrals=True)

    #gals=gals[gals['Type']==0]
    #gals=gals[(gals['Mvir']>1e1)&(gals['Mvir']<10**1.2)]

    gals=gals[gals['BulgeStellarMass']/gals['StellarMass']<0.3]
    gals=gals[gals['StellarDiskScaleLength']>0] 
    #gals=gals[(gals['Sfr']>0)]
  
    AM=np.sqrt(gals['AMstars'][:,0]**2+gals['AMstars'][:,1]**2+gals['AMstars'][:,2]**2)
    AMhalo=gals['Spin']*np.sqrt(2)*gals['Mvir']*gals['Vvir']*gals['Rvir']
    
    #Convergence selection
    selection=(AMhalo>1e-2) #Tiamat unresolved for $j_H\sim10^{1} \textrm{km kpc s}^{-1}$, this is in Mpc
    gals=gals[selection]
    AM=AM[selection]
    AMhalo=AMhalo[selection]


    AMratio=AM/AMhalo
    diskMass=(gals['StellarMass']-gals['BulgeStellarMass'])*1e10
    massRatio=diskMass/(gals['Mvir']*1e10)


    #plot_jm_vs_m(AMratio,massRatio,gals)

    plot_y_vs_x(np.log10((AM*1e3)/(diskMass*1e-10)),diskMass,axes1[j,ii],MassRatio=0,range=([7,9],[1,2.7]),redshiftAxes=axes4[0])
    #plt.savefig('/home/mmarshal/results/plots/AngularMomentumMass.pdf', format='pdf')
  
    plot_y_vs_x(AMratio/massRatio,diskMass,axes2[j,ii],MassRatio=0,Okamura_Fig8=0,range=([7,9],[0.25,1.1]),redshiftAxes=axes4[2])
    #plt.savefig('/home/mmarshal/results/plots/AngularMomentumMass.pdf', format='pdf')
  
    plot_y_vs_x(AMratio/massRatio,gals['Mvir']*1e10,axes3[j,ii],MassRatio=0,Okamura_Fig8=0,range=([9.5,12.3],[0.25,1.1]),redshiftAxes=axes4[3])

    
    if j==1:
      axes1[j,ii].set_xlabel(r'$\log(M_{\ast~\rm{disc}})$')
      axes2[j,ii].set_xlabel(r'$\log(M_{\ast~\rm{disc}})$')
      axes3[j,ii].set_xlabel(r'$\log(M_{\rm{vir}})$')
    else: 
      axes1[j,ii].set_xticklabels([])
      axes2[j,ii].set_xticklabels([])
      axes3[j,ii].set_xticklabels([])
    if ii==0:
      axes1[j,ii].set_ylabel(r'$\log(j_\ast (\rm{kpc km s}^{-1})$')
      axes2[j,ii].set_ylabel(r'$j_\ast/j_{\rm{H}}$')
      axes3[j,ii].set_ylabel(r'$j_\ast/j_{\rm{H}}$')
      #axes[j,ii].set_ylabel(r'$\log(R_e/\mathrm{kpc})$')
      #axes2[j,ii].set_ylabel(r'$\log(R_e/\mathrm{kpc})$')
    else:
      axes1[j,ii].set_yticklabels([])
      axes2[j,ii].set_yticklabels([])
      axes3[j,ii].set_yticklabels([])
    #axes[j,ii].set_title('$z=${}'.format(redshift[snapshot]))
    #axes1[j,ii].text(-1, -3.5, r'$z={}$'.format(redshift[snapshot]))
    axes2[j,ii].text(9, -3.5, r'$z={}$'.format(redshift[snapshot]))
    axes3[j,ii].text(11, -3.5, r'$z={}$'.format(redshift[snapshot]))
    
  ##Okamura+17
  axes4[3].plot([10.9,12.3],[0.77]*2,'--',color=colors[4],label='Okamura et al. (2018)')
  axes4[3].plot([10.9,12.3],[0.70]*2,':',color=colors[4],label='__nolegend')
  axes4[3].plot([10.9,12.3],[0.96]*2,':',color=colors[4],label='__nolegend__')


  axes4[0].set_ylabel(r'$\log(j_\ast (\rm{kpc km s}^{-1})$')
  axes4[2].set_ylabel(r'$j_\ast/j_{\rm{H}}$')
  #axes4[3].set_ylabel(r'$j_\ast/j_{\rm{DM}}$')
  axes4[3].set_yticks([])  
  axes4[0].set_xlabel(r'$\log(M_{\ast~\rm{disc}})$')
  axes4[2].set_xlabel(r'$\log(M_{\ast~\rm{disc}})$')
  axes4[3].set_xlabel(r'$\log(M_{\rm{vir}})$')
  axes4[1].axis('off')
  axes4[4].axis('off')
  axes4[3].legend(loc=(1.03,0.1))
  axes3[2,0].axis('off')
  axes3[2,1].axis('off')
  axes3[2,2].axis('off')
  axes2[2,0].axis('off')
  axes2[2,1].axis('off')
  axes2[2,2].axis('off')
  axes1[2,0].axis('off')
  axes1[2,1].axis('off')
  axes1[2,2].axis('off')

  fig4.subplots_adjust(left=0.125,right=0.9,bottom=0.15)
  fig4.savefig('/home/mmarshal/results/plots/Paper1/AngularMomentumEvolution.pdf', format='pdf')
  plt.show()
