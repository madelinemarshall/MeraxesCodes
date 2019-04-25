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
from _load_data import load_data

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (6.6,3.2)
matplotlib.rcParams['lines.linewidth'] = 2.5
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4

cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
}


def plot_observations(axes,color,mass):
  #Kormendy & Ho (2013):
  #MBH/10^9=(0.49\pm0.6)(Mbulge/10^11)^(1.17\pm0.08), intrinsic scatter 0.28 dex (p571)
  #logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
  #axes.errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)\n- Fit',capsize=3,linewidth=2.5, color=colors[1],zorder=0)
  
  kh=pd.read_csv('/home/mmarshal/simulation_codes/data/KormendyHo_ellipticals_MBHMbulge.csv')
  logMBH=np.log10(kh['log(M_BH) [Msun]'])
  logMbulge=kh['log(M) [Msun]']-2*np.log10(cosmo['h']/0.705)
  axes.plot(logMbulge,logMBH,'o',color=colors[1],label="Kormendy \& Ho (2013)\n- Ellipticals",markersize=3)
  
  kh=pd.read_csv('/home/mmarshal/simulation_codes/data/KormendyHo_classicalbulges_MBHMbulge.csv')
  logMBH=np.log10(kh['log(M_BH) [Msun]'])
  logMbulge=kh['log(M) [Msun]']-2*np.log10(cosmo['h']/0.705)
  axes.plot(logMbulge,logMBH,'^',color=colors[1],label="Kormendy \& Ho (2013)\n- Classical Bulges",markersize=3)

  if mass == 'stellar':
    rv=pd.read_csv('/home/mmarshal/simulation_codes/data/ReinesVolonteri.csv',comment='#')
    logMstar=rv['logM*[Msun]'] -2*np.log10(cosmo['h']/0.70)
    logMBH=rv['logMBH']
    axes.plot(logMstar,logMBH,'o',color=colors[0],label="Reines \& Volonteri (2015)",markersize=3)

  else:
    scott=pd.read_csv('/home/mmarshal/simulation_codes/data/Scott2013_MBHMbulge.csv')
    logMBH=np.log10(scott['M_BH [$10^8$ M$_\odot$]']*1e8)
    logMbulge=np.log10(scott['M_sph [$10^{10}$ M$_\odot$]']*1e10)
    axes.plot(logMbulge,logMBH,'s',label="Scott et al. (2013)",color=colors[-2],markersize=3)

    #jiang=pd.read_csv('/home/mmarshal/simulation_codes/data/Jiang2011_MBHMbulge.csv')
    #logMBH=jiang['logMbh [Msun]']
    #logMbulge=np.log10(jiang['Msph [Msun]']) -1*np.log10(cosmo['h']/0.71)
    #axes.plot(logMbulge,logMBH,'>',color=colors[3],label="Jiang et al. (2011)",markersize=3)

    grah=pd.read_csv('/home/mmarshal/simulation_codes/data/GrahamScott2015_MBHMbulge.csv')
    logMBH=np.log10(grah['Mbh [10+5Msun]']*1e5)
    logMbulge=np.log10(grah['Msph [GMsun]']*1e9)
    axes.plot(logMbulge,logMBH,'o',color=colors[4],label="Graham and Scott (2015)",markersize=3)
    axes.plot(0,0,'o',color=colors[0],label="Reines \& Volonteri (2015)",markersize=3)


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


def plot_MBHMstellar(filename,snapshots,mass_bulge,split,bulge_type,contours,color,axes):
  slope=np.zeros(len(snapshots))
  slope_errs=np.zeros(len(snapshots))
  inter=np.zeros(len(snapshots))
  inter_errs=np.zeros(len(snapshots))
  ii=0
  for snap in snapshots:
    gals=load_data(filename,snap,['GhostFlag','Mvir','StellarMass','BlackHoleMass','CentralGal','BulgeStellarMass','MergerBulgeStellarMass'],centrals=True)
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
    
    cp.contour_plot(logMstel,logMBH,xlab=None,ylab=None,xlims=[9.1,11.6],ylims=[6,9.5],axes=axes,colors='gray',linewidth=0.9)#,levels=np.logspace(-2,1,7))
    
    #logMBH_cut=logMBH[(logMstel>9.5)]#&(logMBH>6.5)]
    #logMstel=logMstel[(logMstel>9.5)]#&(logMBH>6.5)]
    #logMBH=logMBH_cut
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
    
    from scipy.stats import linregress
    sl,inter,rval,stderr_s,sterr_i=lsqfity(logMstel,logMBH)
    axes.plot([9.5,11.6],inter+sl*np.array([9.5,11.6]),'k',lw=2.5,label='$z=0$ (Tiamat-125-HR)')
    print("corr. coeff. = {}, mass_bulge = {}".format(lsqfity(logMstel,logMBH),mass_bulge))
    ii+=1
  #axes.set_xlim([9.51,11.6])
  #axes.set_ylim([6.5,9.5])


def lsqfity(X, Y):
    """
    Calculate a "MODEL-1" least squares fit.

    The line is fit by MINIMIZING the residuals in Y only.

    The equation of the line is:     Y = my * X + by.

    Equations are from Bevington & Robinson (1992)
    Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
    pp: 104, 108-109, 199.

    Data are input and output as follows:

    my, by, ry, smy, sby = lsqfity(X,Y)
    X     =    x data (vector)
    Y     =    y data (vector)
    my    =    slope
    by    =    y-intercept
    ry    =    correlation coefficient
    smy   =    standard deviation of the slope
    sby   =    standard deviation of the y-intercept

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Determine the size of the vector.
    n = len(X)

    # Calculate the sums.

    Sx = np.sum(X)
    Sy = np.sum(Y)
    Sx2 = np.sum(X ** 2)
    Sxy = np.sum(X * Y)
    Sy2 = np.sum(Y ** 2)

    # Calculate re-used expressions.
    num = n * Sxy - Sx * Sy
    den = n * Sx2 - Sx ** 2

    # Calculate my, by, ry, s2, smy and sby.
    my = num / den
    by = (Sx2 * Sy - Sx * Sxy) / den
    ry = num / (np.sqrt(den) * np.sqrt(n * Sy2 - Sy ** 2))

    diff = Y - by - my * X

    s2 = np.sum(diff * diff) / (n - 2)
    smy = np.sqrt(n * s2 / den)
    sby = np.sqrt(Sx2 * s2 / den)

    return my, by, ry, smy, sby   



if __name__=='__main__':
  ##SETUP
  data_folder='/home/mmarshal/data_dragons/'
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,213:0.55,250:0}
  color={52:'C0',63:'C1',78:'C2',100:'C3',116:'C4',134:'aqua',158:'pink',213:'k',250:'k'}
  meraxes_loc='/output/meraxes.hdf5'

  ##OPTIONS	
  plot_type=1 #Plot MBH vs Mstellar and Mbulge (1) or vs Mbulge for different bulge types (2)
  z0=1 #Plot z=0.55 relation?
  filename_125='paper2_T125'
  mass_bulge=0 #Plot the total stellar mass (0) or bulge stellar mass (1)?
  split=True #Is bulge split into MDBH and IDBH?
  bulge_type=0 #If wanting only one bulge component, which type? 1:IDBH, 2:MDBH
  contour=1 #Plot contours of z=2 distribution?
  
  ##PLOT
  fig,ax = plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0},sharex=True,sharey=True)
  #fig,axes=plt.subplots(1,1)
  plot_type
  if plot_type==1: 
    logMBH=np.log10(0.49)+1.17*np.log10(np.array([10**8.17,10**12])/10**11)+9
    #ax[1].errorbar([8.17,12],logMBH,yerr=0.28,linestyle='--',label='Kormendy \& Ho (2013)\n- Fit',capsize=3,linewidth=2.5, color=colors[1],zorder=0)
    plot_observations(ax[0],color,"bulge")
    plot_observations(ax[1],color,"stellar")
    if z0:
      plot_MBHMstellar(filename_125,[250],True,split,bulge_type,contour,color,ax[0])
      plot_MBHMstellar(filename_125,[250],False,split,bulge_type,contour,color,ax[1])
    ax[1].set_xlabel(r'$\log(M_\ast/M_\odot)$')
    ax[0].set_xlabel(r'$\log(M_{\textrm{bulge}}/M_\odot)$')
  elif plot_type==2: 
    plot_observations(ax[0],color,"bulge")
    plot_observations(ax[1],color,"stellar")
    if z0:
      plot_MBHMstellar(filename_125,[250],True,split,1,contour,color,ax[0])
      plot_MBHMstellar(filename_125,[250],True,split,2,contour,color,ax[1])
    ax[1].set_xlabel(r'$\log(M_{\textrm{merger-driven bulge}})$')
    ax[0].set_xlabel(r'$\log(M_{\textrm{instability-driven bulge}})$')


  lgd=ax[0].legend(fontsize='small',loc='upper center', bbox_to_anchor=(2.35, 0.7))  
  ax[0].set_ylabel(r'$\log(M_{\rm{BH}}/M_\odot)$') 
  if plot_type==1:
    plt.savefig('/home/mmarshal/results/plots/Paper2/MBHMStellarRelation_z0.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  else:
    plt.savefig('/home/mmarshal/results/plots/Paper2/MBHMStellarRelation_BulgeType_z0.pdf', format='pdf',bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

  ##PERFORM STATISTICS
  #find_fit()
