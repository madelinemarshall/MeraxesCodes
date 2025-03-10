import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
from astrodatapy.number_density import number_density
from scipy.integrate import quad#rature
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  }
colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                  '#ff7f00','#a65628','#f781bf','#98ff98']*4
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','--','-.',':']*4

def plot_convertedBHMF(ax,z,make_legend,labels_dict):
    obs    = number_density(feature='GSMF',z_target=z,quiet=1,h=cosmo['h'])
    xlim    = (7, 13)
    ylim    = (-7, -0.5)
    xlabel  = r"$\log_{10}[M_*/{\rm M_{\odot}}]$"
    ylabel  = r"$\log_{10}[\rm \phi/Mpc^{-3} dex^{-1}]$"

    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Qin2017_Tiamat") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        if label in labels_dict:
            if make_legend:
              continue
            datatype   = labels_dict[label][0]
            color   = labels_dict[label][1]
            marker   = labels_dict[label][2]
            linestyle   = labels_dict[label][3]
        else:
            datatype   = obs.target_observation['DataType'][ii]
            color      = colors[ii]#'gray'
            marker     = markers[j_data]
            linestyle  = linestyles[k_func]
            labels_dict[label]=[datatype,color,marker,linestyle]

        data[:,1:] = np.log10(data[:,1:])
        if datatype == 'data':
            ax.errorbar(data[:,0]-2.5,  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        elif datatype == 'dataULimit':
            ax.errorbar(data[:,0]-2.5,  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        else:
            ax.plot(data[:,0]-2.5,data[:,1],label=label,color=color,linestyle=linestyle,lw=1)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1


def plot_BHMF(ax,z):
    obs    = number_density(feature='BHMF',z_target=z,quiet=1,h=cosmo['h'])
    xlim    = (6, 10)
    ylim    = (-6, -0.5)
    xlabel  = r"$\log_{10}[M_{\rm BH}/{\rm M_{\odot}}]$"
    ylabel  = r"$\log_{10}[\rm \phi/Mpc^{-3} dex^{-1}]$"

    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Shen2012") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        datatype   = obs.target_observation['DataType'][ii]
        color      = 'k'#obs_color[z]
        marker     = markers[j_data]
        linestyle  = linestyles[k_func]
        data[:,1:] = np.log10(data[:,1:])
        if datatype == 'data':
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker,markersize=4)
            j_data +=1
        elif datatype == 'dataULimit':
            ax.errorbar(data[:,0],  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker)
            j_data +=1
        else:
            ax.plot(data[:,0],data[:,1],label=label,color=color,linestyle=linestyle,lw=3)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1


def plot_obsSMF(ax,z,make_legend,labels_dict,data_combined_M,data_combined_phi,data_combined_phi_err):
    obs    = number_density(feature='GSMF',z_target=z,quiet=1,h=cosmo['h'])
    xlim    = (7, 13)
    ylim    = (-7, -0.5)
    xlabel  = r"$\log_{10}[M_*/{\rm M_{\odot}}]$"
    ylabel  = r"$\log_{10}[\rm \phi/Mpc^{-3} dex^{-1}]$"

    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Qin2017_Tiamat") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR"):
            continue
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        if label in labels_dict:
            if make_legend:
              continue
            datatype   = labels_dict[label][0]
            color   = labels_dict[label][1]
            marker   = labels_dict[label][2]
            linestyle   = labels_dict[label][3]
        else:
            datatype   = obs.target_observation['DataType'][ii]
            color      = colors[ii]#'gray'
            marker     = markers[j_data]
            linestyle  = linestyles[k_func]
            labels_dict[label]=[datatype,color,marker,linestyle]

        data[:,1:] = np.log10(data[:,1:])
        if datatype == 'data':
            #ax.errorbar(data[:,0]-2.5,  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
            #            label=label,color=color,fmt=marker,markersize=3,lw=1)
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        elif datatype == 'dataULimit':
            #ax.errorbar(data[:,0]-2.5,  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
            #            label=label,color=color,fmt=marker,markersize=3,lw=1)
            ax.errorbar(data[:,0],  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker,markersize=3,lw=1)
            j_data +=1
        else:
            ax.plot(data[:,0],data[:,1],label=label,color=color,linestyle=linestyle,lw=1)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1
        #data_combined_phi.append(data[:,1])
        #data_combined_M.append(data[:,0])
        #data_combined_phi_err.append(data[:,3]-data[:,1])
        if obs.target_observation.index[ii] != "Thanjavur2016":
          data_combined_phi=np.concatenate((data_combined_phi,data[:,1]))
          data_combined_M=np.concatenate((data_combined_M,data[:,0]))
          data_combined_phi_err=np.concatenate((data_combined_phi_err,data[:,1]-data[:,3]))
    
    return (data_combined_M,data_combined_phi,data_combined_phi_err)


def best_fit_SMF(xdata,ydata,yerr,ax):
    selection=(xdata>8)&(xdata<11.8)
    xdata=xdata[selection]
    ydata=ydata[selection]
    yerr=yerr[selection]
    
    masses=np.linspace(8,12,100)
    p0=[4*10**-3,10**-3,-0.35,-1.47,11] #give M,MBH in log
    popt, pcov = curve_fit(SMF_form, xdata, ydata, p0=p0, sigma=yerr,absolute_sigma=True)
    ax.plot(masses, (SMF_form(masses, *popt)), 'k',linewidth=2)
    #popt, pcov = curve_fit(SMF_form, xdata, ydata, p0=p0)#, sigma=yerr,absolute_sigma=True)
    #ax.plot(masses, (SMF_form(masses, *popt)), 'b--',linewidth=2.5)
    return popt


def plot_Davis(axes):
    davis=pd.read_csv('/home/mmarshal/simulation_codes/data/DavisVika_BHMF.csv',header=None)
    axes.plot(davis[0],davis[1],'o',color='red',label="Davis et al. (2014)",markersize=3)
    axes.fill_between(davis[0],davis[2],davis[3],color='red',alpha=0.2)


def plot_Shankar(axes):
    shan=pd.read_csv('/home/mmarshal/simulation_codes/data/Shankar_BHMF.csv',header=None)
    axes.plot(shan[0],shan[1],'o',color='blue',label="Shankar et al. (2009)",markersize=3)
    axes.fill_between(shan[0],shan[3],shan[2],color='blue',alpha=0.2)


def SMF_GAMA(M):
  #From https://arxiv.org/pdf/1111.5707.pdf (GAMA)
  p1=3.96*10**-3
  p2=0.79*10**-3
  a1=-0.35
  a2=-1.47
  Ms=10.66 #give M,MBH in log
  #h=0.678
  return np.array(SMF_form(M-0.2833,p1,p2,a1,a2,Ms))+np.log10(0.9086)
  #return (p1*(10**M/10**Ms)**a1+p2*(10**M/10**Ms)**a2) * np.exp(-10**M/10**Ms) / 10**Ms * np.log(10)*10**M
 

def SMF_form(M,p1,p2,a1,a2,Ms):
  return np.log10((p1*(10**M/10**Ms)**a1+p2*(10**M/10**Ms)**a2) * np.exp(-10**M/10**Ms) / 10**Ms * np.log(10)*10**M )


def func(M,MBH,SMF_type,popt=None):
  #Black-hole bulge-mass relation from Kormendy & Ho 2013
  sigma=0.28
  med_BH=10**9*0.49*(10**M/10**11)**1.17 
  if SMF_type=='SMF_GAMA':
    stellarMasses=10**SMF_GAMA(M)
  elif SMF_type=='SMF_best':
    stellarMasses=10**SMF_form(M,*popt)
  return 1/(np.sqrt(2*np.pi)*sigma*10**MBH) * stellarMasses * np.exp(-(np.log(10**MBH/med_BH))**2/(2*sigma**2))*np.log(10)*10**MBH


def calc_BHMF(ax,SMF_type,popt=None):
  MBH=np.linspace(5,10,100)
  phi_BH=np.zeros(100)
  for ii in range(0,100):
    phi_BH[ii],val=quad(func,6,16,args=(MBH[ii],SMF_type,popt))
  ax.plot(MBH,np.log10(phi_BH),'k')#to h=0.678 from 0.7 (y)
  np.save('data/BestFitBHMF_M.npy',MBH)
  np.save('data/BestFitBHMF_phi.npy',np.log10(phi_BH))


if __name__=='__main__':
  #fig,ax=plt.subplots(1,1)
  #plot_convertedBHMF(ax,0,0,{})
  #plot_BHMF(ax,0)
  #plt.show()
  data_combined_phi=np.empty(0)
  data_combined_M=np.empty(0)
  data_combined_phi_err=np.empty(0)

  fig,ax=plt.subplots(1,1)
  M=np.linspace(8,12)

  data_combined_M,data_combined_phi,data_combined_phi_err=plot_obsSMF(ax,0,0,{},data_combined_M,data_combined_phi,data_combined_phi_err)
  popt=best_fit_SMF(data_combined_M,data_combined_phi,data_combined_phi_err,ax)

  #ax.plot(M+0.2833,SMF_GAMA(M),'r') ##Chabrier + h=0.7 to 0.678 conversion (x), to h=0.678 from 0.7 (y)
  ax.plot(M,SMF_GAMA(M),'r') ##Chabrier + h=0.7 to 0.678 conversion (x), to h=0.678 from 0.7 (y)
  ax.set_xlim([8,12])
  ax.set_ylim([-5,0])
  ax.set_xlabel(r'$\log(M_\ast / M_\odot)$')
  ax.set_ylabel(r'$\log(\Phi~\rm{ dex}^{-1}~\rm{ Mpc}^{-3})$')
  plt.show()

  fig,ax=plt.subplots(1,1)


  calc_BHMF(ax,'SMF_GAMA')
  calc_BHMF(ax,'SMF_best',popt=popt)
  plot_Davis(ax)
  plot_Shankar(ax)
  #plot_BHMF(ax,0)
  #plt.yscale('log')
  ax.set_xlabel(r'$\log(M_{BH} / M_\odot)$')
  ax.set_ylabel(r'$\log(\Phi~\rm{ dex}^{-1}~\rm{ Mpc}^{-3})$')
  ax.set_xlim(6,10)
  ax.set_ylim(-5,-2)
  plt.show()

