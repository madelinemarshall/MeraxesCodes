import numpy as np
from dragons import meraxes
import os
#import matplotlib
import matplotlib.pyplot as plt
import sys
import pandas as pd
from astrodatapy.number_density import number_density
import magcalc as mc

#Sets plot defaults
import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.2,3.4)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
color_maps     = ['Reds', 'Blues', 'Greens'] *4
markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
linestyles     = ['-','--','-.',':']*4


def load_data(filename,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)]
  return gals


def load_mags(filename,snapshot):
  redshift={52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  #if (filename=='bulges_update1102_full')&(snapshot==158):
  #  dust=pd.read_csv('/home/mmarshal/PhD/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'_dust.txt',sep=' ')
  #  i775=np.array(dust)[:,0]
  #  m1600=np.array(dust)[:,1]
  #  return m1600
  #else:
  MUV=pd.read_hdf('/home/mmarshal/results/mags_output/'+filename+'/mags_6_'+format(snapshot,'03d')+'.hdf5')['M1600-100']
  AUV=mc.reddening(1600,MUV,redshift[snapshot])
  MUV_dust=MUV+AUV
  return MUV_dust


def plot_LF(lums,axes,**kwargs):
    maxval=-7
    minval=-23.5
    hist, bin_edges = np.histogram(lums,range=(minval,maxval),bins=40)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    Max=Max[hist>5]
    hist=hist[hist>5]
    hist_plus=hist+np.sqrt(hist)
    hist_minus=hist-np.sqrt(hist)
    phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/100**3)
    phi_plus=np.log10(hist_plus/(bin_edges[1]-bin_edges[0])/100**3)
    phi_minus=np.log10(hist_minus/(bin_edges[1]-bin_edges[0])/100**3)
    axes.plot(Max,phi,**kwargs)
    axes.fill_between(Max, phi_plus, phi_minus, alpha=0.5,color=kwargs['color'],label='__nolegend__')
    #axes.plot(Max,np.log10(hist/(bin_edges[1]-bin_edges[0])/100.**3),**kwargs)
    return axes


def plot_obsGLF(ax,z,make_legend,labels_dict):
    colors         = ['#984ea3','#98ff98','#ff7f00','#a65628','#f781bf','#377eb8','#4daf4a']*4
    feature = 'GLF_UV'
    xlim    = (-15,-24)
    ylim    = (-7.5,-1)
    xlabel  = r"$M_{1600}$"
    ylabel  = r"$\log_{10}[\rm \phi/Mpc^{-3} dex^{-1}]$"

    obs    = number_density(feature=feature,z_target=z,quiet=1,h=cosmo['h'])
    j_data = 0
    k_func = 0
    for ii in range(obs.n_target_observation):
        if (obs.target_observation.index[ii] == "Qin2017_Tiamat") or (obs.target_observation.index[ii] == "Qin2017_Tiamat125_HR") or (obs.target_observation.index[ii] == 'Bouwens2015_ulimit'):
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
            ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3],data[:,2]- data[:,1]],\
                        label=label,color=color,fmt=marker,markersize=3)
            j_data +=1
        elif datatype == 'dataULimit':
            ax.errorbar(data[:,0],  data[:,1], yerr = -0.2*data[:,1], uplims=True,\
                        label=label,color=color,fmt=marker,markersize=3)
            j_data +=1
        else:
            ax.plot(data[:,0],data[:,1],label=label,color=color,linestyle=linestyle,lw=3)
            ax.fill_between(data[:,0], data[:,2],data[:,3],color=color,alpha=0.5)
            k_func +=1

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    #leg = ax.legend(loc='lower left')
    #leg.get_frame().set_alpha(0.5)


def make_legend(ax):
    plot_LF(mags_bulges,ax,**{'linestyle':'-','label':'Modified Meraxes','linewidth':1,'color':'k'})
    plot_LF(mags_default,ax,**{'linestyle':'-','label':'Default Meraxes','linewidth':1,'color':[0.35,0.35,0.35]})
    labels_dict={}
    for snapshot in [52,63,78,100,115]:
      plot_obsGLF(ax,redshift[snapshot],1,labels_dict)
    ax.legend(fontsize='small',ncol=3,loc=(-0.1,0.2))
    ax.axis('off')
    ax.set_xlim(8,8.1)
    ax.set_ylim(-6,-5.8)
    #ax.set_visible(False)


if __name__=="__main__":
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
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  redshift={52:8,63:7,78:6,100:5,115:4,134:3,158:2}
  col={63:'turquoise',78:'r',100:'b',115:'g',134:'purple',158:'orange'} 
  filename='tuned_reion'

  fig,axes=plt.subplots(2,5,gridspec_kw = {'wspace':0,'height_ratios':[4,1]})

  ii=-1
  j=0
  labels_dict={}
  for snapshot in [52,63,78,100,115]:
    ii+=1
  #  if ii==3:
  #    j+=1
  #    ii=0
    mags_default=load_mags('default_reion',snapshot)
    mags_bulges=load_mags(filename,snapshot)
    #mags_2=load_mags(filename2,snapshot)

    plot_obsGLF(axes[j,ii],redshift[snapshot],0,labels_dict)
    #plot_LF(mags_2,axes[ii],'Bulge Model - Croton SF',':')
    plot_LF(mags_bulges,axes[j,ii],**{'linestyle':'-','label':'Modified Meraxes','linewidth':1,'color':'k','zorder':100})
    plot_LF(mags_default,axes[j,ii],**{'linestyle':'-','label':'Default Meraxes','linewidth':1,'color':[0.35,0.35,0.35],'zorder':101})

    #plot_obsGLF(axes[ii],redshift[snapshot],hubble_h=cosmo['h'],markersize=3,legend=True,silent=True,color=[0.5,0.5,0.5],alpha=1.0)

    if ii==0:
      axes[j,ii].set_ylabel(r'$\log(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
    else:
      axes[j,ii].set_yticklabels([])
    axes[j,ii].set_xlabel(r'$M_{1600}$')
    axes[j,ii].set_xlim([-23.2,-15.2])
    axes[j,ii].set_ylim([-5.5,-0.5])
    #axes[ii].set_title('z={}'.format(redshift[snapshot]))
    axes[j,ii].text(-18.55, -4.9, r'$z={}$'.format(redshift[snapshot]),weight='bold',size='large')
    axes[j,ii].grid(color=[0.8,0.8,0.8],linestyle='--') 
  plt.tight_layout()
  #fig.subplots_adjust(right=0.9,bottom=0.2)
  make_legend(axes[1,1])
  axes[1,0].axis('off')
  axes[1,1].axis('off')
  axes[1,2].axis('off')
  axes[1,3].axis('off')
  axes[1,4].axis('off')

  #axes[2].set_ylim([-7,-1])

  #axes[4].legend(['Bulge Model','Default Meraxes'])#,'Bulge Model - Croton SF'])
  plt.savefig('/home/mmarshal/results/plots/GLF.pdf', format='pdf')#,bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.show()

