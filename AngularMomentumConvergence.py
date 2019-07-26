#Does a quick comparison with the Okumara+17 j/m for stellar disks.
#Ours isn't crazy (0.5 relative to their 0.77)
import numpy as np
from dragons import meraxes
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/mmarshal/simulation_codes')

#Sets plot defaults
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (7.6,4)
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

def load_data(filename,snapshot,cosmo):
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=['StellarMass','GhostFlag','BlackHoleMass','BulgeStellarMass','AMstars','Spin','Vvir','Rvir','Mvir'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals


def plot_hist2d(xdata,ydata,axes,xlims,ylims):
  H, xedges, yedges, img=axes.hist2d(xdata, ydata, bins=40, range=[xlims,ylims], weights=None, cmin=1, cmax=None, data=None,cmap='Blues',norm=matplotlib.colors.LogNorm())
  axes.set_ylim(ylims)
  axes.set_xlim(xlims)
  plt.tight_layout()


def plot_avg(xdata,ydata,axes,xlims,bin_width,redshiftPlot=False,T125=0):
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
    if np.size(y)>=50:
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
  
  if redshiftPlot:
    #axes.fill_between(bin_centre,pct_r_16,pct_r_84,color=color[snapshot],alpha=0.2)
    if T125:
      axes.plot(bin_centre,avg_r,'--',color=color[snapshot],label=r'$z={}$'.format(redshift[snapshot]))
    else:
      axes.plot(bin_centre,avg_r,color=color[snapshot],label=r'$z={}$'.format(redshift[snapshot]))
  else:
    axes.errorbar(bin_centre,avg_r,yerr=np.array([avg_r-pct_r_16,pct_r_84-avg_r]),color='k',marker='s',markersize=4,label='M19 - Median')
  #return bin_centre,avg_r,np.array([avg_r-pct_r_16,pct_r_84-avg_r])



def plot_y_vs_x(y,x,axes,T125=0,range=None,redshiftAxes=None):
  if not range:
    range=([np.log10(min(x)),np.log10(max(x))],[(min(y)),(max(y))])
  xlims=range[0]
  ylims=range[1]
  axes.plot(ylims,ylims,'--',color=[0.5,0.5,0.5])
  selection=(y>0)&(x>0)
  #plot_hist2d(np.log10(x[selection]),(y[selection]),axes,xlims,ylims)
  plot_avg(np.log10(x[selection]),y[selection],axes,xlims,0.25)
    
  
  if redshiftAxes:
    plot_avg(np.log10(x[selection]),(y[selection]),redshiftAxes,xlims,0.25,redshiftPlot=True,T125=T125)
    redshiftAxes.set_ylim(ylims)
    redshiftAxes.set_xlim(xlims)


def plot_AMF(data,boxwidth,axes,**kwargs):
    maxval=np.nanmax(np.log10(data[data>0])) 
    minval=np.nanmin(np.log10(data[data>0]))
    hist, bin_edges = np.histogram(np.log10(data[data>0]),range=(minval,maxval),bins=60)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    #Max=Max[hist>5]
    #hist=hist[hist>5]
    #hist_plus=hist+np.sqrt(hist)
    #hist_minus=hist-np.sqrt(hist)
    phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    #phi_plus=np.log10(hist_plus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    #phi_minus=np.log10(hist_minus/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    axes.plot(Max,phi,**kwargs)
    #axes.fill_between(Max, phi_plus, phi_minus, alpha=0.5,color=kwargs['color'],label='__nolegend__')
    axes.set_ylabel(r'$\phi$')
    return axes


if __name__=='__main__':
  filename='paper1'#str(sys.argv[1])#'bulges_update0915_ddsf'
  filename_T125='paper1_T125'
  
  fig1,axes1=plt.subplots(3,3,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1]})
  fig2,axes2=plt.subplots(3,3,gridspec_kw = {'wspace':0, 'hspace':0,'height_ratios':[3,3,1]})
  fig4,axes4=plt.subplots(1,3)#,sharey=True,gridspec_kw={'wspace':0})
  fig5,axes5=plt.subplots(1,1)#,sharey=True,gridspec_kw={'wspace':0})
  fig6,axes6=plt.subplots(1,1)#,sharey=True,gridspec_kw={'wspace':0})
  fig7,axes7=plt.subplots(1,1)#,sharey=True,gridspec_kw={'wspace':0})
  
  redshift={37:10,43:9,52:8,63:7,78:6,100:5,116:4,134:3,158:2,192:1,250:0}
  color={63:'#e41a1c',78:'#377eb8',100:'#4daf4a',116:'#984ea3',\
                  134:'#ff7f00',158:'#f781bf',192:'#a65628',250:'black'}
 
 
  ii=-1
  j=0
  #for snapshot in np.flip([63,78,100,116,134,158],0):
  for snapshot in np.flip([100,116,134,158,250],0):
    ii+=1
    if ii==3:
      j+=1
      ii=0
    
    if snapshot<160: 
      gals=load_data(filename,snapshot,cosmo)
      AM=np.sqrt(gals['AMstars'][:,0]**2+gals['AMstars'][:,1]**2+gals['AMstars'][:,2]**2)
      diskMass=gals['StellarMass']-gals['BulgeStellarMass']
      AMhalo=gals['Spin']*np.sqrt(2)*gals['Vvir']*gals['Rvir']
      plot_y_vs_x(np.log10((AMhalo)),gals['Mvir'],axes1[j,ii],T125=0,redshiftAxes=axes4[2],range=([0,3],[0,3]))
      plot_y_vs_x(np.log10((AMhalo)/(gals['Mvir'])),gals['Mvir'],axes1[j,ii],T125=0,redshiftAxes=axes4[0],range=([0,3],[0,0.7]))
      plot_AMF(AMhalo,100,axes5,**{'color':color[snapshot],'label':'T z={}'.format(redshift[snapshot])})
      plot_AMF(AM,100,axes6,**{'color':color[snapshot],'label':'T z={}'.format(redshift[snapshot])})
      plot_AMF(AM[diskMass>0]/diskMass[diskMass>0]*1e3,100,axes7,**{'color':color[snapshot],'label':'T z={}'.format(redshift[snapshot])})
    

    gals_T125=load_data(filename_T125,snapshot,cosmo)
  
    AM_T125=np.sqrt(gals_T125['AMstars'][:,0]**2+gals_T125['AMstars'][:,1]**2+gals_T125['AMstars'][:,2]**2)
    diskMass_T125=gals_T125['StellarMass']-gals_T125['BulgeStellarMass']
    AMhalo_T125=gals_T125['Spin']*np.sqrt(2)*gals_T125['Vvir']*gals_T125['Rvir']

    plot_y_vs_x(np.log10((AMhalo_T125)/(gals_T125['Mvir'])),gals_T125['Mvir'],axes1[j,ii],T125=1,redshiftAxes=axes4[0],range=([0,3],[0,0.7]))
    plot_y_vs_x(np.log10((AMhalo_T125)),gals_T125['Mvir'],axes1[j,ii],T125=1,redshiftAxes=axes4[2],range=([0,3],[0,3]))
    
    plot_AMF(AMhalo_T125,125/cosmo['h'],axes5,**{'color':color[snapshot],'label':'T125 z={}'.format(redshift[snapshot]),'linestyle':'--'})
    plot_AMF(AM_T125,125/cosmo['h'],axes6,**{'color':color[snapshot],'label':'T125 z={}'.format(redshift[snapshot]),'linestyle':'--'})
    plot_AMF(AM_T125[diskMass_T125>0]/diskMass_T125[diskMass_T125>0]*1e3,125/cosmo['h'],axes7,**{'color':color[snapshot],'label':'T125 z={}'.format(redshift[snapshot]),'linestyle':'--'})

  axes5.set_xlabel(r'Halo j')
  axes5.legend()
  axes6.set_xlabel(r'Galaxy J (km/s Mpc)')
  axes6.legend()
  axes7.legend()
  axes7.set_xlabel(r'Galaxy j (km/s kpc)')
  plt.show()
