import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dragons import meraxes


def load_data(filename,meraxes_loc,snapshot,prop,cosmo):
  gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
      snapshot=snapshot,props=[prop,'GhostFlag'],\
      h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  gals=gals[(gals[prop]*1e10>1e4)]
  return gals


if __name__ == '__main__':  
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
  redshift={52:8,63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  #prop='BlackHoleMass'
  prop='StellarMass' 
  #prop='BulgeStellarMass'
  #prop='Mvir'
  
  filename='tuned_reion'
  filename125='tuned_reion_T125'
 
  fig, axes = plt.subplots(2, 4,gridspec_kw = {'wspace':0, 'hspace':0})
  ii=-1
  j=0
  for snapshot in [52, 63,78,100,116,134,158]:
    ii+=1
    if ii==4:
      j+=1
      ii=0
    gals=load_data(filename,meraxes_loc,snapshot,prop,cosmo)
    gals_125=load_data(filename125,meraxes_loc,snapshot,prop,cosmo)

    #axes[j,ii].hist(np.log10(gals[prop]*1e10),density=True,alpha=0.5)
    #axes[j,ii].hist(np.log10(gals_125[prop]*1e10),density=True,alpha=0.5)
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(4,13),bins=80)
    phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/(100)**3) 
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    axes[j,ii].plot(Max,phi)

    hist, bin_edges = np.histogram(np.log10(gals_125[prop][gals_125[prop]>0]*1e10),range=(4,13),bins=80)
    phi=np.log10(hist/(bin_edges[1]-bin_edges[0])/(125/cosmo['h'])**3) 
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    axes[j,ii].plot(Max,phi)
  plt.savefig('/home/mmarshal/results/plots/BulgeHist_lessseed.pdf',format='pdf')
  plt.show()
