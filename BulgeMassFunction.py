## Compares the distributions/functions of various galaxy properties when altering certain
##parameters in the meraxes model.
# Input: - filenames of the directories which contain the simulations of interest (set of str, any number)
#        - the snapshot to calculate & plot the distributions (int) 

import numpy as np
from dragons import meraxes
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

#Setup
cosmo = {'omega_M_0' : 0.308,
'omega_lambda_0' : 0.692,
'omega_b_0' : 0.04839912,
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}
data_folder='/home/mmarshal/data_dragons/'
meraxes_loc='/output/meraxes.hdf5'

filename=[]
nfiles=len(sys.argv)-2
snapshot=int(sys.argv[-1])
#redshift=bool(sys.argv[-1])
gals={}
for ii in range(1,nfiles+1):
  filename.append(str(sys.argv[ii]))

for ii in range(0,nfiles):
  gals["{}".format(filename[ii])]=meraxes.io.read_gals(data_folder+filename[ii]+meraxes_loc,\
                                          snapshot=snapshot,\
                                          h=cosmo['h'],quiet=True)
  gals["{}".format(filename[ii])] = gals["{}".format(filename[ii])][(gals["{}".format(filename[ii])]["GhostFlag"]==0)]#remove ghosts


#Plot histograms of the BH mass distribution for each of the different models
prop='StellarMass'
maxval2=np.nanmax(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]*1e10))
minval2=np.nanmin(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]*1e10))

for ii in range(0,nfiles):
  #Plot histograms of stellar mass distribution for each of the different models
  prop='StellarMass'
  hist, bin_edges = np.histogram(np.log10(gals["{}".format(filename[ii])][prop]*1e10),range=(minval2,maxval2))
  plt.plot(bin_edges[:-1],hist/(bin_edges[1]-bin_edges[0])/100**3,label='Total Stellar Mass')
  prop='BulgeStellarMass'
  hist, bin_edges = np.histogram(np.log10(gals["{}".format(filename[ii])][prop]*1e10),range=(minval2,maxval2))
  plt.plot(bin_edges[:-1],hist/(bin_edges[1]-bin_edges[0])/100**3,':',label='Bulge Mass')
  hist, bin_edges = np.histogram(np.log10(gals["{}".format(filename[ii])]['StellarMass']*1e10-gals["{}".format(filename[ii])]['BulgeStellarMass']*1e10),range=(minval2,maxval2))
  plt.plot(bin_edges[:-1],hist/(bin_edges[1]-bin_edges[0])/100**3,'--',label='Disk Mass')
plt.title('Snapshot {}'.format(snapshot))
plt.xlabel('log(Mass)')
plt.ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
plt.yscale('log')
plt.legend()

plt.savefig('/home/mmarshal/results/plots/BulgeFunction.pdf',format='pdf')
plt.show()
