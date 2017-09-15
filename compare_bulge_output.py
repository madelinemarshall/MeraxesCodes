## Compares the distributions/functions of various galaxy properties when altering certain
##parameters in the meraxes model.
# Input: - filenames of the directories which contain the simulations of interest (set of str, any number)
#        - the snapshot to calculate & plot the distributions (int) 

import numpy as np
from dragons import meraxes
import os
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
maxval1=np.nanmax(np.log10(gals["{}".format(filename[0])]['BlackHoleMass'][gals["{}".format(filename[0])]['BlackHoleMass']>0]*1e10))
minval1=np.nanmin(np.log10(gals["{}".format(filename[0])]['BlackHoleMass'][gals["{}".format(filename[0])]['BlackHoleMass']>0]*1e10))
prop='StellarMass'
maxval2=np.nanmax(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]*1e10))
minval2=np.nanmin(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]*1e10))
prop='ColdGas'
maxval3=np.nanmax(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]*1e10))
minval3=np.nanmin(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]*1e10))
prop='Sfr'
maxval4=np.nanmax(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]))
minval4=np.nanmin(np.log10(gals["{}".format(filename[0])][prop][gals["{}".format(filename[0])][prop]>0]))


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
for ii in range(0,nfiles):
  hist, bin_edges = np.histogram(np.log10(gals["{}".format(filename[ii])]['BlackHoleMass']*1e10),range=(minval1,maxval1))
  ax1.plot(bin_edges[:-1],hist/(bin_edges[1]-bin_edges[0])/100**3,label=filename[ii])
  #Plot histograms of stellar mass distribution for each of the different models
  prop='StellarMass'
  hist, bin_edges = np.histogram(np.log10(gals["{}".format(filename[ii])][prop]*1e10),range=(minval2,maxval2))
  ax2.plot(bin_edges[:-1],hist/(bin_edges[1]-bin_edges[0])/100**3,label=filename[ii])
  #Plot histograms of cold gas mass distribution for each of the different models
  prop='ColdGas'
  hist, bin_edges = np.histogram(np.log10(gals["{}".format(filename[ii])][prop]*1e10),range=(minval3,maxval3))
  ax3.plot(bin_edges[:-1],hist/(bin_edges[1]-bin_edges[0])/100**3,label=filename[ii])
  #Plot histograms of SFR distribution for each of the different models
  prop='Sfr'
  hist, bin_edges = np.histogram(np.log10(gals["{}".format(filename[ii])][prop]),range=(minval4,maxval4))
  ax4.plot(bin_edges[:-1],hist/(bin_edges[1]-bin_edges[0])/100**3,label=filename[ii])


ax1.set_title('BH Mass Function')
ax1.set_xlabel('log(Mass)')
ax1.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax1.set_yscale('log')
ax1.legend()

ax2.set_xlabel('log(Mass)')
ax2.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax2.set_title('Stellar Mass Function')
ax2.set_yscale('log')

ax3.set_xlabel('log(Mass)')
ax3.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax3.set_title('Cold Gas Mass Function')
ax3.set_yscale('log')


ax4.set_xlabel('log(SFR)')
ax4.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax4.set_title('SFR Function')
ax4.set_yscale('log')

f.subplots_adjust(hspace=.5)
plt.legend()
plt.suptitle('Snapshot '+str(snapshot))
plt.show()
