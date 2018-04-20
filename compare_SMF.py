import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
sys.path.append('Yuxiang/')
from _plot_obsGSMF import plot_obsGSMF

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
snapshot=int(sys.argv[2])
redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}


#Load data:
#Default parameters
gals_default=meraxes.io.read_gals(data_folder+'bulges_correctBHMF_tiamat125'+meraxes_loc,\
                                        snapshot=snapshot,\
                                        h=cosmo['h'],quiet=True)
gals_default = gals_default[(gals_default["GhostFlag"]==0)]#remove ghosts
size_def=125/cosmo['h']
size_bulges=125/cosmo['h']
filename1=str(sys.argv[1])
#filename2=str(sys.argv[2])

gals_bulges=meraxes.io.read_gals(data_folder+filename1+meraxes_loc,\
                                         snapshot=snapshot,\
                                         h=cosmo['h'],quiet=True)
gals_bulges = gals_bulges[(gals_bulges["GhostFlag"]==0)]#remove ghosts

#gals_bulges_2=meraxes.io.read_gals(data_folder+filename2+meraxes_loc,\
#                                         snapshot=snapshot,\
#                                         h=cosmo['h'],quiet=True)
#gals_bulges_2 = gals_bulges_2[(gals_bulges_2["GhostFlag"]==0)]#remove ghosts

#Plot StellarM Mass Function
prop='StellarMass'
maxval=12.5 
minval=8
hist_default, bin_edges = np.histogram(np.log10(gals_default[prop][gals_default[prop]>0]*1e10),range=(minval,maxval),bins=30)
hist_bulges, bin_edges = np.histogram(np.log10(gals_bulges[prop][gals_bulges[prop]>0]*1e10),range=(minval,maxval),bins=30)
#hist_bulges_2, bin_edges = np.histogram(np.log10(gals_bulges_2[prop][gals_bulges_2[prop]>0]*1e10),range=(minval,maxval),bins=30)

fig,axes=plt.subplots(1,1)
bin_edges=np.array(bin_edges, dtype=np.float128)
Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.

axes.plot(Max,np.log10(hist_default/(bin_edges[1]-bin_edges[0])/size_def**3),label='Original Bulge Model',linewidth=2)
axes.plot(Max,np.log10(hist_bulges/(bin_edges[1]-bin_edges[0])/size_bulges**3),'--',label=filename1,linewidth=2)
plot_obsGSMF(axes,redshift[snapshot],hubble_h=cosmo['h'],markersize=3,legend=False,silent=False,color=[0.5,0.5,0.5],alpha=1.0)
#plt.plot(Max,hist_bulges_2/(bin_edges[1]-bin_edges[0])/100.**3,'--',label=filename2)

plt.xlabel('log(Mass)')
plt.ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
plt.title('Stellar Mass Function')
plt.xlim([8,12.5])
plt.ylim([-5.8,-1.2])
plt.legend()
plt.show()
