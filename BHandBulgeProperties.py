##Plots BH mass, stellar mass and bulge fraction functions, bulge fraction against BH mass, bulge fraction
##against stellar mass, and histograms of bulge fraction in bins of BH mass. Also performs KS test to see 
##if distribution of bulge fractions is the same for log(MBH)>7 and other mass samples. 
#Input: snapshot

import numpy as np
from dragons import meraxes
import os
import sys
import matplotlib.pyplot as plt
from chi_sq import chi_sq
import scipy.stats as stats

cosmo = {'omega_M_0' : 0.308,
'omega_lambda_0' : 0.692,
'omega_b_0' : 0.04839912,
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}

sim = 'bulges_update1102_full'#0901_ddsf' #str(sys.argv[1])
fmeraxes = '/home/mmarshal/data_dragons/'+sim+'/output/meraxes.hdf5'
snapshot=int(sys.argv[1])
gals=meraxes.io.read_gals(fmeraxes,\
				snapshot=snapshot,\
                                h=cosmo['h'],quiet=True)

##Consider relation between BH mass and bulge fraction at fixed snapshot
SMCut=0 #1e-2 #in units of 1e10 Msol

bulgeFrac=gals['BulgeStellarMass'][gals['StellarMass']>SMCut]/gals['StellarMass'][gals['StellarMass']>SMCut]
BHMass=np.log10(gals['BlackHoleMass'][gals['StellarMass']>SMCut]*1e10)
n_bins=15
med_bulgeFrac=np.zeros(n_bins)
std_bulgeFrac=np.zeros(n_bins)
centres=np.zeros(n_bins)
M_start=3.5
M_end=np.max(BHMass)
bin_width=(M_end-M_start)/n_bins
for ii in range(0,n_bins):
  centres[ii]=M_start+bin_width/2*(2*ii+1)
  med_bulgeFrac[ii]=np.median(bulgeFrac[(BHMass>=centres[ii]-bin_width/2)&(BHMass<centres[ii]+bin_width/2)])
  std_bulgeFrac[ii]=np.std(bulgeFrac[(BHMass>=centres[ii]-bin_width/2)&(BHMass<centres[ii]+bin_width/2)])

plt.plot(BHMass,bulgeFrac,'k.')
plt.plot(centres,med_bulgeFrac,'r-')
plt.plot(centres,med_bulgeFrac+std_bulgeFrac,'r--')
plt.plot(centres,med_bulgeFrac-std_bulgeFrac,'r--')
plt.xlabel("log(BH Mass)")
plt.ylabel("Bulge to Total Stellar Mass Fraction")
plt.ylim(0,1)
plt.xlim(3,9.5)
plt.show()


##Consider relation between stellar mass and bulge fraction at fixed snapshot
SM=np.log10(gals['StellarMass']*1e10)[gals['StellarMass']>SMCut]
n_bins=15
med_bulgeFrac_SM=np.zeros(n_bins)
std_bulgeFrac_SM=np.zeros(n_bins)
centres_SM=np.zeros(n_bins)
SM_start=4
SM_end=np.max(SM)
bin_width_SM=(SM_end-SM_start)/n_bins
for ii in range(0,n_bins):
  centres_SM[ii]=SM_start+bin_width_SM/2*(2*ii+1)
  med_bulgeFrac_SM[ii]=np.median(bulgeFrac[(SM>=centres_SM[ii]-bin_width_SM/2)&(SM<centres_SM[ii]+bin_width_SM/2)])
  
  std_bulgeFrac_SM[ii]=np.std(bulgeFrac[(SM>=centres_SM[ii]-bin_width_SM/2)&(SM<centres_SM[ii]+bin_width_SM/2)])

plt.plot(SM,bulgeFrac,'k.')
plt.plot(centres_SM,med_bulgeFrac_SM,'r-')
plt.plot(centres_SM,med_bulgeFrac_SM+std_bulgeFrac_SM,'r--')
plt.plot(centres_SM,med_bulgeFrac_SM-std_bulgeFrac_SM,'r--')
plt.xlabel("log(Stellar Mass)")
plt.ylabel("Bulge to Total Stellar Mass Fraction")
plt.ylim(0,1)
plt.xlim(4,12)
plt.show()

##Consider distribution of bulge fractions for galaxies with BH in specific bins
#Plot histograms and make KS test to see if distributions are statistically different
#hist_tot,bin_edges=np.histogram(bulgeFrac[BHMass>4],range=(0,1),bins=10)
#hist_45,bin_edges=np.histogram(bulgeFrac[(BHMass>=4)&(BHMass<5)],range=(0,1),bins=10)
#hist_56,bin_edges=np.histogram(bulgeFrac[(BHMass>=5)&(BHMass<6)],range=(0,1),bins=10)
#hist_67,bin_edges=np.histogram(bulgeFrac[(BHMass>=6)&(BHMass<7)],range=(0,1),bins=10)
#hist_g7,bin_edges=np.histogram(bulgeFrac[(BHMass>=7)],range=(0,1),bins=10)

hist_tot,bin_edges=np.histogram(bulgeFrac[BHMass>4],range=(0,1),bins=10)
hist_45,bin_edges=np.histogram(bulgeFrac[(BHMass>=4)&(BHMass<M_end-3)],range=(0,1),bins=10)
hist_56,bin_edges=np.histogram(bulgeFrac[(BHMass>=M_end-3)&(BHMass<M_end-2)],range=(0,1),bins=10)
hist_67,bin_edges=np.histogram(bulgeFrac[(BHMass>=M_end-2)&(BHMass<M_end-1)],range=(0,1),bins=10)
hist_g7,bin_edges=np.histogram(bulgeFrac[(BHMass>=M_end-1)],range=(0,1),bins=10)

bin_centres=bin_edges[1:]-(bin_edges[1]-bin_edges[0])/2

plt.plot(bin_centres,hist_tot/len(bulgeFrac[BHMass>4]),'-',label="All greater than 4")
plt.plot(bin_centres,hist_45/len(bulgeFrac[(BHMass>=4)&(BHMass<M_end-3)]),'-',label="4 - {:.2f}".format(M_end-3))
plt.plot(bin_centres,hist_56/len(bulgeFrac[(BHMass>=M_end-3)&(BHMass<M_end-2)]),'-',label="{:.2f} - {:.2f}".format(M_end-3,M_end-2))
plt.plot(bin_centres,hist_67/len(bulgeFrac[(BHMass>=M_end-2)&(BHMass<M_end-1)]),'-',label="{:.2f} - {:.2f}".format(M_end-2,M_end-1))
plt.plot(bin_centres,hist_g7/len(bulgeFrac[(BHMass>=M_end-1)]),'-',label="Greater than {:.2f}".format(M_end-1))
plt.xlabel("Bulge to Total Stellar Mass Fraction")
plt.ylabel("Proportion of Galaxies in the Mass Range Found in Bin")
plt.title("Histogram of Bulge Fraction for Different BH Mass Bins, Snapshot {}".format(snapshot))
plt.legend()
plt.show()


#KS_g7_47,p_g7_47=stats.ks_2samp(bulgeFrac[(BHMass>4)&(BHMass<7)],bulgeFrac[(BHMass>=7)])
#KS_g7_57,p_g7_57=stats.ks_2samp(bulgeFrac[(BHMass>5)&(BHMass<7)],bulgeFrac[(BHMass>=7)])
#KS_g7_46,p_g7_46=stats.ks_2samp(bulgeFrac[BHMass>=7],bulgeFrac[(BHMass>4)&(BHMass<6)])
KS_g7_47,p_g7_47=stats.ks_2samp(bulgeFrac[(BHMass>4)&(BHMass<M_end-1)],bulgeFrac[(BHMass>=M_end-1)])
KS_g7_57,p_g7_57=stats.ks_2samp(bulgeFrac[(BHMass>5)&(BHMass<M_end-1)],bulgeFrac[(BHMass>=M_end-1)])
KS_g7_46,p_g7_46=stats.ks_2samp(bulgeFrac[BHMass>=M_end-1],bulgeFrac[(BHMass>4)&(BHMass<M_end-2)])
print(p_g7_47,p_g7_57,p_g7_46,'M_end=',M_end)

##Plot BH and Stellar Mass Functions
f, (ax1, ax2, ax3) = plt.subplots(1, 3)
maxval=np.nanmax((BHMass[BHMass>0]))
minval=np.nanmin((BHMass[BHMass>0]))
hist_default, bin_edges = np.histogram((BHMass[BHMass>0]),range=(minval,maxval),bins=30)
bin_edges=np.array(bin_edges, dtype=np.float128)
Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
ax1.plot(Max,hist_default/(bin_edges[1]-bin_edges[0])/100.**3)
ax1.set_xlabel('log(BH / M$_\odot$))')
ax1.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax1.set_title('BH Mass Function')
ax1.set_yscale('log')
ax1.set_xlim(3,9.5)

maxval=np.nanmax(SM[SM>0])
minval=np.nanmin(SM[SM>0])
hist_default, bin_edges = np.histogram(SM[SM>0],range=(minval,maxval),bins=30)
bin_edges=np.array(bin_edges, dtype=np.float128)
Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
ax2.plot(Max,hist_default/(bin_edges[1]-bin_edges[0])/100.**3)
ax2.set_xlabel('log(Stellar Mass / M$_\odot$)')
ax2.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax2.set_title('Stellar Mass Function')
ax2.set_yscale('log')
ax2.set_xlim(0,12)

hist,edges=np.histogram(bulgeFrac,range=(0,1),bins=30)
edges=np.array(edges, dtype=np.float128)
Max=edges[0:-1] + (edges[1]-edges[0])/2.
ax3.set_yscale('log')
ax3.set_xlabel('Bulge to Total Stellar Mass Fraction')
ax3.set_ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
ax3.set_title('Bulge Fraction Function')
ax3.plot(Max,hist/(edges[1]-edges[0])/100.**3)
ax3.set_xlim(0,1)

plt.suptitle('Snapshot {}'.format(snapshot))
plt.show()

