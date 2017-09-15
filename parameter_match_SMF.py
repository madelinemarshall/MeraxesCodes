import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys
import pandas as pd

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
#meraxes_loc='/meraxes.hdf5'
snapshot=int(sys.argv[3])
#scalefactor=5
scalefactor=1

#Load data:
#Default parameters
gals_default=meraxes.io.read_gals(data_folder+'default'+meraxes_loc,\
                                        snapshot=snapshot,\
                                        h=cosmo['h'],quiet=True)
gals_default = gals_default[(gals_default["GhostFlag"]==0)]#remove ghosts

filename1=str(sys.argv[1])
filename2=str(sys.argv[2])

gals_bulges=meraxes.io.read_gals(data_folder+filename1+meraxes_loc,\
                                         snapshot=snapshot,\
                                         h=cosmo['h'],quiet=True)
gals_bulges = gals_bulges[(gals_bulges["GhostFlag"]==0)]#remove ghosts

gals_bulges_2=meraxes.io.read_gals(data_folder+filename2+meraxes_loc,\
                                         snapshot=snapshot,\
                                         h=cosmo['h'],quiet=True)
gals_bulges_2 = gals_bulges_2[(gals_bulges_2["GhostFlag"]==0)]#remove ghosts

#Plot StellarM Mass Function
prop='StellarMass'
maxval=np.nanmax(np.log10(gals_default[prop][gals_default[prop]>0]*1e10)) 
minval=np.nanmin(np.log10(gals_default[prop][gals_default[prop]>0]*1e10))
hist_default, bin_edges = np.histogram(np.log10(gals_default[prop][gals_default[prop]>0]*1e10),range=(minval,maxval),bins=30)
hist_bulges, bin_edges = np.histogram(np.log10(gals_bulges[prop][gals_bulges[prop]>0]*1e10),range=(minval,maxval),bins=30)
hist_bulges_2, bin_edges = np.histogram(np.log10(gals_bulges_2[prop][gals_bulges_2[prop]>0]*1e10),range=(minval,maxval),bins=30)

bin_edges=np.array(bin_edges, dtype=np.float128)

Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
plt.plot(Max,np.log10(hist_default/(bin_edges[1]-bin_edges[0])/100.**3),label='Default')
plt.plot(Max,np.log10(hist_bulges*scalefactor/(bin_edges[1]-bin_edges[0])/100.**3),'--',label=filename1)
plt.plot(Max,np.log10(hist_bulges_2*scalefactor/(bin_edges[1]-bin_edges[0])/100.**3),'--',label=filename2)

if snapshot==63: #redshift 7
  dat=pd.read_json('reduced_data/duncan2014_z7.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Duncan+14')
  dat=pd.read_json('reduced_data/song2015_z7.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Song+15')
  dat=pd.read_json('reduced_data/gonzalez2011_z7.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Gonzalez+11')
  dat=pd.read_json('reduced_data/grazian2015_z7.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Grazian+15')
if snapshot==78: #redshift 6
  dat=pd.read_json('reduced_data/duncan2014_z6.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Duncan+14')
  dat=pd.read_json('reduced_data/song2015_z6.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Song+15')
  dat=pd.read_json('reduced_data/gonzalez2011_z6.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Gonzalez+11')
  dat=pd.read_json('reduced_data/grazian2015_z6.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Grazian+15')
  dat=pd.read_json('reduced_data/spitler2017_5.50z6.50.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Spitler+17')
if snapshot==100: #redshift 5
  dat=pd.read_json('reduced_data/song2015_z5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Song+15')
  dat=pd.read_json('reduced_data/caputi2010_4.25z5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Caputi+10')
  dat=pd.read_json('reduced_data/davidzon2017_4.5z5.5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Davidzon+17')
  dat=pd.read_json('reduced_data/duncan2014_z5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Duncan+14')
  dat=pd.read_json('reduced_data/gonzalez2011_z5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Gonzalez+11')
  dat=pd.read_json('reduced_data/grazian2015_z5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Grazian+15')
  dat=pd.read_json('reduced_data/spitler2017_4.50z5.50.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Spitler+17')
if snapshot==115: #redshift 4
  dat=pd.read_json('reduced_data/song2015_z4.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Song+15')
  dat=pd.read_json('reduced_data/muzzin2013_3.0z4.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Muzzin+13')
  dat=pd.read_json('reduced_data/caputi2010_3.5z4.25.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Caputi+10')
  dat=pd.read_json('reduced_data/davidzon2017_3.5z4.5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Davidzon+17')
  dat=pd.read_json('reduced_data/duncan2014_z4.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Duncan+14')
  dat=pd.read_json('reduced_data/grazian2015_z4.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Grazian+15')
  dat=pd.read_json('reduced_data/ilbert2013_3.0z4.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Ilbert+13')
  dat=pd.read_json('reduced_data/marchesini2009_3.0z4.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Marchesini+09')
  dat=pd.read_json('reduced_data/spitler2017_3.50z4.50.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Spitler+17')
if snapshot==134: #redshift 3
  dat=pd.read_json('reduced_data/mortlock2011_z3.25.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Mortlock+11')
  dat=pd.read_json('reduced_data/muzzin2013_3.0z4.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Muzzin+13')
  dat=pd.read_json('reduced_data/muzzin2013_2.5z3.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Muzzin+13')
  dat=pd.read_json('reduced_data/caputi2010_3z3.5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Caputi+10')
  dat=pd.read_json('reduced_data/davidzon2017_3.0z3.5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Davidzon+17')
  dat=pd.read_json('reduced_data/davidzon2017_2.5z3.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Davidzon+17')
  dat=pd.read_json('reduced_data/ilbert2013_3.0z4.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Ilbert+13')
  dat=pd.read_json('reduced_data/ilbert2013_2.5z3.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Ilbert+13')
  dat=pd.read_json('reduced_data/marchesini2009_3.0z4.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Marchesini+09')
  dat=pd.read_json('reduced_data/marchesini2009_2.0z3.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Marchesini+09')
  dat=pd.read_json('reduced_data/spitler2017_2.50z3.00.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Spitler+17')
  dat=pd.read_json('reduced_data/spitler2017_3.00z3.50.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Spitler+17')
if snapshot==158: #redshift 2
  dat=pd.read_json('reduced_data/mortlock2011_z2.25.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Mortlock+11')
  dat=pd.read_json('reduced_data/muzzin2013_1.5z2.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Muzzin+13')
  dat=pd.read_json('reduced_data/muzzin2013_2.0z2.5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Muzzin+13')
  dat=pd.read_json('reduced_data/davidzon2017_1.5z2.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Davidzon+17')
  dat=pd.read_json('reduced_data/davidzon2017_2.0z2.5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Davidzon+17')
  dat=pd.read_json('reduced_data/ilbert2013_1.5z2.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Ilbert+13')
  dat=pd.read_json('reduced_data/ilbert2013_2.0z2.5.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Ilbert+13')
  dat=pd.read_json('reduced_data/marchesini2009_1.3z2.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Marchesini+09')
  dat=pd.read_json('reduced_data/marchesini2009_2.0z3.0.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Marchesini+09')
  dat=pd.read_json('reduced_data/spitler2017_1.50z2.00.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Spitler+17')
  dat=pd.read_json('reduced_data/spitler2017_2.00z2.50.json',orient='split')
  plt.errorbar(dat['logM'],dat['phi'],yerr=[dat['dphi_lo'],dat['dphi_up']],fmt='o',label='Spitler+17')


plt.xlabel('log(Mass)')
plt.ylabel(r'$\log\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
plt.title('Stellar Mass Function, snapshot {}'.format(snapshot))
plt.xlim([7.5,11])
#plt.ylim([10.**-4,10.**-0.8])
plt.legend()
plt.show()

