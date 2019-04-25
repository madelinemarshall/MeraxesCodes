import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import fsolve
from scipy.optimize import root
from matplotlib import rc
rc('text', usetex=True)

cosmo = {'omega_m' : 0.3,
'omega_R' : 0,
'omega_k' : 1 - 0.3 - 0.7,
'omega_lambda' : 0.7,
'h' : 0.7,
'n' : 1,
'sigma_8' : 0.9,
'c' : 3e8,
'G' : 4.302e-9, #Mpc/Msun (km/s)^2
'mu' : 0.59,
'm_p' : 1.6726219e-27, #kg
'k_b' : 1.38064852e-23 #(m/s)^2 kg / K
}

def scale(z):
    return (1+z)**-1

def Hubble(z):
    return 100*cosmo['h']*np.sqrt(cosmo['omega_m']*scale(z)**-3 + \
    cosmo['omega_R']*scale(z)**-4 + cosmo['omega_lambda'])

def dH(z):
    return cosmo['c']/Hubble(z) #kpc

def calc_ang_diam_dist(z):
    return integrate.quad(dH,0,z)[0]*scale(z) #kpc

def spatial_resolution(D,z):
    return 1.22 * 1600*1e-10*(1+z) * calc_ang_diam_dist(z) / D 

zz=np.linspace(5,10,10)
JWST_spatial_res=np.zeros_like(zz)
HST_spatial_res=np.zeros_like(zz)
for ii in range(0,len(zz)):
  JWST_spatial_res[ii]=spatial_resolution(6.5,zz[ii])/2
  HST_spatial_res[ii]=spatial_resolution(2.4,zz[ii])/2
plt.plot(zz,HST_spatial_res)
plt.plot(zz,JWST_spatial_res)
#plt.yscale('log')
#plt.ylim(10**-1.52,10**1)
plt.ylim(0,1.5)

plt.show()
